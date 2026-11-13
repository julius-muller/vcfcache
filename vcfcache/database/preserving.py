"""Input-preserving annotation with exact biallelic cache matches.

The sample is the immutable backbone. Cache construction may still split
alleles; sample annotation never does. The cache answers only the records it
can answer exactly; every other record reaches the recipe unchanged, as direct
annotation would see it. An unusual input therefore costs hits, not correctness.

The cached path runs coordinate-disjoint shards concurrently and publishes
their outputs by BGZF block copy, so no record is decoded or recompressed after
it has been annotated. Cache access is restricted to windows around the input
when the cache is much larger than the input, and to whole contigs otherwise.
"""

import errno
import gzip
import hashlib
import os
import re
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile
import threading
import time
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HIT = "VCFCACHE_INTERNAL_HIT"
KEY = "%CHROM\t%POS\t%REF\t%ALT\n"
WINDOW = 16384
# Window lookups pay per-window index seeks; they win once the cache holds far
# more records than the input (measured: wins at 250x, loses at 4x).
WINDOW_RATIO = 20
# The input is cut into one shard per core, of equal record count. Measured on
# a genome and an exome (vep-rs core, 8 cores): 8 shards were as fast as or
# faster than 16 or 24, and faster than 4. Small inputs use fewer, larger shards.
SHARDS_PER_THREAD = 1
MIN_SHARD = 5000
# Every SAMPLE-th record's position is kept to place the cuts.
SAMPLE = 1024
UNBOUNDED = 2**31 - 1
IDENTITY = "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%FORMAT\n"


def count(bcftools, path):
    """Fail closed instead of treating an unreadable index as zero misses."""
    value = subprocess.check_output([bcftools, "index", "-n", str(path)], text=True)
    try:
        result = int(value.strip())
    except ValueError as error:
        raise RuntimeError(f"Cannot count indexed records in {path}") from error
    if result < 0:
        raise RuntimeError(f"Invalid record count in {path}")
    return result


def identity(bcftools, path):
    """Hash unchanged columns and all FORMAT data without transferring INFO text.

    Visit contigs in name order so different numeric header IDs cannot create
    false differences. No normalization, allele reordering or field exclusions.
    """
    samples = subprocess.check_output([bcftools, "query", "-l", str(path)])
    index = subprocess.check_output(
        [bcftools, "index", "--stats", str(path)], text=True
    )
    contigs = sorted(line.split("\t")[0] for line in index.splitlines() if line)
    command = [bcftools, "query", "-f", IDENTITY]
    if contigs:
        command += ["-r", ",".join(contigs)]
    command.append(str(path))
    # A file-backed stderr avoids a pipe deadlock on malformed inputs.
    with tempfile.TemporaryFile() as errors:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=errors)
        assert process.stdout is not None
        with process.stdout:
            digest = hashlib.file_digest(process.stdout, "sha256").hexdigest()
        if process.wait():
            errors.seek(0)
            raise RuntimeError(errors.read().decode(errors="replace"))
    return samples, digest


def validate_recipe(bcftools, source, output):
    """Reject truncation and altered keys/genotypes before output is published."""
    expected, actual = count(bcftools, source), count(bcftools, output)
    if actual != expected:
        raise RuntimeError(
            f"Annotation record count changed: {expected} input, {actual} output; "
            "recipes must retain all records, including unannotated records"
        )
    with ThreadPoolExecutor(max_workers=2) as pool:
        before = pool.submit(identity, bcftools, source)
        after = pool.submit(identity, bcftools, output)
        if before.result() != after.result():
            raise RuntimeError(
                "Annotation recipe changed variant identity, record order, "
                "samples or FORMAT data; refusing to publish"
            )
    return actual


def index_stats(bcftools, path):
    """Return (contig, length, records) of occupied contigs from the index alone."""
    listing = subprocess.check_output(
        [bcftools, "index", "--stats", str(path)], text=True
    )
    occupied = []
    for line in listing.splitlines():
        fields = line.split("\t")
        if len(fields) >= 3 and fields[2].isdigit() and int(fields[2]):
            length = int(fields[1]) if fields[1].isdigit() else UNBOUNDED
            occupied.append((fields[0], length, int(fields[2])))
    return occupied


# One awk process per contig does the per-record work at C speed and in
# parallel; Python only reads its few result lines. Within a locus, a key seen
# twice marks the position. Keys compare without case, as bcftools compares
# alleles.
SCAN = r"""
{
  if ($1 "\t" $2 != locus) { delete seen; locus = $1 "\t" $2 }
  key = toupper($3 "\t" $4)
  if ((key in seen) && !(locus in dup)) { dup[locus] = 1; print "D\t" $2 "\t" $1 }
  seen[key] = 1
  if ((NR - 1) % sample == 0) print "S\t" (NR - 1) "\t" $2
  if (windowed) { b = int(($2 - 1) / window); if (!(b in block)) { block[b] = 1; print "B\t" b } }
}
END { print "N\t" NR }
"""


def scan_contig(bcftools, source, contig, windowed):
    """Read one contig's variant keys once.

    bcftools pairs records one-to-one: of two identical input keys only the
    first would take a cached annotation, and the recipe result of the second
    would later be paired with the first. Every record at such a position is
    therefore annotated by the recipe; all other records still use the cache.

    Returns:
        dict with the record count, (index, POS) of every SAMPLE-th record,
        positions holding duplicate keys, and occupied windows when windowed.
    """
    query = shlex.join([bcftools, "query", "-f", KEY, "-r", contig, str(source)])
    awk = shlex.join(
        ["awk", "-F", "\t", "-v", f"sample={SAMPLE}", "-v", f"window={WINDOW}"]
        + ["-v", f"windowed={int(windowed)}", SCAN]
    )
    result = subprocess.run(
        ["bash", "-c", f"set -o pipefail; {query} | {awk}"],
        capture_output=True,
        check=False,
    )
    if result.returncode:
        raise RuntimeError(result.stderr.decode(errors="replace"))
    samples, duplicates, blocks, records = [], [], set(), None
    for line in result.stdout.split(b"\n"):
        fields = line.split(b"\t")
        if fields[0] == b"S":
            samples.append((int(fields[1]), int(fields[2])))
        elif fields[0] == b"D":
            duplicates.append(int(fields[1]))
        elif fields[0] == b"B":
            blocks.add(int(fields[1]))
        elif fields[0] == b"N":
            records = int(fields[1])
    if records is None:
        raise RuntimeError(f"Key scan of {contig} in {source} returned no count")
    return dict(
        records=records,
        samples=samples,
        duplicates=duplicates,
        blocks=blocks if windowed else None,
    )


def plan_ranges(order, scans, parts, span=True):
    """Cut the input into about equal record counts, between positions.

    A cut falls before a sampled record's POS, so records sharing a position
    always stay in one shard. With ``span`` a shard may cover several contigs,
    read in the order listed; without it every contig starts a new shard.

    Returns:
        Shards, each a list of (contig, first POS, last POS) in file order.
    """
    total = sum(scans[c]["records"] for c in order)
    if not total:
        return []
    parts = max(1, min(parts, total))
    cuts = set()
    offset = 0
    targets = [round(k * total / parts) for k in range(1, parts)]
    for contig in order:
        records = scans[contig]["records"]
        for target in targets:
            if offset <= target < offset + records:
                local = target - offset
                chosen = None
                for index, position in scans[contig]["samples"]:
                    if index > local:
                        break
                    chosen = (index, position)
                if chosen is None or chosen[0] == 0:
                    cuts.add((contig, 1))
                else:
                    cuts.add((contig, chosen[1]))
        offset += records
    shards: list = [[]]
    for contig in order:
        starts = sorted(p for c, p in cuts if c == contig)
        # bcftools concat -a, used at duplicate positions, orders contigs by
        # header; a file in another order keeps each contig in its own shard.
        if shards[-1] and (not span or (starts and starts[0] == 1)):
            shards.append([])
        bounds = [p for p in starts if p > 1]
        first = 1
        for position in bounds:
            shards[-1].append((contig, first, position - 1))
            shards.append([])
            first = position
        shards[-1].append((contig, first, UNBOUNDED))
    return [shard for shard in shards if shard]


def shard_rows(shard, scans):
    """BED rows of one shard: whole ranges, or its windows clipped to them."""
    rows = []
    for contig, first, last in shard:
        blocks = scans[contig]["blocks"]
        if blocks is None:
            rows.append(f"{contig}\t{first - 1}\t{last}")
            continue
        previous = None
        for block in sorted(blocks):
            start = max(block * WINDOW, first - 1)
            end = min((block + 1) * WINDOW, last)
            if start >= end:
                continue
            if previous and start <= previous[1]:
                previous[1] = max(previous[1], end)
                continue
            if previous:
                rows.append(f"{contig}\t{previous[0]}\t{previous[1]}")
            previous = [start, end]
        if previous:
            rows.append(f"{contig}\t{previous[0]}\t{previous[1]}")
    return rows


def shard_duplicates(shard, scans):
    """Return CHROM-tab-POS lines of the duplicate positions in one shard."""
    return [
        f"{contig}\t{position}"
        for contig, first, last in shard
        for position in scans[contig]["duplicates"]
        if first <= position <= last
    ]


def contig_order(bcftools, source, stats):
    """Occupied contigs in the order of the header's contig lines.

    This is the order bcftools uses when it merges by position, and the
    fallback output order. The index's own order is not used: tabix numbers
    contigs in file order.
    """
    occupied = {contig for contig, _, _ in stats}
    header = subprocess.check_output([bcftools, "view", "-h", str(source)], text=True)
    declared = [
        line.split("ID=", 1)[1].split(",")[0].rstrip(">")
        for line in header.splitlines()
        if line.startswith("##contig=<")
    ]
    order = [contig for contig in declared if contig in occupied]
    # Contigs used but not declared follow in index order.
    order += [c for c, _, _ in stats if c not in set(order)]
    return order


def read_index(bcftools, source, stats):
    """Contigs in file order with the file offsets of their 16 kb windows.

    Reads the input's CSI or TBI index only, no record. Each reference's
    pseudo-bin gives where its records start in the file and how many there
    are; the deepest bins give where each window's records start. The counts
    must agree with `bcftools index --stats`, otherwise None is returned and
    the caller falls back to header order and a key scan.

    Returns:
        [(contig, records, first offset, end offset, {window: offset})] in
        file order, offsets in compressed bytes; or None.
    """
    path = next(
        (
            Path(f"{source}{s}")
            for s in (".csi", ".tbi")
            if Path(f"{source}{s}").exists()
        ),
        None,
    )
    if path is None:
        return None
    try:
        data = gzip.open(path).read()
        offset = 4

        def take(fmt):
            nonlocal offset
            values = struct.unpack_from(fmt, data, offset)
            offset += struct.calcsize(fmt)
            return values

        names = None
        tbi = data[:4] == b"TBI\x01"
        if tbi:
            (n_ref,) = take("<i")
            *_, l_nm = take("<7i")
            names = data[offset : offset + l_nm].split(b"\0")[:-1]
            offset += l_nm
            min_shift, depth = 14, 5
        elif data[:4] == b"CSI\x01":
            min_shift, depth, l_aux = take("<3i")
            if l_aux >= 28:
                l_nm = struct.unpack_from("<i", data, offset + 24)[0]
                names = data[offset + 28 : offset + 28 + l_nm].split(b"\0")[:-1]
            offset += l_aux
            (n_ref,) = take("<i")
        else:
            return None
        if names is None:
            header = subprocess.check_output(
                [bcftools, "view", "-h", str(source)], text=True
            )
            names = [
                line.split("ID=", 1)[1].split(",")[0].rstrip(">").encode()
                for line in header.splitlines()
                if line.startswith("##contig=<")
            ]
        meta_bin = ((1 << ((depth + 1) * 3)) - 1) // 7 + 1
        first_leaf = ((1 << (depth * 3)) - 1) // 7
        head = struct.Struct("<Ii" if tbi else "<IQi")
        pair = struct.Struct("<QQ")
        refs = []
        for tid in range(n_ref):
            (n_bin,) = take("<i")
            meta, windows = None, {}
            for _ in range(n_bin):
                fields = head.unpack_from(data, offset)
                offset += head.size
                bin_, n_chunk = fields[0], fields[-1]
                if bin_ == meta_bin:
                    meta = [
                        pair.unpack_from(data, offset + 16 * i) for i in range(n_chunk)
                    ]
                elif bin_ >= first_leaf and n_chunk:
                    # Chunks are sorted, so the first starts the window.
                    windows[bin_ - first_leaf] = pair.unpack_from(data, offset)[0] >> 16
                offset += 16 * n_chunk
            if tbi:
                (n_intv,) = take("<i")
                offset += 8 * n_intv
            if meta and len(meta) == 2 and meta[1][0]:
                refs.append(
                    (
                        names[tid].decode(),
                        meta[1][0],
                        meta[0][0] >> 16,
                        meta[0][1] >> 16,
                        windows,
                    )
                )
    except (OSError, struct.error, IndexError, EOFError):
        return None
    expected = {c: n for c, _, n in stats}
    if min_shift != 14 or {c: n for c, n, *_ in refs} != expected:
        return None
    return sorted(refs, key=lambda ref: ref[2])


def plan_from_index(refs, parts, span=True):
    """Cut the input at window starts into parts of about equal file bytes.

    Records are compressed alike, so equal bytes are about equal records. A cut
    falls before a window's first position, never inside a position.
    """
    if not refs:
        return []
    start, end = refs[0][2], max(r[3] for r in refs)
    marks = [
        (offset, contig, window)
        for contig, _, _, _, windows in refs
        for window, offset in sorted(windows.items())
    ]
    cuts = set()
    for k in range(1, parts):
        target = start + k * (end - start) / parts
        hit = next((m for m in marks if m[0] >= target), None)
        if hit:
            cuts.add((hit[1], hit[2] * WINDOW + 1))
    order = [r[0] for r in refs]
    first_window = {c: min(w) * WINDOW + 1 if w else 1 for c, _, _, _, w in refs}
    cuts = {(c, 1 if p <= first_window[c] else p) for c, p in cuts}
    shards: list = [[]]
    for contig in order:
        starts = sorted(p for c, p in cuts if c == contig)
        # bcftools concat -a, used at duplicate positions, orders contigs by
        # header; a file in another order keeps each contig in its own shard.
        if shards[-1] and (not span or (starts and starts[0] == 1)):
            shards.append([])
        first = 1
        for position in (p for p in starts if p > 1):
            shards[-1].append((contig, first, position - 1))
            shards.append([])
            first = position
        shards[-1].append((contig, first, UNBOUNDED))
    return [shard for shard in shards if shard]


def duplicate_scan(bcftools, source, bed):
    """Start a key scan of one shard; its result is read with scan_result."""
    query = shlex.join(
        [bcftools, "query", "-f", KEY, "-R", str(bed), "--regions-overlap", "0"]
        + [str(source)]
    )
    awk = shlex.join(
        ["awk", "-F", "\t", "-v", f"sample={UNBOUNDED}", "-v", f"window={WINDOW}"]
        + ["-v", "windowed=0", SCAN]
    )
    return subprocess.Popen(
        ["bash", "-c", f"set -o pipefail; {query} | {awk}"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def scan_result(process):
    """Return the duplicate positions (CHROM-tab-POS) and record count of a scan."""
    out, err = process.communicate()
    if process.returncode:
        raise RuntimeError(err.decode(errors="replace"))
    rows, records = [], None
    for line in out.split(b"\n"):
        fields = line.split(b"\t")
        if fields[0] == b"D":
            rows.append(f"{fields[2].decode()}\t{fields[1].decode()}")
        elif fields[0] == b"N":
            records = int(fields[1])
    if records is None:
        raise RuntimeError("Shard key scan returned no count")
    return rows, records


BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def bgzf_blocks(data):
    """(offset, compressed size, uncompressed size) of every complete block."""
    blocks, offset = [], 0
    while offset + 18 <= len(data):
        size = struct.unpack_from("<H", data, offset + 16)[0] + 1
        if offset + size > len(data):
            break
        blocks.append(
            (offset, size, struct.unpack_from("<I", data, offset + size - 4)[0])
        )
        offset += size
    return blocks


def bgzf_compress(data, level=6, threads=1):
    """BGZF-compress bytes, ending with the EOF block, as htslib expects."""

    def block(start):
        chunk = data[start : start + 0xFF00]
        deflate = zlib.compressobj(level, zlib.DEFLATED, -15)
        body = deflate.compress(chunk) + deflate.flush()
        return b"".join(
            (
                b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00",
                struct.pack("<H", len(body) + 25),
                body,
                struct.pack("<II", zlib.crc32(chunk), len(chunk)),
            )
        )

    starts = range(0, len(data), 0xFF00)
    # zlib releases the GIL, so blocks compress in parallel.
    with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        return b"".join(pool.map(block, starts)) + BGZF_EOF


# One virtual-offset field moved by one compressed byte, as a little-endian
# 64-bit lane: (1 << 16) is the compressed-offset unit of a virtual offset.
LANE = (1 << 16).to_bytes(8, "little")


class Csi:
    """A CSI index kept as raw bytes, with every virtual-offset field marked.

    refs[tid] = (start, end, bins): the byte range of that reference's bin
    records and {bin: (record offset, n_chunk)}. mask has LANE at every bin
    loffset and chunk offset (never at the meta bin), so moving all offsets
    of a byte range by delta compressed bytes is a single integer addition.
    """

    def __init__(self, path):
        data = gzip.open(path).read()
        if data[:4] != b"CSI\x01":
            raise ValueError(f"Not a CSI index: {path}")
        self.min_shift, self.depth, l_aux = struct.unpack_from("<3i", data, 4)
        self.aux = data[16 : 16 + l_aux]
        self.meta = ((1 << ((self.depth + 1) * 3)) - 1) // 7 + 1
        offset = 16 + l_aux
        (n_ref,) = struct.unpack_from("<i", data, offset)
        offset += 4
        mask = bytearray(len(data))
        unpack = struct.Struct("<IQi").unpack_from
        self.refs = []
        for _ in range(n_ref):
            (n_bin,) = struct.unpack_from("<i", data, offset)
            offset += 4
            start = offset
            bins = {}
            for _ in range(n_bin):
                bin_, _, n_chunk = unpack(data, offset)
                bins[bin_] = (offset, n_chunk)
                if bin_ != self.meta:
                    mask[offset + 4 : offset + 12] = LANE
                    mask[offset + 16 : offset + 16 + 16 * n_chunk] = LANE * (
                        2 * n_chunk
                    )
                offset += 16 + 16 * n_chunk
            self.refs.append((start, offset, bins))
        if offset > len(data):
            raise ValueError(f"Truncated CSI index: {path}")
        self.n_no_coor = (
            struct.unpack_from("<Q", data, offset)[0]
            if offset + 8 <= len(data)
            else None
        )
        self.data, self.mask = data, bytes(mask)

    def moved(self, start, end, delta):
        """Bytes [start, end) with every virtual offset moved by delta bytes."""
        value = int.from_bytes(self.data[start:end], "little")
        lanes = int.from_bytes(self.mask[start:end], "little")
        # Each lane holds a value below 2**64 before and after the move, so
        # no carry or borrow crosses a lane boundary.
        return (value + delta * lanes).to_bytes(end - start, "little")

    def record(self, bin_):
        """(loffset, [(beg, end), ...]) of one bin of a reference's bins."""
        offset, n_chunk = bin_
        loffset = struct.unpack_from("<Q", self.data, offset + 4)[0]
        chunks = list(
            struct.iter_unpack(
                "<QQ", self.data[offset + 16 : offset + 16 + 16 * n_chunk]
            )
        )
        return loffset, chunks


def copy_range(source, target, start, length):
    """Append bytes [start, start + length) of one open file to another."""
    while length > 0:
        try:
            copied = os.copy_file_range(source.fileno(), target.fileno(), length, start)
        except (AttributeError, OSError):
            copied = None
        if copied is None:
            # No kernel copy here: copy through a buffer instead.
            source.seek(start)
            while length:
                block = source.read(min(length, 1 << 24))
                if not block:
                    raise RuntimeError("Shard ended early")
                target.write(block)
                length -= len(block)
            return
        if copied == 0:
            raise RuntimeError("Shard ended early")
        start += copied
        length -= copied


def bcf_header(path):
    """Header blocks of a bcftools BCF and its decompressed text, or None.

    bcftools writes the header into its own BGZF blocks, so records start on
    a block boundary; anything else is reported as None.
    """
    with open(path, "rb") as stream:
        data = stream.read(1 << 20)
        while True:
            text, seen, end, length = b"", 0, 0, None
            complete = False
            for offset, size, usize in bgzf_blocks(data):
                if offset + size > len(data):
                    break
                text += zlib.decompress(data[offset + 18 : offset + size - 8], -15)
                seen += usize
                end = offset + size
                if length is None and len(text) >= 9:
                    length = 9 + struct.unpack_from("<I", text, 5)[0]
                if length is not None and seen >= length:
                    complete = True
                    break
            if complete:
                break
            more = stream.read(1 << 22)
            if not more:
                return None
            data += more
    if text[:5] != b"BCF\x02\x02" or seen != length:
        return None
    return data[:end], text


def gather(parts, staged, threads=1):
    """Join shard BCFs by block copy and join their CSI indexes arithmetically.

    Each shard was written with its own index; the joined index moves every
    file offset by where that shard's records land in the joined file, so the
    output is never read again. Returns False, writing nothing, whenever the
    shards do not have the expected shape; the caller then uses bcftools.
    """
    headers, indexes = [], []
    for part in parts:
        with open(part, "rb") as stream:
            stream.seek(-len(BGZF_EOF), os.SEEK_END)
            if stream.read() != BGZF_EOF:
                return False
        header = bcf_header(part)
        if header is None or (headers and header[1] != headers[0][1]):
            return False
        headers.append(header)
        indexes.append(Csi(Path(f"{part}.csi")))
    first = indexes[0]
    if len({(i.min_shift, i.depth, i.aux, len(i.refs)) for i in indexes}) != 1:
        return False
    spans, deltas = [], []
    base = len(headers[0][0])
    for part, (raw, _) in zip(parts, headers, strict=True):
        end = len(raw)
        length = part.stat().st_size - end - len(BGZF_EOF)
        if length < 0:
            return False
        spans.append((end, length))
        # Moving a block by delta bytes adds delta to the compressed part of
        # every virtual offset in it.
        deltas.append(base - end)
        base += length

    def copy():
        with open(staged, "wb") as out:
            out.write(headers[0][0])
            out.flush()
            for part, (end, length) in zip(parts, spans, strict=True):
                with open(part, "rb") as source:
                    copy_range(source, out, end, length)
            out.write(BGZF_EOF)

    pack = struct.Struct("<IQi").pack
    index = [
        b"CSI\x01",
        struct.pack("<3i", first.min_shift, first.depth, len(first.aux)),
        first.aux,
        struct.pack("<i", len(first.refs)),
    ]
    with ThreadPoolExecutor(max_workers=1) as pool:
        copying = pool.submit(copy)
        for tid in range(len(first.refs)):
            holders = [
                (i, d) for i, d in zip(indexes, deltas, strict=True) if i.refs[tid][2]
            ]
            # A reference split between shards shares only the bins around
            # each cut, one per level at most, plus the meta bin. Those are
            # merged here; all other bin records are moved as raw bytes.
            shared, seen = {first.meta}, set()
            if len(holders) > 1:
                for i, _ in holders:
                    keys = i.refs[tid][2].keys()
                    shared.update(keys & seen)
                    seen.update(keys)
            pieces, merged, n_bin = [], {}, 0
            for i, delta in holders:
                start, end, bins = i.refs[tid]
                cuts = sorted((bins[b], b) for b in shared if b in bins)
                n_bin += len(bins) - len(cuts)
                position = start
                for (offset, n_chunk), bin_ in cuts:
                    if offset > position:
                        pieces.append(i.moved(position, offset, delta))
                    position = offset + 16 + 16 * n_chunk
                    loffset, chunks = i.record((offset, n_chunk))
                    shift = delta << 16
                    target = merged.get(bin_)
                    if bin_ == first.meta:
                        (beg, fin), counts = chunks
                        if target is None:
                            merged[bin_] = [
                                loffset,
                                [(beg + shift, fin + shift), counts],
                            ]
                        else:
                            (begin, _), (mapped, unmapped) = target[1]
                            target[1] = [
                                (begin, fin + shift),
                                (mapped + counts[0], unmapped + counts[1]),
                            ]
                        continue
                    moved = [(a + shift, b + shift) for a, b in chunks]
                    if target is None:
                        merged[bin_] = [loffset + shift, moved]
                    else:
                        target[0] = min(target[0], loffset + shift)
                        target[1].extend(moved)
                if end > position:
                    pieces.append(i.moved(position, end, delta))
            index.append(struct.pack("<i", n_bin + len(merged)))
            index += pieces
            for bin_, (loffset, chunks) in merged.items():
                index.append(pack(bin_, loffset, len(chunks)))
                index += [struct.pack("<QQ", a, b) for a, b in chunks]
        if any(i.n_no_coor is not None for i in indexes):
            index.append(struct.pack("<Q", sum(i.n_no_coor or 0 for i in indexes)))
        compressed = bgzf_compress(b"".join(index), threads=threads)
        copying.result()
    # Written after the data, so the index is never older than the file.
    Path(f"{staged}.csi").write_bytes(compressed)
    return True


# A recipe that reads its input once, as `bcftools view ${INPUT_BCF}`, can
# take the misses from a pipe; it then starts while the cache lookup runs.
STREAM = re.compile(r"view\s+\\?\$\{INPUT_BCF\}")


def streamable(recipe):
    """Whether the recipe reads ${INPUT_BCF} exactly once, through bcftools view."""
    return recipe.count("INPUT_BCF") == 1 and bool(STREAM.search(recipe))


def start_recipe(manager, work, recipe_input, annotated):
    """Start the recipe in the background, reading from ``recipe_input``."""
    auxiliary = manager.output_dir / "auxiliary"
    auxiliary.mkdir(exist_ok=True)
    recipe = manager._substitute_variables(
        manager.nfa_config_content["annotation_cmd"],
        extra_vars={
            "INPUT_BCF": str(recipe_input),
            "OUTPUT_BCF": str(annotated),
            "AUXILIARY_DIR": str(auxiliary),
        },
    )
    (work / "recipe.sh").write_text(f"#!/bin/bash\nset -euo pipefail\n\n{recipe}\n")
    manager.logger.debug(f"Command (started early): {recipe}")
    logs = [(work / "stdout.txt").open("wb"), (work / "stderr.txt").open("wb")]
    process = subprocess.Popen(
        ["bash", "-c", "set -euo pipefail\n" + recipe],
        cwd=work,
        stdout=logs[0],
        stderr=logs[1],
    )
    return dict(
        process=process,
        recipe=recipe,
        logs=logs,
        auxiliary=auxiliary,
        started=time.monotonic(),
    )


def finish_recipe(manager, work, started):
    """Wait for a recipe started by start_recipe; return its exit status."""
    status = started["process"].wait()
    for log in started["logs"]:
        log.close()
    manager.logger.info(
        f"Recipe finished {time.monotonic() - started['started']:.3f}s after its "
        f"early start, exit status {status}"
    )
    for filename, output in (
        ("stdout.txt", "annotation_tool.log"),
        ("stderr.txt", "annotation_tool_err.log"),
    ):
        path = work / filename
        if path.exists() and path.stat().st_size:
            shutil.copyfile(path, manager.output_dir / output)
    auxiliary = started["auxiliary"]
    if auxiliary.exists() and not any(auxiliary.iterdir()):
        auxiliary.rmdir()
    return status


def run_recipe(manager, work, command, recipe_input, annotated):
    """Run the annotation recipe once; its output is checked by check_recipe."""
    auxiliary = manager.output_dir / "auxiliary"
    auxiliary.mkdir(exist_ok=True)
    recipe = manager._substitute_variables(
        manager.nfa_config_content["annotation_cmd"],
        extra_vars={
            "INPUT_BCF": str(recipe_input),
            "OUTPUT_BCF": str(annotated),
            "AUXILIARY_DIR": str(auxiliary),
        },
    )
    result = command(recipe)
    for filename, output in (
        ("stdout.txt", "annotation_tool.log"),
        ("stderr.txt", "annotation_tool_err.log"),
    ):
        path = work / filename
        if path.exists() and path.stat().st_size:
            shutil.copyfile(path, manager.output_dir / output)
    if not any(auxiliary.iterdir()):
        auxiliary.rmdir()
    return result


def check_recipe(manager, work, recipe_input, annotated, tag, bcftools):
    """Refuse recipe output that lacks the tag or altered any record."""
    try:
        validation_start = time.monotonic()
        manager._validate_info_tag(annotated, tag)
        validate_recipe(bcftools, recipe_input, annotated)
        manager.logger.info(
            f"Recipe identity/FORMAT validation completed in {time.monotonic() - validation_start:.3f}s"
        )
    except Exception:
        if annotated.exists():
            shutil.copyfile(annotated, work / "failed_recipe.bcf")
            if Path(str(annotated) + ".csi").exists():
                shutil.copyfile(
                    Path(str(annotated) + ".csi"), work / "failed_recipe.bcf.csi"
                )
        raise


def run(manager, cache=None):
    """Run the recipe once, on the full input or the union of cache misses."""
    from vcfcache.database.workflow_manager import BcftoolsCommand

    work = manager.work_dir / ("annotate" if cache else "annotate-nocache")
    work.mkdir(parents=True, exist_ok=True)
    bcftools = str(manager.params_file_content["bcftools_cmd"])
    threads = manager.params_file_content.get("threads", 1)
    if threads < 1:
        raise ValueError("threads must be positive")
    level = getattr(manager, "output_level", 1)
    if level not in range(10):
        raise ValueError("Output compression level must be between 0 and 9")
    tool = shlex.quote(bcftools)
    source = Path(manager.input_file)
    requested_output = manager._resolve_output_bcf(source.stem)
    if requested_output != "-":
        destination = Path(requested_output)
        if destination.resolve() == source.resolve() or (
            destination.exists() and destination.samefile(source)
        ):
            raise ValueError("Output must not overwrite the input file")
    total = count(bcftools, source)
    tag = manager.nfa_config_content["must_contain_info_tag"]

    def command(text, folder=work):
        return BcftoolsCommand(text, manager.logger, folder).run()

    def quoted(path):
        return shlex.quote(str(path))

    # Stage on the destination filesystem so failed validation leaves an old
    # deliverable untouched. The original input is never opened for writing.
    stdout = requested_output == "-"
    target = work / "stdout.bcf" if stdout else Path(requested_output)
    target.parent.mkdir(parents=True, exist_ok=True)
    stats = {}
    missing = total
    with tempfile.TemporaryDirectory(
        prefix=".vcfcache-publish-", dir=target.parent
    ) as temporary:
        staged = Path(temporary) / "output.bcf"
        if not cache:
            manager.logger.info(f"Annotating {total} unchanged records")
            result = run_recipe(manager, work, command, source, staged)
            check_recipe(manager, work, source, staged, tag, bcftools)
        else:
            headers = [
                subprocess.check_output([bcftools, "view", "-h", str(path)], text=True)
                for path in (source, cache)
            ]
            for path, header in zip((source, cache), headers, strict=True):
                if f"##INFO=<ID={HIT}," in header:
                    raise ValueError(f"Reserved INFO/{HIT} already exists in {path}")
            # An input that already declares the tag may carry annotation the
            # cache did not supply; only then must hits be marked explicitly.
            # Otherwise an empty tag identifies a miss, which saves the marker
            # and its removal on every record.
            marked = f"##INFO=<ID={tag}," in headers[0]
            windowed = bool(total) and count(bcftools, cache) > WINDOW_RATIO * total
            # One parallel pass over the input's keys places equal-record cuts
            # and finds duplicate positions and, for small inputs, windows.
            parts = max(
                1, min(int(SHARDS_PER_THREAD * threads), -(-total // MIN_SHARD))
            )
            plan_start = time.monotonic()
            # Large inputs are cut from the index alone, in file order, and
            # each shard's duplicate scan runs beside its cache lookup. Small
            # inputs need their windows first, so they are scanned up front.
            stats_ = index_stats(bcftools, source)
            refs = read_index(bcftools, source, stats_) if total else None
            header = contig_order(bcftools, source, stats_)
            # File order, from the index; header order when it cannot be read.
            order = [r[0] for r in refs] if refs else header
            span = order == header
            scans = None
            if refs and not windowed and parts > 1:
                shards = plan_from_index(refs, parts, span)
                planning = "index"
            else:
                with ThreadPoolExecutor(max_workers=threads) as pool:
                    scans = dict(
                        zip(
                            order,
                            pool.map(
                                lambda c: scan_contig(bcftools, source, c, windowed),
                                order,
                            ),
                            strict=True,
                        )
                    )
                if sum(s["records"] for s in scans.values()) != total:
                    raise RuntimeError("Key scan read a different number of records")
                shards = plan_ranges(order, scans, parts, span)
                planning = "scan"
            manager.logger.info(
                f"Planned {len(shards)} shard(s) from the {planning} in "
                f"{time.monotonic() - plan_start:.3f}s"
            )
            # Shards share the cores; a lone shard gets all of them.
            shard_threads = max(1, threads // max(1, len(shards)))
            folders = []
            for number in range(max(1, len(shards))):
                folder = work / f"shard-{number}"
                folder.mkdir(exist_ok=True)
                folders.append(folder)
            regions = [""] * len(folders)
            # Positions with duplicate keys, per shard; see scan_contig.
            duplicates: list = (
                [shard_duplicates(shard, scans) for shard in shards]
                if scans is not None
                else [[] for _ in shards]
            ) or [[]]
            scanned = [0] * len(folders)
            mark = f"-m +{HIT} " if marked else ""
            strip = f"-x INFO/{HIT} " if marked else ""
            annotated = work / "recipe.bcf"
            # A recipe that streams its input starts now: it loads its own
            # data while the cache lookup runs, then reads each shard's misses
            # from a pipe as soon as that shard and those before it are done.
            stream = bool(shards) and streamable(
                manager.nfa_config_content["annotation_cmd"]
            )
            ready = [threading.Event() for _ in folders]
            failed = threading.Event()
            feed_errors: list = []
            started = None
            fifo = work / "misses.fifo"

            def feed():
                """Write every shard's misses, in order, into the recipe's pipe."""
                try:
                    descriptor = None
                    while descriptor is None:
                        try:
                            descriptor = os.open(fifo, os.O_WRONLY | os.O_NONBLOCK)
                        except OSError as error:
                            if error.errno != errno.ENXIO:
                                raise
                            if started["process"].poll() is not None or failed.is_set():
                                return
                            time.sleep(0.02)
                    os.set_blocking(descriptor, True)
                    with os.fdopen(descriptor, "wb") as pipe:
                        for number, folder in enumerate(folders):
                            ready[number].wait()
                            if failed.is_set():
                                return
                            subprocess.run(
                                [bcftools, "view", "--no-version"]
                                + (["-H"] if number else [])
                                + [str(folder / "misses.bcf")],
                                stdout=pipe,
                                check=True,
                            )
                except BaseException as error:  # reported after the lookup
                    feed_errors.append(error)

            if stream:
                fifo.unlink(missing_ok=True)
                os.mkfifo(fifo)
                started = start_recipe(manager, work, fifo, annotated)
                feeder = threading.Thread(target=feed, daemon=True)
                feeder.start()

            def lookup(number):
                """Overlay exact cache hits on one shard; signal its misses."""
                try:
                    return select(number)
                except BaseException:
                    failed.set()
                    raise
                finally:
                    ready[number].set()

            def select(number):
                """Overlay exact cache hits on one shard and select its misses."""
                folder = folders[number]
                scan = None
                if shards:
                    bed = folder / "regions.bed"
                    if scans is None:
                        rows = [f"{c}\t{a - 1}\t{b}" for c, a, b in shards[number]]
                    else:
                        rows = shard_rows(shards[number], scans)
                    bed.write_text("".join(f"{row}\n" for row in rows))
                    regions[number] = f"-R {quoted(bed)} --regions-overlap 0 "
                    if scans is None:
                        scan = duplicate_scan(bcftools, source, bed)
                # Exact compares allele SETS, not ALT order. Route every
                # multiallelic record to the recipe rather than reusing
                # allele-indexed annotation values.
                command(
                    f"{tool} annotate --no-version {regions[number]}--pair-logic exact "
                    f"-i 'N_ALT=1' -k {mark}-a {quoted(cache)} "
                    f"{quoted(source)} -c INFO -o overlay.bcf -Ob1 -W "
                    f"--threads {shard_threads}",
                    folder,
                )
                if scan:
                    duplicates[number], scanned[number] = scan_result(scan)
                # Without duplicates the misses are written directly; with
                # them, the misses elsewhere and every record at a duplicate
                # position (taken from the input, never from the overlay) are
                # merged by position.
                dups = folder / "duplicates.tsv"
                if duplicates[number]:
                    dups.write_text("".join(f"{row}\n" for row in duplicates[number]))
                elsewhere = f"-T ^{quoted(dups)} " if duplicates[number] else ""
                target = "elsewhere.bcf -W" if duplicates[number] else "misses.bcf"
                if marked:
                    command(
                        f"{tool} view --no-version {elsewhere}"
                        f"-i 'INFO/{HIT}!=1 || INFO/{tag}=\"\"' overlay.bcf -Ou | "
                        f"{tool} annotate --no-version {strip}-o {target} -Ob",
                        folder,
                    )
                else:
                    command(
                        f"{tool} view --no-version {elsewhere}-i 'INFO/{tag}=\"\"' "
                        f"overlay.bcf -o {target} -Ob --threads {shard_threads}",
                        folder,
                    )
                if duplicates[number]:
                    command(
                        f"{tool} view --no-version -R {quoted(dups)} --regions-overlap 0 "
                        f"{quoted(source)} -o duplicated.bcf -Ob -W --threads 1 && "
                        f"{tool} concat --no-version -a elsewhere.bcf duplicated.bcf "
                        "-o misses.bcf -Ob --threads 1",
                        folder,
                    )
                return count(bcftools, folder / "overlay.bcf")

            manager.logger.info(
                f"Step 1/4: Exact biallelic cache overlay in {len(folders)} shard(s), "
                f"cache access by {'windows' if windowed else 'contigs'}; input unchanged"
            )
            try:
                with ThreadPoolExecutor(max_workers=threads) as pool:
                    if sum(pool.map(lookup, range(len(folders)))) != total:
                        raise RuntimeError("Cache overlay changed input record count")
            except BaseException:
                failed.set()
                if started:
                    started["process"].kill()
                    finish_recipe(manager, work, started)
                raise
            if scans is None and shards and sum(scanned) != total:
                raise RuntimeError("Shard key scans read a different number of records")
            sites = sum(len(d) for d in duplicates)
            if sites:
                manager.logger.warning(
                    f"Duplicate variant keys at {sites} position(s): every record "
                    "there is annotated by the recipe"
                )
            manager.logger.info("Step 2/4: Select actual misses, not pre-existing tags")
            recipe_input = work / "misses.bcf"
            parts = " ".join(quoted(folder / "misses.bcf") for folder in folders)
            command(
                f"{tool} concat --no-version --naive {parts} -o {quoted(recipe_input)} && "
                f"{tool} index -f {quoted(recipe_input)}"
            )
            missing = count(bcftools, recipe_input)
            # The recipe output is checked while the shards are merged; nothing
            # is published before that check has passed.
            checker = ThreadPoolExecutor(max_workers=1)
            checked = None
            if started:
                manager.logger.info(
                    f"Step 3/4: Annotating {missing} unchanged records (started early)"
                )
                feeder.join()
                status = finish_recipe(manager, work, started)
                fifo.unlink(missing_ok=True)
                if feed_errors:
                    raise RuntimeError(
                        f"Streaming misses to the recipe failed: {feed_errors[0]}"
                    )
                if missing and status:
                    raise subprocess.CalledProcessError(status, started["recipe"])
            elif missing:
                manager.logger.info(f"Step 3/4: Annotating {missing} unchanged records")
                run_recipe(manager, work, command, recipe_input, annotated)
            if missing:
                checked = checker.submit(
                    check_recipe, manager, work, recipe_input, annotated, tag, bcftools
                )
            # At a duplicate position the overlay may hold a cached value that
            # the merge only replaces where the recipe wrote the tag. Should
            # the recipe leave such a record without it, merge onto a base
            # that holds the input's own records there instead.
            clean = False
            if sites and missing:
                every = work / "duplicates.tsv"
                every.write_text("".join(f"{row}\n" for d in duplicates for row in d))
                if checked:
                    checked.result()
                untagged = subprocess.run(
                    [bcftools, "view", "-H", "-R", str(every), "--regions-overlap", "0"]
                    + ["-i", f'INFO/{tag}=""', str(annotated)],
                    capture_output=True,
                    check=True,
                ).stdout
                clean = bool(untagged.strip())
                if clean:
                    manager.logger.info(
                        "Recipe left records at duplicate positions without "
                        f"INFO/{tag}; merging onto the input's own records there"
                    )
            manager.logger.info("Step 4/4: Overlay annotations; retain every record")

            def merge(number):
                """Add recipe annotations to one shard and drop any marker."""
                overlay = (
                    f"{regions[number]}--pair-logic exact -a {quoted(annotated)} -c INFO "
                    if missing
                    else ""
                )
                base = "overlay.bcf"
                if clean and duplicates[number]:
                    base = "base.bcf"
                    command(
                        f"{tool} view --no-version -T ^duplicates.tsv overlay.bcf "
                        "-o kept.bcf -Ob1 -W --threads 1 && "
                        f"{tool} concat --no-version -a kept.bcf duplicated.bcf "
                        "-o base.bcf -Ob1 -W --threads 1",
                        folders[number],
                    )
                if overlay or strip:
                    command(
                        f"{tool} annotate --no-version {overlay}{strip}{base} "
                        f"-o final.bcf -Ob{level} -W --threads {shard_threads}",
                        folders[number],
                    )
                else:
                    command(
                        f"{tool} view --no-version overlay.bcf -o final.bcf "
                        f"-Ob{level} -W --threads {shard_threads}",
                        folders[number],
                    )

            try:
                with ThreadPoolExecutor(max_workers=threads) as pool:
                    list(pool.map(merge, range(len(folders))))
            finally:
                # A recipe that altered records explains a failed merge best.
                checker.shutdown(wait=True)
                if checked:
                    checked.result()
            # Gather by BGZF block copy: shards are coordinate-disjoint and in
            # file order, so nothing is decoded or recompressed here. Their
            # indexes, written while each shard was compressed, are joined too.
            gather_start = time.monotonic()
            try:
                joined = gather(
                    [folder / "final.bcf" for folder in folders], staged, threads
                )
            except (
                OSError,
                ValueError,
                OverflowError,
                struct.error,
                zlib.error,
            ) as error:
                manager.logger.warning(
                    f"Joining shard indexes failed ({error}); indexing"
                )
                joined = False
            if joined:
                result = subprocess.CompletedProcess(["gather"], 0)
            else:
                parts = " ".join(quoted(folder / "final.bcf") for folder in folders)
                result = command(
                    f"{tool} concat --no-version --naive {parts} -o {quoted(staged)} && "
                    f"{tool} index -f --threads {threads} {quoted(staged)}"
                )
            manager.logger.info(
                f"Gathered {len(folders)} shard(s) "
                f"{'with joined indexes' if joined else 'and indexed the output'} "
                f"in {time.monotonic() - gather_start:.3f}s"
            )
            stats = {
                "shards": len(folders),
                "duplicate_sites": sites,
                "joined_index": joined,
                "recipe_streamed": stream,
                "planning": planning if shards else None,
                "duplicate_clean_base": clean,
                "cache_access": "windows" if windowed else "contigs",
                "output_level": level,
                "hit_marker": marked,
            }
        if count(bcftools, staged) != total:
            raise RuntimeError("Final output record count differs from input")
        if stdout:
            with staged.open("rb") as stream:
                shutil.copyfileobj(stream, sys.stdout.buffer)
        else:
            os.replace(staged, target)
            os.replace(Path(str(staged) + ".csi"), Path(str(target) + ".csi"))
    manager.last_run_stats = {
        "input_variants": total,
        "output_variants": total,
        "missing_variants": missing if cache else None,
        "missing_annotated": missing if cache else None,
        "dropped_variants": 0,
        "input_preprocessing": "none",
        "cache_matching": "biallelic-exact" if cache else None,
        **stats,
    }
    return result
