# SPDX-License-Identifier: GPL-3.0-only
"""Bounded, advisory resident-buffer prediction; never changes record semantics."""

import time

import pysam


def predict_buffer(input_file, cache, records, *, max_seconds=0.5):
    """Sample at most eight records per contig/128 total; abstain if insufficient.

    The estimate intentionally includes decoded payload/allocator allowance. It
    predicts the native record buffer, not RSS or a guaranteed resource bound.
    Any indexed-read problem abstains and leaves the existing native path intact.
    """
    start = time.monotonic()
    sampled, matches, total = 0, 0, 0
    try:
        with (
            pysam.VariantFile(str(input_file)) as source,
            pysam.VariantFile(str(cache)) as annotations,
        ):
            for contig in source.header.contigs:
                if sampled >= 128 or time.monotonic() - start >= max_seconds:
                    break
                if contig not in annotations.header.contigs:
                    continue
                for index, record in enumerate(source.fetch(contig)):
                    if (
                        index >= 8
                        or sampled >= 128
                        or time.monotonic() - start >= max_seconds
                    ):
                        break
                    # Native reserves 1 KiB plus encoded/decoded record storage.
                    size = 1024 + 4 * len(str(record).encode())
                    for candidate_number, candidate in enumerate(
                        annotations.fetch(contig, record.start, record.start + 1)
                    ):
                        if (
                            candidate_number >= 64
                            or time.monotonic() - start >= max_seconds
                        ):
                            break
                        if (
                            candidate.start == record.start
                            and candidate.alleles == record.alleles
                        ):
                            size += 3 * len(str(candidate).encode())
                            matches += 1
                            break
                    total += size
                    sampled += 1
    except (OSError, ValueError):
        return None
    if sampled < 16 or matches < 8:
        return None
    return dict(
        sampled_records=sampled,
        exact_cache_matches=matches,
        estimated_buffer_bytes=int(records * total / sampled),
        wall_seconds=time.monotonic() - start,
        estimate_note="Heuristic decoded-payload allowance, not a buffer/RSS guarantee; no whole-cache scan",
    )
