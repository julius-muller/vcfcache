//! Experimental store only: no automatic build, discovery or production use.
//!
//! Locus-first candidate retrieval intentionally retains all alleles at a locus.
//! Full variant keys are stored, but bcftools remains the semantic authority:
//! its padded-REF/allele-pair matching is not equivalent to exact-key matching.
use anyhow::{ensure, Context, Result};
use memmap2::Mmap;
use rust_htslib::bcf::{self, Read};
use serde_json::json;
use std::{
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
    process::{Command, Stdio},
    time::Instant,
};

const MAGIC: &[u8] = b"VCFIDX01";
const ENTRY: usize = 32;

fn mapped(path: &Path) -> Result<Mmap> {
    // Immutable store files are opened read-only. Callers must not truncate or
    // modify them while mapped, as with every mmap-backed database.
    let file = File::open(path)?;
    if file.metadata()?.len() == 0 {
        return Ok(memmap2::MmapMut::map_anon(1)?.make_read_only()?);
    }
    Ok(unsafe { Mmap::map(&file)? })
}

fn u32_at(bytes: &[u8], offset: usize) -> Result<u32> {
    Ok(u32::from_le_bytes(
        bytes
            .get(offset..offset + 4)
            .context("truncated store")?
            .try_into()?,
    ))
}

fn u64_at(bytes: &[u8], offset: usize) -> Result<u64> {
    Ok(u64::from_le_bytes(
        bytes
            .get(offset..offset + 8)
            .context("truncated store")?
            .try_into()?,
    ))
}

fn hash(key: &[u8]) -> u64 {
    key.iter().fold(0xcbf29ce484222325_u64, |h, b| {
        (h ^ u64::from(*b)).wrapping_mul(0x100000001b3)
    })
}

fn locus(reader: &bcf::Reader, record: &bcf::Record) -> Result<Vec<u8>> {
    let mut key = reader
        .header()
        .rid2name(record.rid().context("missing contig")?)?
        .to_vec();
    key.push(0);
    key.extend_from_slice(&record.pos().to_le_bytes());
    Ok(key)
}

fn header_end(payload: &[u8]) -> Result<usize> {
    ensure!(
        payload.starts_with(b"BCF\x02\x02"),
        "store payload is not plain BCF2.2"
    );
    let end = 9 + u32_at(payload, 5)? as usize;
    ensure!(end <= payload.len(), "truncated BCF header");
    Ok(end)
}

pub fn build(cache: &Path, output: &Path) -> Result<()> {
    let start = Instant::now();
    fs::create_dir(output).context("store output must be a new directory")?;
    let payload_path = output.join("payload.bcf");
    // -Ou uses uncompressed BGZF *blocks*, not a raw mmap-able BCF stream.
    // Unwrap those blocks once during explicit conversion, never per sample.
    let mut producer = Command::new("bcftools")
        .arg("view")
        .arg("--no-version")
        .arg(cache)
        .arg("-Ou")
        .stdout(Stdio::piped())
        .spawn()?;
    let decoded = Command::new("gzip")
        .arg("-dc")
        .stdin(producer.stdout.take().context("missing BCF pipe")?)
        .stdout(File::create(&payload_path)?)
        .status();
    let produced = producer.wait()?;
    ensure!(
        decoded?.success() && produced.success(),
        "cannot create plain BCF payload"
    );
    let payload = mapped(&payload_path)?;
    let mut offset = header_end(&payload)?;
    let mut reader = bcf::Reader::from_path(&payload_path)?;
    let mut record = reader.empty_record();
    let mut key_writer = BufWriter::new(File::create(output.join("keys.bin"))?);
    let mut key_offset = 0u64;
    let mut entries = Vec::new();
    while let Some(result) = reader.read(&mut record) {
        result?;
        let prefix = locus(&reader, &record)?;
        let mut key = prefix.clone();
        for allele in record.alleles() {
            key.extend_from_slice(&(u32::try_from(allele.len())?).to_le_bytes());
            key.extend_from_slice(allele);
        }
        let length =
            8usize + u32_at(&payload, offset)? as usize + u32_at(&payload, offset + 4)? as usize;
        ensure!(offset + length <= payload.len(), "truncated payload record");
        entries.push((
            hash(&prefix),
            offset as u64,
            key_offset,
            u32::try_from(length)?,
            u32::try_from(key.len())?,
        ));
        key_writer.write_all(&key)?;
        key_offset += key.len() as u64;
        offset += length;
    }
    ensure!(
        offset == payload.len(),
        "payload and record parser disagree"
    );
    key_writer.flush()?;
    entries.sort_unstable_by_key(|entry| (entry.0, entry.1));
    let mut index = BufWriter::new(File::create(output.join("index.bin"))?);
    index.write_all(MAGIC)?;
    for &(h, offset, key_offset, length, key_length) in &entries {
        index.write_all(&h.to_le_bytes())?;
        index.write_all(&offset.to_le_bytes())?;
        index.write_all(&key_offset.to_le_bytes())?;
        index.write_all(&length.to_le_bytes())?;
        index.write_all(&key_length.to_le_bytes())?;
    }
    index.flush()?;
    let checksum = Command::new("sha256sum").arg(cache).output()?;
    ensure!(checksum.status.success(), "source checksum failed");
    let checksum = String::from_utf8(checksum.stdout)?;
    let manifest = json!({"experimental":true,"format":"VCFIDX01","records":entries.len(),
        "source":cache,"source_sha256":checksum.split_whitespace().next(),
        "build_seconds":start.elapsed().as_secs_f64(),
        "payload_bytes":payload.len(),"key_bytes":key_offset,"index_bytes":8+entries.len()*ENTRY,
        "note":"Read-only candidate store, not a replacement for bcftools INFO-overlay semantics. Build time includes decompression, indexing and source hashing."});
    fs::write(
        output.join("manifest.json"),
        serde_json::to_vec_pretty(&manifest)?,
    )?;
    Ok(())
}

pub fn query(store: &Path, input: &Path) -> Result<()> {
    ensure!(
        store.join("manifest.json").is_file(),
        "store build did not finish"
    );
    let index = mapped(&store.join("index.bin"))?;
    let payload = mapped(&store.join("payload.bcf"))?;
    let keys = mapped(&store.join("keys.bin"))?;
    ensure!(
        index.starts_with(MAGIC) && (index.len() - 8) % ENTRY == 0,
        "invalid store index"
    );
    let count = (index.len() - 8) / ENTRY;
    let mut reader = bcf::Reader::from_path(input)?;
    let mut record = reader.empty_record();
    let mut selected = Vec::new();
    let mut previous = None;
    while let Some(result) = reader.read(&mut record) {
        result?;
        let coordinate = (record.rid(), record.pos());
        if previous == Some(coordinate) {
            continue;
        }
        previous = Some(coordinate);
        let prefix = locus(&reader, &record)?;
        let h = hash(&prefix);
        let (mut left, mut right) = (0, count);
        while left < right {
            let middle = (left + right) / 2;
            if u64_at(&index, 8 + middle * ENTRY)? < h {
                left = middle + 1;
            } else {
                right = middle;
            }
        }
        while left < count && u64_at(&index, 8 + left * ENTRY)? == h {
            let base = 8 + left * ENTRY;
            let offset = u64_at(&index, base + 8)? as usize;
            let key_offset = u64_at(&index, base + 16)? as usize;
            let length = u32_at(&index, base + 24)? as usize;
            let key_length = u32_at(&index, base + 28)? as usize;
            let key = keys
                .get(key_offset..key_offset + key_length)
                .context("invalid key offset")?;
            if key.starts_with(&prefix) {
                ensure!(offset + length <= payload.len(), "invalid payload offset");
                selected.push((offset, length));
            }
            left += 1;
        }
    }
    selected.sort_unstable();
    selected.dedup();
    let mut output = BufWriter::new(std::io::stdout().lock());
    output.write_all(&payload[..header_end(&payload)?])?;
    for (offset, length) in selected {
        output.write_all(&payload[offset..offset + length])?;
    }
    output.flush()?;
    Ok(())
}
