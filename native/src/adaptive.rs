//! One-pass normalized spooling and index-only cache-work estimates.
use anyhow::{ensure, Context, Result};
use rust_htslib::{
    bcf::{self, Read},
    htslib,
};
use serde::Serialize;
use std::{collections::HashMap, ffi::CString, fs::File, io::BufWriter, path::Path, ptr};

const WINDOW: i64 = 16 * 1024;

#[derive(Serialize)]
struct Window {
    contig: String,
    start: i64,
    end: i64,
    records: u64,
    estimated_cache_bytes: u64,
}

#[derive(Serialize)]
struct Contig {
    name: String,
    cache_records: u64,
    estimated_cache_bytes: u64,
    requested_cache_bytes: u64,
}

/// The index owns its data; every temporary iterator is destroyed immediately.
struct CacheIndex {
    index: *mut htslib::hts_idx_t,
    header: bcf::header::HeaderView,
    file_size: u64,
}

impl CacheIndex {
    fn records(&self, name: &str) -> u64 {
        let Ok(rid) = self.header.name2rid(name.as_bytes()) else {
            return 0;
        };
        let (mut mapped, mut unmapped) = (0, 0);
        if unsafe { htslib::hts_idx_get_stat(self.index, rid as i32, &mut mapped, &mut unmapped) }
            == 0
        {
            mapped
        } else {
            0
        }
    }
    fn open(path: &Path) -> Result<Self> {
        let reader = bcf::Reader::from_path(path)?;
        let name = CString::new(path.to_str().context("non-UTF8 cache path")?)?;
        let index = unsafe { htslib::bcf_index_load3(name.as_ptr(), ptr::null(), 0) };
        ensure!(!index.is_null(), "cannot load cache CSI index");
        Ok(Self {
            index,
            header: reader.header().clone(),
            file_size: path.metadata()?.len(),
        })
    }

    fn ranges(&self, name: &str, start: i64, end: i64) -> Vec<(u64, u64)> {
        let Ok(rid) = self.header.name2rid(name.as_bytes()) else {
            return vec![];
        };
        // hts_itr_query computes candidate chunks without reading cache records.
        let iterator = unsafe {
            htslib::hts_itr_query(
                self.index,
                rid as i32,
                start,
                end,
                Some(htslib::bcf_readrec),
            )
        };
        if iterator.is_null() {
            return vec![];
        }
        let mut result = Vec::new();
        unsafe {
            for i in 0..(*iterator).n_off {
                let chunk = *(*iterator).off.add(i as usize);
                // Virtual offsets encode compressed block offsets in high bits.
                // Include the final BGZF block conservatively, capped at EOF.
                result.push((chunk.u >> 16, ((chunk.v >> 16) + 65536).min(self.file_size)));
            }
            htslib::hts_itr_destroy(iterator);
        }
        result
    }
}

impl Drop for CacheIndex {
    fn drop(&mut self) {
        unsafe { htslib::hts_idx_destroy(self.index) };
    }
}

fn union_bytes(ranges: &mut [(u64, u64)]) -> u64 {
    ranges.sort_unstable();
    let mut sum = 0;
    let (mut begin, mut end) = (0, 0);
    for &(left, right) in ranges.iter() {
        if left > end {
            sum += end - begin;
            begin = left;
            end = right;
        } else {
            end = end.max(right);
        }
    }
    sum + end - begin
}

/// Minimal safe ownership wrapper for public HTSlib indexed wb0 output.
/// rust-htslib's Writer does not expose BGZF level-zero or index-on-write.
pub(super) struct Spool {
    file: *mut htslib::htsFile,
    index_name: CString,
}

impl Spool {
    pub(super) fn new(path: &Path, header: &bcf::header::HeaderView) -> Result<Self> {
        Self::open(path, header, c"wb0", 0)
    }

    pub(super) fn new_compressed(
        path: &Path,
        header: &bcf::header::HeaderView,
        threads: usize,
    ) -> Result<Self> {
        Self::open(path, header, c"wb", threads)
    }

    fn open(
        path: &Path,
        header: &bcf::header::HeaderView,
        mode: &std::ffi::CStr,
        threads: usize,
    ) -> Result<Self> {
        let name = CString::new(path.to_str().context("non-UTF8 output path")?)?;
        let file = unsafe { htslib::hts_open(name.as_ptr(), mode.as_ptr()) };
        ensure!(!file.is_null(), "cannot create normalized spool");
        let spool = Self {
            file,
            index_name: CString::new(format!("{}.csi", path.display()))?,
        };
        unsafe {
            if threads > 0 {
                ensure!(
                    htslib::hts_set_threads(file, threads as i32) == 0,
                    "cannot set writer threads"
                );
            }
            ensure!(
                htslib::bcf_hdr_write(file, header.as_ptr()) == 0,
                "cannot write spool header"
            );
            // The header may already have dispatched BGZF blocks to the
            // encoder pool. Stabilize its virtual offset before idx_init;
            // otherwise large headers can produce a bad first CSI chunk.
            ensure!(
                htslib::bgzf_flush(htslib::hts_get_bgzfp(file)) == 0,
                "cannot flush BCF header before indexing"
            );
            ensure!(
                htslib::bcf_idx_init(file, header.as_ptr(), 14, spool.index_name.as_ptr()) == 0,
                "cannot initialize spool index"
            );
        }
        Ok(spool)
    }

    pub(super) fn write(
        &self,
        header: &bcf::header::HeaderView,
        record: &bcf::Record,
    ) -> Result<()> {
        unsafe {
            ensure!(
                htslib::bcf_write(self.file, header.as_ptr(), record.inner) == 0,
                "cannot write normalized record"
            );
        }
        Ok(())
    }

    pub(super) fn finish(mut self) -> Result<()> {
        let index_span = super::profiling::Span::new("bcf.index_save_flush");
        ensure!(
            unsafe { htslib::bcf_idx_save(self.file) } == 0,
            "cannot save normalized CSI"
        );
        drop(index_span);
        let _close_span = super::profiling::Span::new("bcf.close_drain_codecs");
        let file = std::mem::replace(&mut self.file, ptr::null_mut());
        ensure!(
            unsafe { htslib::hts_close(file) } == 0,
            "cannot finish normalized spool"
        );
        Ok(())
    }
}

impl Drop for Spool {
    fn drop(&mut self) {
        if !self.file.is_null() {
            unsafe { htslib::hts_close(self.file) };
        }
    }
}

pub fn prepare(output: &Path, metadata: &Path, cache: &Path) -> Result<()> {
    let _span = super::profiling::Span::new("adaptive.prepare_decode_index_spool");
    let mut reader = bcf::Reader::from_stdin()?;
    let spool = Spool::new(output, reader.header())?;
    let index = CacheIndex::open(cache)?;
    let mut record = reader.empty_record();
    let mut windows: Vec<Window> = vec![];
    let mut previous = None;
    while let Some(result) = reader.read(&mut record) {
        result?;
        let coordinate = (record.rid().context("record without contig")?, record.pos());
        ensure!(record.pos() >= 0, "record without position");
        if let Some(last) = previous {
            ensure!(coordinate >= last, "normalized input is not sorted");
        }
        previous = Some(coordinate);
        let name = std::str::from_utf8(reader.header().rid2name(coordinate.0)?)?;
        let start = record.pos() / WINDOW * WINDOW;
        if let Some(last) = windows
            .last_mut()
            .filter(|last| last.contig == name && last.start == start)
        {
            last.records += 1;
        } else {
            windows.push(Window {
                contig: name.into(),
                start,
                end: start + WINDOW,
                records: 1,
                estimated_cache_bytes: 0,
            });
        }
        spool.write(reader.header(), &record)?;
    }
    spool.finish()?;
    let mut all_ranges: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    let mut names = Vec::new();
    for window in &mut windows {
        let mut ranges = index.ranges(&window.contig, window.start, window.end);
        window.estimated_cache_bytes = union_bytes(&mut ranges);
        if !all_ranges.contains_key(&window.contig) {
            names.push(window.contig.clone());
        }
        all_ranges
            .entry(window.contig.clone())
            .or_default()
            .extend(ranges);
    }
    let contigs: Vec<_> = names
        .into_iter()
        .map(|name| {
            let requested_cache_bytes =
                union_bytes(all_ranges.get_mut(&name).expect("collected contig"));
            let estimated_cache_bytes = union_bytes(&mut index.ranges(&name, 0, i64::MAX));
            Contig {
                cache_records: index.records(&name),
                name,
                estimated_cache_bytes,
                requested_cache_bytes,
            }
        })
        .collect();
    serde_json::to_writer(
        BufWriter::new(File::create(metadata)?),
        &serde_json::json!({
            "window_size": WINDOW, "windows": windows, "contigs": contigs,
            "estimate_note": "CSI candidate compressed byte ranges including conservative final-block allowance; not observed I/O"
        }),
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::union_bytes;
    #[test]
    fn overlapping_chunks_are_counted_once() {
        assert_eq!(union_bytes(&mut [(20, 30), (0, 10), (5, 25), (40, 45)]), 35);
        assert_eq!(union_bytes(&mut []), 0);
    }
}
