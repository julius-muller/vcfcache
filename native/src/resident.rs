//! Opt-in bounded-memory BCF overlay. No cache-format or recipe changes.
use anyhow::{ensure, Context, Result};
use rust_htslib::{
    bcf::{
        self,
        header::{HeaderRecord, TagType},
        Read,
    },
    htslib,
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{HashMap, HashSet},
    ffi::CString,
    fs::File,
    io::{BufRead, Write},
    path::{Path, PathBuf},
    ptr,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    },
    time::Instant,
};

#[derive(Deserialize)]
struct Job {
    input: PathBuf,
    cache: PathBuf,
    misses: PathBuf,
    annotated: PathBuf,
    output: PathBuf,
    stats: PathBuf,
    tag: String,
    workers: usize,
    cpus: usize,
    memory_bytes: usize,
    preserve: bool,
    lookup: String,
    #[serde(default = "default_window_bp")]
    window_bp: i64,
}

fn default_window_bp() -> i64 {
    16384
}

#[derive(Default, Serialize)]
struct Stats {
    input_variants: u64,
    missing_variants: u64,
    output_variants: u64,
    annotated_output: u64,
    missing_annotated: u64,
    candidate_records: u64,
    lookup_queries: u64,
    buffer_bytes: usize,
    phase_timings: Vec<serde_json::Value>,
}

// Record buffers dominate memory. This deliberately conservative estimate also
// reserves space for unpacked structures and allocator overhead. It is a buffer
// budget, not a claim to cap whole-process RSS (headers/HTSlib live separately).
fn bytes(record: &bcf::Record) -> usize {
    let r = record.inner();
    let mut size = 1024 + 2 * (r.shared.m + r.indiv.m);
    // Updated INFO strings may live outside shared.s until serialization.
    // Include even borrowed vptr blocks, deliberately overestimating buffers.
    if !r.d.info.is_null() {
        for i in 0..r.n_info() as usize {
            size += unsafe { (*r.d.info.add(i)).vptr_len as usize };
        }
    }
    size
}

fn timed(stats: &mut Stats, phase: &str, start: Instant) {
    stats
        .phase_timings
        .push(serde_json::json!({"phase":phase,"wall_seconds":start.elapsed().as_secs_f64()}));
}

fn next_record(reader: &mut bcf::Reader) -> Result<Option<bcf::Record>> {
    let mut record = reader.empty_record();
    match reader.read(&mut record) {
        Some(result) => {
            result?;
            Ok(Some(record))
        }
        None => Ok(None),
    }
}

fn coordinate(
    record: &bcf::Record,
    header: &bcf::header::HeaderView,
) -> Result<Option<(u32, i64)>> {
    let name = record
        .header()
        .rid2name(record.rid().context("record without contig")?)?;
    Ok(header.name2rid(name).ok().map(|rid| (rid, record.pos())))
}

fn event(value: serde_json::Value) -> Result<()> {
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    serde_json::to_writer(&mut out, &value)?;
    out.write_all(b"\n")?;
    out.flush()?;
    Ok(())
}

fn unsupported(message: &str) -> anyhow::Error {
    anyhow::anyhow!("UNSUPPORTED: {message}")
}

/// Copy INFO declarations verbatim, retaining descriptions and IDX translation.
fn merge_header(
    base: &bcf::header::HeaderView,
    incoming: &bcf::header::HeaderView,
) -> Result<Arc<bcf::header::HeaderView>> {
    let header = bcf::Header::from_template(base);
    unsafe {
        let src = incoming.as_ptr();
        for i in 0..(*src).nhrec {
            let rec = *(*src).hrec.add(i as usize);
            if (*rec).type_ != htslib::BCF_HL_INFO as i32 {
                continue;
            }
            // bcf_hdr_add_hrec takes ownership. Copy rather than borrowing.
            let duplicate = htslib::bcf_hrec_dup(rec);
            ensure!(!duplicate.is_null(), "cannot copy INFO declaration");
            ensure!(
                htslib::bcf_hdr_add_hrec(header.inner_ptr(), duplicate) >= 0,
                "cannot merge INFO header"
            );
        }
        ensure!(
            htslib::bcf_hdr_sync(header.inner_ptr()) == 0,
            "cannot sync INFO header"
        );
        Ok(Arc::new(bcf::header::HeaderView::from_ptr(
            htslib::bcf_hdr_dup(header.inner_ptr()),
        )))
    }
}

fn tags(header: &bcf::header::HeaderView) -> Result<Vec<(Vec<u8>, TagType)>> {
    header
        .header_records()
        .into_iter()
        .filter_map(|row| match row {
            HeaderRecord::Info { values, .. } => Some(values["ID"].as_bytes().to_vec()),
            _ => None,
        })
        .map(|name| Ok((name.clone(), header.info_type(&name)?.0)))
        .collect()
}

fn compatible(base: &bcf::header::HeaderView, incoming: &bcf::header::HeaderView) -> Result<()> {
    for (name, _) in tags(incoming)? {
        if let Ok(left) = base.info_type(&name) {
            ensure!(
                left == incoming.info_type(&name)?,
                "Incompatible INFO header for {}",
                String::from_utf8_lossy(&name)
            );
        }
    }
    Ok(())
}

/// bcftools -c INFO copies declared present values, including missing values. Flags are
/// special: an absent source flag removes an existing target flag on a match.
fn overlay(
    target: &mut bcf::Record,
    source: &bcf::Record,
    fields: &[(Vec<u8>, TagType)],
) -> Result<()> {
    for (name, kind) in fields {
        match kind {
            TagType::Integer => {
                if let Some(v) = source.info(name).integer()? {
                    target.push_info_integer(name, &v)?;
                }
            }
            TagType::Float => {
                if let Some(v) = source.info(name).float()? {
                    target.push_info_float(name, &v)?;
                }
            }
            TagType::String => {
                if let Some(v) = source.info(name).string()? {
                    if v.is_empty() || (v.len() == 1 && v[0].is_empty()) {
                        // Keep the present zero-length (flag-like) BCF entry,
                        // not the text "1" that bcftools query renders for it.
                        target.push_info_flag(name)?;
                    } else {
                        target.push_info_string(name, &v)?;
                    }
                }
            }
            TagType::Flag => {
                if source.info(name).flag()? {
                    target.push_info_flag(name)?;
                } else {
                    target.clear_info_flag(name)?;
                }
            }
        }
    }
    Ok(())
}

/// Direct public HTSlib iterator: only unpack records at requested loci.
/// rust-htslib's synced IndexedReader unpacks and copies every candidate.
struct CacheReader {
    file: *mut htslib::htsFile,
    index: *mut htslib::hts_idx_t,
    header: Arc<bcf::header::HeaderView>,
}

/// HTSlib supports sharing one codec pool across independent file handles.
/// The owning scope joins all readers before destroying this pool.
struct CodecPool(*mut htslib::hts_tpool);
unsafe impl Send for CodecPool {}
unsafe impl Sync for CodecPool {}
impl Drop for CodecPool {
    fn drop(&mut self) {
        unsafe { htslib::hts_tpool_destroy(self.0) };
    }
}

impl CacheReader {
    fn open(path: &Path) -> Result<Self> {
        let name = CString::new(path.to_str().context("non-UTF8 cache path")?)?;
        let file = unsafe { htslib::hts_open(name.as_ptr(), c"rb".as_ptr()) };
        ensure!(!file.is_null(), "cannot open cache");
        let fallback_header = bcf::Reader::from_path(path)?.header().clone();
        let mut reader = Self {
            file,
            index: ptr::null_mut(),
            header: Arc::new(fallback_header),
        };
        let header = unsafe { htslib::bcf_hdr_read(file) };
        ensure!(!header.is_null(), "cannot read cache header");
        reader.header = Arc::new(unsafe { bcf::header::HeaderView::from_ptr(header) });
        reader.index = unsafe { htslib::bcf_index_load3(name.as_ptr(), ptr::null(), 0) };
        ensure!(!reader.index.is_null(), "cannot read cache CSI");
        Ok(reader)
    }

    fn query<F>(&mut self, rid: u32, start: i64, end: i64, mut visit: F) -> Result<u64>
    where
        F: FnMut(&mut bcf::Record) -> Result<()>,
    {
        let itr = unsafe {
            htslib::hts_itr_query(
                self.index,
                rid as i32,
                start,
                end,
                Some(htslib::bcf_readrec),
            )
        };
        if itr.is_null() {
            return Ok(0);
        }
        let mut record = self.header.empty_record();
        let mut count = 0;
        let result = (|| {
            loop {
                let rc = unsafe {
                    htslib::hts_itr_next(
                        htslib::hts_get_bgzfp(self.file),
                        itr,
                        record.inner.cast(),
                        ptr::null_mut(),
                    )
                };
                if rc == -1 {
                    break;
                }
                ensure!(rc >= 0, "cache iterator read failed");
                if record.pos() >= start && record.pos() < end {
                    count += 1;
                    visit(&mut record)?;
                }
            }
            Ok(count)
        })();
        unsafe { htslib::hts_itr_destroy(itr) };
        result
    }
}

impl Drop for CacheReader {
    fn drop(&mut self) {
        unsafe {
            if !self.index.is_null() {
                htslib::hts_idx_destroy(self.index);
            }
            if !self.file.is_null() {
                htslib::hts_close(self.file);
            }
        }
    }
}

struct Entry {
    record: bcf::Record,
    missing: bool,
}

fn possible_padded_match(left: &[&[u8]], right: &[&[u8]]) -> bool {
    if left.len() != 2 || right.len() != 2 {
        return true;
    }
    let (long, short) = if left[0].len() >= right[0].len() {
        (left, right)
    } else {
        (right, left)
    };
    long[0].starts_with(short[0])
        && long[1].starts_with(short[1])
        && long[1].get(short[1].len()..) == long[0].get(short[0].len()..)
}

fn lookup(
    entries: &mut [Entry],
    cache: &Path,
    memory: &AtomicUsize,
    limit: usize,
    mode: &str,
    window_bp: i64,
    pool: Option<&CodecPool>,
) -> Result<(u64, u64)> {
    let mut worker_span = super::profiling::Span::new("native.lookup.worker");
    worker_span.details(serde_json::json!({"input_records":entries.len(), "mode":mode}));
    let setup_span = super::profiling::Span::new("native.lookup.setup");
    let mut reader = CacheReader::open(cache)?;
    if let Some(pool) = pool {
        let mut setting = htslib::htsThreadPool {
            pool: pool.0,
            qsize: 2,
        };
        ensure!(
            unsafe { htslib::hts_set_thread_pool(reader.file, &mut setting) } == 0,
            "cannot attach cache codec pool"
        );
    }
    unsafe {
        htslib::hts_set_cache_size(reader.file, 2 * 1024 * 1024);
    }
    let fields = tags(&reader.header)?;
    drop(setup_span);
    lookup_records(
        entries,
        &mut reader,
        memory,
        limit,
        mode,
        window_bp,
        &fields,
    )
}

fn lookup_records(
    entries: &mut [Entry],
    reader: &mut CacheReader,
    memory: &AtomicUsize,
    limit: usize,
    mode: &str,
    window_bp: i64,
    fields: &[(Vec<u8>, TagType)],
) -> Result<(u64, u64)> {
    let mut loci: HashMap<(u32, i64), Vec<usize>> = HashMap::new();
    let mut windows = Vec::<(u32, i64, i64)>::new();
    for (i, e) in entries.iter().enumerate() {
        let name = e
            .record
            .header()
            .rid2name(e.record.rid().context("missing contig")?)?;
        let Ok(rid) = reader.header.name2rid(name) else {
            continue;
        };
        let pos = e.record.pos();
        loci.entry((rid, pos)).or_default().push(i);
        let start = if mode == "contigs" {
            0
        } else {
            pos / window_bp * window_bp
        };
        let end = if mode == "contigs" {
            i64::MAX
        } else {
            start + window_bp
        };
        if let Some(last) = windows
            .last_mut()
            .filter(|last| last.0 == rid && last.2 >= start)
        {
            last.2 = last.2.max(end);
        } else {
            windows.push((rid, start, end));
        }
    }
    let mut seen = HashSet::new();
    let mut query_span = super::profiling::Span::new("native.lookup.query_decode_overlay");
    let mut candidates = 0;
    for &(rid, start, end) in &windows {
        candidates += reader.query(rid, start, end, |source| {
            let Some(indices) = loci.get(&(rid, source.pos())) else {
                return Ok(());
            };
            source.unpack();
            if source.alleles().len() != 2 {
                return Err(unsupported("multiallelic cache candidate"));
            }
            for &i in indices {
                let target = &mut entries[i].record;
                if source.alleles() != target.alleles() {
                    if source
                        .alleles()
                        .iter()
                        .zip(target.alleles())
                        .all(|(a, b)| a.eq_ignore_ascii_case(b))
                    {
                        return Err(unsupported("case-insensitive allele pairing"));
                    }
                    // Different REF representations can be compatible under
                    // bcftools pairing; never turn a padded match into a miss.
                    if source.alleles()[0] != target.alleles()[0]
                        && possible_padded_match(&source.alleles(), &target.alleles())
                    {
                        return Err(unsupported("non-exact REF pairing"));
                    }
                    continue;
                }
                if !seen.insert(i) {
                    return Err(unsupported("duplicate matching cache records"));
                }
                let before = bytes(target);
                overlay(target, source, fields)?;
                let after = bytes(target);
                if after > before {
                    let total =
                        memory.fetch_add(after - before, Ordering::Relaxed) + after - before;
                    if total > limit {
                        return Err(unsupported("native memory buffer budget exceeded"));
                    }
                }
            }
            Ok(())
        })?;
    }
    query_span.details(serde_json::json!({"candidates":candidates, "queries":windows.len(), "buffer_bytes":memory.load(Ordering::Relaxed)}));
    Ok((candidates, windows.len() as u64))
}

fn execute(job: &Job) -> Result<()> {
    ensure!(
        [1024, 4096, 16384].contains(&job.window_bp),
        "unsupported native query window"
    );
    let _attempt = super::profiling::Span::new("native.execute");
    let prepare_span = super::profiling::Span::new("native.prepare_decode_split_stream");
    let mut stats = Stats::default();
    let start = Instant::now();
    let mut input = bcf::Reader::from_path(&job.input)?;
    let cache_header = bcf::Reader::from_path(&job.cache)?.header().clone();
    compatible(input.header(), &cache_header)?;
    let mut header = merge_header(input.header(), &cache_header)?;
    if !matches!(
        header.info_type(job.tag.as_bytes()),
        Ok((TagType::String, _))
    ) {
        return Err(unsupported(
            "required INFO tag is not a declared String field",
        ));
    }
    let memory = AtomicUsize::new(0);
    let mut entries = Vec::new();
    let mut previous = None;
    let mut keys = HashSet::new();
    for result in input.records() {
        let mut record = result?;
        let coordinate = (record.rid().context("missing contig")?, record.pos());
        ensure!(
            record.pos() >= 0 && previous.is_none_or(|p| coordinate >= p),
            "input is not coordinate sorted"
        );
        if previous != Some(coordinate) {
            keys.clear();
        }
        previous = Some(coordinate);
        if record.alleles().len() != 2 {
            return Err(unsupported("unsplit multiallelic or non-variant input"));
        }
        if record.alleles()[1] == b"*" {
            continue;
        }
        if !keys.insert(
            record
                .alleles()
                .iter()
                .map(|v| v.to_vec())
                .collect::<Vec<_>>(),
        ) {
            return Err(unsupported("duplicate input key"));
        }
        record.translate(&mut header)?;
        let total = memory.fetch_add(bytes(&record), Ordering::Relaxed) + bytes(&record);
        if total > job.memory_bytes {
            return Err(unsupported("native memory buffer budget exceeded"));
        }
        entries.push(Entry {
            record,
            missing: false,
        });
    }
    stats.input_variants = entries.len() as u64;
    drop(prepare_span);
    timed(&mut stats, "prepare", start);
    let start = Instant::now();
    let lookup_span = super::profiling::Span::new("native.lookup");
    let workers = job.workers.max(1).min(entries.len().max(1));
    let size = entries.len().max(1).div_ceil(workers);
    let codec_threads = job.cpus.saturating_sub(workers);
    let pool = if codec_threads > 0 {
        let ptr = unsafe { htslib::hts_tpool_init(codec_threads as i32) };
        ensure!(!ptr.is_null(), "cannot create cache codec pool");
        Some(CodecPool(ptr))
    } else {
        None
    };
    // Keep equal-POS groups together so a cache locus is handled only once.
    let mut chunks = Vec::new();
    let mut rest = entries.as_mut_slice();
    while !rest.is_empty() {
        let mut end = size.min(rest.len());
        while end < rest.len()
            && rest[end - 1].record.rid() == rest[end].record.rid()
            && rest[end - 1].record.pos() == rest[end].record.pos()
        {
            end += 1;
        }
        let (chunk, remaining) = rest.split_at_mut(end);
        chunks.push(chunk);
        rest = remaining;
    }
    std::thread::scope(|scope| -> Result<()> {
        let handles: Vec<_> = chunks
            .into_iter()
            .map(|chunk| {
                let memory = &memory;
                let pool = pool.as_ref();
                scope.spawn(move || {
                    lookup(
                        chunk,
                        &job.cache,
                        memory,
                        job.memory_bytes,
                        &job.lookup,
                        job.window_bp,
                        pool,
                    )
                })
            })
            .collect();
        for handle in handles {
            let (c, q) = handle
                .join()
                .map_err(|_| anyhow::anyhow!("lookup worker panicked"))??;
            stats.candidate_records += c;
            stats.lookup_queries += q;
        }
        Ok(())
    })?;
    drop(pool);
    drop(lookup_span);
    timed(&mut stats, "lookup", start);
    stats.buffer_bytes = memory.load(Ordering::Relaxed);
    let classify_span = super::profiling::Span::new("native.classify_misses");
    let mut buffer = bcf::record::Buffer::new();
    for entry in &mut entries {
        entry.missing = !super::annotated(&entry.record, job.tag.as_bytes(), &mut buffer)?;
        stats.missing_variants += u64::from(entry.missing);
    }
    drop(classify_span);
    if stats.missing_variants > 0 {
        let _span = super::profiling::Span::new("native.spool_misses");
        let writer = super::adaptive::Spool::new(&job.misses, &header)?;
        for entry in entries.iter().filter(|e| e.missing) {
            writer.write(&header, &entry.record)?;
        }
        writer.finish()?;
    }
    event(
        serde_json::json!({"event":"misses-ready","input_variants":stats.input_variants,"missing_variants":stats.missing_variants,"buffer_bytes":stats.buffer_bytes}),
    )?;
    let wait_span = super::profiling::Span::new("native.await_annotation");
    let mut line = String::new();
    ensure!(
        std::io::stdin().lock().read_line(&mut line)? > 0,
        "controller disconnected"
    );
    let response: serde_json::Value = serde_json::from_str(&line)?;
    ensure!(response["action"] == "resume", "controller aborted");
    drop(wait_span);
    let preflight_span = super::profiling::Span::new("native.merge_preflight");
    let start = Instant::now();
    let mut recipe_reader = None;
    let mut recipe_fields = Vec::new();
    if stats.missing_variants > 0 {
        let mut path = response
            .get("annotated")
            .and_then(|v| v.as_str())
            .map(PathBuf::from)
            .unwrap_or_else(|| job.annotated.clone());
        let mut annotations = bcf::Reader::from_path(&path)?;
        compatible(&header, annotations.header())?;
        let merged = merge_header(&header, annotations.header())?;
        header = merged;
        for entry in &mut entries {
            entry.record.translate(&mut header)?;
        }
        // Arbitrary output matching stays with bcftools on the miss-only file;
        // the normal exact-biallelic path below requires no extra process.
        let mut map = HashMap::new();
        for (i, e) in entries.iter().enumerate() {
            map.insert(
                (
                    e.record
                        .header()
                        .rid2name(e.record.rid().unwrap())?
                        .to_vec(),
                    e.record.pos(),
                    e.record
                        .alleles()
                        .iter()
                        .map(|v| v.to_vec())
                        .collect::<Vec<_>>(),
                ),
                i,
            );
        }
        // Validate before touching retained records. Unusual allele pairing or
        // duplicates are resolved by bcftools WITHOUT invoking the recipe twice.
        let mut seen_keys = HashSet::new();
        let positions: HashSet<_> = map
            .keys()
            .map(|(name, pos, _)| (name.clone(), *pos))
            .collect();
        let mut compatibility_merge = false;
        let mut previous_recipe = None;
        for row in annotations.records() {
            let record = row?;
            stats.missing_annotated += 1;
            if let Some(coord) = coordinate(&record, &header)? {
                if previous_recipe.is_some_and(|previous| coord < previous) {
                    compatibility_merge = true;
                }
                previous_recipe = Some(coord);
            }
            let key = (
                record
                    .header()
                    .rid2name(record.rid().context("recipe output without contig")?)?
                    .to_vec(),
                record.pos(),
                record
                    .alleles()
                    .iter()
                    .map(|v| v.to_vec())
                    .collect::<Vec<_>>(),
            );
            if (!map.contains_key(&key) && positions.contains(&(key.0.clone(), key.1)))
                || (map.contains_key(&key) && !seen_keys.insert(key))
            {
                compatibility_merge = true;
            }
        }
        if compatibility_merge {
            let original = job.misses.with_file_name("merge-original.bcf");
            let writer = super::adaptive::Spool::new(&original, &header)?;
            for entry in &entries {
                writer.write(&header, &entry.record)?;
            }
            writer.finish()?;
            event(
                serde_json::json!({"event":"merge-fallback","input":original,"reason":"non-exact or duplicate recipe output"}),
            )?;
            line.clear();
            ensure!(
                std::io::stdin().lock().read_line(&mut line)? > 0,
                "controller disconnected"
            );
            let response: serde_json::Value = serde_json::from_str(&line)?;
            ensure!(response["action"] == "resume", "controller aborted");
            path = PathBuf::from(
                response["annotated"]
                    .as_str()
                    .context("missing resolved recipe path")?,
            );
        }
        annotations = bcf::Reader::from_path(&path)?;
        recipe_fields = tags(annotations.header())?;
        recipe_reader = Some(annotations);
    }
    timed(&mut stats, "merge-preflight", start);
    drop(preflight_span);
    let start = Instant::now();
    let write_span = super::profiling::Span::new("native.construct_encode_output");
    let writer =
        super::adaptive::Spool::new_compressed(&job.output, &header, job.cpus.saturating_sub(1))?;
    let mut first_output = None;
    let mut pending = match &mut recipe_reader {
        Some(reader) => next_record(reader)?,
        None => None,
    };
    let mut active_locus = None;
    let mut locus_annotations = Vec::new();
    // Consume originals, freeing each record after its only final write. Recipe
    // payloads are held for ONE coordinate, not appended to the whole buffer.
    for mut entry in entries {
        let coord = (
            entry.record.rid().context("missing output contig")?,
            entry.record.pos(),
        );
        if active_locus != Some(coord) {
            locus_annotations.clear();
            active_locus = Some(coord);
            while let Some(record) = pending.as_ref() {
                let incoming = coordinate(record, &header)?;
                if incoming.is_some_and(|v| v > coord) {
                    break;
                }
                let record = pending.take().unwrap();
                if incoming == Some(coord) {
                    locus_annotations.push(record);
                }
                pending = next_record(recipe_reader.as_mut().unwrap())?;
            }
        }
        if let Some(annotation) = locus_annotations
            .iter()
            .find(|r| r.alleles() == entry.record.alleles())
        {
            overlay(&mut entry.record, annotation, &recipe_fields)?;
        }
        let annotated = super::annotated(&entry.record, job.tag.as_bytes(), &mut buffer)?;
        if annotated {
            stats.annotated_output += 1;
        }
        if annotated || job.preserve {
            first_output.get_or_insert((
                entry.record.rid().context("missing output contig")?,
                entry.record.pos(),
            ));
            writer.write(&header, &entry.record)?;
            stats.output_variants += 1;
        }
    }
    drop(write_span);
    let finish_span = super::profiling::Span::new("native.finalize_output");
    writer.finish()?;
    drop(finish_span);
    if let Some((rid, pos)) = first_output {
        let _span = super::profiling::Span::new("native.verify_output_index");
        let mut verify = CacheReader::open(&job.output)?;
        ensure!(
            verify.query(rid, pos, pos + 1, |_| Ok(()))? > 0,
            "new CSI cannot retrieve the first output locus"
        );
    }
    timed(&mut stats, "merge-and-write", start);
    let mut file = File::create(&job.stats)?;
    serde_json::to_writer(&mut file, &stats)?;
    file.flush()?;
    event(serde_json::json!({"event":"complete"}))
}

/// Route at most 2K records (a complete equal-POS group) in each 128-MiB buffer.
/// Reuses resident INFO contracts without retaining a whole annotated genome.
#[allow(clippy::too_many_arguments)]
pub fn route_stream(
    cache: &Path,
    hits: &Path,
    misses: &Path,
    order: &Path,
    stats: &Path,
    tag: &str,
    threads: usize,
    indexed_misses: bool,
) -> Result<()> {
    let _span = super::profiling::Span::new("adaptive.native_route");
    let mut counts = super::Counts::default();
    let mut input = bcf::Reader::from_stdin()?;
    let result = (|| -> Result<()> {
        let pool = if threads > 0 {
            let ptr = unsafe { htslib::hts_tpool_init(threads as i32) };
            ensure!(!ptr.is_null(), "cannot create cache codec pool");
            Some(CodecPool(ptr))
        } else {
            None
        };
        let mut reader = CacheReader::open(cache)?;
        compatible(input.header(), &reader.header)?;
        let mut header = merge_header(input.header(), &reader.header)?;
        if !matches!(header.info_type(tag.as_bytes()), Ok((TagType::String, _))) {
            return Err(unsupported(
                "required INFO tag is not a declared String field",
            ));
        }
        if let Some(pool) = &pool {
            let mut setting = htslib::htsThreadPool {
                pool: pool.0,
                qsize: 2,
            };
            ensure!(
                unsafe { htslib::hts_set_thread_pool(reader.file, &mut setting) } == 0,
                "cannot attach route codec pool"
            );
        }
        unsafe {
            htslib::hts_set_cache_size(reader.file, 2 * 1024 * 1024);
        }
        let fields = tags(&reader.header)?;
        let template = bcf::Header::from_template(&header);
        let mut hits = bcf::Writer::from_path(hits, &template, true, bcf::Format::Bcf)?;
        let mut plain = if indexed_misses {
            None
        } else {
            Some(bcf::Writer::from_path(
                misses,
                &template,
                true,
                bcf::Format::Bcf,
            )?)
        };
        let indexed = if indexed_misses {
            Some(super::adaptive::Spool::new(misses, &header)?)
        } else {
            None
        };
        let mut order = std::io::BufWriter::new(File::create(order)?);
        let mut pending = next_record(&mut input)?;
        let mut previous = None;
        let mut keys = HashSet::new();
        let mut info_buffer = bcf::record::Buffer::new();
        const LIMIT: usize = 128 * 1024 * 1024;
        while pending.is_some() {
            let memory = AtomicUsize::new(0);
            let mut entries = Vec::with_capacity(2048);
            while let Some(mut record) = pending.take() {
                let coord = (record.rid().context("missing input contig")?, record.pos());
                if entries.len() >= 2048 && previous != Some(coord) {
                    pending = Some(record);
                    break;
                }
                ensure!(
                    record.pos() >= 0 && previous.is_none_or(|last| coord >= last),
                    "input is not coordinate sorted"
                );
                if previous != Some(coord) {
                    keys.clear();
                }
                previous = Some(coord);
                if record.alleles().len() != 2 {
                    return Err(unsupported("unsplit multiallelic input"));
                }
                if !keys.insert(
                    record
                        .alleles()
                        .iter()
                        .map(|v| v.to_vec())
                        .collect::<Vec<_>>(),
                ) {
                    return Err(unsupported("duplicate input key"));
                }
                record.translate(&mut header)?;
                let retained = memory.fetch_add(bytes(&record), Ordering::Relaxed) + bytes(&record);
                if retained > LIMIT {
                    return Err(unsupported("bounded router buffer budget exceeded"));
                }
                entries.push(Entry {
                    record,
                    missing: false,
                });
                pending = next_record(&mut input)?;
            }
            let (candidates, queries) = lookup_records(
                &mut entries,
                &mut reader,
                &memory,
                LIMIT,
                "regions",
                16384,
                &fields,
            )?;
            counts.native_candidate_records += candidates;
            counts.native_lookup_queries += queries;
            for entry in entries {
                let hit = super::annotated(&entry.record, tag.as_bytes(), &mut info_buffer)?;
                if hit {
                    hits.write(&entry.record)?;
                } else {
                    if let Some(writer) = &indexed {
                        writer.write(&header, &entry.record)?;
                    } else {
                        plain
                            .as_mut()
                            .context("missing plain spool writer")?
                            .write(&entry.record)?;
                    }
                    counts.missing_variants += 1;
                }
                order.write_all(&[u8::from(!hit)])?;
                counts.input_variants += 1;
            }
        }
        order.flush()?;
        drop(hits);
        drop(plain);
        if let Some(writer) = indexed {
            writer.finish()?;
        }
        // Readers must detach their codec pool before its owning scope ends.
        drop(reader);
        drop(pool);
        Ok(())
    })();
    if let Err(error) = result {
        if error.to_string().starts_with("UNSUPPORTED:") {
            // Drain the upstream selector, so a conservative fallback does not
            // turn into SIGPIPE under the controller's pipefail contract.
            for record in input.records() {
                record?;
            }
            counts.fallback_reason = Some(error.to_string());
        } else {
            return Err(error);
        }
    }
    super::save(stats, &counts)
}

pub fn run(path: &Path) -> Result<()> {
    let job: Job = serde_json::from_reader(File::open(path)?)?;
    ensure!(
        job.cpus > 0 && job.workers > 0 && job.memory_bytes > 0,
        "invalid native resource budget"
    );
    match execute(&job) {
        Err(error) if error.to_string().starts_with("UNSUPPORTED:") => {
            event(serde_json::json!({"event":"fallback","reason":error.to_string()}))
        }
        result => result,
    }
}
