// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Julius Müller
//! Stream routing and ordinal reconstruction; bcftools retains matching semantics.

use anyhow::{bail, ensure, Context, Result};
use clap::{Parser, Subcommand};
use rust_htslib::bcf::{self, header::TagType, Read};
use serde::Serialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read as IoRead, Write};
use std::path::{Path, PathBuf};

mod adaptive;
mod profiling;
mod resident;
#[cfg(feature = "experimental-store")]
mod store;

#[derive(Parser)]
#[command(version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Retain a small input in memory, perform native BCF lookup, await recipe.
    Resident {
        #[arg(long)]
        job: PathBuf,
    },
    /// EXPERIMENT: explicitly build an immutable mmap lookup sidecar.
    #[cfg(feature = "experimental-store")]
    StoreBuild {
        #[arg(long)]
        cache: PathBuf,
        #[arg(long)]
        output: PathBuf,
    },
    /// EXPERIMENT: emit candidate annotation BCF, retaining bcftools matching.
    #[cfg(feature = "experimental-store")]
    StoreQuery {
        #[arg(long)]
        store: PathBuf,
        #[arg(long)]
        input: PathBuf,
    },
    /// Spool normalized records and collect sparse lookup metadata in one pass.
    Prepare {
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        metadata: PathBuf,
        #[arg(long)]
        cache: PathBuf,
    },
    /// Split a cache-annotated uncompressed BCF stream, recording original order.
    Route {
        #[arg(long)]
        hits: PathBuf,
        #[arg(long)]
        misses: PathBuf,
        #[arg(long)]
        order: PathBuf,
        #[arg(long)]
        stats: PathBuf,
        #[arg(long)]
        tag: String,
        /// Experimental index-on-write miss spool; legacy plain BCF is the default.
        #[arg(long)]
        indexed_misses: bool,
    },
    /// Restore input order, filter annotation-less records and emit plain BCF.
    Assemble {
        #[arg(long)]
        hits: PathBuf,
        #[arg(long)]
        misses: PathBuf,
        #[arg(long)]
        order: PathBuf,
        #[arg(long)]
        stats: PathBuf,
        #[arg(long)]
        tag: String,
        #[arg(long)]
        preserve_unannotated: bool,
        #[arg(long)]
        output: Option<PathBuf>,
        #[arg(long, default_value_t = 0)]
        threads: usize,
        /// Trust the router's hit classification instead of parsing the same INFO again.
        #[arg(long)]
        trust_routed_hits: bool,
    },
}

#[derive(Default, Serialize)]
struct Counts {
    input_variants: u64,
    missing_variants: u64,
    output_variants: u64,
    annotated_output: u64,
    missing_annotated: u64,
    fallback_reason: Option<String>,
    trusted_hit_records: u64,
}

fn save(path: &Path, counts: &Counts) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    serde_json::to_writer(&mut writer, counts)?;
    writer.flush()?;
    Ok(())
}

fn annotated(record: &bcf::Record, tag: &[u8], buffer: &mut bcf::record::Buffer) -> Result<bool> {
    // bcftools INFO/TAG=="" matches an absent, empty, or dot-valued entry,
    // including an individual missing entry within a comma-separated vector.
    Ok(match record.info_shared_buffer(tag, buffer).string()? {
        None => false,
        Some(values) => reusable(&values),
    })
}

fn reusable(values: &[&[u8]]) -> bool {
    // A literal comma is a nonempty BCF string even when its elements are empty.
    (values.len() > 1 || values.iter().any(|v| !v.is_empty())) && values.iter().all(|v| *v != b".")
}

#[cfg(test)]
mod tests {
    use super::reusable;

    #[test]
    fn missing_string_entries_match_bcftools() {
        assert!(!reusable(&[]));
        assert!(!reusable(&[b""]));
        assert!(!reusable(&[b"."]));
        assert!(!reusable(&[b".", b"a"]));
        assert!(reusable(&[b"a", b""]));
        assert!(reusable(&[b"", b""]));
        assert!(reusable(&[b"C|gene|."]));
    }
}

fn route(
    hits: &Path,
    misses: &Path,
    order: &Path,
    stats: &Path,
    tag: &str,
    indexed_misses: bool,
) -> Result<()> {
    let _span = profiling::Span::new("adaptive.route_split_spool");
    let mut reader = bcf::Reader::from_stdin()?;
    let header = bcf::Header::from_template(reader.header());
    let mut hit_writer = bcf::Writer::from_path(hits, &header, true, bcf::Format::Bcf)?;
    let mut plain_misses = if indexed_misses {
        None
    } else {
        Some(bcf::Writer::from_path(
            misses,
            &header,
            true,
            bcf::Format::Bcf,
        )?)
    };
    let indexed = if indexed_misses {
        Some(adaptive::Spool::new(misses, reader.header())?)
    } else {
        None
    };
    let mut order_writer = BufWriter::new(File::create(order)?);
    let mut counts = Counts::default();
    let supported = matches!(
        reader.header().info_type(tag.as_bytes()),
        Ok((TagType::String, _))
    );
    if !supported {
        counts.fallback_reason = Some("required INFO tag is not a declared String field".into());
    }
    let mut previous = None;
    let mut locus_keys: HashMap<Vec<u8>, bool> = HashMap::new();
    let mut first_key = Vec::with_capacity(32);
    let mut key = Vec::with_capacity(32);
    let mut first_hit = false;
    let mut info_buffer = bcf::record::Buffer::new();
    let mut record = reader.empty_record();
    while let Some(result) = reader.read(&mut record) {
        result?;
        let coordinate = (record.rid().context("record without contig")?, record.pos());
        if let Some(last) = previous {
            ensure!(
                coordinate >= last,
                "input records are not coordinate sorted"
            );
            if coordinate != last {
                locus_keys.clear();
            }
        }
        let same_locus = previous == Some(coordinate);
        previous = Some(coordinate);
        let hit = supported && annotated(&record, tag.as_bytes(), &mut info_buffer)?;
        // Partitioning unsplit multiallelic or ambiguous duplicate keys can
        // change bcftools' second INFO overlay. Request explicit legacy replay.
        let alleles = record.alleles();
        if alleles.len() != 2 {
            counts.fallback_reason = Some("unsplit multiallelic or non-variant record".into());
        }
        key.clear();
        for allele in alleles {
            key.extend_from_slice(allele);
            key.push(0);
        }
        if same_locus {
            if locus_keys.is_empty() {
                locus_keys.insert(first_key.clone(), first_hit);
            }
            if let Some(was_hit) = locus_keys.insert(key.clone(), hit) {
                if was_hit != hit {
                    counts.fallback_reason =
                        Some("duplicate variant key occurs in both partitions".into());
                }
            }
        } else {
            first_key.clear();
            first_key.extend_from_slice(&key);
            first_hit = hit;
        }
        if hit {
            hit_writer.write(&record)?;
        } else {
            if let Some(writer) = &indexed {
                writer.write(reader.header(), &record)?;
            } else {
                plain_misses
                    .as_mut()
                    .context("missing plain spool writer")?
                    .write(&record)?;
            }
            counts.missing_variants += 1;
        }
        order_writer.write_all(&[u8::from(!hit)])?;
        counts.input_variants += 1;
    }
    order_writer.flush()?;
    // Explicitly close before Python compresses/indexes the miss spool.
    drop(hit_writer);
    drop(plain_misses);
    if let Some(writer) = indexed {
        writer.finish()?;
    }
    save(stats, &counts)
}

fn assemble(
    hits: &Path,
    misses: &Path,
    order: &Path,
    stats: &Path,
    tag: &str,
    preserve: bool,
    output: OutputOptions<'_>,
) -> Result<()> {
    let mut hits = bcf::Reader::from_path(hits)?;
    let _span = profiling::Span::new("adaptive.assemble_encode");
    let mut misses = if misses == Path::new("-") {
        bcf::Reader::from_stdin()?
    } else {
        bcf::Reader::from_path(misses)?
    };
    // The miss overlay's header is the original cache-annotated header plus
    // the recipe's new definitions, so it is the authoritative output header.
    let header = bcf::Header::from_template(misses.header());
    ensure!(
        hits.header().samples() == misses.header().samples(),
        "sample headers differ"
    );
    let mut writer = match output.path {
        Some(path) => bcf::Writer::from_path(path, &header, false, bcf::Format::Bcf)?,
        None => bcf::Writer::from_stdout(&header, true, bcf::Format::Bcf)?,
    };
    if output.threads > 0 && output.path.is_some() {
        writer.set_threads(output.threads)?;
    }
    let mut hits_record = hits.empty_record();
    let mut misses_record = misses.empty_record();
    let mut order = BufReader::new(File::open(order)?);
    let mut counts = Counts::default();
    let mut info_buffer = bcf::record::Buffer::new();
    let mut previous = None;
    let mut bytes = [0; 65536];
    loop {
        let length = order.read(&mut bytes)?;
        if length == 0 {
            break;
        }
        for partition in &bytes[..length] {
            let record = match partition {
                0 => {
                    hits.read(&mut hits_record)
                        .context("hit spool ended early")??;
                    &mut hits_record
                }
                1 => {
                    misses
                        .read(&mut misses_record)
                        .context("miss spool ended early")??;
                    counts.missing_variants += 1;
                    &mut misses_record
                }
                _ => bail!("invalid routing manifest byte"),
            };
            writer.translate(record);
            let coordinate = (record.rid().context("record without contig")?, record.pos());
            if let Some(last) = previous {
                ensure!(
                    coordinate >= last,
                    "reconstructed records are not coordinate sorted"
                );
            }
            previous = Some(coordinate);
            let known_hit = output.trust_routed_hits && *partition == 0;
            let has_annotation = known_hit || annotated(record, tag.as_bytes(), &mut info_buffer)?;
            counts.trusted_hit_records += u64::from(known_hit);
            counts.input_variants += 1;
            if has_annotation {
                counts.annotated_output += 1;
                if *partition == 1 {
                    counts.missing_annotated += 1;
                }
            }
            if has_annotation || preserve {
                writer.write(record)?;
                counts.output_variants += 1;
            }
        }
    }
    ensure!(hits.read(&mut hits_record).is_none(), "extra hit records");
    ensure!(
        misses.read(&mut misses_record).is_none(),
        "extra miss records"
    );
    drop(writer);
    save(stats, &counts)
}

struct OutputOptions<'a> {
    path: Option<&'a Path>,
    threads: usize,
    trust_routed_hits: bool,
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Resident { job } => resident::run(&job),
        #[cfg(feature = "experimental-store")]
        Command::StoreBuild { cache, output } => store::build(&cache, &output),
        #[cfg(feature = "experimental-store")]
        Command::StoreQuery { store, input } => store::query(&store, &input),
        Command::Prepare {
            output,
            metadata,
            cache,
        } => adaptive::prepare(&output, &metadata, &cache),
        Command::Route {
            hits,
            misses,
            order,
            stats,
            tag,
            indexed_misses,
        } => route(&hits, &misses, &order, &stats, &tag, indexed_misses),
        Command::Assemble {
            hits,
            misses,
            order,
            stats,
            tag,
            preserve_unannotated,
            output,
            threads,
            trust_routed_hits,
        } => assemble(
            &hits,
            &misses,
            &order,
            &stats,
            &tag,
            preserve_unannotated,
            OutputOptions {
                path: output.as_deref(),
                threads,
                trust_routed_hits,
            },
        ),
    }
}
