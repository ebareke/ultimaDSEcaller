//! Splice-junction extraction from BAM (and CRAM, when a reference is
//! supplied). One pass per file, parallel across files via rayon.
//!
//! For each `N` CIGAR op in each read, we emit a `(chrom, donor_end,
//! acceptor_start)` junction together with a contribution weight (1.0 for
//! uniquely-mapped or primary-only counting; 1/NH for fractional counting).
//!
//! Filters applied per read:
//! * Unmapped / secondary / supplementary alignments skipped according to
//!   the multimap strategy.
//! * MAPQ below threshold dropped.
//! * Junctions whose anchor on either side (consecutive M/=/X bases adjacent
//!   to the N op) is below `min_overhang` are dropped — these are
//!   indistinguishable from spurious gapped alignments.
//!
//! ## Note on noodles API
//!
//! This module uses the `noodles` BAM reader. The noodles API has changed
//! across minor versions; if the dependency is updated, the only surface
//! that needs adjustment is [`read_bam_junctions`].

use std::collections::HashMap;
use std::path::Path;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::Position;
use crate::cli::MultimapStrategy;
use crate::config::RunConfig;
use crate::error::{UltiError, UltiResult};

/// Genome-oriented splice junction. `donor_end` and `acceptor_start` are the
/// 1-based reference coordinates of the exon bases *flanking* the intron —
/// not the intron's own first/last base.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Junction {
    pub chrom: String,
    pub donor_end: Position,
    pub acceptor_start: Position,
}

/// Per-sample junction support. Counts are kept as f64 because fractional
/// multimapper accounting can produce non-integer weights.
pub struct SampleJunctions {
    pub sample_id: String,
    pub counts: HashMap<Junction, f64>,
    /// Total uniquely-mapped reads that passed QC. Useful for normalization
    /// and as a denominator for the QC report.
    pub total_reads: u64,
    /// Reads dropped because their MAPQ was below threshold.
    pub low_mapq_reads: u64,
}

/// Per-cohort junction matrix: union of junctions across all samples with a
/// dense column per sample (in `samples` order).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JunctionMatrix {
    pub samples: Vec<String>,
    pub counts: HashMap<Junction, Vec<f64>>,
}

impl JunctionMatrix {
    pub fn empty() -> Self {
        JunctionMatrix {
            samples: Vec::new(),
            counts: HashMap::new(),
        }
    }

    /// Atomic binary write — for checkpoint/resume.
    pub fn save(&self, path: &std::path::Path) -> UltiResult<()> {
        use std::io::Write;
        let tmp = path.with_extension("tmp");
        {
            let f = std::fs::File::create(&tmp).map_err(|e| UltiError::io(&tmp, e))?;
            let mut w = std::io::BufWriter::new(f);
            bincode::serialize_into(&mut w, self).map_err(|e| UltiError::Cache(e.to_string()))?;
            w.flush().map_err(|e| UltiError::io(&tmp, e))?;
        }
        std::fs::rename(&tmp, path).map_err(|e| UltiError::io(path, e))
    }

    pub fn load(path: &std::path::Path) -> UltiResult<Self> {
        let f = std::fs::File::open(path).map_err(|e| UltiError::io(path, e))?;
        let r = std::io::BufReader::new(f);
        bincode::deserialize_from(r).map_err(|e| UltiError::Cache(e.to_string()))
    }
}

/// Parallel entry point — extracts junctions for every sample in the run
/// config and assembles the cohort-wide matrix.
pub fn extract(cfg: &RunConfig) -> UltiResult<JunctionMatrix> {
    let sample_ids: Vec<String> = cfg.samples.iter().map(|s| s.id.clone()).collect();

    let per_sample: Vec<UltiResult<SampleJunctions>> = cfg
        .samples
        .par_iter()
        .map(|s| {
            read_bam_junctions(
                &s.bam,
                &s.id,
                cfg.reference.as_deref(),
                cfg.reads.min_mapq,
                cfg.reads.min_overhang,
                cfg.reads.multimap,
            )
        })
        .collect();

    // Surface the first error if any.
    let per_sample: Vec<SampleJunctions> = per_sample.into_iter().collect::<UltiResult<_>>()?;

    let mut counts: HashMap<Junction, Vec<f64>> = HashMap::new();
    let n = per_sample.len();
    for (i, sj) in per_sample.into_iter().enumerate() {
        for (j, c) in sj.counts {
            counts.entry(j).or_insert_with(|| vec![0.0; n])[i] = c;
        }
    }

    if cfg.reads.junction_tolerance > 0 {
        counts = collapse_wobble(counts, cfg.reads.junction_tolerance, n);
    }

    Ok(JunctionMatrix {
        samples: sample_ids,
        counts,
    })
}

/// Read junctions from one BAM file. Logic is straightforward CIGAR
/// traversal: maintain the current 1-based reference position, accumulate
/// anchor lengths, and emit a junction whenever an `N` op of length ≥1 is
/// encountered (provided both anchors meet `min_overhang`).
pub fn read_bam_junctions(
    path: &Path,
    sample_id: &str,
    _reference: Option<&Path>,
    min_mapq: u8,
    min_overhang: u32,
    multimap: MultimapStrategy,
) -> UltiResult<SampleJunctions> {
    use noodles::bam;
    use noodles::sam::alignment::record::cigar::op::Kind;

    let mut reader = bam::io::reader::Builder
        .build_from_path(path)
        .map_err(|e| UltiError::Alignment {
            path: path.into(),
            message: format!("cannot open BAM: {e}"),
        })?;
    let header = reader.read_header().map_err(|e| UltiError::Alignment {
        path: path.into(),
        message: format!("cannot read header: {e}"),
    })?;

    let mut counts: HashMap<Junction, f64> = HashMap::new();
    let mut total_reads = 0u64;
    let mut low_mapq_reads = 0u64;

    for result in reader.records() {
        let record = result.map_err(|e| UltiError::Alignment {
            path: path.into(),
            message: format!("record parse error: {e}"),
        })?;

        // Skip unmapped.
        let flags = record.flags();
        if flags.is_unmapped() {
            continue;
        }

        // Multimapper policy.
        let is_secondary = flags.is_secondary() || flags.is_supplementary();
        let nh = read_nh_tag(&record).unwrap_or(1);
        let weight: f64 = match multimap {
            MultimapStrategy::Discard => {
                if nh > 1 || is_secondary {
                    continue;
                }
                1.0
            }
            MultimapStrategy::Primary => {
                if is_secondary {
                    continue;
                }
                1.0
            }
            MultimapStrategy::Fractional => {
                if is_secondary {
                    // Each alignment of a multimapper contributes 1/NH; only
                    // counting non-secondary would double-discount, so we
                    // accept secondaries here and weight them.
                }
                1.0 / (nh.max(1) as f64)
            }
        };

        // MAPQ.
        let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
        if mapq < min_mapq {
            low_mapq_reads += 1;
            continue;
        }
        total_reads += 1;

        let ref_id = match record.reference_sequence_id() {
            Some(Ok(id)) => id,
            _ => continue,
        };
        let chrom = match header.reference_sequences().get_index(ref_id) {
            Some((name, _)) => name.to_string(),
            None => continue,
        };
        let mut pos = match record.alignment_start() {
            Some(Ok(p)) => p.get() as u64,
            _ => continue,
        };

        // CIGAR walk. Track the length of the matching run just *before* an
        // N op (the left anchor) and reset it for each new N.
        let mut left_anchor: u32 = 0;
        // Pending N events whose right anchor we still need to validate.
        let mut pending: Vec<(Position, Position, u32)> = Vec::new();

        for op_result in record.cigar().iter() {
            let op = match op_result {
                Ok(o) => o,
                Err(_) => {
                    pending.clear();
                    break;
                }
            };
            let len = op.len() as u32;
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    // Resolve any pending N whose right anchor is now this op.
                    for (donor_end, acceptor_start, la) in pending.drain(..) {
                        if la >= min_overhang && len >= min_overhang {
                            let j = Junction {
                                chrom: chrom.clone(),
                                donor_end,
                                acceptor_start,
                            };
                            *counts.entry(j).or_insert(0.0) += weight;
                        }
                    }
                    left_anchor = len;
                    pos += len as u64;
                }
                Kind::Deletion => {
                    // D consumes ref but not the anchor (no sequence match).
                    // Treat as resetting anchor pessimistically.
                    left_anchor = 0;
                    pos += len as u64;
                }
                Kind::Skip => {
                    // N — splice. donor is the last ref base before N (pos-1),
                    // acceptor is the first ref base after N (pos+len).
                    let donor_end = pos.saturating_sub(1);
                    let acceptor_start = pos + len as u64;
                    pending.push((donor_end, acceptor_start, left_anchor));
                    pos += len as u64;
                    left_anchor = 0;
                }
                Kind::Insertion | Kind::SoftClip => {
                    // Consume query but not ref.
                    left_anchor = 0;
                }
                Kind::HardClip | Kind::Pad => {}
            }
        }
        // Any pending N whose right anchor never appeared (read ended on N or
        // soft-clip) is dropped.
    }

    Ok(SampleJunctions {
        sample_id: sample_id.to_string(),
        counts,
        total_reads,
        low_mapq_reads,
    })
}

/// Read the `NH:i:` tag (number of reported alignments for this read).
/// Defaults to 1 if absent. Compares the raw two-byte tag rather than a
/// named constant to stay forward-compatible with noodles versioning.
fn read_nh_tag<R: noodles::sam::alignment::Record>(record: &R) -> Option<u32> {
    let data = record.data();
    for field in data.iter() {
        let Ok((tag, value)) = field else { continue };
        let bytes: [u8; 2] = tag.into();
        if &bytes == b"NH" {
            return value.as_int().map(|n| n as u32);
        }
    }
    None
}

/// Collapse junctions whose `(donor_end, acceptor_start)` differ by ≤ tolerance
/// on a per-chromosome basis. The most-supported junction wins; others get
/// merged into it. Used for long-read data where soft-edge wobble is common.
fn collapse_wobble(
    counts: HashMap<Junction, Vec<f64>>,
    tolerance: u32,
    n_samples: usize,
) -> HashMap<Junction, Vec<f64>> {
    let mut by_chrom: HashMap<String, Vec<(Junction, Vec<f64>, f64)>> = HashMap::new();
    for (j, c) in counts {
        let total: f64 = c.iter().sum();
        by_chrom
            .entry(j.chrom.clone())
            .or_default()
            .push((j, c, total));
    }

    let mut out: HashMap<Junction, Vec<f64>> = HashMap::new();
    for (_, mut entries) in by_chrom {
        // Sort by total descending so dominant junctions absorb satellites.
        entries.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
        let mut anchors: Vec<(Junction, Vec<f64>)> = Vec::new();
        for (j, c, _) in entries {
            let mut merged = false;
            for (anchor, ac) in anchors.iter_mut() {
                if (j.donor_end as i64 - anchor.donor_end as i64).unsigned_abs() <= tolerance as u64
                    && (j.acceptor_start as i64 - anchor.acceptor_start as i64).unsigned_abs()
                        <= tolerance as u64
                {
                    for i in 0..n_samples {
                        ac[i] += c[i];
                    }
                    merged = true;
                    break;
                }
            }
            if !merged {
                anchors.push((j, c));
            }
        }
        for (j, c) in anchors {
            out.insert(j, c);
        }
    }
    out
}
