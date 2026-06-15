//! BAM pileup — mean read depth inside each requested interval.
//!
//! Two paths:
//! * **Indexed**: when a `.bai` (or `.csi`) sibling exists *and* the
//!   requested regions cover a small fraction of the genome (default
//!   threshold: 5%), per-region querying is dramatically faster.
//! * **Streaming**: one full-BAM pass, accumulating into all regions at
//!   once via a `Lapper`. Best when many regions or no index.
//!
//! `pileup_regions` dispatches automatically. Use [`pileup_regions_streaming`]
//! or [`pileup_regions_indexed`] explicitly if you want to force a mode.
//!
//! `mean_depth = total_match_overlap / interval_length`.
//!
//! Reads filtered out:
//! * unmapped, secondary, supplementary
//! * MAPQ below `min_mapq`
//! * any CIGAR-parse error (the read is silently skipped)

use std::collections::HashMap;
use std::path::Path;

use rust_lapper::{Interval, Lapper};

use crate::error::{UltiError, UltiResult};

/// `(chrom, start, end)` triple — interval bounds are 1-based inclusive, GTF-style.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Region {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

/// Per-interval results, in the same order the caller supplied `regions`.
pub struct PileupResult {
    pub mean_depth: Vec<f64>,
    pub total_match_bases: Vec<u64>,
    /// When `record_per_position == true` was set on the call, each region
    /// gets a `Vec<u32>` of per-base depths (length = region.end − region.start + 1).
    /// `None` otherwise — sparing the memory cost when callers only need
    /// the scalar summary.
    pub per_position: Option<Vec<Vec<u32>>>,
}

/// Options controlling pileup behavior.
#[derive(Debug, Clone, Copy, Default)]
pub struct PileupOpts {
    pub record_per_position: bool,
}

/// Dispatcher — picks indexed vs streaming based on availability and the
/// fraction of the genome the regions cover. Uses default `PileupOpts`
/// (no per-position output) for backward compatibility.
pub fn pileup_regions(
    path: &Path,
    regions: &[Region],
    min_mapq: u8,
) -> UltiResult<PileupResult> {
    pileup_regions_with_opts(path, regions, min_mapq, PileupOpts::default())
}

/// Like [`pileup_regions`] but exposes the [`PileupOpts`] knob to request
/// per-position depth tracks.
pub fn pileup_regions_with_opts(
    path: &Path,
    regions: &[Region],
    min_mapq: u8,
    opts: PileupOpts,
) -> UltiResult<PileupResult> {
    if regions.is_empty() {
        return Ok(PileupResult {
            mean_depth: Vec::new(),
            total_match_bases: Vec::new(),
            per_position: if opts.record_per_position { Some(Vec::new()) } else { None },
        });
    }
    if has_index(path) && small_query(regions) {
        match pileup_regions_indexed(path, regions, min_mapq, opts) {
            Ok(r) => return Ok(r),
            Err(e) => {
                tracing::warn!(
                    "indexed pileup failed ({e}); falling back to streaming"
                );
            }
        }
    }
    pileup_regions_streaming(path, regions, min_mapq, opts)
}

fn has_index(path: &Path) -> bool {
    let mut bai = path.to_path_buf();
    bai.set_extension(format!(
        "{}.bai",
        path.extension().and_then(|s| s.to_str()).unwrap_or("bam")
    ));
    if bai.exists() {
        return true;
    }
    let mut csi = path.to_path_buf();
    csi.set_extension(format!(
        "{}.csi",
        path.extension().and_then(|s| s.to_str()).unwrap_or("bam")
    ));
    csi.exists()
}

fn small_query(regions: &[Region]) -> bool {
    // Heuristic: indexed access wins when total region length is small.
    // Threshold: 100 Mb of total query span (covers e.g. all introns of a
    // few thousand human genes — exactly the IR / sashimi use case).
    let total: u64 = regions.iter().map(|r| r.end.saturating_sub(r.start) + 1).sum();
    total < 100_000_000
}

/// Query a BAM by region using its `.bai` / `.csi` index. Errors if no
/// index is present.
pub fn pileup_regions_indexed(
    path: &Path,
    regions: &[Region],
    min_mapq: u8,
    opts: PileupOpts,
) -> UltiResult<PileupResult> {
    use noodles::bam;
    use noodles::core::{Position as NPos, Region as NRegion};
    use noodles::sam::alignment::record::cigar::op::Kind;

    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .map_err(|e| UltiError::Alignment {
            path: path.into(),
            message: format!("indexed open failed: {e}"),
        })?;
    let header = reader.read_header().map_err(|e| UltiError::Alignment {
        path: path.into(),
        message: format!("indexed header: {e}"),
    })?;

    let mut total_bases = vec![0u64; regions.len()];
    let mut per_pos: Option<Vec<Vec<u32>>> = if opts.record_per_position {
        Some(
            regions
                .iter()
                .map(|r| vec![0u32; (r.end - r.start + 1) as usize])
                .collect(),
        )
    } else {
        None
    };
    for (region_idx, r) in regions.iter().enumerate() {
        let start = NPos::try_from(r.start.max(1) as usize).map_err(|_| {
            UltiError::Alignment {
                path: path.into(),
                message: "region start out of range".into(),
            }
        })?;
        let end = NPos::try_from(r.end.max(r.start) as usize).map_err(|_| {
            UltiError::Alignment {
                path: path.into(),
                message: "region end out of range".into(),
            }
        })?;
        let region = NRegion::new(r.chrom.as_bytes().to_vec(), start..=end);
        let query = match reader.query(&header, &region) {
            Ok(q) => q,
            Err(_) => continue, // chrom may be absent from header
        };
        for result in query {
            let record = match result {
                Ok(r) => r,
                Err(_) => continue,
            };
            let flags = record.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                continue;
            }
            let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
            if mapq < min_mapq {
                continue;
            }
            let mut pos = match record.alignment_start() {
                Some(Ok(p)) => p.get() as u64,
                _ => continue,
            };
            for op_result in record.cigar().iter() {
                let Ok(op) = op_result else { break };
                let len = op.len() as u64;
                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                        let run_start = pos;
                        let run_end_excl = pos + len;
                        let r_lo = run_start.max(r.start);
                        let r_hi = run_end_excl.min(r.end + 1);
                        if r_hi > r_lo {
                            total_bases[region_idx] += r_hi - r_lo;
                            if let Some(pp) = per_pos.as_mut() {
                                let lo = (r_lo - r.start) as usize;
                                let hi = (r_hi - r.start) as usize;
                                for p in &mut pp[region_idx][lo..hi] {
                                    *p = p.saturating_add(1);
                                }
                            }
                        }
                        pos = run_end_excl;
                    }
                    Kind::Skip => pos += len,
                    Kind::Insertion | Kind::SoftClip => {}
                    Kind::HardClip | Kind::Pad => {}
                }
            }
        }
    }
    let mean_depth: Vec<f64> = regions
        .iter()
        .zip(total_bases.iter())
        .map(|(r, &b)| (b as f64) / ((r.end - r.start + 1).max(1) as f64))
        .collect();
    Ok(PileupResult {
        mean_depth,
        total_match_bases: total_bases,
        per_position: per_pos,
    })
}

/// Stream the whole BAM, accumulating into all regions in one pass.
pub fn pileup_regions_streaming(
    path: &Path,
    regions: &[Region],
    min_mapq: u8,
    opts: PileupOpts,
) -> UltiResult<PileupResult> {
    use noodles::bam;
    use noodles::sam::alignment::record::cigar::op::Kind;

    // Build per-chromosome interval lookup; each interval's `val` is its index
    // in the caller's input vector.
    let mut by_chrom: HashMap<String, Vec<Interval<u64, usize>>> = HashMap::new();
    for (i, r) in regions.iter().enumerate() {
        by_chrom.entry(r.chrom.clone()).or_default().push(Interval {
            start: r.start,
            stop: r.end + 1, // Lapper uses half-open
            val: i,
        });
    }
    let lappers: HashMap<String, Lapper<u64, usize>> = by_chrom
        .into_iter()
        .map(|(c, ivs)| (c, Lapper::new(ivs)))
        .collect();

    let mut total_bases = vec![0u64; regions.len()];
    let mut per_pos: Option<Vec<Vec<u32>>> = if opts.record_per_position {
        Some(
            regions
                .iter()
                .map(|r| vec![0u32; (r.end - r.start + 1) as usize])
                .collect(),
        )
    } else {
        None
    };

    let mut reader = bam::io::reader::Builder
        .build_from_path(path)
        .map_err(|e| UltiError::Alignment {
            path: path.into(),
            message: format!("cannot open BAM for pileup: {e}"),
        })?;
    let header = reader.read_header().map_err(|e| UltiError::Alignment {
        path: path.into(),
        message: format!("cannot read header: {e}"),
    })?;

    for result in reader.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
        if mapq < min_mapq {
            continue;
        }
        let ref_id = match record.reference_sequence_id() {
            Some(Ok(id)) => id,
            _ => continue,
        };
        let chrom = match header.reference_sequences().get_index(ref_id) {
            Some((name, _)) => name.to_string(),
            None => continue,
        };
        let Some(lapper) = lappers.get(&chrom) else {
            continue;
        };
        let mut pos = match record.alignment_start() {
            Some(Ok(p)) => p.get() as u64,
            _ => continue,
        };

        for op_result in record.cigar().iter() {
            let Ok(op) = op_result else {
                break;
            };
            let len = op.len() as u64;
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                    let run_start = pos;
                    let run_end_excl = pos + len;
                    // Intersect with each region this run touches.
                    for hit in lapper.find(run_start, run_end_excl) {
                        let r_lo = hit.start.max(run_start);
                        let r_hi = hit.stop.min(run_end_excl);
                        if r_hi > r_lo {
                            total_bases[hit.val] += r_hi - r_lo;
                            if let Some(pp) = per_pos.as_mut() {
                                let region = &regions[hit.val];
                                let lo = (r_lo - region.start) as usize;
                                let hi = (r_hi - region.start) as usize;
                                for p in &mut pp[hit.val][lo..hi] {
                                    *p = p.saturating_add(1);
                                }
                            }
                        }
                    }
                    pos = run_end_excl;
                }
                Kind::Skip => {
                    pos += len; // intron — no depth contribution
                }
                Kind::Insertion | Kind::SoftClip => {} // query-only
                Kind::HardClip | Kind::Pad => {}
            }
        }
    }

    let mean_depth: Vec<f64> = regions
        .iter()
        .zip(total_bases.iter())
        .map(|(r, &b)| {
            let len = (r.end - r.start + 1).max(1) as f64;
            b as f64 / len
        })
        .collect();

    Ok(PileupResult {
        mean_depth,
        total_match_bases: total_bases,
        per_position: per_pos,
    })
}
