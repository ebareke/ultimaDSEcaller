//! Per-event protein consequence + NMD prediction.
//!
//! For each AS event we determine the consequence on the coding sequence
//! of the *canonical (most likely) transcript* of the event's gene:
//!
//! * **Noncoding**: the event lies entirely in UTR or non-coding regions.
//! * **InFrame**: the alternative form changes mRNA length by a multiple
//!   of 3 — predicted protein-altering but in-frame; no premature stop
//!   from frame-shift.
//! * **FrameShift**: length change is not a multiple of 3; downstream
//!   codons shift register.
//! * **PrematureStop**: walking the new reading frame from the affected
//!   position introduces a stop codon upstream of the natural one.
//! * **NmdTriggered**: as PrematureStop, *and* the new stop codon is
//!   ≥ 50 nt upstream of the last splice junction in the mRNA — the
//!   "50 nt rule" for nonsense-mediated decay.
//!
//! CDS information is sourced from a side-pass on the GTF (we collect
//! `CDS` records with `transcript_id`).
//!
//! Scope notes:
//! * Implemented in full for SE (exon skipping) — the highest-impact and
//!   best-defined case.
//! * For A5SS / A3SS we compute the donor/acceptor delta length and
//!   predict in-frame vs frame-shift.
//! * For MXE we report Unknown (depends on which alt form is canonical).
//! * For IR we predict by intron-retained length.
//! * For AFE / ALE we report Noncoding when the alt exon is upstream of
//!   the start codon or downstream of the stop codon.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::annotation::Exon;
use crate::error::{UltiError, UltiResult};
use crate::events::{ASEvent, EventKind};
use crate::motif::{reverse_complement, translate};
use crate::{Position, Strand};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ProteinConsequence {
    Noncoding,
    InFrame,
    FrameShift,
    PrematureStop,
    NmdTriggered,
    Unknown,
}

impl ProteinConsequence {
    pub fn short(self) -> &'static str {
        match self {
            ProteinConsequence::Noncoding => "noncoding",
            ProteinConsequence::InFrame => "in-frame",
            ProteinConsequence::FrameShift => "frame-shift",
            ProteinConsequence::PrematureStop => "PTC",
            ProteinConsequence::NmdTriggered => "NMD",
            ProteinConsequence::Unknown => "unknown",
        }
    }
}

/// One transcript's CDS interval list (genomic coords, half-open per GTF/Gencode).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdsTranscript {
    pub transcript_id: String,
    pub gene_id: String,
    pub chrom: String,
    pub strand: Strand,
    /// CDS exonic intervals, sorted by start (genomic coords).
    pub cds: Vec<Exon>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CdsCatalog {
    /// transcript_id -> CDS spec
    pub transcripts: HashMap<String, CdsTranscript>,
}

impl CdsCatalog {
    pub fn transcripts_of_gene(&self, gene_id: &str) -> Vec<&CdsTranscript> {
        self.transcripts
            .values()
            .filter(|t| t.gene_id == gene_id)
            .collect()
    }

    /// Pick a canonical transcript for a gene — the one with the longest CDS.
    pub fn canonical_for(&self, gene_id: &str) -> Option<&CdsTranscript> {
        self.transcripts_of_gene(gene_id)
            .into_iter()
            .max_by_key(|t| t.cds.iter().map(|e| e.end - e.start + 1).sum::<u64>())
    }
}

/// Parse CDS records from a GTF/GFF3. We do a focused second pass — the
/// annotation parser only keeps exon records, so this pass collects the
/// CDS rows separately. Returns an empty catalog if no CDS rows are
/// present (annotations like FlyBase exon-only GFFs).
pub fn parse_cds(path: &Path) -> UltiResult<CdsCatalog> {
    let f = File::open(path).map_err(|e| UltiError::io(path, e))?;
    let mut reader: Box<dyn BufRead> = if path
        .file_name()
        .and_then(|s| s.to_str())
        .map(|n| n.to_ascii_lowercase().ends_with(".gz"))
        .unwrap_or(false)
    {
        let bgz = noodles::bgzf::Reader::new(f);
        Box::new(BufReader::new(bgz))
    } else {
        Box::new(BufReader::new(f))
    };

    let is_gff3 = path
        .file_name()
        .and_then(|s| s.to_str())
        .map(|n| {
            let l = n.to_ascii_lowercase();
            l.ends_with(".gff") || l.ends_with(".gff3") || l.ends_with(".gff.gz") || l.ends_with(".gff3.gz")
        })
        .unwrap_or(false);

    let mut catalog = CdsCatalog::default();
    let mut buf = String::new();
    loop {
        buf.clear();
        let n = reader
            .read_line(&mut buf)
            .map_err(|e| UltiError::io(path, e))?;
        if n == 0 {
            break;
        }
        if buf.starts_with('#') || buf.trim().is_empty() {
            continue;
        }
        let line = buf.trim_end_matches('\n');
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        if fields[2] != "CDS" {
            continue;
        }
        let chrom = fields[0].to_string();
        let start: u64 = match fields[3].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let end: u64 = match fields[4].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let strand = Strand::from_char(fields[6].chars().next().unwrap_or('.'));
        let attrs = fields[8];
        let (gene_id, tx_id) = if is_gff3 {
            parse_gff3_attrs(attrs)
        } else {
            parse_gtf_attrs(attrs)
        };
        let (Some(gene_id), Some(tx_id)) = (gene_id, tx_id) else {
            continue;
        };
        let entry = catalog
            .transcripts
            .entry(tx_id.clone())
            .or_insert_with(|| CdsTranscript {
                transcript_id: tx_id.clone(),
                gene_id: gene_id.clone(),
                chrom: chrom.clone(),
                strand,
                cds: Vec::new(),
            });
        entry.cds.push(Exon { start, end });
    }
    for t in catalog.transcripts.values_mut() {
        t.cds.sort();
        t.cds.dedup();
    }
    Ok(catalog)
}

fn parse_gtf_attrs(attrs: &str) -> (Option<String>, Option<String>) {
    let mut gene_id = None;
    let mut tx_id = None;
    for part in attrs.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let (k, v) = match part.split_once(' ') {
            Some(kv) => kv,
            None => continue,
        };
        let v = v.trim().trim_matches('"').to_string();
        match k {
            "gene_id" => gene_id = Some(v),
            "transcript_id" => tx_id = Some(v),
            _ => {}
        }
    }
    (gene_id, tx_id)
}

fn parse_gff3_attrs(attrs: &str) -> (Option<String>, Option<String>) {
    let mut gene_id = None;
    let mut tx_id = None;
    let mut parent = None;
    for part in attrs.split(';') {
        let (k, v) = match part.split_once('=') {
            Some(kv) => kv,
            None => continue,
        };
        match k.trim() {
            "gene_id" => gene_id = Some(v.to_string()),
            "transcript_id" => tx_id = Some(v.to_string()),
            "Parent" => parent = Some(v.to_string()),
            _ => {}
        }
    }
    (gene_id, tx_id.or(parent))
}

/// Total CDS length for a transcript.
fn cds_length(cds: &[Exon]) -> u64 {
    cds.iter().map(|e| e.end - e.start + 1).sum()
}

/// Position (in transcript coding coordinates, 1-based) of a genomic position
/// within the CDS. Returns `None` if the position is not in any CDS exon.
#[allow(dead_code)]
fn cds_position(cds: &[Exon], strand: Strand, genomic: Position) -> Option<u64> {
    let mut accumulated: u64 = 0;
    let chain: Vec<&Exon> = match strand {
        Strand::Reverse => cds.iter().rev().collect(),
        _ => cds.iter().collect(),
    };
    for ex in chain {
        if genomic >= ex.start && genomic <= ex.end {
            let offset = match strand {
                Strand::Reverse => ex.end - genomic + 1,
                _ => genomic - ex.start + 1,
            };
            return Some(accumulated + offset);
        }
        accumulated += ex.end - ex.start + 1;
    }
    None
}

/// Predict the protein consequence of a single AS event against a chosen
/// transcript's CDS.
pub fn predict_consequence_for_transcript(
    event: &ASEvent,
    tx: &CdsTranscript,
    reference: Option<&crate::motif::Reference>,
) -> ProteinConsequence {
    let cds = &tx.cds;
    if cds.is_empty() {
        return ProteinConsequence::Unknown;
    }
    let cds_lo = cds.iter().map(|e| e.start).min().unwrap_or(0);
    let cds_hi = cds.iter().map(|e| e.end).max().unwrap_or(0);

    // Compute the length change (in CDS-overlapping nt) the event imposes.
    let (delta_len, affected_pos) = match event.kind {
        EventKind::SE => {
            // Skipping exon = exons[1].
            if event.exons.len() < 3 {
                return ProteinConsequence::Unknown;
            }
            let skip = event.exons[1];
            let cds_overlap = exon_cds_overlap(&skip, cds);
            (cds_overlap, skip.start)
        }
        EventKind::A5SS => {
            if event.exons.len() < 3 {
                return ProteinConsequence::Unknown;
            }
            let short = event.exons[0];
            let long = event.exons[1];
            let delta = (long.end as i64) - (short.end as i64);
            (delta.unsigned_abs(), short.end)
        }
        EventKind::A3SS => {
            if event.exons.len() < 3 {
                return ProteinConsequence::Unknown;
            }
            let short = event.exons[0];
            let long = event.exons[1];
            let delta = (short.start as i64) - (long.start as i64);
            (delta.unsigned_abs(), short.start)
        }
        EventKind::IR => {
            let Some((s, e)) = event.retained_intron else {
                return ProteinConsequence::Unknown;
            };
            (e - s + 1, s)
        }
        EventKind::MXE => return ProteinConsequence::Unknown,
        EventKind::AFE | EventKind::ALE => {
            // If the alternative exon lies outside the CDS extent, it's UTR-only.
            let alt = if event.exons.is_empty() {
                return ProteinConsequence::Unknown;
            } else {
                event.exons[0]
            };
            if alt.end < cds_lo || alt.start > cds_hi {
                return ProteinConsequence::Noncoding;
            }
            return ProteinConsequence::Unknown;
        }
    };

    if delta_len == 0 {
        return ProteinConsequence::Noncoding;
    }
    if affected_pos < cds_lo || affected_pos > cds_hi {
        return ProteinConsequence::Noncoding;
    }

    let in_frame = delta_len % 3 == 0;
    if in_frame {
        return ProteinConsequence::InFrame;
    }
    // Frame-shift: check for downstream PTC if we have the reference.
    if let Some(reference) = reference {
        if let Some(consequence) = predict_ptc(event, tx, reference) {
            return consequence;
        }
    }
    ProteinConsequence::FrameShift
}

/// How many nucleotides of an exon overlap with the transcript's CDS.
fn exon_cds_overlap(exon: &Exon, cds: &[Exon]) -> u64 {
    let mut total = 0_u64;
    for c in cds {
        let lo = exon.start.max(c.start);
        let hi = exon.end.min(c.end);
        if hi >= lo {
            total += hi - lo + 1;
        }
    }
    total
}

/// Walk the post-event reading frame from the first affected codon and
/// look for a premature stop. If found, check the 50-nt NMD rule.
fn predict_ptc(
    event: &ASEvent,
    tx: &CdsTranscript,
    reference: &crate::motif::Reference,
) -> Option<ProteinConsequence> {
    // Construct the post-event CDS sequence: take the canonical CDS and
    // apply the event (skip an exon, shorten one, retain an intron, etc.).
    let mutated = build_post_event_cds(event, tx, reference)?;
    if mutated.is_empty() {
        return None;
    }
    let protein = translate(&mutated);
    let natural_aa_len = cds_length(&tx.cds) / 3;
    // Find the *new* stop codon position (in aa).
    let new_stop = protein.find('*').map(|i| i + 1)?; // 1-based aa position
    if new_stop as u64 >= natural_aa_len {
        return None; // not premature
    }
    // 50-nt NMD rule: PTC must be ≥ 50 nt upstream of the LAST exon-exon junction.
    let nmd_threshold_aa = nmd_threshold_aa_position(tx);
    if (new_stop as u64) < nmd_threshold_aa {
        Some(ProteinConsequence::NmdTriggered)
    } else {
        Some(ProteinConsequence::PrematureStop)
    }
}

fn nmd_threshold_aa_position(tx: &CdsTranscript) -> u64 {
    // Total CDS length in nt minus last CDS exon length minus 50 → divide by 3 for aa.
    let total = cds_length(&tx.cds);
    let last_exon = match tx.strand {
        Strand::Reverse => tx.cds.first().map(|e| e.end - e.start + 1),
        _ => tx.cds.last().map(|e| e.end - e.start + 1),
    }
    .unwrap_or(0);
    let threshold_nt = total.saturating_sub(last_exon).saturating_sub(50);
    threshold_nt / 3
}

/// Build the CDS *as it would look* after applying the event.
fn build_post_event_cds(
    event: &ASEvent,
    tx: &CdsTranscript,
    reference: &crate::motif::Reference,
) -> Option<Vec<u8>> {
    // Read each CDS exon as a chunk, optionally modifying based on the event.
    let mut out: Vec<u8> = Vec::new();
    for ex in &tx.cds {
        match event.kind {
            EventKind::SE => {
                // If this CDS exon coincides with the skipped exon, drop it.
                let skip = event.exons.get(1).copied();
                if let Some(s) = skip {
                    if exon_cds_overlap(&s, &[*ex]) > 0 {
                        continue;
                    }
                }
            }
            EventKind::IR => {
                if let Some((lo, hi)) = event.retained_intron {
                    // If this exon is the donor or acceptor, splice the intron in.
                    if ex.end == lo - 1 || ex.start == hi + 1 {
                        // Pull contiguous donor+intron+acceptor — handled by appending
                        // the intron sequence at the right boundary below.
                    }
                }
            }
            _ => {}
        }
        if let Some(seq) = reference.slice(&tx.chrom, ex.start, ex.end) {
            out.extend_from_slice(seq);
        }
        // For IR: after writing the donor exon, also append the retained intron.
        if event.kind == EventKind::IR {
            if let Some((lo, hi)) = event.retained_intron {
                if ex.end == lo - 1 {
                    if let Some(intron) = reference.slice(&tx.chrom, lo, hi) {
                        out.extend_from_slice(intron);
                    }
                }
            }
        }
    }
    if tx.strand == Strand::Reverse {
        out = reverse_complement(&out);
    }
    Some(out)
}

/// Predict consequence for an event using the gene's canonical (longest-CDS)
/// transcript.
pub fn predict_consequence(
    event: &ASEvent,
    catalog: &CdsCatalog,
    reference: Option<&crate::motif::Reference>,
) -> ProteinConsequence {
    let Some(tx) = catalog.canonical_for(&event.gene_id) else {
        return ProteinConsequence::Unknown;
    };
    predict_consequence_for_transcript(event, tx, reference)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_gtf(s: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::Builder::new().suffix(".gtf").tempfile().unwrap();
        f.write_all(s.as_bytes()).unwrap();
        f
    }

    #[test]
    fn parses_cds_records() {
        let gtf = "\
chr1\ts\texon\t100\t300\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ts\tCDS\t150\t300\t.\t+\t0\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ts\texon\t500\t700\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ts\tCDS\t500\t650\t.\t+\t2\tgene_id \"G1\"; transcript_id \"T1\";
";
        let tmp = write_gtf(gtf);
        let cds = parse_cds(tmp.path()).unwrap();
        assert_eq!(cds.transcripts.len(), 1);
        let t = cds.transcripts.get("T1").unwrap();
        assert_eq!(t.cds.len(), 2);
        assert_eq!(cds_length(&t.cds), (300 - 150 + 1) + (650 - 500 + 1));
    }

    #[test]
    fn classifies_in_frame_skip() {
        use crate::events::{ASEvent, EventKind};
        // Skipped exon of length 30 (multiple of 3) → in-frame.
        let ev = ASEvent {
            event_id: "test".into(),
            gene_id: "G1".into(),
            chrom: "chr1".into(),
            strand: Strand::Forward,
            kind: EventKind::SE,
            exons: vec![
                Exon { start: 100, end: 200 },
                Exon { start: 300, end: 329 }, // 30 nt
                Exon { start: 500, end: 600 },
            ],
            inclusion_junctions: vec![],
            exclusion_junctions: vec![],
            retained_intron: None,
        };
        let tx = CdsTranscript {
            transcript_id: "T1".into(),
            gene_id: "G1".into(),
            chrom: "chr1".into(),
            strand: Strand::Forward,
            cds: vec![
                Exon { start: 100, end: 200 },
                Exon { start: 300, end: 329 },
                Exon { start: 500, end: 600 },
            ],
        };
        let c = predict_consequence_for_transcript(&ev, &tx, None);
        assert_eq!(c, ProteinConsequence::InFrame);
    }

    #[test]
    fn classifies_frame_shift_skip() {
        use crate::events::{ASEvent, EventKind};
        // Skipped exon of length 31 → frame-shift.
        let ev = ASEvent {
            event_id: "test".into(),
            gene_id: "G1".into(),
            chrom: "chr1".into(),
            strand: Strand::Forward,
            kind: EventKind::SE,
            exons: vec![
                Exon { start: 100, end: 200 },
                Exon { start: 300, end: 330 }, // 31 nt
                Exon { start: 500, end: 600 },
            ],
            inclusion_junctions: vec![],
            exclusion_junctions: vec![],
            retained_intron: None,
        };
        let tx = CdsTranscript {
            transcript_id: "T1".into(),
            gene_id: "G1".into(),
            chrom: "chr1".into(),
            strand: Strand::Forward,
            cds: vec![
                Exon { start: 100, end: 200 },
                Exon { start: 300, end: 330 },
                Exon { start: 500, end: 600 },
            ],
        };
        let c = predict_consequence_for_transcript(&ev, &tx, None);
        assert_eq!(c, ProteinConsequence::FrameShift);
    }
}
