//! Splice-site motif analysis against a reference FASTA.
//!
//! For each junction `(chrom, donor_end, acceptor_start)` we look up the
//! 2 nt immediately *after* the donor exon and the 2 nt immediately
//! *before* the acceptor exon (intron sides). The canonical eukaryotic
//! splice signals are:
//!
//! | Motif      | Donor | Acceptor | Frequency |
//! |------------|-------|----------|-----------|
//! | U2 GT-AG   | `GT`  | `AG`     | ~99%      |
//! | U2 GC-AG   | `GC`  | `AG`     | ~0.5%     |
//! | U12 AT-AC  | `AT`  | `AC`     | ~0.05%    |
//!
//! On the `-` strand a junction with donor `CT` and acceptor `AC` is the
//! reverse complement of GT-AG; we report it as such rather than calling
//! it non-canonical. Motif results feed the cryptic-junction filter and
//! the per-event consequence annotation.

use std::collections::HashMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::{UltiError, UltiResult};
use crate::{Position, Strand};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SpliceMotif {
    GtAg,
    GcAg,
    AtAc,
    GtAgReverse,
    GcAgReverse,
    AtAcReverse,
    NonCanonical,
    Unknown,
}

impl SpliceMotif {
    pub fn is_canonical(self) -> bool {
        !matches!(self, SpliceMotif::NonCanonical | SpliceMotif::Unknown)
    }

    pub fn short(self) -> &'static str {
        match self {
            SpliceMotif::GtAg => "GT-AG",
            SpliceMotif::GcAg => "GC-AG",
            SpliceMotif::AtAc => "AT-AC",
            SpliceMotif::GtAgReverse => "CT-AC",
            SpliceMotif::GcAgReverse => "CT-GC",
            SpliceMotif::AtAcReverse => "GT-AT",
            SpliceMotif::NonCanonical => "non-canonical",
            SpliceMotif::Unknown => "unknown",
        }
    }

    pub fn from_dinucleotides(donor: &[u8], acceptor: &[u8]) -> Self {
        if donor.len() != 2 || acceptor.len() != 2 {
            return SpliceMotif::Unknown;
        }
        let d = [donor[0].to_ascii_uppercase(), donor[1].to_ascii_uppercase()];
        let a = [acceptor[0].to_ascii_uppercase(), acceptor[1].to_ascii_uppercase()];
        if d == [b'G', b'T'] && a == [b'A', b'G'] {
            return SpliceMotif::GtAg;
        }
        if d == [b'G', b'C'] && a == [b'A', b'G'] {
            return SpliceMotif::GcAg;
        }
        if d == [b'A', b'T'] && a == [b'A', b'C'] {
            return SpliceMotif::AtAc;
        }
        if d == [b'C', b'T'] && a == [b'A', b'C'] {
            return SpliceMotif::GtAgReverse;
        }
        if d == [b'C', b'T'] && a == [b'G', b'C'] {
            return SpliceMotif::GcAgReverse;
        }
        if d == [b'G', b'T'] && a == [b'A', b'T'] {
            return SpliceMotif::AtAcReverse;
        }
        SpliceMotif::NonCanonical
    }
}

/// In-memory reference container.
pub struct Reference {
    sequences: HashMap<String, Vec<u8>>,
}

impl Reference {
    pub fn load(path: &Path) -> UltiResult<Self> {
        use noodles::fasta;
        let mut reader = fasta::reader::Builder::default()
            .build_from_path(path)
            .map_err(|e| UltiError::reference(path, format!("cannot open: {e}")))?;
        let mut sequences = HashMap::new();
        for result in reader.records() {
            let record = result.map_err(|e| UltiError::reference(path, e.to_string()))?;
            let name = std::str::from_utf8(record.name())
                .map_err(|e| UltiError::reference(path, format!("non-UTF-8 chrom name: {e}")))?
                .to_string();
            let seq = record.sequence().as_ref().to_vec();
            sequences.insert(name, seq);
        }
        if sequences.is_empty() {
            return Err(UltiError::reference(path, "no sequences found"));
        }
        Ok(Reference { sequences })
    }

    /// 1-based inclusive `[start, end]` slice. Returns `None` if out of range.
    pub fn slice(&self, chrom: &str, start: Position, end: Position) -> Option<&[u8]> {
        let seq = self.sequences.get(chrom)?;
        if start == 0 || end == 0 || end < start {
            return None;
        }
        let lo = (start - 1) as usize;
        let hi = end as usize;
        if hi > seq.len() {
            return None;
        }
        Some(&seq[lo..hi])
    }

    pub fn chromosomes(&self) -> impl Iterator<Item = &String> {
        self.sequences.keys()
    }

    /// Classify a junction by reading the donor and acceptor dinucleotides.
    pub fn classify_junction(
        &self,
        chrom: &str,
        donor_end: Position,
        acceptor_start: Position,
    ) -> SpliceMotif {
        if acceptor_start <= donor_end + 2 {
            return SpliceMotif::Unknown;
        }
        let intron_first = donor_end + 1;
        let intron_last = acceptor_start - 1;
        let donor = self.slice(chrom, intron_first, intron_first + 1);
        let acceptor = self.slice(chrom, intron_last - 1, intron_last);
        match (donor, acceptor) {
            (Some(d), Some(a)) => SpliceMotif::from_dinucleotides(d, a),
            _ => SpliceMotif::Unknown,
        }
    }
}

/// Annotate every junction in a list with its motif classification.
pub fn classify_junctions(
    reference: &Reference,
    junctions: &[(String, Position, Position)],
) -> Vec<SpliceMotif> {
    junctions
        .iter()
        .map(|(c, d, a)| reference.classify_junction(c, *d, *a))
        .collect()
}

/// Reverse-complement a DNA sequence. Used by the protein consequence
/// walk when the transcript is on the `-` strand.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'N' => b'N',
            _ => b'N',
        })
        .collect()
}

/// Translate a coding sequence to a protein string. Standard genetic code;
/// ambiguous codons → `X`; stop codons → `*`.
pub fn translate(cds: &[u8]) -> String {
    let mut out = String::with_capacity(cds.len() / 3);
    for chunk in cds.chunks_exact(3) {
        out.push(codon_to_aa(chunk));
    }
    out
}

fn codon_to_aa(c: &[u8]) -> char {
    let c = [
        c[0].to_ascii_uppercase(),
        c[1].to_ascii_uppercase(),
        c[2].to_ascii_uppercase(),
    ];
    match &c {
        b"TTT" | b"TTC" => 'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
        b"ATT" | b"ATC" | b"ATA" => 'I',
        b"ATG" => 'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
        b"TAT" | b"TAC" => 'Y',
        b"TAA" | b"TAG" | b"TGA" => '*',
        b"CAT" | b"CAC" => 'H',
        b"CAA" | b"CAG" => 'Q',
        b"AAT" | b"AAC" => 'N',
        b"AAA" | b"AAG" => 'K',
        b"GAT" | b"GAC" => 'D',
        b"GAA" | b"GAG" => 'E',
        b"TGT" | b"TGC" => 'C',
        b"TGG" => 'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
        _ => 'X',
    }
}

/// Score bonus for canonical motifs — feeds the cryptic-junction
/// confidence formula.
pub fn motif_bonus(m: SpliceMotif) -> f64 {
    match m {
        SpliceMotif::GtAg | SpliceMotif::GtAgReverse => 1.0,
        SpliceMotif::GcAg | SpliceMotif::GcAgReverse => 0.6,
        SpliceMotif::AtAc | SpliceMotif::AtAcReverse => 0.4,
        SpliceMotif::NonCanonical => -0.5,
        SpliceMotif::Unknown => 0.0,
    }
}

/// Strand-aware donor / acceptor extraction — returns transcript-strand
/// dinucleotides for `-` strand genes.
pub fn transcript_dinucleotides(
    reference: &Reference,
    chrom: &str,
    donor_end: Position,
    acceptor_start: Position,
    strand: Strand,
) -> Option<([u8; 2], [u8; 2])> {
    if acceptor_start <= donor_end + 2 {
        return None;
    }
    let intron_first = donor_end + 1;
    let intron_last = acceptor_start - 1;
    let d = reference.slice(chrom, intron_first, intron_first + 1)?;
    let a = reference.slice(chrom, intron_last - 1, intron_last)?;
    let d2 = [d[0].to_ascii_uppercase(), d[1].to_ascii_uppercase()];
    let a2 = [a[0].to_ascii_uppercase(), a[1].to_ascii_uppercase()];
    if strand == Strand::Reverse {
        let drc = reverse_complement(&a2);
        let arc = reverse_complement(&d2);
        Some(([drc[0], drc[1]], [arc[0], arc[1]]))
    } else {
        Some((d2, a2))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classifies_canonical_motifs() {
        assert_eq!(SpliceMotif::from_dinucleotides(b"GT", b"AG"), SpliceMotif::GtAg);
        assert_eq!(SpliceMotif::from_dinucleotides(b"gc", b"ag"), SpliceMotif::GcAg);
        assert_eq!(SpliceMotif::from_dinucleotides(b"AT", b"AC"), SpliceMotif::AtAc);
        assert_eq!(
            SpliceMotif::from_dinucleotides(b"CT", b"AC"),
            SpliceMotif::GtAgReverse
        );
        assert_eq!(
            SpliceMotif::from_dinucleotides(b"AA", b"AA"),
            SpliceMotif::NonCanonical
        );
    }

    #[test]
    fn translates_standard_codons() {
        assert_eq!(translate(b"ATGGCATAA"), "MA*");
        assert_eq!(translate(b"TTTCCCNNN"), "FPX");
    }

    #[test]
    fn reverse_complement_basics() {
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
    }
}
