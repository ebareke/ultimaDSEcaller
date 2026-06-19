//! Annotation engine: parses GTF / GFF3 (plain or bgzipped, auto-detected
//! by extension) and builds per-gene splice graphs.
//!
//! ## Splice graph model
//!
//! For each gene:
//! * **Nodes** are unique `(start, end)` exonic intervals across all
//!   transcripts of that gene.
//! * **Edges** are annotated introns — a directed edge from exon *A* to exon
//!   *B* indicates that some transcript contains *A* immediately followed by
//!   *B* (i.e. they are consecutive exons of that transcript).
//! * Each edge stores the set of supporting transcript IDs, so downstream
//!   event detection can reason about transcript-level structure.
//! * For fast spatial lookup the gene also keeps a `Lapper` indexed by exon
//!   coordinates.
//!
//! This is the SUPPA/rMATS-canonical form. It deliberately does *not*
//! subdivide overlapping exons into "exonic parts" — that ambiguates AFE/ALE
//! and complicates A5SS/A3SS detection. Partial overlap is instead
//! represented by two distinct exon nodes sharing one endpoint.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use petgraph::graph::{DiGraph, NodeIndex};
use rust_lapper::{Interval, Lapper};
use serde::{Deserialize, Serialize};

use crate::error::{UltiError, UltiResult};
use crate::{Position, Strand};

/// An exon — the unit of a splice-graph node.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Exon {
    pub start: Position,
    pub end: Position,
}

impl Exon {
    #[inline]
    pub fn len(&self) -> Position {
        self.end.saturating_sub(self.start) + 1
    }

    /// True for an invalid (reversed) interval; a real exon spans >= 1 base.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.end < self.start
    }
}

/// Annotated splice junction (intron) between two exons.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotatedIntron {
    pub donor_end: Position,
    pub acceptor_start: Position,
    pub transcripts: HashSet<String>,
}

/// Per-gene splice graph.
pub struct SpliceGraph {
    pub gene_id: String,
    pub gene_name: Option<String>,
    pub chrom: String,
    pub strand: Strand,
    pub graph: DiGraph<Exon, AnnotatedIntron>,
    /// Lookup: exon → node index. Keyed on the exon tuple, so two exons that
    /// share a start coordinate but differ in end live on separate nodes.
    pub exon_index: HashMap<Exon, NodeIndex>,
    /// Spatial index for fast "which exons overlap [a, b]?" queries.
    pub lapper: Lapper<Position, NodeIndex>,
    /// Original per-transcript exon chains (ordered 5' → 3' on the *transcript*
    /// strand). Useful for transcript-aware event walks.
    pub transcripts: BTreeMap<String, Vec<Exon>>,
}

impl SpliceGraph {
    /// Returns all junctions present in the graph as `(donor_end, acceptor_start)`.
    pub fn junctions(&self) -> Vec<(Position, Position)> {
        let mut out = Vec::with_capacity(self.graph.edge_count());
        for e in self.graph.edge_indices() {
            let w = &self.graph[e];
            out.push((w.donor_end, w.acceptor_start));
        }
        out
    }
}

/// Container for the full annotation across all genes.
pub struct Annotation {
    pub genes: HashMap<String, SpliceGraph>,
    pub gene_by_chrom: HashMap<String, Vec<String>>,
}

impl Annotation {
    pub fn gene_count(&self) -> usize {
        self.genes.len()
    }
}

/// Parse a GTF or GFF3 file (extension-detected) and build the per-gene
/// splice graphs.
pub fn parse(path: &Path) -> UltiResult<Annotation> {
    let format = detect_format(path)?;
    let reader = open_reader(path)?;
    let raw = read_exons(reader, format, path)?;
    Ok(build_annotation(raw))
}

#[derive(Debug, Clone, Copy)]
enum Format {
    Gtf,
    Gff3,
}

fn detect_format(path: &Path) -> UltiResult<Format> {
    let mut name = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    // Strip .gz so "foo.gtf.gz" still detects as GTF.
    if let Some(stripped) = name.strip_suffix(".gz") {
        name = stripped.to_string();
    }
    if name.ends_with(".gtf") {
        Ok(Format::Gtf)
    } else if name.ends_with(".gff") || name.ends_with(".gff3") {
        Ok(Format::Gff3)
    } else {
        Err(UltiError::Annotation {
            path: path.into(),
            line: 0,
            message: format!(
                "cannot infer annotation format from extension; expected .gtf, .gff or .gff3 \
                 (got `{name}`)"
            ),
        })
    }
}

fn open_reader(path: &Path) -> UltiResult<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|e| UltiError::io(path, e))?;
    let gz = path
        .file_name()
        .and_then(|s| s.to_str())
        .map(|n| n.to_ascii_lowercase().ends_with(".gz"))
        .unwrap_or(false);
    if gz {
        // bgzipped is more common in bioinformatics than vanilla gzip, but noodles' bgzf
        // reader handles both bgzf and plain gzip transparently.
        let bgz = noodles::bgzf::Reader::new(file);
        Ok(Box::new(BufReader::new(bgz)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Intermediate per-transcript exon record harvested from the file before we
/// reorganize into the graph.
#[derive(Debug, Clone)]
struct RawExon {
    chrom: String,
    start: Position,
    end: Position,
    strand: Strand,
    gene_id: String,
    gene_name: Option<String>,
    transcript_id: String,
}

fn read_exons<R: BufRead>(reader: R, format: Format, path: &Path) -> UltiResult<Vec<RawExon>> {
    let mut out = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| UltiError::io(path, e))?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t');
        let chrom = next_field(&mut fields, i, "seqname", path)?;
        let _source = fields.next();
        let feat = next_field(&mut fields, i, "feature", path)?;
        if feat != "exon" {
            continue;
        }
        let start: Position = next_field(&mut fields, i, "start", path)?
            .parse()
            .map_err(|_| ann_err(path, i, "start is not an integer"))?;
        let end: Position = next_field(&mut fields, i, "end", path)?
            .parse()
            .map_err(|_| ann_err(path, i, "end is not an integer"))?;
        let _score = fields.next();
        let strand_str = next_field(&mut fields, i, "strand", path)?;
        let strand = Strand::from_char(strand_str.chars().next().unwrap_or('.'));
        let _frame = fields.next();
        let attrs = fields.next().unwrap_or("");
        let (gene_id, gene_name, tx_id) = match format {
            Format::Gtf => parse_gtf_attrs(attrs),
            Format::Gff3 => parse_gff3_attrs(attrs),
        };
        let gene_id =
            gene_id.ok_or_else(|| ann_err(path, i, "exon record missing gene_id / Parent"))?;
        let tx_id =
            tx_id.ok_or_else(|| ann_err(path, i, "exon record missing transcript_id / Parent"))?;
        out.push(RawExon {
            chrom: chrom.to_string(),
            start,
            end,
            strand,
            gene_id,
            gene_name,
            transcript_id: tx_id,
        });
    }
    Ok(out)
}

fn next_field<'a>(
    it: &mut std::str::Split<'a, char>,
    line_no: usize,
    name: &str,
    path: &Path,
) -> UltiResult<&'a str> {
    it.next()
        .ok_or_else(|| ann_err(path, line_no, &format!("missing field `{name}`")))
}

fn ann_err(path: &Path, line: usize, message: &str) -> UltiError {
    UltiError::Annotation {
        path: path.into(),
        line,
        message: message.into(),
    }
}

/// Parse GTF-style attributes: `gene_id "ENSG..."; transcript_id "ENST..."; ...`.
fn parse_gtf_attrs(attrs: &str) -> (Option<String>, Option<String>, Option<String>) {
    let mut gene_id = None;
    let mut gene_name = None;
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
            "gene_name" => gene_name = Some(v),
            "transcript_id" => tx_id = Some(v),
            _ => {}
        }
    }
    (gene_id, gene_name, tx_id)
}

/// Parse GFF3-style attributes: `ID=...;Parent=...;gene_name=...`.
fn parse_gff3_attrs(attrs: &str) -> (Option<String>, Option<String>, Option<String>) {
    let mut gene_id = None;
    let mut gene_name = None;
    let mut tx_id = None;
    let mut parent = None;
    for part in attrs.split(';') {
        let (k, v) = match part.split_once('=') {
            Some(kv) => kv,
            None => continue,
        };
        match k.trim() {
            "gene_id" => gene_id = Some(v.to_string()),
            "gene_name" | "Name" => gene_name = Some(v.to_string()),
            "transcript_id" => tx_id = Some(v.to_string()),
            "Parent" => parent = Some(v.to_string()),
            _ => {}
        }
    }
    // GFF3 exon records typically have `Parent=transcript:ENST...` rather than
    // explicit gene_id/transcript_id. Fall back to deriving them.
    let tx_id = tx_id.or_else(|| parent.clone());
    // We don't have gene_id directly in a GFF3 exon line — that requires a
    // second pass linking Parent transcripts to their genes. For now, when
    // gene_id is unset, we use the transcript_id as a degenerate proxy. A
    // future improvement is to build the transcript→gene map from the file's
    // mRNA records.
    let gene_id = gene_id.or_else(|| tx_id.clone());
    (gene_id, gene_name, tx_id)
}

fn build_annotation(raw: Vec<RawExon>) -> Annotation {
    // Restructure RawExon -> per-gene "GeneSeed" so build_from_seeds can be
    // shared with the cache path.
    let mut by_gene: HashMap<String, Vec<RawExon>> = HashMap::new();
    for r in raw {
        by_gene.entry(r.gene_id.clone()).or_default().push(r);
    }
    let mut seeds: Vec<GeneSeed> = Vec::with_capacity(by_gene.len());
    for (gene_id, exons) in by_gene {
        let chrom = exons[0].chrom.clone();
        let strand = exons[0].strand;
        let gene_name = exons.iter().find_map(|e| e.gene_name.clone());
        let mut transcripts: BTreeMap<String, Vec<Exon>> = BTreeMap::new();
        for e in &exons {
            transcripts
                .entry(e.transcript_id.clone())
                .or_default()
                .push(Exon {
                    start: e.start,
                    end: e.end,
                });
        }
        for chain in transcripts.values_mut() {
            chain.sort();
            chain.dedup();
        }
        seeds.push(GeneSeed {
            gene_id,
            gene_name,
            chrom,
            strand,
            transcripts,
        });
    }
    build_from_seeds(seeds)
}

/// Shared seed used by both the GTF/GFF3 parser and the cache loader.
pub(crate) struct GeneSeed {
    pub gene_id: String,
    pub gene_name: Option<String>,
    pub chrom: String,
    pub strand: Strand,
    pub transcripts: BTreeMap<String, Vec<Exon>>,
}

/// Entry point used by the cache loader to rebuild a live Annotation.
pub fn build_from_cached(genes: Vec<crate::cache::CachedGene>) -> Annotation {
    let seeds = genes
        .into_iter()
        .map(|g| GeneSeed {
            gene_id: g.gene_id,
            gene_name: g.gene_name,
            chrom: g.chrom,
            strand: g.strand,
            transcripts: g.transcripts,
        })
        .collect();
    build_from_seeds(seeds)
}

fn build_from_seeds(seeds: Vec<GeneSeed>) -> Annotation {
    let mut genes = HashMap::with_capacity(seeds.len());
    let mut gene_by_chrom: HashMap<String, Vec<String>> = HashMap::new();

    for GeneSeed {
        gene_id,
        gene_name,
        chrom,
        strand,
        transcripts,
    } in seeds
    {
        // 2) Collect unique exons (by tuple) — these become nodes.
        let unique_exons: HashSet<Exon> = transcripts
            .values()
            .flat_map(|chain| chain.iter().copied())
            .collect();

        // 3) Build the directed splice graph.
        let mut graph: DiGraph<Exon, AnnotatedIntron> = DiGraph::new();
        let mut exon_index: HashMap<Exon, NodeIndex> = HashMap::with_capacity(unique_exons.len());
        for ex in &unique_exons {
            let idx = graph.add_node(*ex);
            exon_index.insert(*ex, idx);
        }

        // 4) Add intron edges from consecutive exons of each transcript.
        // Genomic-order coords are kept; for "-" strand transcripts the
        // *transcript* order is reversed, but the donor/acceptor labels are
        // still correct because donor = upstream-in-transcript = downstream-in-
        // genome on the minus strand. To keep the graph genome-oriented we
        // always add edges in genomic-coordinate order.
        let mut edge_pool: HashMap<(NodeIndex, NodeIndex, Position, Position), HashSet<String>> =
            HashMap::new();
        for (tx_id, chain) in &transcripts {
            for window in chain.windows(2) {
                let (a, b) = (window[0], window[1]);
                let (ai, bi) = (exon_index[&a], exon_index[&b]);
                let donor_end = a.end;
                let acceptor_start = b.start;
                edge_pool
                    .entry((ai, bi, donor_end, acceptor_start))
                    .or_default()
                    .insert(tx_id.clone());
            }
        }
        for ((ai, bi, de, as_), txs) in edge_pool {
            graph.add_edge(
                ai,
                bi,
                AnnotatedIntron {
                    donor_end: de,
                    acceptor_start: as_,
                    transcripts: txs,
                },
            );
        }

        // 5) Spatial index.
        let intervals: Vec<Interval<Position, NodeIndex>> = unique_exons
            .iter()
            .map(|ex| Interval {
                start: ex.start,
                stop: ex.end + 1,
                val: exon_index[ex],
            })
            .collect();
        let lapper = Lapper::new(intervals);

        gene_by_chrom
            .entry(chrom.clone())
            .or_default()
            .push(gene_id.clone());

        genes.insert(
            gene_id.clone(),
            SpliceGraph {
                gene_id,
                gene_name,
                chrom,
                strand,
                graph,
                exon_index,
                lapper,
                transcripts,
            },
        );
    }

    Annotation {
        genes,
        gene_by_chrom,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_gtf(s: &str) -> tempfile::NamedTempFile {
        use std::io::Write;
        let mut f = tempfile::Builder::new().suffix(".gtf").tempfile().unwrap();
        f.write_all(s.as_bytes()).unwrap();
        f
    }

    #[test]
    fn shared_start_does_not_collide() {
        // Two exons with the same start but different ends — the prototype
        // would have lost one of them. Here they should be separate nodes.
        let gtf = "\
chr1\ts\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ts\texon\t100\t250\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";
chr1\ts\texon\t300\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ts\texon\t300\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";
";
        let tmp = write_gtf(gtf);
        let ann = parse(tmp.path()).unwrap();
        let g = ann.genes.get("G1").unwrap();
        assert_eq!(g.graph.node_count(), 3);
        assert_eq!(g.graph.edge_count(), 2);
    }
}
