//! End-to-end-ish tests that do not require a BAM (which would need
//! crafting a binary BAM, out of scope here). They exercise the
//! annotation → event-detection → stats → output pipeline against a
//! synthetic in-memory junction matrix.

use std::collections::HashMap;
use std::io::Write;

use ultimadsecaller::annotation;
use ultimadsecaller::events::{self, EventKind};
use ultimadsecaller::junctions::JunctionMatrix;

fn write_gtf(s: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new().suffix(".gtf").tempfile().unwrap();
    f.write_all(s.as_bytes()).unwrap();
    f
}

#[test]
fn detects_canonical_skipped_exon() {
    // Two transcripts: T1 includes the middle exon, T2 skips it.
    let gtf = "\
chr1\ts\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ts\texon\t300\t350\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ts\texon\t500\t600\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ts\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";
chr1\ts\texon\t500\t600\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";
";
    let tmp = write_gtf(gtf);
    let ann = annotation::parse(tmp.path()).unwrap();
    let jm = JunctionMatrix {
        samples: vec!["s1".into()],
        counts: HashMap::new(),
    };
    let evs = events::detect_all(&ann, &jm);
    let se: Vec<_> = evs.iter().filter(|e| e.kind == EventKind::SE).collect();
    assert!(!se.is_empty(), "should detect at least one SE event");
    let ev = se[0];
    assert_eq!(ev.exons.len(), 3);
    assert_eq!(ev.exons[0].start, 100);
    assert_eq!(ev.exons[1].start, 300);
    assert_eq!(ev.exons[2].start, 500);
    assert_eq!(ev.inclusion_junctions, vec![(200, 300), (350, 500)]);
    assert_eq!(ev.exclusion_junctions, vec![(200, 500)]);
}

#[test]
fn detects_a5ss() {
    // Two donor variants of the same exon body (start = 100; ends 200 and 250),
    // both joined to a common downstream exon at 500-600.
    let gtf = "\
chr1\ts\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ts\texon\t500\t600\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ts\texon\t100\t250\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";
chr1\ts\texon\t500\t600\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";
";
    let tmp = write_gtf(gtf);
    let ann = annotation::parse(tmp.path()).unwrap();
    let jm = JunctionMatrix {
        samples: vec!["s1".into()],
        counts: HashMap::new(),
    };
    let evs = events::detect_all(&ann, &jm);
    let a5: Vec<_> = evs.iter().filter(|e| e.kind == EventKind::A5SS).collect();
    assert!(!a5.is_empty(), "should detect at least one A5SS event");
}

#[test]
fn fdr_monotonic_in_p() {
    use ultimadsecaller::stats::bh_fdr;
    let p = [0.001, 0.005, 0.02, 0.5];
    let adj = bh_fdr(&p);
    // BH-adjusted p-values must be non-decreasing in the original sort order.
    for w in adj.windows(2) {
        assert!(w[0] <= w[1] + 1e-12, "BH should be monotonic, got {adj:?}");
    }
}
