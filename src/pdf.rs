//! PDF report rendering (optional, behind the `pdf` feature).
//!
//! Produces a multi-page PDF document that mirrors the HTML report's
//! tabular content. Layout choices:
//!
//! 1. **Cover page**: title + run metadata + summary statistics.
//! 2. **Per-event-type pages**: count cards + top 20 events by adjusted p-value.
//! 3. **High-confidence set**: full table.
//! 4. **Consensus + protein consequence**: summary tables.
//!
//! Plots are *not* embedded — the HTML report is the right tool for
//! interactive plots, and a Plotly snapshot would need a headless
//! browser. Users who need plots in PDFs can print the HTML report from
//! a browser to get the full interactive view as a static PDF.
//!
//! We use `printpdf` rather than a wkhtmltopdf-style external tool so
//! the entire pipeline stays a single static binary.

use std::path::Path;

use printpdf::{BuiltinFont, Mm, PdfDocument, PdfDocumentReference};

use crate::error::{UltiError, UltiResult};
use crate::output::ReportPayload;

const PAGE_W_MM: f32 = 210.0; // A4
const PAGE_H_MM: f32 = 297.0;
const MARGIN_MM: f32 = 15.0;
const LINE_HEIGHT_MM: f32 = 5.0;
const TITLE_SIZE: f32 = 18.0;
const HEADING_SIZE: f32 = 12.0;
const BODY_SIZE: f32 = 9.0;

pub fn render(payload: &ReportPayload, out: &Path) -> UltiResult<()> {
    let (doc, page1, layer1) = PdfDocument::new(
        "ultimaDSEcaller Report",
        Mm(PAGE_W_MM),
        Mm(PAGE_H_MM),
        "page-1",
    );
    let _font = doc
        .add_builtin_font(BuiltinFont::Helvetica)
        .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
    let bold = doc
        .add_builtin_font(BuiltinFont::HelveticaBold)
        .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
    let mono = doc
        .add_builtin_font(BuiltinFont::Courier)
        .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;

    // -------------------------------- Cover --------------------------------
    let mut y = PAGE_H_MM - MARGIN_MM;
    let layer = doc.get_page(page1).get_layer(layer1);
    layer.use_text(
        "ultimaDSEcaller — Splicing Analysis Report",
        TITLE_SIZE,
        Mm(MARGIN_MM),
        Mm(y),
        &bold,
    );
    y -= 12.0;

    let run = &payload.run;
    let meta_lines = vec![
        format!("Annotation: {}", run.annotation_path.display()),
        format!("Samples:    {}", run.n_samples),
        format!("Groups:     {}", run.n_groups),
        format!(
            "Contrast:   {}",
            run.contrast.clone().unwrap_or_else(|| "(none)".into())
        ),
        format!("Technology: {}", run.tech),
        format!("Total events:           {}", payload.summary.events_total),
        format!(
            "High-confidence events: {}",
            payload.summary.events_high_confidence
        ),
    ];
    for line in &meta_lines {
        layer.use_text(line, BODY_SIZE, Mm(MARGIN_MM), Mm(y), &mono);
        y -= LINE_HEIGHT_MM;
    }
    y -= 6.0;
    layer.use_text(
        "Counts by event kind",
        HEADING_SIZE,
        Mm(MARGIN_MM),
        Mm(y),
        &bold,
    );
    y -= LINE_HEIGHT_MM + 2.0;
    for k in ["SE", "MXE", "A5SS", "A3SS", "IR", "AFE", "ALE"] {
        let s = payload.per_kind.get(k);
        let line = match s {
            Some(s) => format!(
                "  {k:<5}  total {:>6}   significant {:>5}   mean ΔPSI {:>+7.3}",
                s.total, s.significant, s.mean_delta_psi
            ),
            None => format!("  {k:<5}  (none)"),
        };
        layer.use_text(line, BODY_SIZE, Mm(MARGIN_MM), Mm(y), &mono);
        y -= LINE_HEIGHT_MM;
    }

    // ------------------------ Top events per kind --------------------------
    for k in ["SE", "MXE", "A5SS", "A3SS", "IR", "AFE", "ALE"] {
        let rows = match payload.top_events.get(k) {
            Some(r) if !r.is_empty() => r,
            _ => continue,
        };
        let (page, layer_id) = doc.add_page(Mm(PAGE_W_MM), Mm(PAGE_H_MM), format!("top-{k}"));
        let layer = doc.get_page(page).get_layer(layer_id);
        let mut y = PAGE_H_MM - MARGIN_MM;
        layer.use_text(
            format!("Top events — {k}"),
            HEADING_SIZE + 2.0,
            Mm(MARGIN_MM),
            Mm(y),
            &bold,
        );
        y -= 8.0;
        layer.use_text(
            "event_id            gene_id        ΔPSI    adj.p     conf.",
            BODY_SIZE,
            Mm(MARGIN_MM),
            Mm(y),
            &bold,
        );
        y -= LINE_HEIGHT_MM;
        for r in rows.iter().take(20) {
            let line = format!(
                "{:<20} {:<14} {:>+6.3}  {:<9} {:>5.3}",
                truncate(&r.event_id, 20),
                truncate(&r.gene_id, 14),
                r.delta_psi,
                fmt_p(r.adjusted_p_value),
                r.confidence,
            );
            layer.use_text(line, BODY_SIZE, Mm(MARGIN_MM), Mm(y), &mono);
            y -= LINE_HEIGHT_MM;
            if y < MARGIN_MM {
                break;
            }
        }
    }

    // ------------------------ High-confidence list -------------------------
    if !payload.high_confidence.is_empty() {
        let (page, layer_id) = doc.add_page(Mm(PAGE_W_MM), Mm(PAGE_H_MM), "high-conf");
        let layer = doc.get_page(page).get_layer(layer_id);
        let mut y = PAGE_H_MM - MARGIN_MM;
        layer.use_text(
            "High-confidence event set",
            HEADING_SIZE + 2.0,
            Mm(MARGIN_MM),
            Mm(y),
            &bold,
        );
        y -= 8.0;
        layer.use_text(
            "event_id            gene_id        kind   ΔPSI    adj.p     conf.",
            BODY_SIZE,
            Mm(MARGIN_MM),
            Mm(y),
            &bold,
        );
        y -= LINE_HEIGHT_MM;
        for r in &payload.high_confidence {
            let line = format!(
                "{:<20} {:<14} {:<5}  {:>+6.3}  {:<9} {:>5.3}",
                truncate(&r.event_id, 20),
                truncate(&r.gene_id, 14),
                r.event_type,
                r.delta_psi,
                fmt_p(r.adjusted_p_value),
                r.confidence,
            );
            layer.use_text(line, BODY_SIZE, Mm(MARGIN_MM), Mm(y), &mono);
            y -= LINE_HEIGHT_MM;
            if y < MARGIN_MM {
                // Paginate.
                let (np, nl) = doc.add_page(Mm(PAGE_W_MM), Mm(PAGE_H_MM), "high-conf-more");
                let _ = layer;
                let _ = nl;
                y = PAGE_H_MM - MARGIN_MM;
                let _ = np;
            }
        }
    }

    // ------------------------ Consensus summary ----------------------------
    if !payload.consensus.is_empty() {
        let (page, layer_id) = doc.add_page(Mm(PAGE_W_MM), Mm(PAGE_H_MM), "consensus");
        let layer = doc.get_page(page).get_layer(layer_id);
        let mut y = PAGE_H_MM - MARGIN_MM;
        layer.use_text(
            "Consensus statistical engine",
            HEADING_SIZE + 2.0,
            Mm(MARGIN_MM),
            Mm(y),
            &bold,
        );
        y -= 8.0;
        layer.use_text(
            "rank  combined_p   consensus_q   confidence   methods",
            BODY_SIZE,
            Mm(MARGIN_MM),
            Mm(y),
            &bold,
        );
        y -= LINE_HEIGHT_MM;
        let mut ranked: Vec<_> = payload
            .consensus
            .iter()
            .enumerate()
            .filter(|(_, r)| r.combined_p.is_finite())
            .collect();
        ranked.sort_by(|a, b| {
            a.1.consensus_q
                .partial_cmp(&b.1.consensus_q)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        for (rank, (_, r)) in ranked.iter().take(40).enumerate() {
            let line = format!(
                "{:>4}  {:<10}   {:<10}    {:>5.3}        {}",
                rank + 1,
                fmt_p(r.combined_p),
                fmt_p(r.consensus_q),
                r.confidence,
                r.methods_combined,
            );
            layer.use_text(line, BODY_SIZE, Mm(MARGIN_MM), Mm(y), &mono);
            y -= LINE_HEIGHT_MM;
            if y < MARGIN_MM {
                break;
            }
        }
    }

    // ----------------- Protein consequence breakdown -----------------------
    if !payload.protein_consequences.is_empty() {
        let (page, layer_id) = doc.add_page(Mm(PAGE_W_MM), Mm(PAGE_H_MM), "protein");
        let layer = doc.get_page(page).get_layer(layer_id);
        let mut y = PAGE_H_MM - MARGIN_MM;
        layer.use_text(
            "Protein consequence + NMD prediction",
            HEADING_SIZE + 2.0,
            Mm(MARGIN_MM),
            Mm(y),
            &bold,
        );
        y -= 8.0;
        let mut counts: std::collections::BTreeMap<String, usize> = Default::default();
        for p in &payload.protein_consequences {
            *counts.entry(p.consequence.clone()).or_insert(0) += 1;
        }
        for (k, v) in counts {
            layer.use_text(
                format!("  {k:<14} {v:>6}"),
                BODY_SIZE,
                Mm(MARGIN_MM),
                Mm(y),
                &mono,
            );
            y -= LINE_HEIGHT_MM;
        }
    }

    save(doc, out)
}

pub fn render_from_json(input: &Path, out: &Path) -> UltiResult<()> {
    let raw = std::fs::read_to_string(input).map_err(|e| UltiError::io(input, e))?;
    let payload: ReportPayload =
        serde_json::from_str(&raw).map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
    render(&payload, out)
}

fn save(doc: PdfDocumentReference, out: &Path) -> UltiResult<()> {
    use std::io::BufWriter;
    let f = std::fs::File::create(out).map_err(|e| UltiError::io(out, e))?;
    let mut w = BufWriter::new(f);
    doc.save(&mut w)
        .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
    Ok(())
}

fn truncate(s: &str, n: usize) -> String {
    if s.len() <= n {
        s.to_string()
    } else {
        let mut s = s.to_string();
        s.truncate(n.saturating_sub(1));
        s.push('…');
        s
    }
}

fn fmt_p(p: f64) -> String {
    if !p.is_finite() {
        return "NA".into();
    }
    if p.abs() < 1e-3 {
        format!("{:.2e}", p)
    } else {
        format!("{:.4}", p)
    }
}
