//! HTML reporter. Produces a single self-contained `.html` file with the
//! results JSON embedded inline. Plotly.js is loaded from CDN so report
//! files stay small (~50 KB) and inherit a maintained, modern plotting
//! library.
//!
//! Layout:
//! 1. Run metadata + executive summary
//! 2. Volcano plot (−log10 q-value vs ΔPSI)
//! 3. Per-event-type PSI distribution histograms
//! 4. Top-20 significant events tables (one per AS kind)
//!
//! Currently *not* rendered yet (clear TODOs in the template):
//! * Sashimi plots — require per-event read pileup (deferred)
//! * MA plot — needs library-size normalization
//! * PCA / UMAP — needs an event-by-sample PSI matrix; small wrinkle on top
//!   of `EventQuant.psi` and a PCA implementation
//! * Sankey isoform plot — needs long-read isoform reconstruction

use std::fs;
use std::path::Path;

use crate::error::{UltiError, UltiResult};
use crate::output::ReportPayload;

const TEMPLATE: &str = include_str!("report_template.html");

pub fn render(payload: &ReportPayload, out: &Path) -> UltiResult<()> {
    let json = serde_json::to_string(payload)
        .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
    let html = TEMPLATE.replace("/*__ULTIMADSE_DATA__*/null", &json);
    fs::write(out, html).map_err(|e| UltiError::io(out, e))
}

pub fn render_from_json(input: &Path, out: &Path) -> UltiResult<()> {
    let raw = fs::read_to_string(input).map_err(|e| UltiError::io(input, e))?;
    let payload: ReportPayload =
        serde_json::from_str(&raw).map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
    render(&payload, out)
}
