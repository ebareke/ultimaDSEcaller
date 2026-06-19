//! Thin progress-bar facade over `indicatif`.
//!
//! The facade serves two purposes:
//! 1. Centralize the bar style so the whole pipeline looks consistent.
//! 2. Make it trivial to silence progress in CI or tests (`UltiASCALLER_NO_PROGRESS=1`)
//!    without touching call sites.
//!
//! Use [`bar`] for a counted bar, [`spinner`] for a length-unknown task,
//! and [`MultiBar`] to display several bars at once.

use std::io::IsTerminal;
use std::sync::Arc;
use std::time::Duration;

use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

const TICK_INTERVAL_MS: u64 = 100;

/// Returns true when progress bars are suppressed (env, non-TTY, or CI).
pub fn suppressed() -> bool {
    if std::env::var("ULTIMADSE_NO_PROGRESS")
        .map(|v| v == "1")
        .unwrap_or(false)
    {
        return true;
    }
    !std::io::stderr().is_terminal()
}

/// Build a deterministic-length progress bar.
pub fn bar(label: &str, total: u64) -> ProgressBar {
    if suppressed() {
        return ProgressBar::hidden();
    }
    let pb = ProgressBar::new(total);
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} {prefix:>14.bold.dim} [{bar:30.cyan/blue}] {pos}/{len} ({elapsed_precise}) {msg}",
        )
        .expect("static template")
        .progress_chars("=>-"),
    );
    pb.set_prefix(label.to_string());
    pb.enable_steady_tick(Duration::from_millis(TICK_INTERVAL_MS));
    pb
}

/// Build a spinner for an unknown-length task.
pub fn spinner(label: &str) -> ProgressBar {
    if suppressed() {
        return ProgressBar::hidden();
    }
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} {prefix:>14.bold.dim} {msg} ({elapsed_precise})",
        )
        .expect("static template"),
    );
    pb.set_prefix(label.to_string());
    pb.enable_steady_tick(Duration::from_millis(TICK_INTERVAL_MS));
    pb
}

/// Coordinator for multiple concurrent bars.
#[derive(Clone)]
pub struct MultiBar {
    inner: Arc<MultiProgress>,
}

impl MultiBar {
    pub fn new() -> Self {
        let mp = MultiProgress::new();
        if suppressed() {
            mp.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }
        MultiBar {
            inner: Arc::new(mp),
        }
    }

    pub fn add_bar(&self, label: &str, total: u64) -> ProgressBar {
        let pb = bar(label, total);
        self.inner.add(pb)
    }

    pub fn add_spinner(&self, label: &str) -> ProgressBar {
        let pb = spinner(label);
        self.inner.add(pb)
    }
}

impl Default for MultiBar {
    fn default() -> Self {
        Self::new()
    }
}
