//! Crate-wide error type with **stable, machine-readable error codes**.
//!
//! Every `UltiError` variant maps to a short `E####` code returned by
//! [`UltiError::code`]. Downstream tooling (CI scripts, wrappers, GUI
//! frontends) should pattern-match on the code rather than the human
//! message — codes are part of the crate's public API contract and
//! follow semver.

use std::path::PathBuf;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum UltiError {
    #[error("[{code}] I/O error at {path:?}: {source}", code = self.code())]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("[{code}] annotation parse error in {path:?} at line {line}: {message}", code = self.code())]
    Annotation {
        path: PathBuf,
        line: usize,
        message: String,
    },

    #[error("[{code}] BAM/CRAM parse error in {path:?}: {message}", code = self.code())]
    Alignment { path: PathBuf, message: String },

    #[error("[{code}] reference FASTA error in {path:?}: {message}", code = self.code())]
    Reference { path: PathBuf, message: String },

    #[error("[{code}] configuration error: {0}", code = self.code())]
    Config(String),

    #[error("[{code}] design matrix error: {0}", code = self.code())]
    Design(String),

    #[error("[{code}] design formula parse error: {0}", code = self.code())]
    Formula(String),

    #[error("[{code}] statistical computation failed: {0}", code = self.code())]
    Stats(String),

    #[error("[{code}] cache error: {0}", code = self.code())]
    Cache(String),

    #[error("[{code}] unsupported feature: {0}", code = self.code())]
    Unsupported(String),

    #[error("[{code}] {0}", code = self.code())]
    Other(#[from] anyhow::Error),
}

pub type UltiResult<T> = Result<T, UltiError>;

impl UltiError {
    /// Stable error code. **Part of the public API** — do not change codes
    /// across minor versions; add new ones for new variants.
    pub fn code(&self) -> &'static str {
        match self {
            UltiError::Io { .. } => "E0001",
            UltiError::Annotation { .. } => "E0010",
            UltiError::Alignment { .. } => "E0020",
            UltiError::Reference { .. } => "E0030",
            UltiError::Config(_) => "E0040",
            UltiError::Design(_) => "E0050",
            UltiError::Formula(_) => "E0051",
            UltiError::Stats(_) => "E0060",
            UltiError::Cache(_) => "E0070",
            UltiError::Unsupported(_) => "E0090",
            UltiError::Other(_) => "E0099",
        }
    }

    pub fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        UltiError::Io {
            path: path.into(),
            source,
        }
    }

    pub fn reference(path: impl Into<PathBuf>, message: impl Into<String>) -> Self {
        UltiError::Reference {
            path: path.into(),
            message: message.into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_are_stable() {
        // Pin known codes — failing this test means a backwards-incompatible
        // change to the public error-code surface.
        assert_eq!(UltiError::Config("x".into()).code(), "E0040");
        assert_eq!(UltiError::Stats("x".into()).code(), "E0060");
        assert_eq!(UltiError::Formula("x".into()).code(), "E0051");
        assert_eq!(UltiError::Cache("x".into()).code(), "E0070");
    }

    #[test]
    fn error_message_includes_code() {
        let e = UltiError::Config("missing thing".into());
        let msg = format!("{e}");
        assert!(msg.contains("E0040"), "got: {msg}");
        assert!(msg.contains("missing thing"));
    }
}
