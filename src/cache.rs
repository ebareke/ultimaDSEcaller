//! Annotation binary cache.
//!
//! Parsing Gencode v45 from a 1.4 GB GTF takes ~30 s. The pipeline runs
//! the parse on every invocation, even when the input hasn't changed.
//! This module persists a binary snapshot (`.ultidse` via `bincode`) and
//! loads it when fresher than the source.
//!
//! Cache contract:
//! * **Versioned**: `CACHE_VERSION` bumped on every schema change. Older
//!   caches are rejected (callers fall back to a full parse).
//! * **mtime-checked**: cache must be ≥ source mtime, else discarded.
//! * **Atomic write**: cache is written to a temp file and renamed.
//!
//! Only the *raw structured* annotation is cached — splice graphs and
//! Lapper indexes are rebuilt on load because they contain non-portable
//! indexes (`NodeIndex` values change between rebuilds).

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use serde::{Deserialize, Serialize};

use crate::Strand;
use crate::annotation::{Annotation, Exon};
use crate::error::{UltiError, UltiResult};

const CACHE_VERSION: u32 = 1;
const CACHE_MAGIC: &[u8; 8] = b"ULTDSE01";

#[derive(Debug, Serialize, Deserialize)]
pub struct CachedAnnotation {
    pub version: u32,
    pub genes: Vec<CachedGene>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CachedGene {
    pub gene_id: String,
    pub gene_name: Option<String>,
    pub chrom: String,
    pub strand: Strand,
    pub transcripts: BTreeMap<String, Vec<Exon>>,
}

/// Convert a live [`Annotation`] into its cacheable form.
pub fn to_cached(ann: &Annotation) -> CachedAnnotation {
    let mut genes = Vec::with_capacity(ann.genes.len());
    for g in ann.genes.values() {
        genes.push(CachedGene {
            gene_id: g.gene_id.clone(),
            gene_name: g.gene_name.clone(),
            chrom: g.chrom.clone(),
            strand: g.strand,
            transcripts: g.transcripts.clone(),
        });
    }
    CachedAnnotation {
        version: CACHE_VERSION,
        genes,
    }
}

/// Rehydrate a [`CachedAnnotation`] into a full [`Annotation`] (rebuilds
/// splice graphs and spatial indexes).
pub fn from_cached(cached: CachedAnnotation) -> UltiResult<Annotation> {
    if cached.version != CACHE_VERSION {
        return Err(UltiError::Cache(format!(
            "cache version {} != expected {}",
            cached.version, CACHE_VERSION
        )));
    }
    Ok(crate::annotation::build_from_cached(cached.genes))
}

/// Default cache path derived from the source: `<source>.ultidse`.
pub fn default_cache_path(source: &Path) -> PathBuf {
    let mut p = source.to_path_buf();
    let mut name = p
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("annotation")
        .to_string();
    name.push_str(".ultidse");
    p.set_file_name(name);
    p
}

/// True iff the cache exists and is at least as new as the source.
pub fn is_fresh(cache: &Path, source: &Path) -> bool {
    let cache_mtime = match File::open(cache)
        .and_then(|f| f.metadata())
        .and_then(|m| m.modified())
    {
        Ok(t) => t,
        Err(_) => return false,
    };
    let src_mtime = match File::open(source)
        .and_then(|f| f.metadata())
        .and_then(|m| m.modified())
    {
        Ok(t) => t,
        Err(_) => return false,
    };
    cache_mtime >= src_mtime || cache_mtime == SystemTime::UNIX_EPOCH
}

/// Load a cached annotation. Returns an error if the file is missing,
/// truncated, has a bad magic header, or has an incompatible version.
pub fn load(path: &Path) -> UltiResult<Annotation> {
    let mut file = File::open(path).map_err(|e| UltiError::io(path, e))?;
    let mut magic = [0u8; 8];
    use std::io::Read;
    file.read_exact(&mut magic)
        .map_err(|e| UltiError::Cache(format!("cannot read magic: {e}")))?;
    if &magic != CACHE_MAGIC {
        return Err(UltiError::Cache("bad magic header".into()));
    }
    let reader = BufReader::new(file);
    let cached: CachedAnnotation =
        bincode::deserialize_from(reader).map_err(|e| UltiError::Cache(e.to_string()))?;
    from_cached(cached)
}

/// Atomically write a binary cache for `ann` to `path`.
pub fn save(path: &Path, ann: &Annotation) -> UltiResult<()> {
    let tmp = {
        let mut t = path.to_path_buf();
        let mut name = t
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("cache")
            .to_string();
        name.push_str(".tmp");
        t.set_file_name(name);
        t
    };
    {
        let f = File::create(&tmp).map_err(|e| UltiError::io(&tmp, e))?;
        let mut w = BufWriter::new(f);
        use std::io::Write;
        w.write_all(CACHE_MAGIC)
            .map_err(|e| UltiError::io(&tmp, e))?;
        let cached = to_cached(ann);
        bincode::serialize_into(&mut w, &cached).map_err(|e| UltiError::Cache(e.to_string()))?;
        w.flush().map_err(|e| UltiError::io(&tmp, e))?;
    }
    std::fs::rename(&tmp, path).map_err(|e| UltiError::io(path, e))?;
    Ok(())
}

/// High-level helper: parse + cache the annotation. Loads from cache if
/// it's fresh, otherwise parses the source and writes a new cache file.
pub fn parse_cached(source: &Path, cache_path: Option<&Path>) -> UltiResult<Annotation> {
    let default_path;
    let cache = match cache_path {
        Some(p) => p,
        None => {
            default_path = default_cache_path(source);
            default_path.as_path()
        }
    };
    if is_fresh(cache, source) {
        match load(cache) {
            Ok(a) => {
                tracing::info!("loaded annotation from cache {:?}", cache);
                return Ok(a);
            }
            Err(e) => {
                tracing::warn!("cache load failed ({e}); falling back to full parse");
            }
        }
    }
    let ann = crate::annotation::parse(source)?;
    if let Err(e) = save(cache, &ann) {
        tracing::warn!("could not write cache {:?}: {e}", cache);
    } else {
        tracing::info!("wrote annotation cache to {:?}", cache);
    }
    Ok(ann)
}
