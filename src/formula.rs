//! Wilkinson-style design formula parser.
//!
//! Supports:
//! * `~ a` — single main effect
//! * `~ a + b` — additive main effects
//! * `~ a + b + a:b` — main effects with explicit interaction
//! * `~ a * b` — sugar for `a + b + a:b`
//! * `~ 0 + a` — no intercept

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::config::Sample;
use crate::error::{UltiError, UltiResult};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Formula {
    pub intercept: bool,
    pub terms: Vec<FormulaTerm>,
    raw: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum FormulaTerm {
    Main(String),
    Interaction(String, String),
}

impl Formula {
    pub fn parse(s: &str) -> UltiResult<Self> {
        let raw = s.trim().to_string();
        let body = raw
            .strip_prefix('~')
            .ok_or_else(|| UltiError::Formula(format!("formula must start with `~`: {raw}")))?
            .trim();
        let mut intercept = true;
        let mut terms: Vec<FormulaTerm> = Vec::new();

        for chunk in body.split('+') {
            let token = chunk.trim();
            if token.is_empty() {
                return Err(UltiError::Formula("empty additive term".into()));
            }
            if token == "0" {
                intercept = false;
                continue;
            }
            if token == "1" {
                intercept = true;
                continue;
            }
            if let Some((l, r)) = token.split_once('*') {
                let (l, r) = (l.trim().to_string(), r.trim().to_string());
                push_unique(&mut terms, FormulaTerm::Main(l.clone()));
                push_unique(&mut terms, FormulaTerm::Main(r.clone()));
                push_unique(&mut terms, FormulaTerm::Interaction(l, r));
                continue;
            }
            if let Some((l, r)) = token.split_once(':') {
                push_unique(
                    &mut terms,
                    FormulaTerm::Interaction(l.trim().into(), r.trim().into()),
                );
                continue;
            }
            push_unique(&mut terms, FormulaTerm::Main(token.into()));
        }
        Ok(Formula {
            intercept,
            terms,
            raw,
        })
    }

    pub fn raw(&self) -> &str {
        &self.raw
    }

    pub fn variables(&self) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        for t in &self.terms {
            match t {
                FormulaTerm::Main(v) => {
                    if !out.contains(v) {
                        out.push(v.clone());
                    }
                }
                FormulaTerm::Interaction(a, b) => {
                    for v in [a, b] {
                        if !out.contains(v) {
                            out.push(v.clone());
                        }
                    }
                }
            }
        }
        out
    }

    /// Build the n×p design matrix from a slice of samples. Categorical
    /// variables are treatment-coded (first level as reference). Numeric
    /// variables pass through as-is.
    pub fn design_matrix(&self, samples: &[Sample]) -> UltiResult<(Vec<String>, Vec<Vec<f64>>)> {
        let vars = self.variables();
        let mut levels: BTreeMap<String, Vec<String>> = BTreeMap::new();
        let mut numeric: BTreeMap<String, Vec<f64>> = BTreeMap::new();

        for v in &vars {
            let mut all_str: Vec<String> = Vec::new();
            let mut all_num: Vec<f64> = Vec::new();
            let mut is_num = true;
            for s in samples {
                let value = sample_value(s, v);
                if let Ok(n) = value.parse::<f64>() {
                    all_num.push(n);
                } else {
                    is_num = false;
                }
                all_str.push(value);
            }
            if is_num && all_num.len() == samples.len() {
                numeric.insert(v.clone(), all_num);
            } else {
                let mut uniq: Vec<String> = Vec::new();
                for x in &all_str {
                    if !uniq.contains(x) {
                        uniq.push(x.clone());
                    }
                }
                levels.insert(v.clone(), uniq);
            }
        }

        let mut col_names: Vec<String> = Vec::new();
        let mut col_values: Vec<Vec<f64>> = Vec::new();
        if self.intercept {
            col_names.push("(Intercept)".into());
            col_values.push(vec![1.0; samples.len()]);
        }

        for term in &self.terms {
            match term {
                FormulaTerm::Main(v) => {
                    if let Some(nums) = numeric.get(v) {
                        col_names.push(v.clone());
                        col_values.push(nums.clone());
                    } else if let Some(lv) = levels.get(v) {
                        for level in lv.iter().skip(1) {
                            col_names.push(format!("{v}={level}"));
                            let vals: Vec<f64> = samples
                                .iter()
                                .map(|s| {
                                    if sample_value(s, v) == *level {
                                        1.0
                                    } else {
                                        0.0
                                    }
                                })
                                .collect();
                            col_values.push(vals);
                        }
                    }
                }
                FormulaTerm::Interaction(a, b) => {
                    let a_cols = expand_term(a, &numeric, &levels, samples);
                    let b_cols = expand_term(b, &numeric, &levels, samples);
                    for (na, va) in &a_cols {
                        for (nb, vb) in &b_cols {
                            col_names.push(format!("{na}:{nb}"));
                            col_values.push(va.iter().zip(vb.iter()).map(|(x, y)| x * y).collect());
                        }
                    }
                }
            }
        }

        let n = samples.len();
        let p = col_names.len();
        if p == 0 {
            return Err(UltiError::Formula(
                "empty design — formula produced no columns".into(),
            ));
        }
        let mut rows: Vec<Vec<f64>> = vec![vec![0.0; p]; n];
        for (j, col) in col_values.iter().enumerate() {
            for (i, v) in col.iter().enumerate() {
                rows[i][j] = *v;
            }
        }
        Ok((col_names, rows))
    }

    /// Contrast vector that selects the named column (or all dummies
    /// derived from a categorical of that name).
    pub fn contrast_for(col_names: &[String], term: &str) -> Vec<f64> {
        let mut c = vec![0.0; col_names.len()];
        for (i, name) in col_names.iter().enumerate() {
            if name == term || name.starts_with(&format!("{term}=")) {
                c[i] = 1.0;
            }
        }
        c
    }
}

fn push_unique(terms: &mut Vec<FormulaTerm>, term: FormulaTerm) {
    let dup = terms.iter().any(|t| match (t, &term) {
        (FormulaTerm::Main(a), FormulaTerm::Main(b)) => a == b,
        (FormulaTerm::Interaction(a, b), FormulaTerm::Interaction(c, d)) => {
            (a == c && b == d) || (a == d && b == c)
        }
        _ => false,
    });
    if !dup {
        terms.push(term);
    }
}

fn sample_value(s: &Sample, v: &str) -> String {
    if v.eq_ignore_ascii_case("group") {
        s.group.clone()
    } else if v.eq_ignore_ascii_case("sample") {
        s.id.clone()
    } else {
        s.covariates.get(v).cloned().unwrap_or_default()
    }
}

fn expand_term(
    name: &str,
    numeric: &BTreeMap<String, Vec<f64>>,
    levels: &BTreeMap<String, Vec<String>>,
    samples: &[Sample],
) -> Vec<(String, Vec<f64>)> {
    if let Some(v) = numeric.get(name) {
        return vec![(name.to_string(), v.clone())];
    }
    if let Some(lv) = levels.get(name) {
        return lv
            .iter()
            .skip(1)
            .map(|level| {
                let col: Vec<f64> = samples
                    .iter()
                    .map(|s| {
                        if sample_value(s, name) == *level {
                            1.0
                        } else {
                            0.0
                        }
                    })
                    .collect();
                (format!("{name}={level}"), col)
            })
            .collect();
    }
    Vec::new()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;
    use std::path::PathBuf;

    fn s(id: &str, group: &str, batch: &str) -> Sample {
        let mut cov = BTreeMap::new();
        cov.insert("batch".to_string(), batch.to_string());
        Sample {
            id: id.into(),
            bam: PathBuf::from(format!("{id}.bam")),
            group: group.into(),
            covariates: cov,
        }
    }

    #[test]
    fn parses_main_and_interaction() {
        let f = Formula::parse("~ group + batch + group:batch").unwrap();
        assert!(f.intercept);
        assert_eq!(f.terms.len(), 3);
    }

    #[test]
    fn star_sugar_expands() {
        let f = Formula::parse("~ group * batch").unwrap();
        assert_eq!(f.terms.len(), 3);
    }

    #[test]
    fn no_intercept() {
        let f = Formula::parse("~ 0 + group").unwrap();
        assert!(!f.intercept);
    }

    #[test]
    fn design_matrix_for_2x2() {
        let samples = vec![
            s("s1", "ctrl", "b1"),
            s("s2", "ctrl", "b2"),
            s("s3", "treat", "b1"),
            s("s4", "treat", "b2"),
        ];
        let f = Formula::parse("~ group + batch").unwrap();
        let (names, rows) = f.design_matrix(&samples).unwrap();
        assert_eq!(names, vec!["(Intercept)", "group=treat", "batch=b2"]);
        assert_eq!(rows[0], vec![1.0, 0.0, 0.0]);
        assert_eq!(rows[3], vec![1.0, 1.0, 1.0]);
    }

    #[test]
    fn interaction_columns() {
        let samples = vec![
            s("s1", "ctrl", "b1"),
            s("s2", "treat", "b1"),
            s("s3", "ctrl", "b2"),
            s("s4", "treat", "b2"),
        ];
        let f = Formula::parse("~ group * batch").unwrap();
        let (names, _) = f.design_matrix(&samples).unwrap();
        assert!(names.iter().any(|n| n.contains(':')));
    }
}
