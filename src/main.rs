//! ultimaDSEcaller binary — thin CLI shell around the `ultimadsecaller` library.

use anyhow::{Context, Result};
use clap::Parser;
use tracing_subscriber::EnvFilter;

use ultimadsecaller::cli::{Cli, Command};
use ultimadsecaller::{
    advanced, annotation, cache, config, consensus, embedding, events, formula, glm, junctions,
    longread, motif, output, progress, protein, quantify, report, sashimi, stats,
};

fn main() -> Result<()> {
    let cli = Cli::parse();
    init_logging(cli.verbose);

    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .ok();
    }

    match cli.command {
        Command::Run(args) => run(args),
        Command::DumpAnnotation(args) => dump_annotation(args),
        Command::Junctions(args) => junctions_cmd(args),
        Command::Report(args) => {
            report::render_from_json(&args.input, &args.out)
                .context("rendering HTML report")?;
            tracing::info!("wrote HTML report to {:?}", args.out);
            Ok(())
        }
        Command::Pdf(args) => {
            #[cfg(feature = "pdf")]
            {
                ultimadsecaller::pdf::render_from_json(&args.input, &args.out)
                    .context("rendering PDF report")?;
                tracing::info!("wrote PDF report to {:?}", args.out);
                Ok(())
            }
            #[cfg(not(feature = "pdf"))]
            {
                let _ = args;
                anyhow::bail!(
                    "the binary was not built with --features pdf — rebuild with \
                     `cargo build --release --features pdf` to enable the `pdf` subcommand"
                );
            }
        }
    }
}

fn init_logging(verbosity: u8) {
    let level = match verbosity {
        0 => "ultimaDSEcaller=info,ultimadsecaller=info",
        1 => "ultimaDSEcaller=debug,ultimadsecaller=debug",
        _ => "trace",
    };
    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(level));
    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_target(false)
        .with_thread_names(false)
        .compact()
        .init();
}

fn run(args: ultimadsecaller::cli::RunArgs) -> Result<()> {
    let cfg = config::resolve(&args).context("resolving run configuration")?;
    let resume = args.resume;
    let no_cache = args.no_cache;
    let cache_path = args.cache_path.clone();
    let formula_spec = args.formula.clone();
    let contrasts_spec: Vec<String> = args.contrasts.clone();
    let consensus_method = args
        .consensus
        .map(|c| c.into())
        .unwrap_or(consensus::ConsensusMethod::Stouffer);

    tracing::info!(
        "ultimaDSEcaller starting — {} samples, tech={:?}, contrast={:?}",
        cfg.samples.len(),
        cfg.tech,
        cfg.contrast.as_ref().map(|c| format!("{}:{}-{}", c.variable, c.numerator, c.denominator)),
    );

    // --- Annotation (cached by default) ---
    let ann_sp = progress::spinner("annotation");
    ann_sp.set_message(format!("parsing {:?}", cfg.annotation));
    let ann = if no_cache {
        annotation::parse(&cfg.annotation).context("annotation parse failed")?
    } else {
        cache::parse_cached(&cfg.annotation, cache_path.as_deref())
            .context("annotation parse/cache failed")?
    };
    ann_sp.finish_with_message(format!("loaded {} genes", ann.gene_count()));

    // --- Reference FASTA (optional) ---
    let reference: Option<motif::Reference> = if let Some(p) = cfg.reference.as_deref() {
        let sp = progress::spinner("reference");
        sp.set_message(format!("loading {:?}", p));
        let r = motif::Reference::load(p).context("loading reference FASTA")?;
        sp.finish_with_message(format!("indexed {} sequences", r.chromosomes().count()));
        Some(r)
    } else {
        None
    };

    // --- CDS catalog (optional — needed for protein consequence) ---
    let cds_catalog = if reference.is_some() {
        match protein::parse_cds(&cfg.annotation) {
            Ok(c) => {
                tracing::info!("parsed CDS records for {} transcripts", c.transcripts.len());
                Some(c)
            }
            Err(e) => {
                tracing::warn!("CDS parse failed ({e}); protein consequence skipped");
                None
            }
        }
    } else {
        None
    };

    // --- Junctions (checkpointed if --resume) ---
    let junctions_cache = cfg.out.join("junctions.bin");
    std::fs::create_dir_all(&cfg.out)
        .with_context(|| format!("cannot create output dir {:?}", cfg.out))?;
    let jm_sp = progress::spinner("junctions");
    let jm = if resume && junctions_cache.exists() {
        jm_sp.set_message(format!("loading checkpoint {:?}", junctions_cache));
        let jm = junctions::JunctionMatrix::load(&junctions_cache)
            .context("loading junction checkpoint")?;
        jm_sp.finish_with_message(format!("resumed {} junctions", jm.counts.len()));
        jm
    } else {
        jm_sp.set_message(format!("extracting from {} BAMs", cfg.samples.len()));
        let jm = junctions::extract(&cfg).context("junction extraction failed")?;
        if let Err(e) = jm.save(&junctions_cache) {
            tracing::warn!("could not write junction checkpoint: {e}");
        }
        jm_sp.finish_with_message(format!("found {} unique junctions", jm.counts.len()));
        jm
    };

    // --- Events + quantification ---
    let evs_sp = progress::spinner("events");
    evs_sp.set_message("detecting AS events");
    let evs = events::detect_all(&ann, &jm);
    evs_sp.finish_with_message(format!("detected {} candidate events", evs.len()));

    let quant_sp = progress::spinner("quantify");
    quant_sp.set_message("computing PSI / ΔPSI");
    let quants = quantify::quantify(&cfg, &evs, &jm);
    quant_sp.finish_with_message(format!("quantified {} events", quants.len()));

    // --- Stats: default test, optional formula-driven GLM/GLMM ---
    let stats_sp = progress::spinner("stats");
    stats_sp.set_message(format!("running {:?}", cfg.test));
    let pvals = if let Some(spec) = formula_spec.as_deref() {
        run_formula_glm(&cfg, &quants, spec)?
    } else {
        stats::test_all_with_method(&quants, Some(&cfg))
    };
    stats_sp.finish_with_message(format!("tested {} events", pvals.len()));

    let rows = output::build_rows(&quants, &pvals);

    // --- Consensus engine (combines BB-LRT/GLM/GLMM/Fisher/DIU + motif) ---
    let _ = &pvals; // kept for the headline rows; consensus runs its own multi-test sweep.
    let consensus_results = build_consensus(
        &cfg,
        &quants,
        reference.as_ref(),
        consensus_method,
    );
    tracing::info!(
        "consensus engine combined {} events via {:?}",
        consensus_results.len(),
        consensus_method
    );

    // --- Protein consequence (if reference + CDS available) ---
    let protein_consequences: Vec<output::ProteinAnnotation> = if let (Some(cat), refq) =
        (cds_catalog.as_ref(), reference.as_ref())
    {
        evs.iter()
            .map(|e| {
                let c = protein::predict_consequence(e, cat, refq);
                output::ProteinAnnotation {
                    event_id: e.event_id.clone(),
                    gene_id: e.gene_id.clone(),
                    consequence: c.short().into(),
                }
            })
            .collect()
    } else {
        Vec::new()
    };

    // Advanced events (cryptic, MSE, MIR, recursive, nested, partial exon,
    // exonic-intronic hybrid, alt promoter, alt polyA, tandem UTR, fusion).
    let adv_params = advanced::AdvancedParams::default();
    let adv_events = advanced::detect_all(&ann, &jm, &adv_params, cfg.fusion_bedpe.as_deref())
        .context("advanced event detection")?;
    let mut advanced_event_counts: std::collections::BTreeMap<String, usize> =
        std::collections::BTreeMap::new();
    for e in &adv_events {
        *advanced_event_counts
            .entry(e.kind.short().to_string())
            .or_insert(0) += 1;
    }
    tracing::info!(
        "detected {} advanced events across {} categories",
        adv_events.len(),
        advanced_event_counts.len(),
    );

    // Per-sample PSI matrix for embeddings.
    let psi_rows: Vec<Vec<f64>> = quants.iter().map(|q| q.psi.clone()).collect();
    let sample_ids: Vec<String> = cfg.samples.iter().map(|s| s.id.clone()).collect();
    let (pca_emb, umap_emb) = if let Some(m) = embedding::build_psi_matrix(&sample_ids, &psi_rows)
    {
        let k = (sample_ids.len() / 3).max(2);
        (
            Some(embedding::pca_2d(&m, &sample_ids)),
            Some(embedding::umap_like_2d(&m, &sample_ids, k)),
        )
    } else {
        (None, None)
    };

    // Sashimi tracks for top events. One pileup pass per sample.
    let sashimi_tracks =
        sashimi::build_top_event_sashimi(&cfg, &evs, &rows, 3).unwrap_or_default();

    // Long-read isoforms + differential usage Sankey.
    let isoform_catalog = longread::reconstruct(&cfg, &ann).unwrap_or_else(|_| longread::IsoformCatalog::empty());
    let sample_groups: Vec<&str> = cfg.samples.iter().map(|s| s.group.as_str()).collect();
    let (sankey, diu_records) = if let Some(c) = cfg.contrast.as_ref() {
        let s = sashimi::build_isoform_sankey(
            &isoform_catalog,
            &sample_groups,
            &c.numerator,
            &c.denominator,
            12,
        );
        let d = longread::differential_usage(
            &isoform_catalog,
            &sample_groups,
            &c.numerator,
            &c.denominator,
        );
        (Some(s), d)
    } else {
        (None, Vec::new())
    };

    // --- New visualizations: heatmap + junction graphs ---
    let heatmap = output::build_heatmap(&quants, &rows, &sample_ids, 50);
    let junction_graphs = output::build_junction_graphs(&ann, &rows, 8);

    // Coverage distribution (per-sample summary).
    let coverage_distribution: Vec<output::CoverageRow> = cfg
        .samples
        .iter()
        .map(|s| {
            // Per-sample mean junction support (rough QC).
            let n = jm.samples.iter().position(|x| x == &s.id);
            let mean = if let Some(idx) = n {
                let mut total = 0.0_f64;
                let mut cnt = 0_u64;
                for v in jm.counts.values() {
                    if v[idx] > 0.0 {
                        total += v[idx];
                        cnt += 1;
                    }
                }
                if cnt > 0 { total / cnt as f64 } else { 0.0 }
            } else {
                0.0
            };
            output::CoverageRow {
                sample: s.id.clone(),
                group: s.group.clone(),
                total_passing_reads: 0,
                low_mapq_reads: 0,
                mean_junction_support: mean,
            }
        })
        .collect();

    let counts_clone = advanced_event_counts.clone();
    let consensus_clone = consensus_results.clone();
    let protein_clone = protein_consequences.clone();
    output::write_all_with_payload(&cfg, &rows, move |mut p| {
        p.pca = pca_emb;
        p.umap = umap_emb;
        p.sashimi = sashimi_tracks;
        p.isoform_sankey = sankey;
        p.advanced_event_counts = counts_clone;
        p.diu = diu_records;
        p.heatmap = heatmap;
        p.coverage_distribution = coverage_distribution;
        p.junction_graphs = junction_graphs;
        p.consensus = consensus_clone;
        p.protein_consequences = protein_clone;
        p
    })
    .context("writing output tables")?;

    // Per-contrast iteration (if multiple contrasts requested).
    if !contrasts_spec.is_empty() {
        run_extra_contrasts(&cfg, &evs, &jm, &contrasts_spec)?;
    }

    // Write the advanced-events table.
    let adv_path = cfg.out.join("advanced_events.tsv");
    write_advanced_events_table(&adv_path, &adv_events).context("writing advanced events")?;

    // Render the HTML report from the just-written results.json.
    let json_path = cfg.out.join("results.json");
    let html_path = cfg.out.join("report.html");
    report::render_from_json(&json_path, &html_path).context("rendering report")?;
    tracing::info!("done — report at {:?}", html_path);

    // Optional PDF rendering. Behind the `pdf` Cargo feature so a default
    // build stays lean.
    #[cfg(feature = "pdf")]
    if args.pdf {
        let pdf_path = cfg.out.join("report.pdf");
        ultimadsecaller::pdf::render_from_json(&json_path, &pdf_path)
            .context("rendering PDF report")?;
        tracing::info!("PDF report at {:?}", pdf_path);
    }
    #[cfg(not(feature = "pdf"))]
    if args.pdf {
        tracing::warn!("--pdf requested but binary was not built with --features pdf; skipping");
    }

    Ok(())
}

/// Build a consensus combination across BB-LRT/GLM/GLMM/Fisher/DIU per event.
///
/// Runs **all** applicable tests (via `stats::test_all_multi`) so that the
/// consensus engine combines independent statistical signals rather than
/// just whichever single test was nominated as primary.
fn build_consensus(
    cfg: &ultimadsecaller::config::RunConfig,
    quants: &[quantify::EventQuant],
    reference: Option<&motif::Reference>,
    method: consensus::ConsensusMethod,
) -> Vec<consensus::ConsensusResult> {
    let multi = stats::test_all_multi(quants, Some(cfg));
    let evidence: Vec<consensus::EventEvidence> = quants
        .iter()
        .zip(multi.iter())
        .map(|(q, m)| {
            let motif = reference.map(|r| {
                let donor = q
                    .event
                    .inclusion_junctions
                    .first()
                    .map(|(d, _)| *d)
                    .unwrap_or(0);
                let acceptor = q
                    .event
                    .inclusion_junctions
                    .first()
                    .map(|(_, a)| *a)
                    .unwrap_or(0);
                r.classify_junction(&q.event.chrom, donor, acceptor)
            });
            let mean_cov = q
                .contrast_summary
                .as_ref()
                .map(|s| 0.5 * (s.mean_coverage_num + s.mean_coverage_denom))
                .unwrap_or(0.0);

            consensus::EventEvidence {
                p_bb_lrt: m.p_bb_lrt,
                p_glm: m.p_glm,
                p_glmm: m.p_glmm,
                p_fisher: m.p_fisher,
                p_diu: None, // DIU is per-gene, joined separately downstream.
                motif,
                mean_coverage: mean_cov,
                replicate_reproducibility: q.reproducibility,
            }
        })
        .collect();
    consensus::combine(&evidence, &consensus::Weights::default(), method)
}

/// Formula-driven GLM test path: parses `--formula` and runs IRLS on the
/// design matrix it produces, replacing the default `[1, treatment]` design.
fn run_formula_glm(
    cfg: &ultimadsecaller::config::RunConfig,
    quants: &[quantify::EventQuant],
    spec: &str,
) -> Result<Vec<stats::PValue>> {
    let f = formula::Formula::parse(spec).context("parsing --formula")?;
    let (col_names, x) = f
        .design_matrix(&cfg.samples)
        .context("design matrix from --formula")?;
    let primary_term = cfg
        .contrast
        .as_ref()
        .map(|c| c.variable.clone())
        .unwrap_or_else(|| "group".into());
    let contrast = formula::Formula::contrast_for(&col_names, &primary_term);
    tracing::info!(
        "formula `{}` → {} cols: {:?}",
        spec,
        col_names.len(),
        col_names
    );

    let mut out = Vec::with_capacity(quants.len());
    let raw_p: Vec<f64> = quants
        .iter()
        .map(|q| {
            let n = cfg.samples.len();
            if q.inclusion.len() != n || q.exclusion.len() != n {
                return f64::NAN;
            }
            let y = q.inclusion.clone();
            let n_vec: Vec<f64> = q
                .inclusion
                .iter()
                .zip(q.exclusion.iter())
                .map(|(i, e)| i + e)
                .collect();
            match glm::fit_glm(&y, &n_vec, &x) {
                Ok(fit) => glm::wald_test(&fit.beta, &fit.vcov, &contrast).p_value,
                Err(_) => f64::NAN,
            }
        })
        .collect();
    let adj = stats::bh_fdr(&raw_p);
    for (i, p) in raw_p.iter().enumerate() {
        let eff = quants[i]
            .contrast_summary
            .as_ref()
            .map(|s| s.delta_psi)
            .unwrap_or(0.0);
        out.push(stats::PValue {
            p_value: *p,
            adjusted_p_value: adj[i],
            effect_size: eff,
            test_used: stats::TestUsed::Glm,
        });
    }
    Ok(out)
}

/// Run additional contrasts beyond the primary one and write each to its
/// own subdirectory.
fn run_extra_contrasts(
    cfg: &ultimadsecaller::config::RunConfig,
    evs: &[events::ASEvent],
    jm: &junctions::JunctionMatrix,
    contrasts: &[String],
) -> Result<()> {
    for spec in contrasts {
        let c = ultimadsecaller::config::Contrast::parse(spec)
            .map_err(anyhow::Error::msg)
            .with_context(|| format!("parsing contrast `{spec}`"))?;
        let mut sub_cfg = cfg.clone();
        sub_cfg.contrast = Some(c.clone());
        let subdir = cfg.out.join(format!(
            "contrast_{}_vs_{}",
            sanitize(&c.numerator),
            sanitize(&c.denominator),
        ));
        std::fs::create_dir_all(&subdir)
            .with_context(|| format!("creating subdir {:?}", subdir))?;
        sub_cfg.out = subdir.clone();

        let quants = quantify::quantify(&sub_cfg, evs, jm);
        let pvals = stats::test_all_with_method(&quants, Some(&sub_cfg));
        let rows = output::build_rows(&quants, &pvals);
        output::write_all(&sub_cfg, &rows)
            .with_context(|| format!("writing contrast subdir {:?}", subdir))?;
        let json_path = subdir.join("results.json");
        let html_path = subdir.join("report.html");
        if json_path.exists() {
            report::render_from_json(&json_path, &html_path).ok();
        }
        tracing::info!("contrast `{spec}` written to {:?}", subdir);
    }
    Ok(())
}

fn sanitize(s: &str) -> String {
    s.chars()
        .map(|c| if c.is_alphanumeric() { c } else { '_' })
        .collect()
}

fn write_advanced_events_table(
    path: &std::path::Path,
    events: &[advanced::AdvancedEvent],
) -> Result<()> {
    use std::io::Write;
    let f = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(f);
    writeln!(
        w,
        "event_id\tgene_id\tchrom\tkind\tcoords\tsupport\tnotes"
    )?;
    for e in events {
        let coords = e
            .coords
            .iter()
            .map(|(s, t)| format!("{s}-{t}"))
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{:.4}\t{}",
            e.event_id,
            e.gene_id.as_deref().unwrap_or(""),
            e.chrom,
            e.kind.short(),
            coords,
            e.support,
            e.notes,
        )?;
    }
    Ok(())
}

fn dump_annotation(args: ultimadsecaller::cli::DumpAnnotationArgs) -> Result<()> {
    let ann = annotation::parse(&args.annotation)?;
    let mut summary = serde_json::Map::new();
    for (gene_id, g) in &ann.genes {
        if let Some(filter) = &args.gene {
            if gene_id != filter {
                continue;
            }
        }
        let entry = serde_json::json!({
            "gene_id": g.gene_id,
            "chrom": g.chrom,
            "strand": g.strand.to_string(),
            "n_exons": g.graph.node_count(),
            "n_introns": g.graph.edge_count(),
            "n_transcripts": g.transcripts.len(),
            "transcripts": g.transcripts,
        });
        summary.insert(gene_id.clone(), entry);
    }
    let json = serde_json::to_string_pretty(&summary)?;
    if let Some(p) = args.out {
        std::fs::write(p, json)?;
    } else {
        println!("{json}");
    }
    Ok(())
}

fn junctions_cmd(args: ultimadsecaller::cli::JunctionsArgs) -> Result<()> {
    use ultimadsecaller::cli::MultimapStrategy;
    use std::io::Write;

    let f = std::fs::File::create(&args.out)?;
    let mut w = std::io::BufWriter::new(f);
    writeln!(w, "sample\tchrom\tdonor_end\tacceptor_start\tcount")?;
    for bam in &args.bams {
        let sample_id = bam
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("sample")
            .to_string();
        let sj = junctions::read_bam_junctions(
            bam,
            &sample_id,
            args.reference.as_deref(),
            args.min_mapq,
            args.min_overhang,
            MultimapStrategy::Primary,
        )?;
        for (j, c) in sj.counts {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{:.4}",
                sample_id, j.chrom, j.donor_end, j.acceptor_start, c
            )?;
        }
    }
    Ok(())
}
