# Security

## Model

ultimaDSEcaller is an offline command-line tool. It reads local genomic
files (BAM/CRAM, GTF/GFF3, FASTA, sample sheets), writes local results, and
makes **no network calls**. There is no server, no account, and no telemetry.
The published HTML report loads the Plotly.js library from a CDN for
interactivity; the analysis itself never leaves your machine, and the report
can be viewed offline if the CDN asset is cached or vendored.

What this leaves, and how it is handled:

| Surface | Handling |
|---|---|
| Untrusted input files (BAM/GTF/FASTA) | Parsed in safe Rust with no `unsafe` blocks; malformed records yield typed errors (`E0010`/`E0020`/`E0030`), not crashes or memory corruption. |
| HTML report rendering | All values are inserted via DOM `textContent`, never `innerHTML`, so a gene/event identifier crafted to contain markup cannot inject script into the report. |
| Memory safety | 100% safe Rust; exercised under `cargo miri` in CI and verifiable under valgrind on the release binary. |
| Dependencies | Pinned in `Cargo.lock`; CI builds with `-D warnings`. |
| Containers | The runtime image bundles pinned tool versions; the caller is a static musl binary with no dynamic-library attack surface. |

Known, documented limitations:

- The interactive report depends on a CDN-hosted Plotly.js. For fully
  air-gapped environments, vendor the library and adjust the report template.
- The tool trusts that input BAMs correspond to the provided annotation and
  reference; it validates structure, not biological provenance.

## Reporting a vulnerability

Email **bareke.eric@gmail.com** with a description and reproduction steps.
Please do not open public issues for exploitable problems before a fix is
available. You can expect an acknowledgement within a few days; fixes are
best-effort but security reports get priority.
