# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Breaking changes

- **Command line interface redesigned** with conventional `--long-options`. Single-dash long
  options were renamed, for example `-db` → `--db`, `-dir` → `--intermediate-dir`,
  `-mm2_options` → `--mm2-options`, `-exclude_partial` → `--exclude-partial`,
  `-infer_genes` → `--infer-genes`, `-sc` → `--copy-identity`, `-overlap` → `--max-overlap`,
  `-gap_open` → `--gap-open` and `-m` → `--minimap2`. Single-letter options are unchanged and
  now also have long forms (e.g. `-g/--gff`, `-a/--min-coverage`, `-p/--threads`). See the
  migration table in the README.
- **Python 3.12 or newer is required** (tested on 3.12, 3.13 and 3.14).
- Progress messages are written to standard error via the `logging` module; errors exit with
  status 1 and a message instead of a traceback.
- GFF3 annotations in which different features share an `ID` are rejected with an error listing
  the conflicting IDs and line numbers. Previously gffutils silently renamed one of the features,
  which could detach child features from their parents. Lines of one multi-line feature (same
  sequence, type, strand and `Parent`) may still share an `ID`.
- Intermediate files of the `--copies` stage are named `reference_all_copies_*` so that every
  stage's SAM file has a unique name.

### Added

- Pluggable alignment backends (`liftoff.align`): `Minimap2Aligner` (default),
  `PrecomputedSamAligner` (`--alignments DIR`, use SAM files produced elsewhere) and
  `CallableAligner` (delegate alignment to a Python function).
- Library API: `liftoff.pipeline.run_liftoff(config, aligner=None)` with a typed
  `LiftoffConfig` and a `LiftoffResult` return value.
- Pyodide/WebAssembly support: every runtime dependency is pure Python or shipped by Pyodide,
  worker processes are only used where available, and the test suite runs under Pyodide in CI.
- `--all-feature-types` to lift every top-level feature type in the annotation instead of only
  genes (and types listed with `-f`).
- `--exclude-feature-types TYPES` (comma-separated, repeatable) to skip top-level feature types
  such as `region`, with a warning for types that are not present as top-level features.
- `--no-cds` to skip CDS status annotation, `-v/--verbose` and `-q/--quiet`.
- `python -m liftoff` entry point and `py.typed` marker.
- `environment.yml` conda environment, GitHub Actions CI (ruff, mypy, pytest on Linux and macOS
  with Python 3.12–3.14, Pyodide tests), pre-commit hooks and this changelog.
- Test suite organised into `tests/unit` and `tests/integration` (auto-applied `unit` /
  `integration` markers) with shared factory fixtures, module-scoped pipeline runs, integration
  scenarios parametrized across alignment backends, unit tests for every module, strict
  resource-warning checks and a 90% branch-coverage threshold. Includes whole-genome and
  chromosome I integration tests against recorded reference outputs and parasail comparison
  corpora.

### Changed

- Packaging migrated from `setup.py` to `pyproject.toml`, built with hatchling; the version is
  derived from git tags with hatch-vcs. The package now uses a `src/` layout, and tests moved to
  the repository root.
- Logo images moved to `docs/images/`; the README references them by relative path.
- README: tested quick-start examples, corrected installation guidance, complete output attribute
  list, notes on running the SQLite feature database under Pyodide, and a future-development note
  on the performance of the pure-Python parasail implementation.
- Code reorganised into `align`, `io` and `mapping` subpackages with typed dataclass models,
  full type hints (mypy strict) and numpydoc docstrings.
- **Removed the `parasail` dependency**: polishing uses a NumPy port of
  `sg_dx_trace_scan_sat` that reproduces parasail's scores and tracebacks exactly.
- **Removed the `pysam` dependency** in favour of a small pure-Python SAM reader with identical
  semantics.
- **Removed the `interlap` and `ujson` dependencies** (replaced by an internal interval index
  and the standard library `json`).
- Dependency minimums raised to biopython 1.87, gffutils 0.14, networkx 3.6, numpy 2.4 and
  pyfaidx 0.9.
- Feature database queries use explicit ordering instead of string-built `IN (...)` lists.
- Alignment jobs run with `--threads` > 1 are collected in a deterministic order, and overlapping
  features are re-mapped in sorted order, making results independent of process scheduling and
  hash randomisation.
- minimap2 options required by Liftoff are merged by whole tokens rather than substring
  matching, and minimap2 failures are reported as errors.
- Faster overall (whole-genome yeast lift-over in about 8 s instead of about 20 s).

### Fixed

- Broken CDSs are now trimmed to the longest intact ORF (at least 180 bases, searched in all
  three frames), with CDS phases recomputed; `valid_ORF` flags describe the trimmed CDS.
  Previously the ORF was found but the CDS was never adjusted.
- Polishing no longer shifts reference exon coordinates: splice-site extensions are computed
  without modifying the exons, instead of being added and imperfectly removed.
- Selecting a polished gene compares valid ORF counts, sequence identity and coverage
  numerically; previously strings were compared (so `'9' > '10'`), and a polished gene whose CDS
  failed to lift always replaced the original.
- `--flank` is applied from the annotated coordinates at every stage, so it no longer widens
  features repeatedly, and single-level features are written at their annotated coordinates
  instead of including the flank.
- `-cds` could not be disabled.
- A missing space when appending `--end-bonus 5` to user-supplied minimap2 options.
- Output files and SQLite connections are now closed properly.
- Piping output into a command that exits early (e.g. `liftoff ... | head`) no longer prints a
  `BrokenPipeError` traceback.
- Plain gzip-compressed FASTA input is reported as a clear error (only BGZF is supported) instead
  of an unhandled pyfaidx exception.
- `--help` no longer shows meaningless `(default: None)` / `(default: False)` annotations.
- The previous tests never compared GFF output lines; the expected outputs were stale and have
  been regenerated.
