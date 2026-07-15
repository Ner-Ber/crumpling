# Data and repository layout

The importable package is isolated under `src/crumpling`. Runtime artifacts do
not live inside the package:

- `data/raw/trial32/` contains the immutable 24-column full-format input.
- `data/raw/trial33/` contains the six-column reduced exploration input.
- `configs/trial32/` contains authored run parameters and the optional manifest.
- `outputs/<run>/` contains disposable generated catalogs and structure tables.
- `notebooks/` contains interactive analysis source.
- `legacy/` contains the original C++ implementation and the former monolithic
  Python port for reference.
- `tests/fixtures/` contains small synthetic inputs suitable for automated tests.

Manifest input paths are resolved relative to the manifest itself, not the
process working directory. Generated output is ignored by Git; raw sample inputs
and configurations remain versionable.
