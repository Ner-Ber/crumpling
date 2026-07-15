# crumpling

`crumpling` detects cascades of bistable-bond snaps in crumpling-sheet
simulations. It turns a chronological 24-column swap table into three related
catalogs:

1. individual snaps;
2. space-time subavalanches (cascades);
3. thermal avalanches assembled from subavalanches.

The original C++ program is retained in `legacy/`. The maintained implementation
is the installable Python package in `src/crumpling`.

## Install

```bash
python -m pip install -e ".[analysis,test]"
```

Core execution only requires NumPy. The `analysis` extra provides pandas,
Jupyter, Matplotlib, Seaborn, and SciPy.
Editable installation requires a current pip; upgrade it first with
`python -m pip install --upgrade pip` if an older environment rejects the
`pyproject.toml`-based editable build.

## Run trial32

```bash
crumpling-swapstat \
  --manifest configs/trial32/trial32.names.txt \
  --config configs/trial32/swapstat.params.ini \
  --output-dir outputs/trial32/current
```

Direct inputs are also supported:

```bash
python -m crumpling \
  --input data/raw/trial32/*__swaps.txt \
  --config configs/trial32/swapstat.params.ini \
  --output-dir outputs/my-run
```

See `docs/running-analysis.md`, `docs/data-layout.md`, and
`docs/output-schemas.md` for details.
