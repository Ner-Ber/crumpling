# Running an analysis

## Command line

Choose either direct input files:

```bash
crumpling-swapstat \
  --input data/raw/trial32/example__swaps.txt \
  --config configs/trial32/swapstat.params.ini \
  --output-dir outputs/example
```

or a multi-file manifest:

```bash
crumpling-swapstat \
  --manifest configs/trial32/trial32.names.txt \
  --config configs/trial32/swapstat.params.ini \
  --output-dir outputs/trial32/current
```

Input tables must contain exactly 24 tab-separated columns, positive
chronologically ordered times, and one header row. The six-column
`__swaps-short.txt` format is intended only for notebook exploration.

## Cursor debugger

Select **crumpling-swapstat: trial32** in Run and Debug and press F5. The launch
configuration runs `python -m crumpling` from the repository root with
repository-absolute arguments; it does not depend on a data-directory working
directory or a custom `PYTHONPATH`.

## Configuration provenance

The trial32 `dt=10` and `dr=400` values were inferred from the source filename.
Other values started from defaults exposed by the legacy C++ program. They are
an example configuration, not a recovered authoritative simulation parameter
file.

The avalanche mechanism values are:

0. relative median-time gap;
1. relative time plus closest spatial candidate;
2. relative time plus fixed distance, largest candidate wins;
3. hybrid absolute/relative time plus fixed distance;
4. hybrid time plus propagation velocity;
5. greedy maximal velocity-compatible grouping.

Every run writes `run_metadata.json` with source paths and the effective parsed
configuration.
