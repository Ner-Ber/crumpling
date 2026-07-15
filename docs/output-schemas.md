# Output schemas

Each run writes a stable directory:

```text
catalogs/snaps.tsv
catalogs/subavalanches.tsv
catalogs/avalanches.tsv
structure/avalanche_structure.tsv
run_metadata.json
```

External event IDs are consistently one-based. Internal Python indices remain
zero-based.

## Snap catalog

One row per bond snap. It includes time, 3D bond-center location, direction,
surface area, source number, and the IDs and sizes of its parent subavalanche and
avalanche. Snap magnitude is treated as a unit event (`0` on a logarithmic
magnitude scale) by the DataFrame API.

## Subavalanche catalog

One row per space-time cascade. Its earthquake-like origin time is the median
(`time_p50`), location is the centroid, and magnitude is
`log10(number of swaps)`. Geometry includes radius of gyration, diameter, and
planar principal inertia moments. `simultaneous_count` includes the event itself.

## Avalanche catalog

One row per thermal avalanche. It includes start/end/median/average time,
centroid, swap and subavalanche counts, logarithmic size, temporal and spatial
gyration, diameter, and direction-corrected swap count.

## Avalanche structure

This is a long-format growth-curve table. Each row identifies an avalanche,
coordinate system (`linear_time`, `log_time`, or `normalized_log_time`), sample
coordinate, cumulative event count, and normalized cumulative count. Counts use
actual counts (`1..N`), correcting the legacy zero-based `0..N-1` output; every
curve therefore reaches normalized count `1.0`.

## Intentional differences from legacy output

- malformed rows fail instead of duplicating the previous event;
- paths are explicit and output is never written beside immutable raw data;
- the malformed double-tab markup header is removed;
- all external IDs use one-based indexing;
- avalanche metrics always use chronological member order;
- avalanche structure uses corrected cumulative counts and a long schema.
