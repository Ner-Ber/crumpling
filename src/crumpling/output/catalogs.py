"""Stable TSV schemas for analyzed event catalogs."""

import csv
import math
import pathlib
from typing import Iterable, Mapping, Sequence

from crumpling import models


SNAP_COLUMNS = (
    "simulation_id",
    "source_file",
    "event_id",
    "source_number",
    "time",
    "x",
    "y",
    "z",
    "surface_area",
    "direction",
    "subavalanche_id",
    "subavalanche_size",
    "avalanche_id",
    "avalanche_size",
)

SUBAVALANCHE_COLUMNS = (
    "simulation_id",
    "source_file",
    "event_id",
    "time_p50",
    "time_start",
    "time_end",
    "duration",
    "x",
    "y",
    "z",
    "n_swaps",
    "magnitude",
    "radius_of_gyration",
    "diameter",
    "inertia_1",
    "inertia_2",
    "simultaneous_count",
    "avalanche_id",
)

AVALANCHE_COLUMNS = (
    "simulation_id",
    "source_file",
    "event_id",
    "time_p50",
    "time_start",
    "time_end",
    "time_average",
    "duration",
    "x",
    "y",
    "z",
    "n_swaps",
    "n_subavalanches",
    "magnitude",
    "direction_corrected_swaps",
    "temporal_gyration",
    "log_temporal_gyration",
    "radius_of_gyration",
    "diameter",
)


def _write(
    path: pathlib.Path,
    columns: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_snaps(
    path: pathlib.Path, simulations: Sequence[models.SimulationResult]
) -> None:
    """Write one row per snap with explicit cluster memberships."""

    rows = []
    for simulation_id, result in enumerate(simulations, start=1):
        for index, event in enumerate(result.events):
            subavalanche_index = result.assignments.swap_to_subavalanche[index]
            avalanche_index = result.assignments.swap_to_avalanche[index]
            rows.append(
                {
                    "simulation_id": simulation_id,
                    "source_file": result.source.name,
                    "event_id": index + 1,
                    "source_number": event.source_number,
                    "time": event.time,
                    "x": event.position[0],
                    "y": event.position[1],
                    "z": event.position[2],
                    "surface_area": event.surface_area,
                    "direction": event.direction,
                    "subavalanche_id": subavalanche_index + 1,
                    "subavalanche_size": result.subavalanches[
                        subavalanche_index
                    ].n_swaps,
                    "avalanche_id": avalanche_index + 1,
                    "avalanche_size": result.avalanches[avalanche_index].n_swaps,
                }
            )
    _write(path, SNAP_COLUMNS, rows)


def write_subavalanches(
    path: pathlib.Path, simulations: Sequence[models.SimulationResult]
) -> None:
    """Write one row per subavalanche."""

    rows = []
    for simulation_id, result in enumerate(simulations, start=1):
        for value in result.subavalanches:
            metric = value.metrics
            rows.append(
                {
                    "simulation_id": simulation_id,
                    "source_file": result.source.name,
                    "event_id": value.index + 1,
                    "time_p50": metric.time_p50,
                    "time_start": metric.time_start,
                    "time_end": metric.time_end,
                    "duration": metric.duration,
                    "x": metric.center[0],
                    "y": metric.center[1],
                    "z": metric.center[2],
                    "n_swaps": value.n_swaps,
                    "magnitude": math.log10(value.n_swaps),
                    "radius_of_gyration": metric.radius_of_gyration,
                    "diameter": metric.diameter,
                    "inertia_1": metric.inertia_1,
                    "inertia_2": metric.inertia_2,
                    "simultaneous_count": metric.simultaneous_count,
                    "avalanche_id": (
                        result.assignments.subavalanche_to_avalanche[value.index]
                        + 1
                    ),
                }
            )
    _write(path, SUBAVALANCHE_COLUMNS, rows)


def write_avalanches(
    path: pathlib.Path, simulations: Sequence[models.SimulationResult]
) -> None:
    """Write one row per thermal avalanche."""

    rows = []
    for simulation_id, result in enumerate(simulations, start=1):
        for value in result.avalanches:
            metric = value.metrics
            rows.append(
                {
                    "simulation_id": simulation_id,
                    "source_file": result.source.name,
                    "event_id": value.index + 1,
                    "time_p50": metric.time_p50,
                    "time_start": metric.time_start,
                    "time_end": metric.time_end,
                    "time_average": metric.time_average,
                    "duration": metric.duration,
                    "x": metric.center[0],
                    "y": metric.center[1],
                    "z": metric.center[2],
                    "n_swaps": value.n_swaps,
                    "n_subavalanches": value.n_subavalanches,
                    "magnitude": math.log10(value.n_swaps),
                    "direction_corrected_swaps": metric.direction_corrected_swaps,
                    "temporal_gyration": metric.temporal_gyration,
                    "log_temporal_gyration": metric.log_temporal_gyration,
                    "radius_of_gyration": metric.spatial_gyration,
                    "diameter": metric.diameter,
                }
            )
    _write(path, AVALANCHE_COLUMNS, rows)
