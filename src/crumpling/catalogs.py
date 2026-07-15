"""Optional pandas catalog builders for interactive analysis."""

from typing import Any

import numpy as np

from crumpling import models

try:
    import pandas as pd
except ImportError:
    pd = None


def _pandas() -> Any:
    if pd is None:
        raise ImportError(
            "Catalog DataFrames require the analysis extra: "
            "pip install 'crumpling[analysis]'"
        )
    return pd


def _temporal_columns(frame: Any) -> Any:
    frame = frame.sort_values("time", kind="mergesort").reset_index(drop=True)
    frame["dt"] = frame["time"].diff()
    frame["log10_dt"] = np.nan
    positive = frame["dt"] > 0
    frame.loc[positive, "log10_dt"] = np.log10(frame.loc[positive, "dt"])
    frame["cumulative_count"] = np.arange(1, len(frame) + 1)
    return frame


def snaps(result: models.SimulationResult) -> Any:
    """Build an earthquake-style snap catalog."""

    pd = _pandas()
    rows = []
    for index, event in enumerate(result.events):
        rows.append(
            {
                "event_id": index + 1,
                "source_number": event.source_number,
                "time": event.time,
                "x": event.position[0],
                "y": event.position[1],
                "z": event.position[2],
                "magnitude": 0.0,
                "surface_area": event.surface_area,
                "subavalanche_id": (
                    result.assignments.swap_to_subavalanche[index] + 1
                ),
                "avalanche_id": result.assignments.swap_to_avalanche[index] + 1,
            }
        )
    return _temporal_columns(pd.DataFrame(rows))


def subavalanches(result: models.SimulationResult) -> Any:
    """Build a cascade catalog with log10(number of swaps) magnitude."""

    pd = _pandas()
    rows = [
        {
            "event_id": value.index + 1,
            "time": value.metrics.time_p50,
            "time_start": value.metrics.time_start,
            "time_end": value.metrics.time_end,
            "duration": value.metrics.duration,
            "x": value.metrics.center[0],
            "y": value.metrics.center[1],
            "z": value.metrics.center[2],
            "n_swaps": value.n_swaps,
            "magnitude": float(np.log10(value.n_swaps)),
            "radius_of_gyration": value.metrics.radius_of_gyration,
            "diameter": value.metrics.diameter,
            "simultaneous_count": value.metrics.simultaneous_count,
            "avalanche_id": (
                result.assignments.subavalanche_to_avalanche[value.index] + 1
            ),
        }
        for value in result.subavalanches
    ]
    return _temporal_columns(pd.DataFrame(rows))


def avalanches(result: models.SimulationResult) -> Any:
    """Build a thermal-avalanche catalog with log10(number of swaps) magnitude."""

    pd = _pandas()
    rows = [
        {
            "event_id": value.index + 1,
            "time": value.metrics.time_p50,
            "time_start": value.metrics.time_start,
            "time_end": value.metrics.time_end,
            "duration": value.metrics.duration,
            "x": value.metrics.center[0],
            "y": value.metrics.center[1],
            "z": value.metrics.center[2],
            "n_swaps": value.n_swaps,
            "n_subavalanches": value.n_subavalanches,
            "magnitude": float(np.log10(value.n_swaps)),
            "radius_of_gyration": value.metrics.spatial_gyration,
            "diameter": value.metrics.diameter,
            "temporal_gyration": value.metrics.temporal_gyration,
        }
        for value in result.avalanches
    ]
    return _temporal_columns(pd.DataFrame(rows))
