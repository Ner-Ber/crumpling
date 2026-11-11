"""Optional pandas catalog builders for interactive analysis.

Clustering itself is unchanged: pass a ``SimulationResult`` from
``crumpling.pipeline`` into ``build_catalogs`` to get the three catalogs
plus parent/child membership links.
"""

from __future__ import annotations

import dataclasses
import pathlib
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


def _parent_id(index: int) -> int | Any:
    """Map a 0-based assignment to a 1-based catalog ID; unset → NA."""

    if index < 0:
        return _pandas().NA
    return index + 1


def snaps(result: models.SimulationResult) -> Any:
    """Build an earthquake-style snap catalog.

    Parent links: ``cascade_id`` and ``avalanche_id``.
    Unassigned parents are pandas ``NA``.
    """

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
                "cascade_id": _parent_id(
                    result.assignments.swap_to_cascade[index]
                ),
                "avalanche_id": _parent_id(
                    result.assignments.swap_to_avalanche[index]
                ),
            }
        )
    return _temporal_columns(pd.DataFrame(rows))


def cascades(result: models.SimulationResult) -> Any:
    """Build a cascade catalog with log10(n_swaps) magnitude.

    Parent link: ``avalanche_id``. Unassigned parents are pandas ``NA``.
    """

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
            "avalanche_id": _parent_id(
                result.assignments.cascade_to_avalanche[value.index]
            ),
        }
        for value in result.cascades
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
            "n_cascades": value.n_cascades,
            "magnitude": float(np.log10(value.n_swaps)),
            "radius_of_gyration": value.metrics.spatial_gyration,
            "diameter": value.metrics.diameter,
            "temporal_gyration": value.metrics.temporal_gyration,
        }
        for value in result.avalanches
    ]
    return _temporal_columns(pd.DataFrame(rows))


@dataclasses.dataclass(frozen=True, slots=True)
class EventCatalogs:
    """Snap / cascade / avalanche catalogs from one simulation analysis.

    Membership (1-based IDs, ``NA`` when unset):
    - ``snaps.cascade_id``, ``snaps.avalanche_id``
    - ``cascades.avalanche_id``
    - ``cascades.member_snap_ids``
    - ``avalanches.member_cascade_ids``
    """

    snaps: Any
    cascades: Any
    avalanches: Any
    source: pathlib.Path


def build_catalogs(result: models.SimulationResult) -> EventCatalogs:
    """Wrap a ``SimulationResult`` as the three catalogs plus membership lists.

    Does not re-run clustering; only reshapes pipeline output. Child-member
    ID columns are attached here so the standalone ``snaps`` / ``cascades`` /
    ``avalanches`` helpers stay schema-stable.
    """

    snap_frame = snaps(result)
    cascade_frame = cascades(result).copy()
    avalanche_frame = avalanches(result).copy()

    snap_members = {
        value.index + 1: tuple(index + 1 for index in value.swap_indices)
        for value in result.cascades
    }
    cascade_members = {
        value.index + 1: tuple(index + 1 for index in value.cascade_indices)
        for value in result.avalanches
    }
    cascade_frame["member_snap_ids"] = cascade_frame["event_id"].map(snap_members)
    avalanche_frame["member_cascade_ids"] = avalanche_frame["event_id"].map(
        cascade_members
    )

    return EventCatalogs(
        snaps=snap_frame,
        cascades=cascade_frame,
        avalanches=avalanche_frame,
        source=result.source,
    )
