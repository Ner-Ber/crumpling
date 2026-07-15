"""Typed domain models for crumpling event analysis."""

import dataclasses
import pathlib
from typing import Sequence


@dataclasses.dataclass(frozen=True, slots=True)
class SwapEvent:
    """One bistable-bond snap read from a 24-column simulation table."""

    source_number: int
    timestep: float
    time: float
    direction: int
    x_size: float
    y_size: float
    surface_area: float
    position: tuple[float, float, float]
    atom1: int
    atom2: int


@dataclasses.dataclass(frozen=True, slots=True)
class SubavalancheMetrics:
    time_p50: float
    time_start: float
    time_end: float
    duration: float
    center: tuple[float, float, float]
    radius_of_gyration: float
    inertia_1: float
    inertia_2: float
    diameter: float
    simultaneous_count: int = 0


@dataclasses.dataclass(frozen=True, slots=True)
class Subavalanche:
    """A space-time cluster of swaps."""

    index: int
    swap_indices: tuple[int, ...]
    metrics: SubavalancheMetrics

    @property
    def n_swaps(self) -> int:
        return len(self.swap_indices)


@dataclasses.dataclass(frozen=True, slots=True)
class AvalancheMetrics:
    time_start: float
    time_end: float
    time_average: float
    log_time_average: float
    time_p50: float
    duration: float
    log_duration: float
    temporal_gyration: float
    log_temporal_gyration: float
    spatial_gyration: float
    diameter: float
    direction_corrected_swaps: int
    center: tuple[float, float, float]


@dataclasses.dataclass(frozen=True, slots=True)
class Avalanche:
    """A thermal avalanche composed of one or more subavalanches."""

    index: int
    subavalanche_indices: tuple[int, ...]
    swap_indices: tuple[int, ...]
    metrics: AvalancheMetrics

    @property
    def n_subavalanches(self) -> int:
        return len(self.subavalanche_indices)

    @property
    def n_swaps(self) -> int:
        return len(self.swap_indices)


@dataclasses.dataclass(frozen=True, slots=True)
class AssignmentMaps:
    swap_to_subavalanche: tuple[int, ...]
    swap_to_avalanche: tuple[int, ...]
    subavalanche_to_avalanche: tuple[int, ...]


@dataclasses.dataclass(frozen=True, slots=True)
class SimulationSpec:
    path: pathlib.Path
    max_simulation_time: float | None = None


@dataclasses.dataclass(frozen=True, slots=True)
class SimulationResult:
    source: pathlib.Path
    events: tuple[SwapEvent, ...]
    subavalanches: tuple[Subavalanche, ...]
    avalanches: tuple[Avalanche, ...]
    assignments: AssignmentMaps


@dataclasses.dataclass(frozen=True, slots=True)
class AnalysisResult:
    simulations: tuple[SimulationResult, ...]


def event_positions(events: Sequence[SwapEvent]) -> tuple[tuple[float, float, float], ...]:
    """Return event positions without exposing a NumPy type in the model layer."""

    return tuple(event.position for event in events)
