"""Space-time clustering of swaps into subavalanches."""

import dataclasses
import math
from typing import Sequence

import numpy as np

from crumpling import config
from crumpling import geometry
from crumpling import models


@dataclasses.dataclass(frozen=True, slots=True)
class DetectionResult:
    subavalanches: tuple[models.Subavalanche, ...]
    swap_to_subavalanche: tuple[int, ...]


def _metrics(
    swap_indices: Sequence[int],
    events: Sequence[models.SwapEvent],
    average_bond_length: float,
) -> models.SubavalancheMetrics:
    ordered = tuple(sorted(swap_indices, key=lambda index: (events[index].time, index)))
    times = np.array([events[index].time for index in ordered], dtype=float)
    points = np.array([events[index].position for index in ordered], dtype=float)
    center = points.mean(axis=0)
    inertia_1, inertia_2 = geometry.planar_inertia_eigenvalues(points)
    return models.SubavalancheMetrics(
        time_p50=float(np.median(times)),
        time_start=float(times[0]),
        time_end=float(times[-1]),
        duration=float(times[-1] - times[0]),
        center=(float(center[0]), float(center[1]), float(center[2])),
        radius_of_gyration=geometry.radius_of_gyration(points),
        inertia_1=inertia_1,
        inertia_2=inertia_2,
        diameter=geometry.maximum_pairwise_distance(points) + average_bond_length,
    )


def _overlap_count(
    current: models.Subavalanche, values: Sequence[models.Subavalanche]
) -> int:
    start = current.metrics.time_start
    finish = current.metrics.time_end
    return sum(
        other.metrics.time_start <= finish and other.metrics.time_end >= start
        for other in values
    )


def detect_subavalanches(
    events: Sequence[models.SwapEvent],
    settings: config.SubavalancheDetectionConfig,
) -> DetectionResult:
    """Detect active space-time clusters using the original closest-group rule."""

    active: list[list[int]] = []
    completed: list[list[int]] = []
    for event_index, event in enumerate(events):
        closest_group: int | None = None
        closest_distance = math.inf
        stopped = [False] * len(active)

        for group_index, group in enumerate(active):
            if event.time - events[group[-1]].time > settings.max_time_gap:
                stopped[group_index] = True
            for previous_index in reversed(group):
                previous = events[previous_index]
                if event.time - previous.time > settings.max_time_gap:
                    break
                distance = math.dist(event.position, previous.position)
                if (
                    distance < settings.max_distance
                    and distance < closest_distance
                ):
                    closest_group = group_index
                    closest_distance = distance

        if closest_group is None:
            active.append([event_index])
            stopped.append(False)
        else:
            active[closest_group].append(event_index)

        completed.extend(
            group for group_index, group in enumerate(active) if stopped[group_index]
        )
        active = [
            group for group_index, group in enumerate(active) if not stopped[group_index]
        ]

    completed.extend(active)
    subavalanches = [
        models.Subavalanche(
            index=index,
            swap_indices=tuple(
                sorted(group, key=lambda swap: (events[swap].time, swap))
            ),
            metrics=_metrics(group, events, settings.average_bond_length),
        )
        for index, group in enumerate(completed)
    ]
    subavalanches = [
        dataclasses.replace(
            value,
            metrics=dataclasses.replace(
                value.metrics,
                simultaneous_count=_overlap_count(value, subavalanches),
            ),
        )
        for value in subavalanches
    ]
    assignments = [-1] * len(events)
    for value in subavalanches:
        for swap_index in value.swap_indices:
            assignments[swap_index] = value.index
    if any(index < 0 for index in assignments):
        raise RuntimeError("subavalanche detection left an event unassigned")
    return DetectionResult(tuple(subavalanches), tuple(assignments))
