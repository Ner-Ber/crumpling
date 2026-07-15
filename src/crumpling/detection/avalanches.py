"""Strategies for combining subavalanches into thermal avalanches."""

import dataclasses
import enum
import math
from typing import Sequence

import numpy as np

from crumpling import config
from crumpling import geometry
from crumpling import models


class AvalancheMechanism(enum.IntEnum):
    RELATIVE_TIME = 0
    CLOSEST_SPATIAL = 1
    LARGEST_FIXED_DISTANCE = 2
    LARGEST_HYBRID_TIME = 3
    LARGEST_VELOCITY = 4
    GREEDY_MAXIMAL_VELOCITY = 5


@dataclasses.dataclass(frozen=True, slots=True)
class DetectionResult:
    avalanches: tuple[models.Avalanche, ...]
    swap_to_avalanche: tuple[int, ...]
    subavalanche_to_avalanche: tuple[int, ...]


def _distance(left: models.Subavalanche, right: models.Subavalanche) -> float:
    return math.dist(left.metrics.center[:2], right.metrics.center[:2])


def _relative_gap(current: models.Subavalanche, previous: models.Subavalanche) -> float:
    return (
        current.metrics.time_p50 - previous.metrics.time_p50
    ) / current.metrics.time_p50


def _time_only_groups(
    ordered: Sequence[int],
    subavalanches: Sequence[models.Subavalanche],
    settings: config.AvalancheDetectionConfig,
) -> list[list[int]]:
    if not ordered:
        return []
    groups = [[ordered[0]]]
    for index in ordered[1:]:
        previous = groups[-1][-1]
        if _relative_gap(subavalanches[index], subavalanches[previous]) >= (
            settings.max_relative_time_gap
        ):
            groups.append([index])
        else:
            groups[-1].append(index)
    return groups


def _closest_groups(
    ordered: Sequence[int],
    subavalanches: Sequence[models.Subavalanche],
    settings: config.AvalancheDetectionConfig,
) -> list[list[int]]:
    active: list[list[int]] = []
    completed: list[list[int]] = []
    for index in ordered:
        current = subavalanches[index]
        stopped = [False] * len(active)
        closest: int | None = None
        closest_distance = math.inf
        for group_index, group in enumerate(active):
            if _relative_gap(current, subavalanches[group[-1]]) > (
                settings.max_relative_time_gap
            ):
                stopped[group_index] = True
            for previous_index in reversed(group):
                previous = subavalanches[previous_index]
                if _relative_gap(current, previous) > settings.max_relative_time_gap:
                    break
                distance = _distance(current, previous)
                if distance < settings.max_distance and distance < closest_distance:
                    closest = group_index
                    closest_distance = distance
        if closest is None:
            active.append([index])
            stopped.append(False)
        else:
            active[closest].append(index)
        for group_index in reversed(range(len(stopped))):
            if stopped[group_index]:
                completed.append(active.pop(group_index))
    completed.extend(active)
    return completed


def _can_join(
    current: models.Subavalanche,
    group: Sequence[int],
    subavalanches: Sequence[models.Subavalanche],
    settings: config.AvalancheDetectionConfig,
    mechanism: AvalancheMechanism,
) -> bool:
    if current.n_swaps < settings.minimum_subavalanche_swaps:
        return False
    for previous_index in reversed(group):
        previous = subavalanches[previous_index]
        absolute_gap = current.metrics.time_p50 - previous.metrics.time_p50
        relative_gap = absolute_gap / current.metrics.time_p50
        if mechanism == AvalancheMechanism.LARGEST_FIXED_DISTANCE:
            if relative_gap > settings.max_relative_time_gap:
                break
        elif (
            absolute_gap > settings.minimum_absolute_time_gap
            and relative_gap > settings.max_relative_time_gap
        ):
            break
        distance = _distance(current, previous)
        if mechanism == AvalancheMechanism.LARGEST_VELOCITY:
            if distance <= absolute_gap * settings.max_velocity:
                return True
        elif distance <= settings.max_distance:
            return True
    return False


def _largest_groups(
    ordered: Sequence[int],
    subavalanches: Sequence[models.Subavalanche],
    settings: config.AvalancheDetectionConfig,
    mechanism: AvalancheMechanism,
) -> list[list[int]]:
    groups: list[list[int]] = []
    sizes: list[int] = []
    for index in ordered:
        eligible = [
            group_index
            for group_index, group in enumerate(groups)
            if _can_join(
                subavalanches[index], group, subavalanches, settings, mechanism
            )
        ]
        if eligible:
            chosen = max(eligible, key=lambda group_index: sizes[group_index])
            groups[chosen].append(index)
            sizes[chosen] += subavalanches[index].n_swaps
        else:
            groups.append([index])
            sizes.append(subavalanches[index].n_swaps)
    return groups


def _velocity_join(
    current: models.Subavalanche,
    group: Sequence[int],
    subavalanches: Sequence[models.Subavalanche],
    settings: config.AvalancheDetectionConfig,
) -> bool:
    for previous_index in reversed(group):
        previous = subavalanches[previous_index]
        gap = current.metrics.time_p50 - previous.metrics.time_p50
        if gap / current.metrics.time_p50 > settings.max_relative_time_gap:
            break
        if _distance(current, previous) <= gap * settings.max_velocity:
            return True
    return False


def _greedy_maximal_groups(
    ordered: Sequence[int],
    subavalanches: Sequence[models.Subavalanche],
    settings: config.AvalancheDetectionConfig,
) -> list[list[int]]:
    candidates: list[list[int]] = []
    for index in ordered:
        eligible = [
            group
            for group in candidates
            if _velocity_join(subavalanches[index], group, subavalanches, settings)
        ]
        if eligible:
            for group in eligible:
                group.append(index)
        else:
            candidates.append([index])

    selected: list[list[int]] = []
    while candidates:
        chosen = max(
            candidates,
            key=lambda group: (
                sum(subavalanches[index].n_swaps for index in group),
                group,
            ),
        )
        chosen_copy = list(chosen)
        selected.append(chosen_copy)
        chosen_set = set(chosen_copy)
        candidates.remove(chosen)
        for group in candidates:
            group[:] = [index for index in group if index not in chosen_set]
        candidates[:] = [group for group in candidates if group]
    return selected


def _metrics(
    swap_indices: Sequence[int], events: Sequence[models.SwapEvent]
) -> models.AvalancheMetrics:
    ordered = tuple(sorted(swap_indices, key=lambda index: (events[index].time, index)))
    times = np.array([events[index].time for index in ordered], dtype=float)
    log_times = np.log10(times)
    points = np.array([events[index].position for index in ordered], dtype=float)
    center = points.mean(axis=0)
    return models.AvalancheMetrics(
        time_start=float(times[0]),
        time_end=float(times[-1]),
        time_average=float(times.mean()),
        log_time_average=float(log_times.mean()),
        time_p50=float(np.median(times)),
        duration=float(times[-1] - times[0]),
        log_duration=float(log_times[-1] - log_times[0]),
        temporal_gyration=float(np.sqrt(np.mean((times - times.mean()) ** 2))),
        log_temporal_gyration=float(
            np.sqrt(np.mean((log_times - log_times.mean()) ** 2))
        ),
        spatial_gyration=geometry.radius_of_gyration(points),
        diameter=geometry.maximum_pairwise_distance(points),
        direction_corrected_swaps=-sum(events[index].direction for index in ordered),
        center=(float(center[0]), float(center[1]), float(center[2])),
    )


def detect_avalanches(
    subavalanches: Sequence[models.Subavalanche],
    events: Sequence[models.SwapEvent],
    settings: config.AvalancheDetectionConfig,
) -> DetectionResult:
    """Run the configured grouping strategy and materialize assignments."""

    mechanism = AvalancheMechanism(settings.mechanism)
    ordered = sorted(
        range(len(subavalanches)),
        key=lambda index: (subavalanches[index].metrics.time_p50, index),
    )
    if mechanism == AvalancheMechanism.RELATIVE_TIME:
        groups = _time_only_groups(ordered, subavalanches, settings)
    elif mechanism == AvalancheMechanism.CLOSEST_SPATIAL:
        groups = _closest_groups(ordered, subavalanches, settings)
    elif mechanism in (
        AvalancheMechanism.LARGEST_FIXED_DISTANCE,
        AvalancheMechanism.LARGEST_HYBRID_TIME,
        AvalancheMechanism.LARGEST_VELOCITY,
    ):
        groups = _largest_groups(ordered, subavalanches, settings, mechanism)
    else:
        groups = _greedy_maximal_groups(ordered, subavalanches, settings)

    avalanches: list[models.Avalanche] = []
    swap_assignments = [-1] * len(events)
    subavalanche_assignments = [-1] * len(subavalanches)
    for avalanche_index, group in enumerate(groups):
        swaps = sorted(
            (
                swap
                for subavalanche_index in group
                for swap in subavalanches[subavalanche_index].swap_indices
            ),
            key=lambda index: (events[index].time, index),
        )
        value = models.Avalanche(
            index=avalanche_index,
            subavalanche_indices=tuple(group),
            swap_indices=tuple(swaps),
            metrics=_metrics(swaps, events),
        )
        avalanches.append(value)
        for subavalanche_index in group:
            subavalanche_assignments[subavalanche_index] = avalanche_index
        for swap_index in swaps:
            swap_assignments[swap_index] = avalanche_index

    if any(index < 0 for index in swap_assignments + subavalanche_assignments):
        raise RuntimeError("avalanche detection left an event unassigned")
    return DetectionResult(
        tuple(avalanches),
        tuple(swap_assignments),
        tuple(subavalanche_assignments),
    )
