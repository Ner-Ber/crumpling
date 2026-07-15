"""Avalanche growth-curve sampling and serialization."""

import csv
import math
import pathlib
from typing import Sequence

import numpy as np

from crumpling import config
from crumpling import models


STRUCTURE_COLUMNS = (
    "simulation_id",
    "source_file",
    "avalanche_id",
    "coordinate",
    "coordinate_value",
    "event_count",
    "normalized_event_count",
)


def _eligible(
    value: models.Avalanche,
    settings: config.StructureConfig,
    window: config.AnalysisWindow,
) -> bool:
    return (
        value.n_swaps >= settings.minimum_swaps
        and value.metrics.time_end >= window.minimum_time
        and (
            window.maximum_time is None
            or value.metrics.time_start <= window.maximum_time
        )
    )


def _sample(
    coordinates: np.ndarray, step: float, endpoint: float
) -> list[tuple[float, int]]:
    samples = np.arange(0.0, endpoint + step * 0.5, step)
    if len(samples) == 0 or samples[-1] < endpoint:
        samples = np.append(samples, endpoint)
    elif samples[-1] > endpoint:
        samples[-1] = endpoint
    counts = np.searchsorted(coordinates, samples, side="right")
    return [(float(sample), int(count)) for sample, count in zip(samples, counts)]


def write_structure(
    path: pathlib.Path,
    simulations: Sequence[models.SimulationResult],
    settings: config.StructureConfig,
    window: config.AnalysisWindow,
) -> None:
    """Write corrected cumulative avalanche growth curves in long format."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=STRUCTURE_COLUMNS, delimiter="\t")
        writer.writeheader()
        for simulation_id, result in enumerate(simulations, start=1):
            for value in result.avalanches:
                if not _eligible(value, settings, window):
                    continue
                times = np.array(
                    [result.events[index].time for index in value.swap_indices],
                    dtype=float,
                )
                elapsed = times - value.metrics.time_start
                modes = [
                    (
                        "linear_time",
                        elapsed,
                        settings.linear_time_step,
                        float(elapsed[-1]),
                    )
                ]
                if value.metrics.time_start > 0:
                    logarithmic = np.log10(times) - math.log10(
                        value.metrics.time_start
                    )
                    modes.append(
                        (
                            "log_time",
                            logarithmic,
                            settings.logarithmic_time_step,
                            float(logarithmic[-1]),
                        )
                    )
                    if value.n_swaps >= 2 and logarithmic[-1] > 0:
                        normalized = logarithmic / logarithmic[-1]
                        modes.append(
                            (
                                "normalized_log_time",
                                normalized,
                                settings.normalized_log_time_step,
                                1.0,
                            )
                        )
                for mode, coordinates, step, endpoint in modes:
                    for coordinate, count in _sample(coordinates, step, endpoint):
                        writer.writerow(
                            {
                                "simulation_id": simulation_id,
                                "source_file": result.source.name,
                                "avalanche_id": value.index + 1,
                                "coordinate": mode,
                                "coordinate_value": coordinate,
                                "event_count": count,
                                "normalized_event_count": count / value.n_swaps,
                            }
                        )
