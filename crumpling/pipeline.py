"""Programmatic orchestration for crumpling analyses."""

import pathlib
from typing import Sequence

from crumpling import config
from crumpling import models
from crumpling import parsing
from crumpling.detection import avalanches
from crumpling.detection import subavalanches


def analyze_events(
    events: Sequence[models.SwapEvent],
    source: pathlib.Path,
    settings: config.AnalysisConfig,
) -> models.SimulationResult:
    """Analyze already-parsed events without performing output I/O."""

    subavalanche_result = subavalanches.detect_subavalanches(
        events, settings.subavalanches
    )
    avalanche_result = avalanches.detect_avalanches(
        subavalanche_result.subavalanches,
        events,
        settings.avalanches,
    )
    return models.SimulationResult(
        source=source,
        events=tuple(events),
        subavalanches=subavalanche_result.subavalanches,
        avalanches=avalanche_result.avalanches,
        assignments=models.AssignmentMaps(
            swap_to_subavalanche=subavalanche_result.swap_to_subavalanche,
            swap_to_avalanche=avalanche_result.swap_to_avalanche,
            subavalanche_to_avalanche=(
                avalanche_result.subavalanche_to_avalanche
            ),
        ),
    )


def analyze_file(
    path: pathlib.Path | str, settings: config.AnalysisConfig
) -> models.SimulationResult:
    """Parse and analyze one full-format swap file."""

    source = pathlib.Path(path).resolve()
    return analyze_events(parsing.read_swap_file(source), source, settings)


def run_analysis(
    simulations: Sequence[models.SimulationSpec],
    settings: config.AnalysisConfig,
) -> models.AnalysisResult:
    """Analyze all simulations in a manifest or direct CLI invocation."""

    return models.AnalysisResult(
        simulations=tuple(analyze_file(spec.path, settings) for spec in simulations)
    )
