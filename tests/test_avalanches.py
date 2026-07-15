"""Avalanche strategy tests."""

import dataclasses
import pathlib

import pytest

from crumpling import config
from crumpling import parsing
from crumpling.detection import avalanches
from crumpling.detection import subavalanches


FIXTURES = pathlib.Path(__file__).parent / "fixtures"


@pytest.mark.parametrize("mechanism", range(6))
def test_all_mechanisms_assign_every_event(mechanism: int) -> None:
    events = parsing.read_swap_file(FIXTURES / "swaps.tsv")
    subavalanche_result = subavalanches.detect_subavalanches(
        events,
        config.SubavalancheDetectionConfig(
            average_bond_length=1.0,
            max_time_gap=2.0,
            max_distance=5.0,
        ),
    )
    settings = dataclasses.replace(
        config.AvalancheDetectionConfig(),
        mechanism=mechanism,
        max_relative_time_gap=0.5,
        minimum_absolute_time_gap=1.0,
        max_distance=25.0,
        max_velocity=25.0,
    )
    result = avalanches.detect_avalanches(
        subavalanche_result.subavalanches, events, settings
    )
    assert sorted(
        swap for value in result.avalanches for swap in value.swap_indices
    ) == list(range(len(events)))
    assert all(index >= 0 for index in result.swap_to_avalanche)
    assert all(index >= 0 for index in result.subavalanche_to_avalanche)


def test_time_only_strategy_splits_large_relative_gap() -> None:
    events = parsing.read_swap_file(FIXTURES / "swaps.tsv")
    subavalanche_result = subavalanches.detect_subavalanches(
        events,
        config.SubavalancheDetectionConfig(
            average_bond_length=1.0,
            max_time_gap=2.0,
            max_distance=5.0,
        ),
    )
    result = avalanches.detect_avalanches(
        subavalanche_result.subavalanches,
        events,
        dataclasses.replace(
            config.AvalancheDetectionConfig(),
            max_relative_time_gap=0.6,
        ),
    )
    assert len(result.avalanches) == 2
