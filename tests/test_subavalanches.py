"""Subavalanche clustering tests."""

import pathlib

from crumpling import config
from crumpling import parsing
from crumpling.detection import subavalanches


FIXTURES = pathlib.Path(__file__).parent / "fixtures"


def test_space_time_boundaries_and_assignments() -> None:
    events = parsing.read_swap_file(FIXTURES / "swaps.tsv")
    settings = config.SubavalancheDetectionConfig(
        average_bond_length=1.0,
        max_time_gap=2.0,
        max_distance=5.0,
    )
    result = subavalanches.detect_subavalanches(events, settings)
    assert [value.swap_indices for value in result.subavalanches] == [
        (0, 1),
        (2,),
        (3,),
    ]
    assert result.swap_to_subavalanche == (0, 0, 1, 2)
    assert result.subavalanches[0].metrics.diameter == 2.0
    assert result.subavalanches[0].metrics.simultaneous_count == 1
