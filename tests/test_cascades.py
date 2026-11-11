"""Cascade clustering tests."""

import pathlib

from crumpling import config
from crumpling import parsing
from crumpling.detection import cascades


FIXTURES = pathlib.Path(__file__).parent / "fixtures"


def test_space_time_boundaries_and_assignments() -> None:
    events = parsing.read_swap_file(FIXTURES / "swaps.tsv")
    settings = config.CascadeDetectionConfig(
        average_bond_length=1.0,
        max_time_gap=2.0,
        max_distance=5.0,
    )
    result = cascades.detect_cascades(events, settings)
    assert [value.swap_indices for value in result.cascades] == [
        (0, 1),
        (2,),
        (3,),
    ]
    assert result.swap_to_cascade == (0, 0, 1, 2)
    assert result.cascades[0].metrics.diameter == 2.0
    assert result.cascades[0].metrics.simultaneous_count == 1
