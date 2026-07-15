"""Input parsing tests."""

import pathlib

import pytest

from crumpling import parsing


FIXTURES = pathlib.Path(__file__).parent / "fixtures"


def test_reads_24_column_swap_file() -> None:
    events = parsing.read_swap_file(FIXTURES / "swaps.tsv")
    assert len(events) == 4
    assert events[0].position == (0.0, 0.0, 0.0)
    assert events[-1].atom2 == 17


def test_rejects_malformed_row() -> None:
    with pytest.raises(parsing.InputFormatError, match="expected 24"):
        parsing.parse_swap_row(["1", "2"], 7)


def test_manifest_paths_are_relative_to_manifest(tmp_path: pathlib.Path) -> None:
    manifest = tmp_path / "inputs.txt"
    manifest.write_text("files 1\n../input.tsv 100\n")
    specs = parsing.read_manifest(manifest)
    assert specs[0].path == (tmp_path / "../input.tsv").resolve()
