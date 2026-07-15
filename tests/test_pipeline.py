"""End-to-end pipeline and output tests."""

import csv
import pathlib

from crumpling import cli
from crumpling.output import catalogs


FIXTURES = pathlib.Path(__file__).parent / "fixtures"


def test_cli_writes_stable_catalogs_and_corrected_structure(
    tmp_path: pathlib.Path,
) -> None:
    output = tmp_path / "run"
    exit_code = cli.main(
        [
            "--input",
            str(FIXTURES / "swaps.tsv"),
            "--config",
            str(FIXTURES / "config.ini"),
            "--output-dir",
            str(output),
        ]
    )
    assert exit_code == 0
    for name in ("snaps", "subavalanches", "avalanches"):
        assert (output / "catalogs" / f"{name}.tsv").is_file()
    assert (output / "run_metadata.json").is_file()

    with (output / "catalogs" / "snaps.tsv").open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        assert tuple(reader.fieldnames or ()) == catalogs.SNAP_COLUMNS
        first = next(reader)
    assert first["event_id"] == "1"
    assert first["subavalanche_id"] == "1"
    assert first["avalanche_id"] == "1"

    with (output / "structure" / "avalanche_structure.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    by_avalanche: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_avalanche.setdefault(row["avalanche_id"], []).append(row)
    for avalanche_rows in by_avalanche.values():
        assert max(float(row["normalized_event_count"]) for row in avalanche_rows) == 1.0
