"""Tests for interactive EventCatalogs wrapping."""

import pathlib

from crumpling import catalogs
from crumpling import config
from crumpling import pipeline


FIXTURES = pathlib.Path(__file__).parent / "fixtures"


def test_build_catalogs_links_snaps_cascades_avalanches() -> None:
    settings = config.load_config(FIXTURES / "config.ini")
    result = pipeline.analyze_file(FIXTURES / "swaps.tsv", settings)
    bundle = catalogs.build_catalogs(result)

    assert len(bundle.snaps) == len(result.events)
    assert len(bundle.cascades) == len(result.cascades)
    assert len(bundle.avalanches) == len(result.avalanches)

    first_snap = bundle.snaps.iloc[0]
    cascade_id = int(first_snap["cascade_id"])
    avalanche_id = int(first_snap["avalanche_id"])
    cascade = bundle.cascades.set_index("event_id").loc[cascade_id]
    assert int(cascade["avalanche_id"]) == avalanche_id
    assert int(first_snap["event_id"]) in cascade["member_snap_ids"]

    avalanche = bundle.avalanches.set_index("event_id").loc[avalanche_id]
    assert cascade_id in avalanche["member_cascade_ids"]
