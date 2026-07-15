"""Command-line interface for swap-stat analysis."""

import argparse
import dataclasses
import json
import pathlib
import sys
from typing import Sequence

from crumpling import config
from crumpling import models
from crumpling import parsing
from crumpling import pipeline
from crumpling.output import catalogs
from crumpling.output import structure


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="crumpling-swapstat",
        description="Cluster full-format bond snaps into subavalanches and avalanches.",
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--input",
        type=pathlib.Path,
        nargs="+",
        help="one or more 24-column full-format swap tables",
    )
    source.add_argument(
        "--manifest",
        type=pathlib.Path,
        help="legacy multi-file input manifest",
    )
    parser.add_argument(
        "--config", type=pathlib.Path, required=True, help="analysis INI file"
    )
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        required=True,
        help="directory for catalogs and structure tables",
    )
    return parser


def _metadata(
    settings: config.AnalysisConfig,
    specs: Sequence[models.SimulationSpec],
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "input_files": [str(spec.path) for spec in specs],
        "effective_config": dataclasses.asdict(settings),
        "intentional_modernizations": [
            "strict 24-column parsing",
            "consistent one-based output IDs",
            "stable TSV schemas",
            "cumulative structure counts use n + 1 semantics",
        ],
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        settings = config.load_config(arguments.config)
        if arguments.manifest:
            specs = parsing.read_manifest(arguments.manifest)
        else:
            specs = tuple(
                models.SimulationSpec(path=path.resolve())
                for path in arguments.input
            )
        result = pipeline.run_analysis(specs, settings)
        output_dir = arguments.output_dir.resolve()
        catalog_dir = output_dir / "catalogs"
        structure_dir = output_dir / "structure"
        catalogs.write_snaps(catalog_dir / "snaps.tsv", result.simulations)
        catalogs.write_subavalanches(
            catalog_dir / "subavalanches.tsv", result.simulations
        )
        catalogs.write_avalanches(
            catalog_dir / "avalanches.tsv", result.simulations
        )
        structure.write_structure(
            structure_dir / "avalanche_structure.tsv",
            result.simulations,
            settings.structure,
            settings.window,
        )
        (output_dir / "run_metadata.json").write_text(
            json.dumps(_metadata(settings, specs), indent=2) + "\n"
        )
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(
        f"Analyzed {sum(len(item.events) for item in result.simulations)} snaps "
        f"from {len(result.simulations)} simulation(s)."
    )
    print(f"Outputs: {output_dir}")
    return 0
