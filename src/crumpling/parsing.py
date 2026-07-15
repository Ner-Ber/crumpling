"""Input parsers for manifests and full-format swap tables."""

import csv
import pathlib
from typing import Sequence

from crumpling import models


EXPECTED_SWAP_COLUMNS = 24


class InputFormatError(ValueError):
    """Raised when an analysis input violates its documented schema."""


def parse_swap_row(fields: Sequence[str], line_number: int) -> models.SwapEvent:
    """Parse one 24-column full-format swap row."""

    if len(fields) != EXPECTED_SWAP_COLUMNS:
        raise InputFormatError(
            f"line {line_number}: expected {EXPECTED_SWAP_COLUMNS} tab-separated "
            f"columns, found {len(fields)}"
        )
    try:
        return models.SwapEvent(
            source_number=int(fields[0]),
            timestep=float(fields[2]),
            time=float(fields[3]),
            direction=int(fields[6]),
            x_size=float(fields[7]),
            y_size=float(fields[8]),
            surface_area=float(fields[9]),
            position=(float(fields[17]), float(fields[18]), float(fields[19])),
            atom1=int(fields[22]),
            atom2=int(fields[23]),
        )
    except ValueError as exc:
        raise InputFormatError(f"line {line_number}: invalid numeric value") from exc


def read_swap_file(path: pathlib.Path | str) -> tuple[models.SwapEvent, ...]:
    """Read a full-format swap table and validate chronological ordering."""

    source = pathlib.Path(path)
    events: list[models.SwapEvent] = []
    with source.open(newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        next(reader, None)
        for line_number, fields in enumerate(reader, start=2):
            events.append(parse_swap_row(fields, line_number))
    if not events:
        raise InputFormatError(f"{source}: no swap records")
    for previous, current in zip(events, events[1:]):
        if current.time < previous.time:
            raise InputFormatError(
                f"{source}: events are not ordered by time "
                f"({previous.time} followed by {current.time})"
            )
    if any(event.time <= 0 for event in events):
        raise InputFormatError(f"{source}: event times must be positive")
    return tuple(events)


def read_manifest(path: pathlib.Path | str) -> tuple[models.SimulationSpec, ...]:
    """Read the legacy whitespace manifest, resolving inputs beside the manifest."""

    manifest_path = pathlib.Path(path).resolve()
    tokens = manifest_path.read_text().split()
    if len(tokens) < 2:
        raise InputFormatError(f"{manifest_path}: incomplete manifest header")
    try:
        count = int(tokens[1])
    except ValueError as exc:
        raise InputFormatError(f"{manifest_path}: invalid simulation count") from exc
    expected = 2 + 2 * count
    if len(tokens) != expected:
        raise InputFormatError(
            f"{manifest_path}: expected {expected} tokens, found {len(tokens)}"
        )
    specs: list[models.SimulationSpec] = []
    for index in range(count):
        offset = 2 + 2 * index
        input_path = pathlib.Path(tokens[offset])
        if not input_path.is_absolute():
            input_path = manifest_path.parent / input_path
        try:
            maximum_time = float(tokens[offset + 1])
        except ValueError as exc:
            raise InputFormatError(
                f"{manifest_path}: invalid maximum time for {input_path}"
            ) from exc
        specs.append(
            models.SimulationSpec(
                path=input_path.resolve(), max_simulation_time=maximum_time
            )
        )
    return tuple(specs)
