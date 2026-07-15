"""Configuration models and INI loading."""

import configparser
import dataclasses
import pathlib


@dataclasses.dataclass(frozen=True, slots=True)
class SubavalancheDetectionConfig:
    average_bond_length: float = 0.0
    max_time_gap: float = 10.0
    max_distance: float = 100.0


@dataclasses.dataclass(frozen=True, slots=True)
class AvalancheDetectionConfig:
    mechanism: int = 0
    minimum_subavalanche_swaps: int = 1
    max_relative_time_gap: float = 0.01
    minimum_absolute_time_gap: float = 0.01
    max_distance: float = 100.0
    max_velocity: float = 20.0


@dataclasses.dataclass(frozen=True, slots=True)
class AnalysisWindow:
    minimum_time: float = 0.0
    maximum_time: float | None = None


@dataclasses.dataclass(frozen=True, slots=True)
class StructureConfig:
    minimum_swaps: int = 0
    linear_time_step: float = 1.0
    logarithmic_time_step: float = 1.0
    normalized_log_time_step: float = 0.1


@dataclasses.dataclass(frozen=True, slots=True)
class AnalysisConfig:
    subavalanches: SubavalancheDetectionConfig = dataclasses.field(
        default_factory=SubavalancheDetectionConfig
    )
    avalanches: AvalancheDetectionConfig = dataclasses.field(
        default_factory=AvalancheDetectionConfig
    )
    window: AnalysisWindow = dataclasses.field(default_factory=AnalysisWindow)
    structure: StructureConfig = dataclasses.field(default_factory=StructureConfig)


def _real(parser: configparser.ConfigParser, section: str, name: str, default: float) -> float:
    return parser.getfloat(section, name, fallback=default)


def _integer(parser: configparser.ConfigParser, section: str, name: str, default: int) -> int:
    return parser.getint(section, name, fallback=default)


def load_config(path: pathlib.Path | str) -> AnalysisConfig:
    """Load and validate the core analysis configuration."""

    config_path = pathlib.Path(path)
    parser = configparser.ConfigParser()
    if not parser.read(config_path):
        raise FileNotFoundError(f"Could not read configuration file: {config_path}")

    maximum_time = _real(parser, "Statistics", "analysis_time_max", -1.0)
    result = AnalysisConfig(
        subavalanches=SubavalancheDetectionConfig(
            average_bond_length=_real(
                parser, "SubavalancheDetection", "avgBondLength", 0.0
            ),
            max_time_gap=_real(
                parser,
                "SubavalancheDetection",
                "subavalanche_dt_threshold",
                10.0,
            ),
            max_distance=_real(
                parser,
                "SubavalancheDetection",
                "subavalanche_dr_threshold",
                100.0,
            ),
        ),
        avalanches=AvalancheDetectionConfig(
            mechanism=_integer(
                parser, "AvalancheDetection", "avalancheDetectionMechanism", 0
            ),
            minimum_subavalanche_swaps=_integer(
                parser,
                "AvalancheDetection",
                "avDetection_subavNswaps_threshold",
                1,
            ),
            max_relative_time_gap=_real(
                parser, "AvalancheDetection", "avDetection_dtt_threshold", 0.01
            ),
            minimum_absolute_time_gap=_real(
                parser, "AvalancheDetection", "avDetection_dt_minimum", 0.01
            ),
            max_distance=_real(
                parser, "AvalancheDetection", "avDetection_dr_threshold", 100.0
            ),
            max_velocity=_real(
                parser,
                "AvalancheDetection",
                "avDetection_velocity_threshold",
                20.0,
            ),
        ),
        window=AnalysisWindow(
            minimum_time=_real(parser, "Statistics", "analysis_time_min", 0.0),
            maximum_time=maximum_time if maximum_time > 0 else None,
        ),
        structure=StructureConfig(
            minimum_swaps=_integer(
                parser, "AvalancheStat", "avStructure_nswaps_threshold", 0
            ),
            linear_time_step=_real(
                parser, "AvalancheStat", "avStructure_dt", 1.0
            ),
            logarithmic_time_step=_real(
                parser, "AvalancheStat", "avStructure_dlogt", 1.0
            ),
            normalized_log_time_step=_real(
                parser, "AvalancheStat", "avStructure_dlogtnorm", 0.1
            ),
        ),
    )
    validate_config(result)
    return result


def validate_config(value: AnalysisConfig) -> None:
    """Raise ``ValueError`` for configurations that cannot be analyzed safely."""

    if value.avalanches.mechanism not in range(6):
        raise ValueError("avalanche detection mechanism must be between 0 and 5")
    if value.subavalanches.max_time_gap < 0 or value.subavalanches.max_distance < 0:
        raise ValueError("subavalanche thresholds cannot be negative")
    if value.avalanches.max_relative_time_gap < 0:
        raise ValueError("relative avalanche time threshold cannot be negative")
    if value.avalanches.minimum_absolute_time_gap < 0:
        raise ValueError("absolute avalanche time threshold cannot be negative")
    if value.avalanches.max_distance < 0 or value.avalanches.max_velocity < 0:
        raise ValueError("avalanche distance and velocity thresholds cannot be negative")
    if value.avalanches.minimum_subavalanche_swaps < 1:
        raise ValueError("minimum subavalanche size must be at least one")
    if value.window.maximum_time is not None:
        if value.window.maximum_time < value.window.minimum_time:
            raise ValueError("analysis maximum time precedes minimum time")
    steps = (
        value.structure.linear_time_step,
        value.structure.logarithmic_time_step,
        value.structure.normalized_log_time_step,
    )
    if any(step <= 0 for step in steps):
        raise ValueError("all avalanche structure sampling steps must be positive")
