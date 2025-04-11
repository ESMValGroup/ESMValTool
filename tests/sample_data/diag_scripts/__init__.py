"""Test diagnostics scripts with sample data."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import yaml

SAMPLE_DATA_DIR = Path(__file__).resolve().parent / "data"
with (SAMPLE_DATA_DIR / "metadata.yml").open("r", encoding="utf-8") as file:
    METADATA = yaml.safe_load(file)


def get_cfg(tmp_dir: Path, input_data: Iterable[str], **kwargs: str) -> dict:
    """Get diagnostic configuration."""
    # Input directories
    in_dir = tmp_dir / "input"
    aux_dir = in_dir / "aux"

    # Output directories
    out_dir = tmp_dir / "output"
    plot_dir = out_dir / "plots"
    run_dir = out_dir / "run"
    work_dir = out_dir / "work"

    for dir_ in (aux_dir, plot_dir, run_dir, work_dir):
        dir_.mkdir(parents=True, exist_ok=True)

    # Setup input data
    metadata: dict[str, dict] = {}
    for filename in input_data:
        filepath = str(SAMPLE_DATA_DIR / filename)
        assert filename in METADATA
        metadata[filepath] = METADATA[filename]
        metadata[filepath]["filename"] = filepath
    metadata_file = in_dir / "metadata.yml"
    with metadata_file.open("w", encoding="utf-8") as file:
        yaml.safe_dump(metadata, file)

    cfg = {
        "auxiliary_data_dir": str(aux_dir),
        "input_data": metadata,
        "input_files": [str(metadata_file)],
        "log_level": "info",
        "output_file_type": "png",
        "plot_dir": str(plot_dir),
        "recipe": "recipe.yml",
        "run_dir": str(run_dir),
        "work_dir": str(work_dir),
        **kwargs,
    }

    return cfg


def load_settings(path: str | Path) -> list[tuple]:
    """Load list of diagnostics settings to parameterize tests."""
    path = Path(path)

    # Try relative path if absolute one does not exist
    if not path.is_file():
        path = Path(__file__).resolve().parent / path

    with path.open("r", encoding="utf-8") as file:
        settings = yaml.safe_load(file)

    parametrize_input: list[tuple] = [
        (s.pop("input_data"), s, s.pop("expected_pngs")) for s in settings
    ]

    return parametrize_input
