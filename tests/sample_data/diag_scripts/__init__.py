"""Test diagnostics scripts with sample data."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
import yaml
from imagehash import hex_to_hash, phash
from PIL import Image

SAMPLE_DATA_DIR = Path(__file__).resolve().parent / "data"
with (SAMPLE_DATA_DIR / "metadata.yml").open("r", encoding="utf-8") as file:
    METADATA = yaml.safe_load(file)

IMAGEHASHES_PATH = Path(__file__).resolve().parent / "imagehashes.yml"
with IMAGEHASHES_PATH.open("r", encoding="utf-8") as file:
    IMAGEHASHES = yaml.safe_load(file)
HASH_SIZE = 16
MAX_PHASH_DISTANCE = 2


def assert_phash(image_key: str, actual_phash: np.ndarray) -> None:
    """Compare imagehashes."""
    assert_msg = (
        f"No expected output for image '{image_key}' in {IMAGEHASHES_PATH}"
    )
    assert image_key in IMAGEHASHES, assert_msg

    expected_phash = hex_to_hash(IMAGEHASHES[image_key])
    distance = expected_phash - actual_phash

    assert_msg = f"Image '{image_key}' changed (Hamming distance: {distance})"
    assert distance < MAX_PHASH_DISTANCE, assert_msg


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
    for file_idx, filename in enumerate(input_data):
        filepath = str(SAMPLE_DATA_DIR / filename)
        assert filename in METADATA, f"{filename} not in metadata.yml file"
        metadata[filepath] = METADATA[filename]
        metadata[filepath]["filename"] = filepath
        metadata[filepath]["recipe_dataset_index"] = file_idx
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


def get_phash(image_path: Path) -> str:
    """Get phash of image."""
    with Image.open(image_path) as img:
        return phash(img, hash_size=HASH_SIZE)


def load_test_setups(path: str | Path) -> list[tuple]:
    """Load test setups (used as input to :func:`pytest.mark.parametrize`)."""
    path = Path(path)

    # Try relative path if absolute one does not exist
    if not path.is_file():
        path = Path(__file__).resolve().parent / path

    with path.open("r", encoding="utf-8") as file:
        setups = yaml.safe_load(file)

    parametrize_input: list[tuple] = [
        (s["input_data"], s["settings"], s["expected_pngs"]) for s in setups
    ]

    return parametrize_input


def write_imagehashes(
    imagehashes: dict[str, str],
    imagehashes_path: Path,
) -> None:
    with imagehashes_path.open("a", encoding="utf-8") as file:
        yaml.safe_dump(imagehashes, file)
