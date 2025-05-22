#!/usr/bin/env python
import tempfile
from pathlib import Path

import pytest


def save_imagehashes() -> None:
    """Write imagehashes."""
    temp_file_kwargs = {
        "prefix": "imagehashes_",
        "suffix": ".yml",
        "dir": Path(__file__).parent,
        "delete": False,
        "delete_on_close": False,
    }
    with tempfile.NamedTemporaryFile(**temp_file_kwargs) as file:
        imagehashes_path = file.name

    pytest_args = [
        "-n0",
        "-m",
        "diagnostic_image_output",
        "--save_imagehashes",
        imagehashes_path,
        str(Path(__file__).parent),
    ]
    pytest.main(pytest_args)

    print("Saved imagehashes to", imagehashes_path)


if __name__ == "__main__":
    save_imagehashes()
