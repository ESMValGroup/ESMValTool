"""Test :mod:`esmvaltool.diag_scripts.monitor.multi_datasets`."""

from pathlib import Path

import matplotlib
import pytest

from esmvaltool.diag_scripts.monitor.multi_datasets import main
from tests.sample_data.diag_scripts import (
    MATPLOTLIB_RC,
    assert_phash,
    get_cfg,
    get_phash,
    load_test_setups,
    write_imagehashes,
)

TEST_SETUPS = load_test_setups("monitor/multi_datasets_setups.yml")
DEFAULT_SETTINGS = {
    "figure_kwargs": {},
    "plot_filename": "{plot_type}_{real_name}_{dataset}_{mip}",
    "plot_folder": "{plot_dir}",
    "savefig_kwargs": {},
}


@pytest.mark.diagnostic_image_output
@pytest.mark.parametrize(
    ("input_data", "settings", "expected_pngs"), TEST_SETUPS
)
def test_diagnostic_image_output(
    pytestconfig,
    tmp_path: Path,
    input_data: list[str],
    settings: dict,
    expected_pngs: list[str],
) -> None:
    """Test if diagnostic image output matches expected output."""
    save_imagehashes = pytestconfig.getoption("save_imagehashes")
    show_images = pytestconfig.getoption("show_images")

    cfg_settings = {**DEFAULT_SETTINGS, **settings}
    cfg = get_cfg(tmp_path, input_data, **cfg_settings)

    # Use fixed matplotlib settings to ensure consistent test results
    with matplotlib.rc_context(fname=MATPLOTLIB_RC):
        main(cfg)

    actual_pngs: list[Path] = []
    imagehashes: dict[str, str] = {}
    for png in expected_pngs:
        image_key = f"monitor.multi_datasets.{png}"
        actual_png = tmp_path / "output" / "plots" / png
        assert actual_png.is_file(), f"{actual_png} does not exist"

        # Skip actual comparison if special options are set
        if save_imagehashes is None and show_images is None:
            assert_phash(image_key, actual_png)

        actual_pngs.append(actual_png)
        imagehashes[image_key] = str(get_phash(actual_png))

    # Save imagehashes if desired
    if save_imagehashes is not None:
        write_imagehashes(imagehashes, Path(save_imagehashes))

    # Show image paths if desired
    if show_images is not None:
        msg = f"Produced images: {','.join(str(p) for p in actual_pngs)}"
        pytest.fail(msg)
