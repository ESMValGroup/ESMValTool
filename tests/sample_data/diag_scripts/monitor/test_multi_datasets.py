"""Test :mod:`esmvaltool.diag_scripts.monitor.multi_datasets`."""

from pathlib import Path

import pytest

from esmvaltool.diag_scripts.monitor.multi_datasets import main
from tests.sample_data.diag_scripts import (
    assert_phash,
    get_cfg,
    get_phash,
    load_test_setups,
    write_imagehashes,
)

TEST_SETUPS = load_test_setups("monitor/multi_datasets_setups.yml")


@pytest.mark.diagnostic_image_output
@pytest.mark.parametrize("input_data,settings,expected_pngs", TEST_SETUPS)
def test_diagnostic_image_output(
    pytestconfig,
    tmp_path: Path,
    input_data: list[str],
    settings: dict,
    expected_pngs: list[str],
) -> None:
    """Test if diagnostic image output matches expected output."""
    # tmp_path = Path.home() / "aaa"
    # tmp_path.mkdir(exist_ok=True)
    save_imagehashes = pytestconfig.getoption("save_imagehashes")

    cfg = get_cfg(tmp_path, input_data, **settings)

    # Force some settings to avoid test fails on different machines
    cfg["plot_filename"] = "{plot_type}_{real_name}_{dataset}_{mip}"
    cfg["plot_folder"] = "{plot_dir}"
    cfg.setdefault("figure_kwargs", {})
    cfg.setdefault("matplotlib_rc_params", {})
    cfg.setdefault("savefig_kwargs", {})
    cfg["figure_kwargs"]["figsize"] = [8.0, 6.0]
    cfg["matplotlib_rc_params"]["backend"] = "agg"
    cfg["savefig_kwargs"]["dpi"] = 100

    main(cfg)

    imagehashes: dict[str, str] = {}
    for png in expected_pngs:
        image_key = f"monitor.multi_datasets.{png}"
        actual_png = tmp_path / "output" / "plots" / png
        assert actual_png.is_file()

        # Skip actual comparison if imagehashes are written
        if save_imagehashes is None:
            assert_phash(image_key, actual_png)

        imagehashes[image_key] = str(get_phash(actual_png))

    # Save imagehashes if desired
    if save_imagehashes is not None:
        write_imagehashes(imagehashes, Path(save_imagehashes))
