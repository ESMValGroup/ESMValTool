"""Test :mod:`esmvaltool.diag_scripts.monitor.multi_datasets`."""

from pathlib import Path

import pytest

from esmvaltool.diag_scripts.monitor.multi_datasets import main
from esmvaltool.utils.testing.regression.compare import compare_png
from tests.sample_data.diag_scripts import get_cfg, load_test_setups

TEST_SETUPS = load_test_setups("monitor/multi_datasets_setups.yml")
DEFAULT_SETTINGS = {
    "figure_kwargs": {"figsize": [5, 5]},
    "plot_filename": "{plot_type}_{real_name}_{dataset}_{mip}",
    "plot_folder": "{plot_dir}",
    "savefig_kwargs": {"dpi": 30, "orientation": "landscape"},
}


@pytest.mark.parametrize("input_data,settings,expected_pngs", TEST_SETUPS)
def test_diagnostic(
    tmp_path: Path,
    input_data: list[str],
    settings: dict,
    expected_pngs: list[str],
) -> None:
    """Test diagnostic with various setups."""
    # tmp_path = Path("/home/b/b309141/tmp/aaaaaaaaaaaaa")
    cfg_settings = {**DEFAULT_SETTINGS, **settings}
    cfg = get_cfg(tmp_path, input_data, **cfg_settings)

    # Explicitly set backend (otherwise, this might lead to different image
    # sizes on different machines, leading to failing tests)
    cfg.setdefault("matplotlib_rc_params", {})
    cfg["matplotlib_rc_params"]["backend"] = "agg"

    main(cfg)

    for png in expected_pngs:
        actual_png = tmp_path / "output" / "plots" / png
        expected_png = (
            Path(__file__).resolve().parent / "expected_output" / png
        )
        assert actual_png.is_file()
        assert expected_png.is_file()
        assert compare_png(expected_png, actual_png)
