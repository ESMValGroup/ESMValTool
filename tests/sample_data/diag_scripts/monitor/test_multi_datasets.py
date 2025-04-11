"""Test :mod:`esmvaltool.diag_scripts.monitor.multi_datasets`."""

from pathlib import Path

import pytest

from esmvaltool.diag_scripts.monitor.multi_datasets import main
from esmvaltool.utils.testing.regression.compare import compare_png
from tests.sample_data.diag_scripts import get_cfg, load_settings

DIAG_SETTINGS = load_settings("monitor/multi_datasets_settings.yml")


@pytest.mark.parametrize("input_data,settings,expected_pngs", DIAG_SETTINGS)
def test_diagnostic(
    input_data: list[str],
    settings: dict,
    expected_pngs: list[str],
) -> None:
    """Test diagnostic with various settings."""
    path = Path.home() / "tmp" / "aaaaaaaaaaaaa"
    cfg = get_cfg(path, input_data, **settings)

    main(cfg)

    for png in expected_pngs:
        actual_png = path / "output" / "plots" / png
        expected_png = Path(__file__).resolve().parent / "output" / png
        assert actual_png.is_file()
        assert expected_png.is_file()
        compare_png(expected_png, actual_png)
