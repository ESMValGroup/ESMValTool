"""Tests for the configuration of the Recipe Test Workflow."""

import re
import shutil
from pathlib import Path

import esmvalcore.config
import esmvalcore.io
import pytest
from esmvalcore.config import Config

import esmvaltool

RTW_DIR = Path(esmvaltool.__file__).parent / "utils" / "recipe_test_workflow"
EXAMPLE_CONFIG_DIR = Path(esmvalcore.config.__file__).parent / "configurations"
SITES = sorted(path.name for path in (RTW_DIR / "site").iterdir())


@pytest.mark.parametrize("site", SITES)
def test_site_configuration(site, tmp_path, monkeypatch):
    """Check that the configuration for each site can be loaded."""
    # Create the configuration directory like the `configure` task does.
    for file in (RTW_DIR / "site" / site / "config").glob("*.yml"):
        shutil.copy(file, tmp_path)
    opt_file = RTW_DIR / "opt" / f"rose-suite-{site}.conf"
    data_sources_file = re.search(
        r'^DATA_SOURCES_FILE="(.*)"$',
        opt_file.read_text(encoding="utf-8"),
        flags=re.MULTILINE,
    )[1]
    shutil.copy(EXAMPLE_CONFIG_DIR / data_sources_file, tmp_path)

    monkeypatch.setenv("OUTPUT_DIR", str(tmp_path / "output"))
    cfg = Config()
    cfg.load_from_dirs([tmp_path])

    assert cfg["output_dir"] == tmp_path / "output"
    session = cfg.start_session("test")
    for project in ("CMIP5", "CMIP6", "OBS", "obs4MIPs"):
        assert esmvalcore.io.load_data_sources(session, project=project)
