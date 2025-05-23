"""Test diagnostic script runs."""

import contextlib
import shutil
import sys
from pathlib import Path
from textwrap import dedent

import esmvalcore
import pytest
import yaml
from esmvalcore._main import run
from packaging import version


def write_config_file(dirname):
    config_file = dirname / "config-user.yml"
    cfg = {
        "output_dir": str(dirname / "output_dir"),
        "rootpath": {
            "default": str(dirname / "input_dir"),
        },
        "drs": {
            "CMIP5": "BADC",
        },
        "log_level": "debug",
    }
    config_file.write_text(yaml.safe_dump(cfg, encoding=None))
    return str(config_file)


@contextlib.contextmanager
def arguments(*args):
    backup = sys.argv
    sys.argv = list(args)
    yield
    sys.argv = backup


def check(result_file):
    """Check the results."""
    result = yaml.safe_load(result_file.read_text())

    print(result)

    required_keys = {
        "input_files",
        "log_level",
        "plot_dir",
        "run_dir",
        "work_dir",
    }
    missing = required_keys - set(result)
    assert not missing


SCRIPTS = [
    "diagnostic.py",
    pytest.param(
        "diagnostic.ncl",
        marks=pytest.mark.skipif(
            sys.platform == "darwin",
            reason="ESMValTool ncl not supported on OSX",
        ),
    ),
    pytest.param(
        "diagnostic.R",
        marks=pytest.mark.skipif(
            sys.platform == "darwin",
            reason="ESMValTool R not supported on OSX",
        ),
    ),
    pytest.param(
        "diagnostic.jl",
        marks=pytest.mark.skip(
            reason="Julia is not installed by default.",
        ),
    ),
]


@pytest.mark.skipif(
    version.parse(esmvalcore.__version__) >= version.parse("2.14.0"),
    reason="ESMValCore >= v2.14.0",
)
@pytest.mark.installation
@pytest.mark.parametrize("script_file", SCRIPTS)
def test_diagnostic_run_config_file(tmp_path, script_file):
    local_script_file = Path(__file__).parent / script_file

    recipe_file = tmp_path / "recipe_test.yml"
    script_file = tmp_path / script_file
    result_file = tmp_path / "result.yml"

    shutil.copy(local_script_file, script_file)

    # Create recipe
    recipe = dedent(f"""
        documentation:
          title: Test recipe
          description: Recipe with no data.
          authors: [andela_bouwe]

        diagnostics:
          diagnostic_name:
            scripts:
              script_name:
                script: {script_file}
                setting_name: {result_file}
        """)
    recipe_file.write_text(str(recipe))

    config_file = write_config_file(tmp_path)
    with arguments(
        "esmvaltool",
        "run",
        "--config_file",
        config_file,
        str(recipe_file),
    ):
        run()

    check(result_file)


@pytest.mark.skipif(
    version.parse(esmvalcore.__version__) < version.parse("2.12.0"),
    reason="ESMValCore < v2.12.0",
)
@pytest.mark.installation
@pytest.mark.parametrize("script_file", SCRIPTS)
def test_diagnostic_run(tmp_path, script_file):
    local_script_file = Path(__file__).parent / script_file

    recipe_file = tmp_path / "recipe_test.yml"
    script_file = tmp_path / script_file
    result_file = tmp_path / "result.yml"
    config_dir = tmp_path / "config"
    config_dir.mkdir(exist_ok=True, parents=True)

    shutil.copy(local_script_file, script_file)

    # Create recipe
    recipe = dedent(f"""
        documentation:
          title: Test recipe
          description: Recipe with no data.
          authors: [andela_bouwe]

        diagnostics:
          diagnostic_name:
            scripts:
              script_name:
                script: {script_file}
                setting_name: {result_file}
        """)
    recipe_file.write_text(str(recipe))

    write_config_file(config_dir)
    with arguments(
        "esmvaltool",
        "run",
        "--config_dir",
        str(config_dir),
        str(recipe_file),
    ):
        run()

    check(result_file)
