"""Test diagnostic script runs."""
import contextlib
import os
import sys
from textwrap import dedent

import pytest
import yaml

from esmvaltool._main import run


def write_config_user_file(dirname):
    config_file = os.path.join(dirname, 'config-user.yml')
    cfg = {
        'output_dir': os.path.join(dirname, 'output_dir'),
        'rootpath': {
            'default': os.path.join(dirname, 'input_dir'),
        },
        'drs': {
            'CMIP5': 'BADC',
        },
        'log_level': 'debug',
    }
    with open(config_file, 'w') as file:
        yaml.safe_dump(cfg, file)
    return config_file


@contextlib.contextmanager
def arguments(*args):
    backup = sys.argv
    sys.argv = list(args)
    yield
    sys.argv = backup


def check(result_file):
    """Check the results."""
    with open(result_file, 'r') as file:
        result = yaml.safe_load(file)

    print(result)

    required_keys = {
        'input_files',
        'log_level',
        'plot_dir',
        'run_dir',
        'work_dir',
    }
    missing = required_keys - set(result)
    assert not missing


SCRIPTS = {
    'diagnostic.py':
    dedent("""
        import yaml
        from esmvaltool.diag_scripts.shared import run_diagnostic

        def main(cfg):
            with open(cfg['setting_name'], 'w') as file:
                yaml.safe_dump(cfg, file)

        if __name__ == '__main__':
            with run_diagnostic() as config:
                main(config)
        """),
    'diagnostic.ncl':
    dedent("""
        begin
            print("INFO    Loading settings from " + getenv("settings"))
            loadscript("$settings")
        end
        print("INFO Writing " + diag_script_info@setting_name)
        n = str_get_nl()
        result = "run_dir: " + config_user_info@run_dir + n +\
                 "work_dir: " + config_user_info@work_dir + n +\
                 "plot_dir: " + config_user_info@plot_dir + n +\
                 "log_level: " + config_user_info@log_level + n +\
                 "input_files: []" + n

        system("echo '" + result + "' > " + diag_script_info@setting_name)
        """),
    'diagnostic.R':
    dedent("""
        library(yaml)

        args <- commandArgs(trailingOnly = TRUE)
        print(paste0("INFO    Loading settings from ", args[1]))
        settings <- yaml::read_yaml(args[1])

        print(paste0("INFO    Writing settings to ", settings$setting_name))
        yaml::write_yaml(settings, settings$setting_name)
        """),
}


@pytest.mark.install
@pytest.mark.parametrize('script_file, script', SCRIPTS.items())
def test_diagnostic_run(tmpdir, script_file, script):

    recipe_file = os.path.join(tmpdir, 'recipe_test.yml')
    script_file = os.path.join(tmpdir, script_file)
    result_file = os.path.join(tmpdir, 'result.yml')

    # Write script to file
    with open(script_file, 'w') as file:
        file.write(script)

    # Create recipe
    recipe = dedent("""
        documentation:
          description: Recipe with no data.
          authors: [ande_bo]

        diagnostics:
          diagnostic_name:
            scripts:
              script_name:
                script: {}
                setting_name: {}
        """.format(script_file, result_file))
    with open(recipe_file, 'w') as file:
        file.write(recipe)

    config_user_file = write_config_user_file(tmpdir)
    with arguments('esmvaltool', '-c', config_user_file, recipe_file):
        run()

    check(result_file)
