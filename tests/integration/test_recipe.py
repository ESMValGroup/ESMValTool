from pprint import pformat
from textwrap import dedent

import pytest
import yaml
from netCDF4 import Dataset

import esmvaltool
from esmvaltool._recipe import read_recipe_file
from esmvaltool._task import DiagnosticTask
from esmvaltool.preprocessor import DEFAULT_ORDER, PreprocessingTask
from .test_diagnostic_run import write_config_user_file

MANDATORY_DATASET_KEYS = (
    'cmor_table',
    'dataset',
    'diagnostic',
    'end_year',
    'field',
    'filename',
    'frequency',
    'institute',
    'long_name',
    'mip',
    'modeling_realm',
    'preprocessor',
    'project',
    'short_name',
    'standard_name',
    'start_year',
    'units',
)

MANDATORY_SCRIPT_SETTINGS_KEYS = (
    'log_level',
    'script',
    'plot_dir',
    'run_dir',
    'work_dir',
)

DEFAULT_PREPROCESSOR_STEPS = (
    'cleanup',
    'cmor_check_data',
    'cmor_check_metadata',
    'concatenate',
    'extract_time',
    'fix_data',
    'fix_file',
    'fix_metadata',
    'load',
    'save',
)


@pytest.fixture
def config_user(tmpdir):
    filename = write_config_user_file(str(tmpdir))
    cfg = esmvaltool._config.read_config_user_file(filename, 'recipe_test')
    cfg['synda_download'] = False
    return cfg


def create_test_file(filename):
    attributes = {
        'tracking_id': 'xyz',
    }
    with Dataset(filename, 'w') as dataset:
        for key, value in attributes.items():
            setattr(dataset, key, value)


def test_recipe(tmpdir, config_user, monkeypatch):
    def find_files(_, filename):
        # Any occurrence of [something] in filename should have
        # been replaced when this function is called.
        assert '[' not in filename

        if filename.endswith('*'):
            filename = filename.rstrip('*') + '1950_2010.nc'
        filename = str(tmpdir / filename)
        create_test_file(filename)
        return [filename]

    monkeypatch.setattr(esmvaltool._data_finder, 'find_files', find_files)

    content = dedent("""
        datasets:
          - dataset: bcc-csm1-1

        preprocessors:
          preprocessor_name:
            extract_levels:
              levels: 85000
              scheme: nearest

        diagnostics:
          diagnostic_name:
            additional_datasets:
              - dataset: GFDL-ESM2G
            variables:
              ta:
                preprocessor: preprocessor_name
                field: T3M
                project: CMIP5
                mip: Amon
                exp: historical
                ensemble: r1i1p1
                start_year: 2000
                end_year: 2002
                additional_datasets:
                  - dataset: MPI-ESM-LR
            scripts:
              script_name:
                script: examples/diagnostic.py
                custom_setting: 1
        """)

    recipe_file = str(tmpdir / 'recipe_test.yml')
    with open(recipe_file, 'w') as file:
        file.write(content)

    recipe = read_recipe_file(str(recipe_file), config_user)
    raw = yaml.safe_load(content)
    # Perform some sanity checks on recipe expansion/normalization
    print("Expanded recipe:")
    assert len(recipe.diagnostics) == len(raw['diagnostics'])
    for diagnostic_name, diagnostic in recipe.diagnostics.items():
        print(pformat(diagnostic))
        source = raw['diagnostics'][diagnostic_name]

        # Check that 'variables' have been read and updated
        assert len(diagnostic['preprocessor_output']) == len(
            source['variables'])
        for variable_name, variables in diagnostic[
                'preprocessor_output'].items():
            assert len(variables) == 3
            for variable in variables:
                for key in MANDATORY_DATASET_KEYS:
                    assert key in variable and variable[key]
                assert variable_name == variable['short_name']

    # Check that the correct tasks have been created
    variables = recipe.diagnostics['diagnostic_name']['preprocessor_output'][
        'ta']
    tasks = {t for task in recipe.tasks for t in task.flatten()}
    preproc_tasks = {t for t in tasks if isinstance(t, PreprocessingTask)}
    diagnostic_tasks = {t for t in tasks if isinstance(t, DiagnosticTask)}

    assert len(preproc_tasks) == 1
    for task in preproc_tasks:
        print("Task", task.name)
        assert task.order == list(DEFAULT_ORDER)
        for product in task.products:
            variable = [
                v for v in variables if v['filename'] == product.filename
            ][0]
            assert product.attributes == variable
            for step in DEFAULT_PREPROCESSOR_STEPS:
                assert step in product.settings
            assert product.files

    assert len(diagnostic_tasks) == 1
    for task in diagnostic_tasks:
        print("Task", task.name)
        assert task.ancestors == list(preproc_tasks)
        assert task.script == 'examples/diagnostic.py'
        for key in MANDATORY_SCRIPT_SETTINGS_KEYS:
            assert key in task.settings and task.settings[key]
        assert task.settings['custom_setting'] == 1
