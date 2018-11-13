import os
from pprint import pformat
from textwrap import dedent
from unittest.mock import create_autospec

import pytest
import yaml
from netCDF4 import Dataset

import esmvaltool
from esmvaltool._recipe import TASKSEP, read_recipe_file
from esmvaltool._task import DiagnosticTask
from esmvaltool.preprocessor import DEFAULT_ORDER, PreprocessingTask
from esmvaltool.preprocessor._io import concatenate_callback

from .test_diagnostic_run import write_config_user_file
from .test_provenance import check_provenance

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
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    attributes = {
        'tracking_id': 'xyz',
    }

    with Dataset(filename, 'w') as dataset:
        for key, value in attributes.items():
            setattr(dataset, key, value)


@pytest.fixture
def patched_datafinder(tmpdir, monkeypatch):
    def find_files(_, filename):
        # Any occurrence of [something] in filename should have
        # been replaced before this function is called.
        assert '[' not in filename

        filename = str(tmpdir / 'input' / filename)
        filenames = []
        if filename.endswith('*'):
            filename = filename.rstrip('*')
            intervals = [
                '1990_1999',
                '2000_2009',
                '2010_2019',
            ]
            for interval in intervals:
                filenames.append(filename + interval + '.nc')
        else:
            filenames.append(filename)

        for file in filenames:
            create_test_file(file)
        return filenames

    monkeypatch.setattr(esmvaltool._data_finder, 'find_files', find_files)


def test_simple_recipe(tmpdir, patched_datafinder, config_user):

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
                start_year: 1999
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

    recipe = read_recipe_file(recipe_file, config_user)
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
            assert len(product.files) == 2
            check_provenance(product)

    assert len(diagnostic_tasks) == 1
    for task in diagnostic_tasks:
        print("Task", task.name)
        assert task.ancestors == list(preproc_tasks)
        assert task.script == 'examples/diagnostic.py'
        for key in MANDATORY_SCRIPT_SETTINGS_KEYS:
            assert key in task.settings and task.settings[key]
        assert task.settings['custom_setting'] == 1


def test_default_preprocessor(tmpdir, patched_datafinder, config_user):

    content = dedent("""
        diagnostics:
          diagnostic_name:
            variables:
              chl:
                project: CMIP5
                mip: Oyr
                exp: historical
                start_year: 2000
                end_year: 2005
                field: TO3Y
                ensemble: r1i1p1
                additional_datasets:
                  - {dataset: CanESM2}
            scripts: null
        """)

    recipe_file = str(tmpdir / 'recipe_test.yml')
    with open(recipe_file, 'w') as file:
        file.write(content)

    recipe = read_recipe_file(recipe_file, config_user)
    assert len(recipe.tasks) == 1
    task = recipe.tasks.pop()
    assert len(task.products) == 1
    product = task.products.pop()
    preproc_dir = os.path.dirname(product.filename)
    assert preproc_dir.startswith(str(tmpdir))

    fix_dir = os.path.join(
        preproc_dir,
        'CMIP5_CanESM2_Oyr_historical_r1i1p1_TO3Y_chl_2000-2005_fixed')
    defaults = {
        'load': {
            'callback':
            concatenate_callback,
            'constraints': ('mass_concentration_of_phytoplankton_expressed_'
                            'as_chlorophyll_in_sea_water'),
        },
        'concatenate': {},
        'fix_file': {
            'project': 'CMIP5',
            'dataset': 'CanESM2',
            'short_name': 'chl',
            'output_dir': fix_dir,
        },
        'fix_data': {
            'project': 'CMIP5',
            'dataset': 'CanESM2',
            'short_name': 'chl',
            'cmor_table': 'CMIP5',
            'mip': 'Oyr',
        },
        'fix_metadata': {
            'project': 'CMIP5',
            'dataset': 'CanESM2',
            'short_name': 'chl',
            'cmor_table': 'CMIP5',
            'mip': 'Oyr',
        },
        'extract_time': {
            'start_year': 2000,
            'end_year': 2006,
            'start_month': 1,
            'end_month': 1,
            'start_day': 1,
            'end_day': 1,
        },
        'cmor_check_metadata': {
            'cmor_table': 'CMIP5',
            'mip': 'Oyr',
            'short_name': 'chl',
        },
        'cmor_check_data': {
            'cmor_table': 'CMIP5',
            'mip': 'Oyr',
            'short_name': 'chl',
        },
        'cleanup': {
            'remove': [fix_dir]
        },
        'save': {
            'compress': False,
            'filename': product.filename,
        }
    }
    assert product.settings == defaults


def test_reference_dataset(tmpdir, patched_datafinder, config_user):

    levels = create_autospec(
        esmvaltool._recipe.get_reference_levels, return_value=[100])
    esmvaltool._recipe.get_reference_levels = levels

    content = dedent("""
        preprocessors:
          test_from_reference:
            regrid:
              target_grid: reference_dataset
              scheme: linear
            extract_levels:
              levels: reference_dataset
              scheme: linear
          test_from_cmor_table:
            extract_levels:
              levels:
                cmor_table: CMIP6
                coordinate: alt16
              scheme: nearest

        diagnostics:
          diagnostic_name:
            variables:
              ta: &var
                preprocessor: test_from_reference
                project: CMIP5
                mip: Amon
                exp: historical
                start_year: 2000
                end_year: 2005
                field: T3M
                ensemble: r1i1p1
                additional_datasets:
                  - {dataset: GFDL-CM3}
                  - {dataset: MPI-ESM-LR}
                reference_dataset: MPI-ESM-LR
              ch4:
                <<: *var
                preprocessor: test_from_cmor_table
                additional_datasets:
                  - {dataset: GFDL-CM3}

            scripts: null
        """)

    recipe_file = str(tmpdir / 'recipe_test.yml')
    with open(recipe_file, 'w') as file:
        file.write(content)

    recipe = read_recipe_file(recipe_file, config_user)

    assert len(recipe.tasks) == 2

    # Check that the reference dataset has been used
    task = next(t for t in recipe.tasks
                if t.name == 'diagnostic_name' + TASKSEP + 'ta')
    assert len(task.products) == 2
    product = next(
        p for p in task.products if p.attributes['dataset'] == 'GFDL-CM3')
    reference = next(
        p for p in task.products if p.attributes['dataset'] == 'MPI-ESM-LR')

    assert product.settings['regrid']['target_grid'] == reference.files[0]
    assert product.settings['extract_levels']['levels'] == levels.return_value

    fix_dir = os.path.splitext(reference.filename)[0] + '_fixed'
    levels.assert_called_once_with(reference.files[0], 'CMIP5', 'MPI-ESM-LR',
                                   'ta', fix_dir, 'air_pressure')

    assert 'regrid' not in reference.settings
    assert 'extract_levels' not in reference.settings

    # Check that levels have been read from CMOR table
    task = next(t for t in recipe.tasks
                if t.name == 'diagnostic_name' + TASKSEP + 'ch4')
    assert len(task.products) == 1
    product = next(iter(task.products))
    assert product.settings['extract_levels']['levels'] == [
        0,
        250,
        750,
        1250,
        1750,
        2250,
        2750,
        3500,
        4500,
        6000,
        8000,
        10000,
        12000,
        14500,
        16000,
        18000,
    ]


def test_custom_preproc_order(tmpdir, patched_datafinder, config_user):
    pass


def test_derive(tmpdir, patched_datafinder, config_user):

    content = dedent("""
        diagnostics:
          diagnostic_name:
            variables:
              toz:
                project: CMIP5
                mip: Amon
                exp: historical
                start_year: 2000
                end_year: 2005
                field: T2Ms
                derive: true
                force_derivation: true
                additional_datasets:
                  - {dataset: GFDL-CM3,  ensemble: r1i1p1}
            scripts: null
        """)

    recipe_file = str(tmpdir / 'recipe_test.yml')
    with open(recipe_file, 'w') as file:
        file.write(content)

    recipe = read_recipe_file(recipe_file, config_user)

    # Check generated tasks
    assert len(recipe.tasks) == 1
    task = recipe.tasks.pop()

    assert task.name == 'diagnostic_name' + TASKSEP + 'toz'
    assert len(task.ancestors) == 2
    assert 'diagnostic_name' + TASKSEP + 'toz_derive_input_ps' in [
        t.name for t in task.ancestors
    ]
    assert 'diagnostic_name' + TASKSEP + 'toz_derive_input_tro3' in [
        t.name for t in task.ancestors
    ]

    # Check product content of tasks
    assert len(task.products) == 1
    product = task.products.pop()
    assert 'derive' in product.settings
    assert product.attributes['short_name'] == 'toz'
    check_provenance(product)
    assert product.files

    ps_product = next(p for a in task.ancestors for p in a.products
                      if p.attributes['short_name'] == 'ps')
    tro3_product = next(p for a in task.ancestors for p in a.products
                        if p.attributes['short_name'] == 'tro3')
    assert ps_product.filename in product.files
    assert tro3_product.filename in product.files
    check_provenance(ps_product)
    check_provenance(tro3_product)


def test_derive_not_needed(tmpdir, patched_datafinder, config_user):

    content = dedent("""
        diagnostics:
          diagnostic_name:
            variables:
              toz:
                project: CMIP5
                mip: Amon
                exp: historical
                start_year: 2000
                end_year: 2005
                field: T2Ms
                derive: true
                force_derivation: false
                additional_datasets:
                  - {dataset: GFDL-CM3,  ensemble: r1i1p1}
            scripts: null
        """)

    recipe_file = str(tmpdir / 'recipe_test.yml')
    with open(recipe_file, 'w') as file:
        file.write(content)

    recipe = read_recipe_file(recipe_file, config_user)

    # Check generated tasks
    assert len(recipe.tasks) == 1
    task = recipe.tasks.pop()

    assert task.name == 'diagnostic_name/toz'
    assert len(task.ancestors) == 1
    ancestor = [t for t in task.ancestors][0]
    assert ancestor.name == 'diagnostic_name/toz_derive_input_toz'

    # Check product content of tasks
    assert len(task.products) == 1
    product = task.products.pop()
    assert 'derive' in product.settings
    assert product.attributes['short_name'] == 'toz'
    check_provenance(product)

    assert len(ancestor.products) == 1
    ancestor_product = ancestor.products.pop()
    assert ancestor_product.filename in product.files
    assert ancestor_product.attributes['short_name'] == 'toz'
    check_provenance(ancestor_product)
