import os
from pprint import pformat
from textwrap import dedent

import iris
import pytest
import yaml
from mock import create_autospec
from six import text_type

import esmvaltool
from esmvaltool._recipe import TASKSEP, read_recipe_file
from esmvaltool._task import DiagnosticTask
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger, get_diagnostic_filename, get_plot_filename)
from esmvaltool.preprocessor import DEFAULT_ORDER, PreprocessingTask
from esmvaltool.preprocessor._io import concatenate_callback

from .test_diagnostic_run import write_config_user_file
from .test_provenance import check_provenance

MANDATORY_DATASET_KEYS = (
    'cmor_table',
    'dataset',
    'diagnostic',
    'end_year',
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
def config_user(tmp_path):
    filename = write_config_user_file(tmp_path)
    cfg = esmvaltool._config.read_config_user_file(filename, 'recipe_test')
    cfg['synda_download'] = False
    return cfg


def create_test_file(filename, tracking_id=None):
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    attributes = {}
    if tracking_id is not None:
        attributes['tracking_id'] = tracking_id
    cube = iris.cube.Cube([], attributes=attributes)

    iris.save(cube, filename)


@pytest.fixture
def patched_datafinder(tmp_path, monkeypatch):
    def tracking_ids(i=0):
        while True:
            yield i
            i += 1

    tracking_id = tracking_ids()

    def find_files(_, filenames):
        # Any occurrence of [something] in filename should have
        # been replaced before this function is called.
        for filename in filenames:
            assert '[' not in filename

        filename = filenames[0]
        filename = str(tmp_path / 'input' / filename)
        filenames = []
        if filename.endswith('*.nc'):
            filename = filename[:-len('*.nc')]
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
            create_test_file(file, next(tracking_id))
        return filenames

    monkeypatch.setattr(esmvaltool._data_finder, 'find_files', find_files)


DEFAULT_DOCUMENTATION = dedent("""
    documentation:
      description: This is a test recipe.
      authors:
        - ande_bo
      references:
        - contact_authors
        - acknow_project
      projects:
        - c3s-magic
    """)


def get_recipe(tempdir, content, cfg):
    """Save and load recipe content."""
    recipe_file = tempdir / 'recipe_test.yml'
    # Add mandatory documentation section
    content = text_type(DEFAULT_DOCUMENTATION + content)
    recipe_file.write_text(content)

    recipe = read_recipe_file(str(recipe_file), cfg)

    return recipe


def test_simple_recipe(tmp_path, patched_datafinder, config_user):

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

    recipe = get_recipe(tmp_path, content, config_user)
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

    assert len(diagnostic_tasks) == 1
    for task in diagnostic_tasks:
        print("Task", task.name)
        assert task.ancestors == list(preproc_tasks)
        assert task.script == 'examples/diagnostic.py'
        for key in MANDATORY_SCRIPT_SETTINGS_KEYS:
            assert key in task.settings and task.settings[key]
        assert task.settings['custom_setting'] == 1


def test_default_preprocessor(tmp_path, patched_datafinder, config_user):

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
                ensemble: r1i1p1
                additional_datasets:
                  - {dataset: CanESM2}
            scripts: null
        """)

    recipe = get_recipe(tmp_path, content, config_user)

    assert len(recipe.tasks) == 1
    task = recipe.tasks.pop()
    assert len(task.products) == 1
    product = task.products.pop()
    preproc_dir = os.path.dirname(product.filename)
    assert preproc_dir.startswith(str(tmp_path))

    fix_dir = os.path.join(
        preproc_dir, 'CMIP5_CanESM2_Oyr_historical_r1i1p1_chl_2000-2005_fixed')
    defaults = {
        'load': {
            'callback': concatenate_callback,
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
            'frequency': 'yr',
        },
        'fix_metadata': {
            'project': 'CMIP5',
            'dataset': 'CanESM2',
            'short_name': 'chl',
            'cmor_table': 'CMIP5',
            'mip': 'Oyr',
            'frequency': 'yr',
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
            'frequency': 'yr',
        },
        'cmor_check_data': {
            'cmor_table': 'CMIP5',
            'mip': 'Oyr',
            'short_name': 'chl',
            'frequency': 'yr',
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


def test_empty_variable(tmp_path, patched_datafinder, config_user):
    """Test that it is possible to specify all information in the dataset."""
    content = dedent("""
        diagnostics:
          diagnostic_name:
            additional_datasets:
              - dataset: CanESM2
                project: CMIP5
                mip: Amon
                exp: historical
                start_year: 2000
                end_year: 2005
                ensemble: r1i1p1
            variables:
              pr:
            scripts: null
        """)

    recipe = get_recipe(tmp_path, content, config_user)
    assert len(recipe.tasks) == 1
    task = recipe.tasks.pop()
    assert len(task.products) == 1
    product = task.products.pop()
    assert product.attributes['short_name'] == 'pr'
    assert product.attributes['dataset'] == 'CanESM2'


def test_reference_dataset(tmp_path, patched_datafinder, config_user,
                           monkeypatch):

    levels = [100]
    get_reference_levels = create_autospec(
        esmvaltool._recipe.get_reference_levels, return_value=levels)
    monkeypatch.setattr(esmvaltool._recipe, 'get_reference_levels',
                        get_reference_levels)

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

    recipe = get_recipe(tmp_path, content, config_user)

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
    assert product.settings['extract_levels']['levels'] == levels

    fix_dir = os.path.splitext(reference.filename)[0] + '_fixed'
    get_reference_levels.assert_called_once_with(
        reference.files[0],
        'CMIP5',
        'MPI-ESM-LR',
        'ta',
        fix_dir
    )

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


def test_custom_preproc_order(tmp_path, patched_datafinder, config_user):

    content = dedent("""
        preprocessors:
          default: &default
            average_region:
              coord1: longitude
              coord2: latitude
            multi_model_statistics:
              span: overlap
              statistics: [mean ]
          custom:
            custom_order: true
            <<: *default

        diagnostics:
          diagnostic_name:
            variables:
              chl_default: &chl
                short_name: chl
                preprocessor: default
                project: CMIP5
                mip: Oyr
                exp: historical
                start_year: 2000
                end_year: 2005
                ensemble: r1i1p1
                additional_datasets:
                  - {dataset: CanESM2}
              chl_custom:
                <<: *chl
                preprocessor: custom
            scripts: null
        """)

    recipe = get_recipe(tmp_path, content, config_user)

    assert len(recipe.tasks) == 2

    default = next(t for t in recipe.tasks if tuple(t.order) == DEFAULT_ORDER)
    custom = next(t for t in recipe.tasks if tuple(t.order) != DEFAULT_ORDER)

    assert custom.order.index('average_region') < custom.order.index(
        'multi_model_statistics')
    assert default.order.index('average_region') > default.order.index(
        'multi_model_statistics')


def test_derive(tmp_path, patched_datafinder, config_user):

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
                derive: true
                force_derivation: true
                additional_datasets:
                  - {dataset: GFDL-CM3,  ensemble: r1i1p1}
            scripts: null
        """)

    recipe = get_recipe(tmp_path, content, config_user)

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
    assert product.files

    ps_product = next(p for a in task.ancestors for p in a.products
                      if p.attributes['short_name'] == 'ps')
    tro3_product = next(p for a in task.ancestors for p in a.products
                        if p.attributes['short_name'] == 'tro3')
    assert ps_product.filename in product.files
    assert tro3_product.filename in product.files


def test_derive_not_needed(tmp_path, patched_datafinder, config_user):

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
                derive: true
                force_derivation: false
                additional_datasets:
                  - {dataset: GFDL-CM3,  ensemble: r1i1p1}
            scripts: null
        """)

    recipe = get_recipe(tmp_path, content, config_user)

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
    assert product.attributes['short_name'] == 'toz'
    assert 'derive' in product.settings

    assert len(ancestor.products) == 1
    ancestor_product = ancestor.products.pop()
    assert ancestor_product.filename in product.files
    assert ancestor_product.attributes['short_name'] == 'toz'
    assert 'derive' not in ancestor_product.settings

    # Check that fixes are applied just once
    fixes = ('fix_file', 'fix_metadata', 'fix_data')
    for fix in fixes:
        assert fix in ancestor_product.settings
        assert fix not in product.settings


def test_derive_with_fx(tmp_path, patched_datafinder, config_user):

    content = dedent("""
        diagnostics:
          diagnostic_name:
            variables:
              nbp_grid:
                project: CMIP5
                mip: Lmon
                exp: historical
                start_year: 2000
                end_year: 2005
                derive: true
                force_derivation: true
                additional_datasets:
                  - {dataset: GFDL-CM3,  ensemble: r1i1p1}
            scripts: null
        """)

    recipe = get_recipe(tmp_path, content, config_user)

    # Check generated tasks
    assert len(recipe.tasks) == 1
    task = recipe.tasks.pop()

    assert task.name == 'diagnostic_name' + TASKSEP + 'nbp_grid'
    assert len(task.ancestors) == 1
    ancestor = [t for t in task.ancestors][0]
    assert ancestor.name == 'diagnostic_name/nbp_grid_derive_input_nbp'

    # Check product content of tasks
    assert len(task.products) == 1
    product = task.products.pop()
    assert 'derive' in product.settings
    assert product.attributes['short_name'] == 'nbp_grid'
    assert 'fx_files' in product.settings['derive']
    assert 'sftlf' in product.settings['derive']['fx_files']
    assert product.settings['derive']['fx_files']['sftlf'] is not None

    assert len(ancestor.products) == 1
    ancestor_product = ancestor.products.pop()
    assert ancestor_product.filename in product.files
    assert ancestor_product.attributes['short_name'] == 'nbp'


def simulate_diagnostic_run(diagnostic_task):
    """Simulate Python diagnostic run."""
    cfg = diagnostic_task.settings
    input_files = [
        p.filename for a in diagnostic_task.ancestors for p in a.products
    ]
    record = {
        'caption': 'Test plot',
        'plot_file': get_plot_filename('test', cfg),
        'statistics': ['mean', 'var'],
        'domains': ['trop', 'et'],
        'plot_type': 'zonal',
        'authors': ['ande_bo'],
        'references': ['acknow_project'],
        'ancestors': input_files,
    }

    diagnostic_file = get_diagnostic_filename('test', cfg)
    create_test_file(diagnostic_file)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, record)

    diagnostic_task._collect_provenance()
    return record


def test_diagnostic_task_provenance(tmp_path, patched_datafinder, config_user):

    script = tmp_path / 'diagnostic.py'
    with script.open('w'):
        pass

    content = dedent("""
        diagnostics:
          diagnostic_name:
            themes:
              - phys
            realms:
              - atmos
            variables:
              chl:
                project: CMIP5
                mip: Oyr
                exp: historical
                start_year: 2000
                end_year: 2005
                ensemble: r1i1p1
                additional_datasets:
                  - dataset: CanESM2
            scripts:
              script_name:
                script: {script}
              script_name2:
                script: {script}
                ancestors: [script_name]
        """.format(script=script))

    recipe = get_recipe(tmp_path, content, config_user)
    diagnostic_task = recipe.tasks.pop()

    simulate_diagnostic_run(next(iter(diagnostic_task.ancestors)))
    record = simulate_diagnostic_run(diagnostic_task)

    # Check resulting product
    product = diagnostic_task.products.pop()
    check_provenance(product)
    for key in ('caption', 'plot_file'):
        assert product.attributes[key] == record[key]
        assert product.entity.get_attribute('attribute:' +
                                            key).pop() == record[key]

    # Check that diagnostic script tags have been added
    with open(
            os.path.join(
                os.path.dirname(esmvaltool.__file__),
                'config-references.yml')) as file:
        tags = yaml.safe_load(file)
    for key in ('statistics', 'domains', 'authors', 'references'):
        assert product.attributes[key] == tuple(
            tags[key][k] for k in record[key])

    # Check that recipe diagnostic tags have been added
    src = yaml.safe_load(DEFAULT_DOCUMENTATION + content)
    for key in ('realms', 'themes'):
        value = src['diagnostics']['diagnostic_name'][key]
        assert product.attributes[key] == tuple(tags[key][k] for k in value)

    # Check that recipe tags have been added
    recipe_record = product.provenance.get_record('recipe:recipe_test.yml')
    assert len(recipe_record) == 1
    for key in ('description', 'references'):
        value = src['documentation'][key]
        if key == 'references':
            value = ', '.join(tags[key][k] for k in value)
        assert recipe_record[0].get_attribute('attribute:' +
                                              key).pop() == value

    # Test that provenance was saved to netcdf, xml and svg plot
    cube = iris.load(product.filename)[0]
    assert 'provenance' in cube.attributes
    prefix = os.path.splitext(product.filename)[0] + '_provenance'
    assert os.path.exists(prefix + '.xml')
    assert os.path.exists(prefix + '.svg')
