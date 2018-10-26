from pprint import pformat

import pytest
import yaml
from netCDF4 import Dataset

import esmvaltool
from esmvaltool._recipe import read_recipe_file
from tests.integration.test_diagnostic_run import write_config_user_file


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
    def find_files(dirnames, filename):
        if filename.endswith('*'):
            filename = filename.rstrip('*') + '1950_2010.nc'
        filename = str(tmpdir / filename)
        create_test_file(filename)
        return [filename]

    monkeypatch.setattr(esmvaltool._data_finder, 'find_files', find_files)

    recipe_content = {
        'preprocessors': {
            'preprocessor_name': {
                'extract_levels': {
                    'levels': 85000,
                    'scheme': 'nearest',
                },
            },
        },
        'diagnostics': {
            'diagnostic_name': {
                'variables': {
                    'ta': {
                        'preprocessor':
                        'preprocessor_name',
                        'field':
                        'T3M',
                        'additional_datasets': [
                            {
                                'dataset': 'bcc-csm1-1',
                                'project': 'CMIP5',
                                'mip': 'Amon',
                                'exp': 'historical',
                                'ensemble': 'r1i1p1',
                                'start_year': 2000,
                                'end_year': 2002
                            },
                        ],
                    },
                },
                'scripts': None,
            }
        }
    }

    recipe_file = str(tmpdir / 'recipe_test.yml')
    with open(recipe_file, 'w') as file:
        yaml.safe_dump(recipe_content, file)

    recipe = read_recipe_file(str(recipe_file), config_user)
    for diagnostic in recipe.diagnostics.values():
        print(pformat(diagnostic))
    for task in recipe.tasks:
        print(task.order)
        for product in task.products:
            print(product)
