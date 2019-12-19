"""Test recipes are well formed."""

import os
from pathlib import Path

import iris
import numpy as np
import pytest

from esmvalcore import _data_finder, _recipe_checks
from esmvalcore._config import read_config_user_file
from esmvalcore._recipe import read_recipe_file

import esmvaltool

from .test_diagnostic_run import write_config_user_file


def _get_recipes():
    recipes_path = Path(esmvaltool.__file__).absolute().parent / 'recipes'
    recipes = recipes_path.glob("**/recipe*.yml")
    return recipes


def _tracking_ids(i=0):
    while True:
        yield i
        i += 1


def _create_test_file(filename, tracking_id=None):
    """Generate dummy data file for recipe checker."""
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    attributes = {}
    if tracking_id is not None:
        attributes['tracking_id'] = tracking_id

    xcoord = iris.coords.DimCoord(np.linspace(0, 5, 5),
                                  standard_name="longitude")
    ycoord = iris.coords.DimCoord(np.linspace(0, 5, 12),
                                  standard_name="latitude")
    zcoord = iris.coords.DimCoord(np.linspace(0, 5, 17),
                                  standard_name="height",
                                  attributes={'positive': 'up'})
    cube = iris.cube.Cube(np.zeros((5, 12, 17), np.float32),
                          dim_coords_and_dims=[(xcoord, 0), (ycoord, 1),
                                               (zcoord, 2)],
                          attributes=attributes)
    iris.save(cube, filename)


def _get_dummy_filenames(drs):
    """Generate list of realistic dummy filename(s) according to drs.

    drs is the directory structure used to find input files in ESMValTool
    """
    dummy_filenames = []

    # Time-invariant (fx) variables don't have years in their filename
    if 'fx' in drs:
        if drs.endswith('[_.]*nc'):
            dummy_filename = drs.replace('[_.]*', '.')
        elif drs.endswith('*.nc'):
            dummy_filename = drs.replace('*', '')
        dummy_filenames.append(dummy_filename)
    # For other variables, add custom (large) intervals in dummy filename
    elif '*' in drs:
        if drs.endswith('[_.]*nc'):
            dummy_filename = drs[:-len('[_.]*nc')]
        elif drs.endswith('*.nc'):
            dummy_filename = drs[:-len('*.nc')]
        # Spread dummy data over multiple files for realistic test
        # Note: adding too many intervals here makes the tests really slow!
        for interval in ['0000_1849', '1850_9999']:
            dummy_filenames.append(dummy_filename + '_' + interval + '.nc')
    # Provide for the possibility of filename drss without *.
    else:
        dummy_filename = drs
        dummy_filenames.append(dummy_filename)
    return dummy_filenames


@pytest.fixture
def config_user(tmp_path):
    """Generate dummy config-user file for testing purposes."""
    filename = write_config_user_file(tmp_path)
    cfg = read_config_user_file(filename, 'recipe_test')
    cfg['synda_download'] = False
    return cfg


@pytest.fixture
def patched_datafinder(tmp_path, monkeypatch):
    """Replace `_datafinder.find_files()`.

    Creates and points to dummy data input files instead of searching for
    existing data.
    """
    def find_files(_, filenames):
        drs = filenames[0]
        dummyfiles = str(tmp_path / 'input' / drs)
        filenames = _get_dummy_filenames(dummyfiles)

        for file in filenames:
            _create_test_file(file, next(tracking_id))

        return filenames

    tracking_id = _tracking_ids()
    monkeypatch.setattr(_data_finder, 'find_files', find_files)


@pytest.fixture
def patched_extract_shape(monkeypatch):
    """Replace `_recipe_checks.extract_shape`.

    Skips check that shapefile exists.
    """
    def extract_shape(settings):
        valid = {
            'method': {'contains', 'representative'},
            'crop': {True, False},
        }

        for key in valid:
            value = settings.get(key)
            if not (value is None or value in valid[key]):
                raise _recipe_checks.RecipeError(
                    "In preprocessor function `extract_shape`: Invalid value"
                    f"'{value}' for argument '{key}', choose from "
                    "{}".format(', '.join(f"'{k}'".lower()
                                          for k in valid[key])))

    monkeypatch.setattr(_recipe_checks, 'extract_shape', extract_shape)


@pytest.mark.parametrize('recipe_file', _get_recipes())
def test_diagnostic_run(recipe_file, config_user, patched_datafinder,
                        patched_extract_shape):
    """Check that recipe files are valid ESMValTool recipes."""
    read_recipe_file(recipe_file, config_user)
