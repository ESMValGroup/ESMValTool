"""Test recipes are well formed."""

import os
from pathlib import Path

import iris
import numpy as np
import pytest

from esmvalcore import _data_finder
from esmvalcore._config import read_config_user_file
from esmvalcore._recipe import read_recipe_file

import esmvaltool

from .test_diagnostic_run import write_config_user_file


def _get_recipes():
    recipes_path = Path(esmvaltool.__file__).absolute().parent / 'recipes'
    recipes = recipes_path.glob("**/recipe*.yml")
    return recipes


@pytest.fixture
def config_user(tmp_path):
    """Generate dummy config-user file for testing purposes."""
    filename = write_config_user_file(tmp_path)
    cfg = read_config_user_file(filename, 'recipe_test')
    cfg['synda_download'] = False
    return cfg


def create_test_file(filename, tracking_id=None):
    """Generate dummy data file for recipe checker."""
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    attributes = {}
    if tracking_id is not None:
        attributes['tracking_id'] = tracking_id
    cube = iris.cube.Cube([], attributes=attributes)
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


def get_dummy_filenames(pattern):
    """Generate list of realistic dummy filename(s) according to pattern."""
    dummy_filenames = []
    # Time-invariant (fx) variables don't have years in their filename
    if 'fx' in pattern:
        if pattern.endswith('[_.]*nc'):
            dummy_filename = pattern.replace('[_.]*', '.')
        elif pattern.endswith('*.nc'):
            dummy_filename = pattern.replace('*', '')
        dummy_filenames.append(dummy_filename)
    # For other variables, add custom (large) intervals in dummy filename
    elif '*' in pattern:
        if pattern.endswith('[_.]*nc'):
            dummy_filename = pattern[:-len('[_.]*nc')]
        elif pattern.endswith('*.nc'):
            dummy_filename = pattern[:-len('*.nc')]
        # Spread dummy data over multiple files for realistic test
        # Note: adding too many intervals here makes the tests really slow!
        for interval in ['1800_1949', '1950_2199']:
            dummy_filenames.append(dummy_filename + '_' + interval + '.nc')
    # Provide for the possibility of filename patterns without *.
    else:
        dummy_filename = pattern
        dummy_filenames.append(dummy_filename)
    return dummy_filenames


@pytest.fixture
def patched_datafinder(tmp_path, monkeypatch):
    """Alternative to `_datafinder.find_files()` for dummy files in tests."""
    def tracking_ids(i=0):
        while True:
            yield i
            i += 1

    tracking_id = tracking_ids()

    def find_files(_, filenames):
        filename_pattern = filenames[0]
        dummyfiles = str(tmp_path / 'input' / filename_pattern)
        filenames = get_dummy_filenames(dummyfiles)

        # Any occurrence of [something] in filename should have
        # been replaced before this function is called.
        for filename in filenames:
            assert '[' not in filename

        for file in filenames:
            create_test_file(file, next(tracking_id))

        return filenames

    monkeypatch.setattr(_data_finder, 'find_files', find_files)


@pytest.mark.parametrize('recipe_file', _get_recipes())
def test_diagnostic_run(recipe_file, config_user, patched_datafinder):
    """Check that recipe files are well formed."""
    read_recipe_file(recipe_file, config_user)
