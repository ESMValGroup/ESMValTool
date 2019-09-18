"""Test recipes are well formed."""

import glob
import os

import iris
import pytest

import esmvalcore
from esmvalcore._recipe import read_recipe_file

from .test_diagnostic_run import write_config_user_file


def _get_recipes():
    recipes_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '..',
        '..',
        'esmvaltool',
        'recipes',
    )
    recipes_path = os.path.abspath(recipes_path)
    recipes = glob.glob(os.path.join(recipes_path, '*.yml'))
    recipes += glob.glob(os.path.join(recipes_path, 'examples', '*.yml'))
    return recipes


@pytest.fixture
def config_user(tmp_path):
    filename = write_config_user_file(tmp_path)
    cfg = esmvalcore._config.read_config_user_file(filename, 'recipe_test')
    cfg['synda_download'] = False
    return cfg


def create_test_file(filename, tracking_id=None):
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    attributes = {}
    if tracking_id is not None:
        attributes['tracking_id'] = tracking_id
    print(attributes)
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
        print(filenames)
        for filename in filenames:
            assert '[' not in filename

        filename = filenames[0]
        filename = str(tmp_path / 'input' / filename)
        filenames = []
        if filename.endswith('*.nc'):
            filename = filename[:-len('*.nc')]
            intervals = [
                '1980_1989',
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

    monkeypatch.setattr(esmvalcore._data_finder, 'find_files', find_files)


@pytest.mark.parametrize('recipe_file', _get_recipes())
def test_diagnostic_run(recipe_file, config_user, patched_datafinder):
    """Check that recipe files are well formed."""
    read_recipe_file(recipe_file, config_user)
