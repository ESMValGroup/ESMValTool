"""Test recipes are well formed."""
from pathlib import Path
from unittest.mock import create_autospec

import pytest
import yaml

import esmvalcore
import esmvalcore._config
import esmvalcore._data_finder
import esmvalcore._recipe
import esmvalcore.cmor.check
import esmvaltool

from .test_diagnostic_run import write_config_user_file


@pytest.fixture(scope='module')
def config_user(tmp_path_factory):
    """Generate dummy config-user dict for testing purposes."""
    path = tmp_path_factory.mktemp('recipe-test')
    filename = write_config_user_file(path)
    # The fixture scope is set to module to avoid very slow
    # test runs, as the following line also reads the CMOR tables
    cfg = esmvalcore._config.read_config_user_file(filename, 'recipe_test', {})
    cfg['synda_download'] = False
    cfg['auxiliary_data_dir'] = str(path / 'auxiliary_data_dir')
    cfg['check_level'] = esmvalcore.cmor.check.CheckLevels['DEFAULT']
    return cfg


def _get_recipes():
    recipes_path = Path(esmvaltool.__file__).absolute().parent / 'recipes'
    recipes = sorted(recipes_path.glob("**/recipe*.yml"))
    ids = tuple(str(p.relative_to(recipes_path)) for p in recipes)
    return recipes, ids


RECIPES, IDS = _get_recipes()


@pytest.mark.parametrize('recipe_file', RECIPES, ids=IDS)
def test_recipe_valid(recipe_file, config_user, monkeypatch):
    """Check that recipe files are valid ESMValTool recipes."""
    # Mock input files
    find_files = create_autospec(esmvalcore._data_finder.find_files,
                                 spec_set=True)
    find_files.side_effect = lambda *_, **__: [
        'test_0000-1849.nc',
        'test_1850-9999.nc',
    ]
    monkeypatch.setattr(esmvalcore._data_finder, 'find_files', find_files)

    # Mock vertical levels
    levels = create_autospec(esmvalcore._recipe.get_reference_levels,
                             spec_set=True)
    levels.side_effect = lambda *_, **__: [1, 2]
    monkeypatch.setattr(esmvalcore._recipe, 'get_reference_levels', levels)

    # Mock valid NCL version
    ncl_version = create_autospec(esmvalcore._recipe_checks.ncl_version,
                                  spec_set=True)
    monkeypatch.setattr(esmvalcore._recipe_checks, 'ncl_version', ncl_version)

    # Mock interpreters installed
    def which(executable):
        if executable in ('julia', 'ncl', 'python', 'Rscript'):
            path = '/path/to/' + executable
        else:
            path = None
        return path

    monkeypatch.setattr(esmvalcore._task, 'which', which)

    # Create a shapefile for extract_shape preprocessor if needed
    recipe = yaml.safe_load(recipe_file.read_text())
    for preproc in recipe.get('preprocessors', {}).values():
        extract_shape = preproc.get('extract_shape')
        if extract_shape and 'shapefile' in extract_shape:
            filename = Path(
                config_user['auxiliary_data_dir']) / extract_shape['shapefile']
            filename.parent.mkdir(parents=True, exist_ok=True)
            filename.touch()

    esmvalcore._recipe.read_recipe_file(recipe_file, config_user)
