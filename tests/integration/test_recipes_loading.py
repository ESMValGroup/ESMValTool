"""Test recipes are well formed."""
from pathlib import Path

import esmvalcore
import esmvalcore._config
import esmvalcore.cmor.check
import pytest
import yaml

import esmvaltool

from esmvalcore import __version__ as core_ver
from packaging import version
from .test_diagnostic_run import write_config_user_file


@pytest.fixture(scope='module')
def config_user(tmp_path_factory):
    """Generate dummy config-user dict for testing purposes."""
    path = tmp_path_factory.mktemp('recipe-test')
    filename = write_config_user_file(path)
    # The fixture scope is set to module to avoid very slow
    # test runs, as the following line also reads the CMOR tables
    cfg = esmvalcore._config.read_config_user_file(filename, 'recipe_test', {})
    cfg['offline'] = True
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
def test_recipe_valid(recipe_file, config_user, mocker):
    """Check that recipe files are valid ESMValTool recipes."""
    # Mock input files
    try:
        # Since ESValCore v2.8.0
        import esmvalcore.local
        module = esmvalcore.local
        method = 'glob'
        # The patched_datafinder fixture does not return the correct input
        # directory structure, so make sure it is set to flat for every project
        from esmvalcore.config import CFG, _config
        mocker.patch.dict(CFG, drs={})
        for project in _config.CFG:
            mocker.patch.dict(_config.CFG[project]['input_dir'], default='/')
    except ImportError:
        # Prior to ESMValCore v2.8.0
        import esmvalcore._data_finder
        module = esmvalcore._data_finder
        method = 'find_files'

    mocker.patch.object(
        module,
        method,
        autospec=True,
        side_effect=lambda *_, **__: [
            'test_0001-1849.nc',
            'test_1850-9999.nc',
        ],
    )

    # Mock vertical levels
    # Account for module change after esmvalcore=2.7
    if version.parse(core_ver) <= version.parse('2.7.1'):
        import esmvalcore._recipe
        mocker.patch.object(
            esmvalcore._recipe,
            'get_reference_levels',
            autospec=True,
            spec_set=True,
            side_effect=lambda *_, **__: [1, 2],
        )
    else:
        import esmvalcore._recipe.recipe
        mocker.patch.object(
            esmvalcore._recipe.recipe,
            'get_reference_levels',
            autospec=True,
            spec_set=True,
            side_effect=lambda *_, **__: [1, 2],
        )

    # Mock valid NCL version
    # Account for module change after esmvalcore=2.7
    if version.parse(core_ver) <= version.parse('2.7.1'):
        import esmvalcore._recipe_checks
        mocker.patch.object(
            esmvalcore._recipe_checks,
            'ncl_version',
            autospec=True,
            spec_set=True,
        )
    else:
        import esmvalcore._recipe.check
        mocker.patch.object(
            esmvalcore._recipe.check,
            'ncl_version',
            autospec=True,
            spec_set=True,
        )

    # Mock interpreters installed
    def which(executable):
        if executable in ('julia', 'ncl', 'python', 'Rscript'):
            path = '/path/to/' + executable
        else:
            path = None
        return path

    mocker.patch.object(
        esmvalcore._task,
        'which',
        autospec=True,
        side_effect=which,
    )

    # Create a shapefile for extract_shape preprocessor if needed
    recipe = yaml.safe_load(recipe_file.read_text())
    for preproc in recipe.get('preprocessors', {}).values():
        extract_shape = preproc.get('extract_shape')
        if extract_shape and 'shapefile' in extract_shape:
            filename = Path(
                config_user['auxiliary_data_dir']) / extract_shape['shapefile']
            filename.parent.mkdir(parents=True, exist_ok=True)
            filename.touch()

    # Account for module change after esmvalcore=2.7
    if version.parse(core_ver) <= version.parse('2.7.1'):
        esmvalcore._recipe.read_recipe_file(recipe_file, config_user)
    else:
        esmvalcore._recipe.recipe.read_recipe_file(recipe_file, config_user)
