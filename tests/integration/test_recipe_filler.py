"""Tests for _data_finder.py."""
import contextlib
import os
import shutil
import sys
import tempfile

import pytest
import yaml

import esmvaltool.utils.recipe_filler as recipe_filler
from esmvalcore.cmor.table import read_cmor_tables
from esmvaltool.utils.recipe_filler import run

# Initialize with standard config developer file
std_config = recipe_filler.read_config_developer_file()
# Initialize CMOR tables
read_cmor_tables(std_config)

# Load test configuration
with open(os.path.join(os.path.dirname(__file__),
                       'recipe_filler.yml')) as file:
    CONFIG = yaml.safe_load(file)


@contextlib.contextmanager
def arguments(*args):
    backup = sys.argv
    sys.argv = list(args)
    yield
    sys.argv = backup


def print_path(path):
    """Print path."""
    txt = path
    if os.path.isdir(path):
        txt += '/'
    if os.path.islink(path):
        txt += ' -> ' + os.readlink(path)
    print(txt)


def tree(path):
    """Print path, similar to the the `tree` command."""
    print_path(path)
    for dirpath, dirnames, filenames in os.walk(path):
        for dirname in dirnames:
            print_path(os.path.join(dirpath, dirname))
        for filename in filenames:
            print_path(os.path.join(dirpath, filename))


def create_file(filename):
    """Create an empty file."""
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    with open(filename, 'a'):
        pass


def create_tree(path, filenames=None, symlinks=None):
    """Create directory structure and files."""
    for filename in filenames or []:
        create_file(os.path.join(path, filename))

    for symlink in symlinks or []:
        link_name = os.path.join(path, symlink['link_name'])
        os.symlink(symlink['target'], link_name)


def write_config_user_file(dirname, file_path, drs):
    config_file = dirname / 'config-user.yml'
    cfg = {
        'log_level': 'info',
        'output_dir': str(dirname / 'recipe_filler_output'),
        'rootpath': {
            'CMIP5': str(dirname / file_path),
            'CMIP6': str(dirname / file_path),
        },
        'drs': {
            'CMIP5': drs,
            'CMIP6': drs,
        },
    }
    config_file.write_text(yaml.safe_dump(cfg, encoding=None))
    return str(config_file)


def write_recipe(dirname, recipe_dict):
    recipe_file = dirname / 'recipe.yml'
    diags = {'diagnostics': recipe_dict}
    recipe_file.write_text(yaml.safe_dump(diags, encoding=None))
    return str(recipe_file)


@pytest.fixture
def root():
    """Root function for tests."""
    dirname = tempfile.mkdtemp()
    yield os.path.join(dirname, 'output1')
    print("Directory structure was:")
    tree(dirname)
    shutil.rmtree(dirname)


def setup_files(tmp_path, root, cfg):
    """Create config, recipe ,output recipe etc."""
    user_config_file = write_config_user_file(tmp_path, root, cfg['drs'])
    diagnostics = {}
    diagnostics["test_diagnostic"] = {}
    diagnostics["test_diagnostic"]["variables"] = {}
    diagnostics["test_diagnostic"]["variables"]["test_var"] = cfg["variable"]
    recipe = write_recipe(tmp_path, diagnostics)
    output_recipe = str(tmp_path / "recipe_auto.yml")

    return user_config_file, recipe, output_recipe


@pytest.mark.parametrize('cfg', CONFIG['has_additional_datasets'])
def test_adding_datasets(tmp_path, root, cfg):
    """Test retrieving additional datasets."""
    create_tree(root, cfg.get('available_files'),
                cfg.get('available_symlinks'))

    user_config_file, recipe, output_recipe = setup_files(tmp_path, root, cfg)

    with arguments(
            'recipe_filler',
            recipe,
            '-c',
            user_config_file,
            '-o',
            output_recipe,
    ):
        run()

    with open(output_recipe, 'r') as file:
        autofilled_recipe = yaml.safe_load(file)
        diag = autofilled_recipe["diagnostics"]["test_diagnostic"]
        var = diag["variables"]["test_var"]
        assert "additional_datasets" in var


@pytest.mark.parametrize('cfg', CONFIG['no_additional_datasets'])
def test_not_adding_datasets(tmp_path, root, cfg):
    """Test retrieving no additional datasets."""
    create_tree(root, cfg.get('available_files'),
                cfg.get('available_symlinks'))

    user_config_file, recipe, output_recipe = setup_files(tmp_path, root, cfg)

    with arguments(
            'recipe_filler',
            recipe,
            '-c',
            user_config_file,
            '-o',
            output_recipe,
    ):
        run()

    with open(output_recipe, 'r') as file:
        autofilled_recipe = yaml.safe_load(file)
        diag = autofilled_recipe["diagnostics"]["test_diagnostic"]
        var = diag["variables"]["test_var"]
        assert "additional_datasets" not in var


def test_bad_var(tmp_path, root):
    """Test a bad variable in the works."""
    cfg = CONFIG['bad_variable'][0]
    user_config_file, recipe, output_recipe = setup_files(tmp_path, root, cfg)

    # this doesn't fail and it shouldn't since it can go on
    # and look for data for other valid variables
    with arguments(
            'recipe_filler',
            recipe,
            '-c',
            user_config_file,
            '-o',
            output_recipe,
    ):
        run()

    with open(output_recipe, 'r') as file:
        autofilled_recipe = yaml.safe_load(file)
        diag = autofilled_recipe["diagnostics"]["test_diagnostic"]
        var = diag["variables"]["test_var"]
        assert "additional_datasets" not in var


def test_no_short_name(tmp_path, root):
    """Test a bad variable in the works."""
    cfg = CONFIG['no_short_name'][0]
    user_config_file, recipe, output_recipe = setup_files(tmp_path, root, cfg)

    # this doesn't fail and it shouldn't since it can go on
    # and look for data for other valid variables
    with arguments(
            'recipe_filler',
            recipe,
            '-c',
            user_config_file,
            '-o',
            output_recipe,
    ):
        run()

    with open(output_recipe, 'r') as file:
        autofilled_recipe = yaml.safe_load(file)
        diag = autofilled_recipe["diagnostics"]["test_diagnostic"]
        var = diag["variables"]["test_var"]
        assert "additional_datasets" not in var
