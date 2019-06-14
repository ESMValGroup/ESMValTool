"""Tests for _data_finder.py."""
import os
import shutil
import tempfile

import pytest
import yaml

import esmvalcore._config
from esmvalcore._data_finder import (get_input_filelist, get_input_fx_filelist,
                                     get_output_file)
from esmvalcore.cmor.table import read_cmor_tables

# Initialize with standard config developer file
esmvalcore._config.CFG = esmvalcore._config.read_config_developer_file()
# Initialize CMOR tables
read_cmor_tables(esmvalcore._config.CFG)

# Load test configuration
with open(os.path.join(os.path.dirname(__file__), 'data_finder.yml')) as file:
    CONFIG = yaml.safe_load(file)


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


@pytest.mark.parametrize('cfg', CONFIG['get_output_file'])
def test_get_output_file(cfg):
    """Test getting output name for preprocessed files."""
    output_file = get_output_file(cfg['variable'], cfg['preproc_dir'])
    assert output_file == cfg['output_file']


@pytest.fixture
def root():
    """Root function for tests."""
    dirname = tempfile.mkdtemp()
    yield os.path.join(dirname, 'output1')
    print("Directory structure was:")
    tree(dirname)
    shutil.rmtree(dirname)


@pytest.mark.parametrize('cfg', CONFIG['get_input_filelist'])
def test_get_input_filelist(root, cfg):
    """Test retrieving input filelist."""
    create_tree(root, cfg.get('available_files'),
                cfg.get('available_symlinks'))

    # Find files
    rootpath = {cfg['variable']['project']: [root]}
    drs = {cfg['variable']['project']: cfg['drs']}
    input_filelist = get_input_filelist(cfg['variable'], rootpath, drs)

    # Test result
    reference = [os.path.join(root, file) for file in cfg['found_files']]
    assert sorted(input_filelist) == sorted(reference)


@pytest.mark.parametrize('cfg', CONFIG['get_input_fx_filelist'])
def test_get_input_fx_filelist(root, cfg):
    """Test retrieving fx filelist."""
    create_tree(root, cfg.get('available_files'),
                cfg.get('available_symlinks'))

    # Find files
    rootpath = {cfg['variable']['project']: [root]}
    drs = {cfg['variable']['project']: cfg['drs']}
    fx_files = get_input_fx_filelist(cfg['variable'], rootpath, drs)

    # Test result
    reference = {
        fx_var: os.path.join(root, filename) if filename else None
        for fx_var, filename in cfg['found_files'].items()
    }
    assert fx_files == reference
