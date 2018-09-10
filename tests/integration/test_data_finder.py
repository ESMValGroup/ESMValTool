import os
import shutil
import tempfile

import pytest
import yaml

import esmvaltool._config
from esmvaltool._data_finder import get_input_filelist

# Initialize with standard config developer file
esmvaltool._config.CFG = esmvaltool._config.read_config_developer_file()

with open(os.path.join(os.path.dirname(__file__), 'data_finder.yml')) as file:
    CONFIG = yaml.safe_load(file)


def print_path(path):
    txt = path
    if os.path.isdir(path):
        txt += '/'
    if os.path.islink(path):
        txt += ' -> ' + os.readlink(path)
    print(txt)


def tree(path):
    """Print path, similar to the the `tree` command."""
    print_path(path)
    for root, dirnames, filenames in os.walk(path):
        for dirname in dirnames:
            print_path(os.path.join(root, dirname))
        for filename in filenames:
            print_path(os.path.join(root, filename))


def create_file(filename):
    """Create an empty file."""
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    with open(filename, 'a'):
        pass


@pytest.fixture
def root():
    dirname = tempfile.mkdtemp()
    yield os.path.join(dirname, 'output1')
    print("Directory structure was:")
    tree(dirname)
    shutil.rmtree(dirname)


@pytest.mark.parametrize('cfg', CONFIG['get_input_filelist'])
def test_get_input_filelist(root, cfg):

    # Create directory structure and files
    for filename in cfg.get('available_files', []):
        create_file(os.path.join(root, filename))

    for symlink in cfg.get('available_symlinks', {}):
        link_name = os.path.join(root, symlink['link_name'])
        os.symlink(symlink['target'], link_name)

    # Find files
    rootpath = {cfg['variable']['project']: root}
    drs = {cfg['variable']['project']: cfg['drs']}
    input_filelist = get_input_filelist(cfg['variable'], rootpath, drs)

    # Test result
    reference = [os.path.join(root, file) for file in cfg['found_files']]
    assert sorted(input_filelist) == sorted(reference)
