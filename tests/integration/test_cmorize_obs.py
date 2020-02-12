"""Tests for the module :mod:`esmvaltool.cmorizers.obs.cmorize_obs`."""

import contextlib
import os
import sys
import yaml
from esmvaltool.cmorizers.obs.cmorize_obs import main as run


def write_config_user_file(dirname):
    config_file = dirname / 'config-user.yml'
    cfg = {
        'output_dir': str(dirname / 'output_dir'),
        'rootpath': {
            'RAWOBS': str(dirname / 'raw_stuff'),
        },
        'log_level': 'debug',
    }
    config_file.write_text(yaml.safe_dump(cfg, encoding=None))
    return str(config_file)


@contextlib.contextmanager
def arguments(*args):
    backup = sys.argv
    sys.argv = list(args)
    yield
    sys.argv = backup


def test_cmorize_obs_woa(tmp_path):
    """Test for example run of cmorize_obs command."""

    config_user_file = write_config_user_file(tmp_path)
    os.makedirs(os.path.join(tmp_path, 'raw_stuff'))
    os.makedirs(os.path.join(tmp_path, 'raw_stuff', 'Tier2'))
    with arguments(
            'cmorize_obs',
            '-c',
            config_user_file,
            '-o',
            'WOA',
    ):
        run()


def test_cmorize_obs_cru(tmp_path):
    """Test for example run of cmorize_obs command."""

    config_user_file = write_config_user_file(tmp_path)
    os.makedirs(os.path.join(tmp_path, 'raw_stuff'))
    os.makedirs(os.path.join(tmp_path, 'raw_stuff', 'Tier2'))
    with arguments(
            'cmorize_obs',
            '-c',
            config_user_file,
            '-o',
            'CRU',
    ):
        run()
