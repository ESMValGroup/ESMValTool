"""Tests for the module :mod:`esmvaltool.cmorizers.data.cmorize_obs`."""

import contextlib
import os
import sys

import iris
import iris.coord_systems
import iris.coords
import iris.cube
import iris.fileformats
import numpy as np
import pytest
import yaml
from cf_units import Unit

from esmvaltool.cmorizers.data.cmorizer import DataCommand


@contextlib.contextmanager
def keep_cwd():
    """Use a context manager since the cmorizer enters and stays in the
    cmorization dir, risking to write test outputs away from test-reports."""
    curr_path = os.getcwd()
    try:
        yield
    finally:
        os.chdir(curr_path)


def write_config_user_file(dirname):
    """Replace config_user file values for testing."""
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


def _create_sample_cube(time_step):
    """Create a quick CMOR-compliant sample cube."""
    coord_sys = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
    cube_data = np.ones((1, 3, 2, 2))
    cube_data[0, 1, 1, 1] = 22.
    time = iris.coords.DimCoord([
        time_step,
    ],
                                standard_name='time',
                                bounds=[[time_step - 0.5, time_step + 0.5]],
                                units=Unit('days since 0000-01-01',
                                           calendar='gregorian'))
    zcoord = iris.coords.DimCoord([0.5, 5., 50.],
                                  var_name='depth',
                                  standard_name='depth',
                                  bounds=[[0., 2.5], [2.5, 25.], [25., 250.]],
                                  units='m',
                                  attributes={'positive': 'down'})
    lons = iris.coords.DimCoord([1.5, 2.5],
                                standard_name='longitude',
                                bounds=[[1., 2.], [2., 3.]],
                                units='K',
                                coord_system=coord_sys)
    lats = iris.coords.DimCoord([1.5, 2.5],
                                standard_name='latitude',
                                bounds=[[1., 2.], [2., 3.]],
                                units='K',
                                coord_system=coord_sys)
    coords_spec = [(time, 0), (zcoord, 1), (lats, 2), (lons, 3)]
    cube = iris.cube.Cube(cube_data, dim_coords_and_dims=coords_spec)
    return cube


def put_dummy_data(data_path):
    """Create a small dummy netCDF file to be cmorized."""
    data_info = [
        # dir_name, file_name_prefix, var_name
        ("temperature", "woa18_decav81B0_t", "t_an"),
        ("salinity", "woa18_decav81B0_s", "s_an"),
        ("oxygen", "woa18_all_o", "o_an"),
        ("nitrate", "woa18_all_n", "n_an"),
        ("phosphate", "woa18_all_p", "p_an"),
        ("silicate", "woa18_all_i", "i_an"),
    ]

    for (dir_name, file_name_prefix, var_name) in data_info:
        file_dir = os.path.join(data_path, dir_name)
        os.makedirs(file_dir)
        for month, step in enumerate(np.arange(0.5, 12.5)):
            gen_cube = _create_sample_cube(step)
            file_name = f"{file_name_prefix}{month:02d}_01.nc"
            file_path = os.path.join(file_dir, file_name)
            gen_cube.var_name = var_name
            iris.save(gen_cube, file_path)


def check_log_file(log_file, no_data=False):
    """Check the cmorization log file."""
    with open(log_file, 'r') as log:
        if no_data:
            msg = "Data for WOA not found"
        else:
            msg = "Fixing data"
        assert any(msg in line for line in log)


def check_output_exists(output_path):
    """Check if cmorizer outputted."""
    # eg Tier2/WOA/OBS6_WOA_clim_2018_Omon_thetao_200001-200012.nc
    output_files = os.listdir(output_path)
    assert len(output_files) == 8
    assert 'OBS6_WOA_clim' in output_files[0]
    out_files = [s.split("_")[5] for s in output_files]
    assert 'thetao' in out_files
    assert 'so' in out_files
    assert 'no3' in out_files
    assert 'po4' in out_files
    assert 'o2' in out_files
    assert 'si' in out_files
    assert 'sos' in out_files
    assert 'tos' in out_files


def check_conversion(output_path):
    """Check basic cmorization."""
    cube = iris.load_cube(os.path.join(output_path,
                                       os.listdir(output_path)[0]))
    assert cube.coord("time").units == Unit('days since 1950-1-1 00:00:00',
                                            calendar='gregorian')
    assert cube.coord("latitude").units == 'degrees'


@contextlib.contextmanager
def arguments(*args):
    """Arrange contextmanager."""
    backup = sys.argv
    sys.argv = list(args)
    yield
    sys.argv = backup


def test_cmorize_obs_woa_no_data(tmp_path):
    """Test for example run of cmorize_obs command."""

    config_user_file = write_config_user_file(tmp_path)
    os.makedirs(os.path.join(tmp_path, 'raw_stuff', 'Tier2'))
    with keep_cwd():
        with pytest.raises(Exception):
            DataCommand().format('WOA', config_user_file)

    log_dir = os.path.join(tmp_path, 'output_dir')
    log_file = os.path.join(log_dir,
                            os.listdir(log_dir)[0], 'run', 'main_log.txt')
    check_log_file(log_file, no_data=True)


def test_cmorize_obs_woa_data(tmp_path):
    """Test for example run of cmorize_obs command."""

    config_user_file = write_config_user_file(tmp_path)
    data_path = os.path.join(tmp_path, 'raw_stuff', 'Tier2', 'WOA')
    put_dummy_data(data_path)
    with keep_cwd():
        DataCommand().format('WOA', config_user_file)

    log_dir = os.path.join(tmp_path, 'output_dir')
    log_file = os.path.join(log_dir,
                            os.listdir(log_dir)[0], 'run', 'main_log.txt')
    check_log_file(log_file, no_data=False)
    output_path = os.path.join(log_dir, os.listdir(log_dir)[0], 'Tier2', 'WOA')
    check_output_exists(output_path)
    check_conversion(output_path)
