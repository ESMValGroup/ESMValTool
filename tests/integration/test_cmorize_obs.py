"""Tests for the module :mod:`esmvaltool.cmorizers.obs.cmorize_obs`."""

import contextlib
import os
import sys

import iris
import numpy as np
import yaml
from cf_units import Unit

from esmvaltool.cmorizers.obs.cmorize_obs import main as run


@contextlib.contextmanager
def keep_cwd():
    """
    Use a context manager since the cmorizer enters
    and stays in the cmorization dir, risking to write
    test outputs away from test-reports.
    """
    curr_path = os.getcwd()
    try:
        yield
    finally:
        os.chdir(curr_path)


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


def _create_sample_cube():
    """Create a quick CMOR-compliant sample cube."""
    coord_sys = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
    cube_data = np.ones((2, 3, 2, 2))
    cube_data[1, 1, 1, 1] = 22.
    time = iris.coords.DimCoord([15, 45],
                                standard_name='time',
                                bounds=[[1., 30.], [30., 60.]],
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
    gen_cube = _create_sample_cube()
    t_path = os.path.join(data_path, "woa13_decav81B0_t00_01.nc")
    # correct var names
    gen_cube.var_name = "t_an"
    iris.save(gen_cube, t_path)
    s_path = os.path.join(data_path, "woa13_decav81B0_s00_01.nc")
    gen_cube.var_name = "s_an"
    iris.save(gen_cube, s_path)
    # incorrect var names
    o2_path = os.path.join(data_path, "woa13_all_o00_01.nc")
    gen_cube.var_name = "o2"
    iris.save(gen_cube, o2_path)
    no3_path = os.path.join(data_path, "woa13_all_n00_01.nc")
    gen_cube.var_name = "no3"
    iris.save(gen_cube, no3_path)
    po4_path = os.path.join(data_path, "woa13_all_p00_01.nc")
    gen_cube.var_name = "po4"
    iris.save(gen_cube, po4_path)
    si_path = os.path.join(data_path, "woa13_all_i00_01.nc")
    gen_cube.var_name = "si"
    iris.save(gen_cube, si_path)


def check_log_file(log_file, no_data=False):
    """Check the cmorization log file."""
    with open(log_file, 'r') as log:
        if no_data:
            msg = "Could not find raw data WOA"
        else:
            msg = "Fixing data"
        assert any(msg in line for line in log)


def check_output_exists(output_path):
    """Check if cmorizer outputted."""
    # eg Tier2/WOA/OBS_WOA_clim_2013v2_Omon_thetao_200001-200002.nc
    output_files = os.listdir(output_path)
    # ['OBS_WOA_clim_2013v2_Omon_thetao_200001-200002.nc',
    # 'OBS_WOA_clim_2013v2_Omon_so_200001-200002.nc']
    assert len(output_files) == 2
    assert 'OBS_WOA_clim' in output_files[0]
    assert 'thetao' in [s.split("_")[5] for s in output_files]
    assert 'so' in [s.split("_")[5] for s in output_files]


def check_conversion(output_path):
    """Check basic cmorization."""
    cube = iris.load_cube(os.path.join(output_path,
                                       os.listdir(output_path)[0]))
    assert cube.coord("time").units == Unit('days since 1950-1-1 00:00:00',
                                            calendar='gregorian')
    assert cube.coord("latitude").units == 'degrees'


@contextlib.contextmanager
def arguments(*args):
    backup = sys.argv
    sys.argv = list(args)
    yield
    sys.argv = backup


def test_cmorize_obs_woa_no_data(tmp_path):
    """Test for example run of cmorize_obs command."""

    config_user_file = write_config_user_file(tmp_path)
    os.makedirs(os.path.join(tmp_path, 'raw_stuff', 'Tier2'))
    with keep_cwd():
        with arguments(
                'cmorize_obs',
                '-c',
                config_user_file,
                '-o',
                'WOA',
        ):
            run()

    log_dir = os.path.join(tmp_path, 'output_dir')
    log_file = os.path.join(log_dir,
                            os.listdir(log_dir)[0], 'run', 'main_log.txt')
    check_log_file(log_file, no_data=True)


def test_cmorize_obs_woa_data(tmp_path):
    """Test for example run of cmorize_obs command."""

    config_user_file = write_config_user_file(tmp_path)
    os.makedirs(os.path.join(tmp_path, 'raw_stuff'))
    data_path = os.path.join(tmp_path, 'raw_stuff', 'Tier2', 'WOA')
    os.makedirs(data_path)
    put_dummy_data(data_path)
    with keep_cwd():
        with arguments(
                'cmorize_obs',
                '-c',
                config_user_file,
                '-o',
                'WOA',
        ):
            run()

    log_dir = os.path.join(tmp_path, 'output_dir')
    log_file = os.path.join(log_dir,
                            os.listdir(log_dir)[0], 'run', 'main_log.txt')
    check_log_file(log_file, no_data=False)
    output_path = os.path.join(log_dir, os.listdir(log_dir)[0], 'Tier2', 'WOA')
    check_output_exists(output_path)
    check_conversion(output_path)
