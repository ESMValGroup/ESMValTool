"""autoassess land-surface diagnostic."""
import os
import logging
import inspect
import sys
import subprocess
import iris
from esmvaltool.diag_scripts.shared import run_diagnostic

logger = logging.getLogger(__name__)

# Diagnostic that takes two datasets (control_model and exp_model
# and observational data


def _make_tmp_dir(cfg):
    """Make the tmp and ancil dirs"""
    tmp_dir = os.path.join(cfg['work_dir'], 'tmp')
    ancil_dir = os.path.join(cfg['work_dir'], 'ancil')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    if not os.path.exists(ancil_dir):
        os.makedirs(ancil_dir)
    return tmp_dir, ancil_dir


def _make_main_dirs(cfg):
    """Create main dirs to hold analysis"""
    suite_loc_m1 = os.path.join(cfg['work_dir'], cfg['control_model'])
    if not os.path.exists(suite_loc_m1):
        os.makedirs(suite_loc_m1)
    suite_loc_m2 = os.path.join(cfg['work_dir'], cfg['exp_model'])
    if not os.path.exists(suite_loc_m2):
        os.makedirs(suite_loc_m2)
    obs_loc = os.path.join(cfg['work_dir'], 'OBS')
    if not os.path.exists(obs_loc):
        os.makedirs(obs_loc)
    return suite_loc_m1, suite_loc_m2, obs_loc


def _make_concatenated_data_dirs(suite_loc_m1, suite_loc_m2):
    """Create dirs to hold cubeList files"""
    suite_data_m1 = os.path.join(suite_loc_m1, 'land_surface')
    if not os.path.exists(suite_data_m1):
        os.makedirs(suite_data_m1)
    suite_data_m2 = os.path.join(suite_loc_m2, 'land_surface')
    if not os.path.exists(suite_data_m2):
        os.makedirs(suite_data_m2)

    # create supermeans directory: land_surface_supermeans
    sup_data_m1 = os.path.join(suite_loc_m1, 'land_surface_supermeans')
    if not os.path.exists(sup_data_m1):
        os.makedirs(sup_data_m1)
    sup_data_m2 = os.path.join(suite_loc_m2, 'land_surface_supermeans')
    if not os.path.exists(sup_data_m2):
        os.makedirs(sup_data_m2)
    return suite_data_m1, suite_data_m2, sup_data_m1, sup_data_m2


def _get_filelists(cfg, obs_type):
    """Put files in lists and return them"""
    files_list_m1 = []
    files_list_m2 = []
    obs_list = []
    for filename, attributes in cfg['input_data'].items():
        base_file = os.path.basename(filename)
        fullpath_file = filename
        if base_file.split('_')[1] == cfg['control_model']:
            files_list_m1.append(fullpath_file)
            if 'fx_files' in attributes:
                files_list_m1.append(attributes['fx_files'][cfg['fx']])
        elif base_file.split('_')[1] == cfg['exp_model']:
            files_list_m2.append(fullpath_file)
            if 'fx_files' in attributes:
                files_list_m2.append(attributes['fx_files'][cfg['fx']])
        elif base_file.split('_')[0] == obs_type:
            obs_list.append(fullpath_file)
    return files_list_m1, files_list_m2, obs_list


def _process_obs(cfg, obs_list, obs_loc):
    """Gather obs files and save them applying specific cases"""
    group_files = [[
        ofile for ofile in obs_list
        if os.path.basename(ofile).split('_')[1] == obs
    ] for obs in cfg['obs_models']]
    for obs_file_group in group_files:
        cubes_list_obs = iris.load(obs_file_group)
        # force add a long_name; supermeans uses extract_strict
        # and for derived vars there is only
        # invalid_standard_name which is an attribute and shit
        for cube in cubes_list_obs:
            if 'invalid_standard_name' in cube.attributes:
                cube.long_name = cube.attributes['invalid_standard_name']
        iris.save(cubes_list_obs, os.path.join(obs_loc, 'cubeList.nc'))


def _process_ctrl_exp_data(files_list_m1, files_list_m2, suite_data_m1,
                           suite_data_m2, sup_data_m1, sup_data_m2):
    """Create and save concatenated cubes for ctrl and exp"""
    cubelist_m1 = iris.load(files_list_m1)
    cubelist_m2 = iris.load(files_list_m2)

    # save to congragated files; save twice for supermeans as well
    cubes_list_path_m1 = os.path.join(suite_data_m1, 'cubeList.nc')
    cubes_list_path_sup_m1 = os.path.join(sup_data_m1, 'cubeList.nc')
    # force add a long_name; supermeans uses extract_strict
    # and for derived vars there is only
    # invalid_standard_name which is an attribute and shit
    for cube in cubelist_m1:
        if 'invalid_standard_name' in cube.attributes:
            cube.long_name = cube.attributes['invalid_standard_name']
    iris.save(cubelist_m1, cubes_list_path_m1)
    iris.save(cubelist_m1, cubes_list_path_sup_m1)
    cubes_list_path_m2 = os.path.join(suite_data_m2, 'cubeList.nc')
    cubes_list_path_sup_m2 = os.path.join(sup_data_m2, 'cubeList.nc')
    for cube in cubelist_m2:
        if 'invalid_standard_name' in cube.attributes:
            cube.long_name = cube.attributes['invalid_standard_name']
    iris.save(cubelist_m2, cubes_list_path_m2)
    iris.save(cubelist_m2, cubes_list_path_sup_m2)
    return (cubes_list_path_m1, cubes_list_path_m2, cubes_list_path_sup_m1,
            cubes_list_path_sup_m2)


def main(cfg):
    """Execute the stratosphere area"""
    logger.setLevel(cfg['log_level'].upper())

    # set the main dirs
    suite_loc_m1, suite_loc_m2, obs_loc = _make_main_dirs(cfg)

    # set the concatenated data dirs
    suite_data_m1, suite_data_m2, sup_data_m1, sup_data_m2 = \
        _make_concatenated_data_dirs(suite_loc_m1, suite_loc_m2)

    # create the ancil and tp dirs
    tmp_dir, ancil_dir = _make_tmp_dir(cfg)

    # get files lists
    files_list_m1, files_list_m2, obs_list = _get_filelists(
        cfg, cfg['obs_type'])

    # spell out the files used
    logger.info("Files for control model: %s", files_list_m1)
    logger.info("Files for exp model: %s", files_list_m2)
    logger.info("Files for obs model: %s", obs_list)

    # load and save control and exp cubelists
    cubes_ctrl, cubes_exp, cubes_sup_ctrl, cubes_sup_exp = \
        _process_ctrl_exp_data(files_list_m1, files_list_m2,
                               suite_data_m1, suite_data_m2,
                               sup_data_m1, sup_data_m2)

    # print the paths
    logger.info("Saved control data cube: %s", cubes_ctrl)
    logger.info("Saved exp data cube: %s", cubes_exp)
    logger.info("Saved control data cube for supermeans: %s",
                cubes_sup_ctrl)
    logger.info("Saved exp data cube for supermeans: %s", cubes_sup_exp)

    cwd = os.path.dirname(
        os.path.abspath(inspect.getfile(inspect.currentframe())))
    command_call = 'python ' + os.path.join(cwd,
                                            'autoassess_source/run_area.py')
    args = {}
    args['--area'] = cfg['area']
    args['--suite-id1'] = cfg['control_model']
    args['--suite-id2'] = cfg['exp_model']
    args['--start-date'] = cfg['start']
    args['--end-date'] = cfg['end']
    args['--obs-dir'] = obs_loc

    if cfg['obs_models']:
        _process_obs(cfg, obs_list, obs_loc)

    args['--out-dir'] = cfg['plot_dir']
    args['--data-dir'] = cfg['work_dir']
    args['--tmp-dir'] = tmp_dir
    args['--ancil-dir'] = ancil_dir
    args_collection = [key + ' ' + args[key] for key in args.keys()]
    sys_call = command_call + ' ' + ' '.join(args_collection)
    logger.info(sys_call)
    # run the thing
    proc = subprocess.Popen(sys_call, stdout=subprocess.PIPE, shell=True)
    out_proc = proc.communicate()
    ret_code = proc.returncode
    logger.info("Diagnostic output: %s", out_proc[0])
    if int(ret_code) != 0:
        logger.info("Diagnostic has failed!")
        sys.exit(1)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
