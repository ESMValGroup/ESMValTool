"""autoassess conservation diagnostic."""
import os
import logging
import inspect
import sys
import subprocess
import iris

from esmvaltool.diag_scripts.shared import run_diagnostic

logger = logging.getLogger(__name__)

# Diagnostic that takes two models (control_model and exp_model


def main(cfg):
    """Execute the conservation area"""
    logger.setLevel(cfg['log_level'].upper())
    if not os.path.exists(cfg['plot_dir']):
        os.makedirs(cfg['plot_dir'])
    if not os.path.exists(cfg['work_dir']):
        os.makedirs(cfg['work_dir'])
    suite_loc_m1 = os.path.join(cfg['work_dir'], cfg['control_model'])
    if not os.path.exists(suite_loc_m1):
        os.makedirs(suite_loc_m1)
    suite_loc_m2 = os.path.join(cfg['work_dir'], cfg['exp_model'])
    if not os.path.exists(suite_loc_m2):
        os.makedirs(suite_loc_m2)
    suite_data_m1 = os.path.join(suite_loc_m1, 'conservation')
    if not os.path.exists(suite_data_m1):
        os.makedirs(suite_data_m1)
    suite_data_m2 = os.path.join(suite_loc_m2, 'conservation')
    if not os.path.exists(suite_data_m2):
        os.makedirs(suite_data_m2)
    tmp_dir = os.path.join(cfg['work_dir'], 'tmp')
    ancil_dir = os.path.join(cfg['work_dir'], 'ancil')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    if not os.path.exists(ancil_dir):
        os.makedirs(ancil_dir)
    # make the masks and conservation dirs
    masks_dir = os.path.join(ancil_dir, 'masks')
    cons_dir = os.path.join(ancil_dir, 'conservation')
    if not os.path.exists(masks_dir):
        os.makedirs(masks_dir)
    if not os.path.exists(cons_dir):
        os.makedirs(cons_dir)
    # grab the masks and place them in the called dirs
    # /masks/qrparm.landfrac_n240.0.pp
    # /conservation/glacialmask_n240.0.pp
    # NEED TO FIGURE OUT WHERE TO STORE
    # THESE TYPE OF FILES IN BASE DIR(S)
    # for now I am moving some dummy files
    base_mask = cfg['landocean_masks']
    cp1 = 'cp ' + os.path.join(base_mask, 'qrparm.landfrac_n240.0.pp') \
          + ' ' + masks_dir
    cp2 = 'cp ' + os.path.join(base_mask, 'glacialmask_n240.0.pp') \
          + ' ' + cons_dir
    proc = subprocess.Popen(cp1, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    proc = subprocess.Popen(cp2, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()

    files_list_m1 = []
    files_list_m2 = []
    obs_list = []
    for filename, attributes in cfg['input_data'].items():
        base_file = os.path.basename(filename)
        fullpath_file = filename
        if base_file.split('_')[1] == cfg[
                'control_model']:
            files_list_m1.append(fullpath_file)
        elif base_file.split('_')[1] == cfg['exp_model']:
            files_list_m2.append(fullpath_file)
        elif base_file.split('_')[0] == 'OBS':
            obs_list.append(fullpath_file)

    # spell out the files used
    logger.info("Files for control model: %s", files_list_m1)
    logger.info("Files for exp model: %s", files_list_m2)
    logger.info("Files for obs model: %s", obs_list)

    # load cubelists
    cubelist_m1 = iris.load(files_list_m1)
    cubelist_m2 = iris.load(files_list_m2)

    # save to congragated files
    cubes_list_path_m1 = os.path.join(suite_data_m1, 'cubeList.nc')
    iris.save(cubelist_m1, cubes_list_path_m1)
    cubes_list_path_m2 = os.path.join(suite_data_m2, 'cubeList.nc')
    iris.save(cubelist_m2, cubes_list_path_m2)
    logger.info("Saved control data cube: %s", cubes_list_path_m1)
    logger.info("Saved exp data cube: %s", cubes_list_path_m2)

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
    args['--obs-dir'] = os.path.dirname(files_list_m1[0])
    if cfg['obs_models'] is not None:
        group_files = [[
            ofile for ofile in obs_list
            if os.path.basename(ofile).split('_')[1] == obs
        ] for obs in cfg['obs_models']]
        for obs_file_group in group_files:
            if obs_file_group:
                cubes_list_obs = iris.load(obs_file_group)
                cubes_list_obs_path = os.path.join(
                    os.path.dirname(obs_file_group[0]),
                    os.path.basename(obs_file_group[0]).split('_')[1] +
                    '_tropical_area_avg.nc')
                iris.save(cubes_list_obs, cubes_list_obs_path)
    args['--out-dir'] = cfg['plot_dir']
    args['--data-dir'] = cfg['work_dir']
    args['--tmp-dir'] = tmp_dir
    args['--ancil-dir'] = ancil_dir
    args_collection = [key + ' ' + args[key] for key in args.keys()]
    sys_call = command_call + ' ' + ' '.join(args_collection)
    logger.info(sys_call)
    # run the thing
    proc = subprocess.Popen(sys_call, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    rc = proc.returncode
    logger.info("Diagnostic output: %s", out)
    if int(rc) != 0:
        logger.info("Diagnostic has failed!")
        sys.exit(1)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
