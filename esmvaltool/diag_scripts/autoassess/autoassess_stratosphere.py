"""autoassess stratosphere diagnostic."""
import logging
import inspect, os
import sys

import iris
import yaml

logger = logging.getLogger(__name__)

# Diagnostic that takes two models (control_model and exp_model
# and observational data (ERA-Interim and MERRA)

def get_cfg():
    """Read diagnostic script configuration from settings.yml."""
    settings_file = sys.argv[1]
    with open(settings_file) as file:
        cfg = yaml.safe_load(file)
    return cfg


def get_input_files(cfg, index=0):
    """Get a dictionary with input files from metadata.yml files."""
    metadata_file = cfg['input_files'][index]
    with open(metadata_file) as file:
        metadata = yaml.safe_load(file)
    return metadata


def main():

    cfg = get_cfg()
    logger.setLevel(cfg['log_level'].upper())

    input_files = get_input_files(cfg)
    os.makedirs(cfg['plot_dir'])
    os.makedirs(cfg['work_dir'])
    suite_loc_m1 = os.path.join(cfg['work_dir'], cfg['control_model'])
    os.makedirs(suite_loc_m1)
    suite_loc_m2 = os.path.join(cfg['work_dir'], cfg['exp_model'])
    os.makedirs(suite_loc_m2)
    suite_data_m1 = os.path.join(suite_loc_m1, 'stratosphere')
    os.makedirs(suite_data_m1)
    suite_data_m2 = os.path.join(suite_loc_m2, 'stratosphere')
    os.makedirs(suite_data_m2)
    tmp_dir = os.path.join(cfg['work_dir'], 'tmp')
    ancil_dir = os.path.join(cfg['work_dir'], 'ancil')
    os.makedirs(tmp_dir)
    os.makedirs(ancil_dir)

    files_list_m1 = []
    files_list_m2 = []
    obs_list = []
    for variable_name, filenames in input_files.items():
        logger.info("Processing variable %s", variable_name)
        # get model data files
        for filename, attributes in filenames.items():
            if os.path.basename(filename).split('_')[1] == cfg['control_model']:
                files_list_m1.append(filename)
            elif os.path.basename(filename).split('_')[1] == cfg['exp_model']:
                files_list_m2.append(filename)
        # get obs files
        for filename, attributes in filenames.items():
            if os.path.basename(filename).split('_')[0] == 'OBS':
                obs_list.append(filename)

    # load cubelists
    cubelist_m1 = iris.load(files_list_m1)
    cubelist_m2 = iris.load(files_list_m2)

    # save to congragated files
    cubes_list_path_m1 = os.path.join(suite_data_m1, 'cubeList.nc')
    iris.save(cubelist_m1, cubes_list_path_m1)
    cubes_list_path_m2 = os.path.join(suite_data_m2, 'cubeList.nc')
    iris.save(cubelist_m2, cubes_list_path_m2)
 
    cwd = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    command_call = 'python ' + os.path.join(cwd, 'autoassess_source/run_area.py')
    args = {}
    args['--area'] = cfg['area']
    args['--suite-id1'] = cfg['control_model']
    args['--suite-id2'] = cfg['exp_model']
    args['--start-date'] = cfg['start']
    args['--end-date'] = cfg['end']
    args['--obs-dir'] = os.path.dirname(files_list_m1[0])
    if cfg['obs_models'] is not None:
        group_files = [[ofile for ofile in obs_list if os.path.basename(ofile).split('_')[1] == obs] for obs in cfg['obs_models']]
        for obs_file_group in group_files:
            cubes_list_obs = iris.load(obs_file_group)
            cubes_list_obs_path = os.path.join(os.path.dirname(obs_file_group[0]), os.path.basename(obs_file_group[0]).split('_')[1] + '_tropical_area_avg.nc')
            iris.save(cubes_list_obs, cubes_list_obs_path)
    args['--out-dir'] = cfg['plot_dir']
    args['--data-dir'] = cfg['work_dir']
    args['--tmp-dir'] = tmp_dir
    args['--ancil-dir'] = ancil_dir
    args_collection = [key + ' ' + args[key] for key in args.keys()]
    sys_call = command_call + ' ' + ' '.join(args_collection)
    logger.info(sys_call)
    os.system(sys_call)

if __name__ == '__main__':
    iris.FUTURE.netcdf_promote = True
    logging.basicConfig(format="%(asctime)s [%(process)d] %(levelname)-8s "
                        "%(name)s,%(lineno)s\t%(message)s")
    main()
