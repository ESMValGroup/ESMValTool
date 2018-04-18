"""Python example diagnostic."""
import logging
import inspect, os
import sys

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import yaml

logger = logging.getLogger(__name__)


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
    suite_loc = os.path.join(cfg['work_dir'], 'u-ab123')
    os.makedirs(suite_loc)
    suite_data = os.path.join(suite_loc, 'stratosphere')
    os.makedirs(suite_data)
    tmp_dir = os.path.join(cfg['work_dir'], 'tmp')
    ancil_dir = os.path.join(cfg['work_dir'], 'ancil')
    os.makedirs(tmp_dir)
    os.makedirs(ancil_dir)

    files_list = []
    for variable_name, filenames in input_files.items():
        logger.info("Processing variable %s", variable_name)
        for filename, attributes in filenames.items():
            files_list.append(filename)

    cubelist = iris.load(files_list)

    cubes_list_path = os.path.join(suite_data, 'cubeList.nc')
    iris.save(cubelist, cubes_list_path)
 
    cwd = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    command_call = 'python ' + os.path.join(cwd, 'autoassess_source/autoassess/run_area.py')
    args = {}
    args['--area'] = cfg['area']
    args['--suite-id1'] = cfg['suite1']
    args['--suite-id2'] = cfg['suite2']
    args['--start-date'] = cfg['start']
    args['--end-date'] = cfg['end']
    args['--obs-dir'] = cfg['obs_data_dir']
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
