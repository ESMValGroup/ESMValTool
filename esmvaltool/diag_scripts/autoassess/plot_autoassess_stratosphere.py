"""autoassess stratosphere diagnostic."""
import logging
import inspect, os
import sys

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
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
    plot_dir = cfg['plot_dir']
    control_model = cfg['control_model']
    exp_model = cfg['exp_model']
    #aa_strato/autoassess_strato_test_1/MPI-ESM-MR_vs_MPI-ESM-LR/stratosphere/MPI-ESM-LR/metrics.csv
    vs = exp_model + '_vs_' + control_model
    file_exp = os.path.join(cfg['plot_dir'], cfg['diag_title'], cfg['diag_name'], vs, 'stratosphere', exp_model, 'metrics.csv')
    file_ref = os.path.join(cfg['plot_dir'], cfg['diag_title'], cfg['diag_name'], vs, 'stratosphere', ref_model, 'metrics.csv') 

    cwd = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    command_call = 'python ' + os.path.join(cwd, 'autoassess_source/autoassess/plot_norm_ac.py')
    args = {}
    args['--exp'] = exp_model
    args['--ref'] = control_model
    args['--plot'] = os.path.join(cfg['plot_dir'], cfg['plot_name'])
    args['--title'] = cfg['title']
    args['--file-exp'] = file_exp
    args['--file-ref'] = file_ref
    args_collection = [key + ' ' + args[key] for key in args.keys()]
    sys_call = command_call + ' ' + ' '.join(args_collection)
    logger.info(sys_call)
    os.system(sys_call)

if __name__ == '__main__':
    iris.FUTURE.netcdf_promote = True
    logging.basicConfig(format="%(asctime)s [%(process)d] %(levelname)-8s "
                        "%(name)s,%(lineno)s\t%(message)s")
    main()
