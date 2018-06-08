"""autoassess stratosphere diagnostic."""
import os
import logging
import inspect
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


def main():

    cfg = get_cfg()
    logger.setLevel(cfg['log_level'].upper())

    control_model = cfg['control_model']
    exp_model = cfg['exp_model']

    vsloc = exp_model + '_vs_' + control_model
    file_exp = os.path.join(
        os.path.dirname(os.path.dirname(cfg['plot_dir'])), cfg['diag_tag'],
        cfg['diag_name'], vsloc, 'stratosphere', exp_model, 'metrics.csv')
    file_ref = os.path.join(
        os.path.dirname(os.path.dirname(cfg['plot_dir'])), cfg['diag_tag'],
        cfg['diag_name'], vsloc, 'stratosphere', control_model, 'metrics.csv')

    cwd = os.path.dirname(
        os.path.abspath(inspect.getfile(inspect.currentframe())))
    plotter_script = os.path.join(cwd, 'autoassess_source/plot_norm_ac.py')
    os.system('chmod +x ' + plotter_script)
    command_call = plotter_script
    args = {}
    args['--exp'] = exp_model
    args['--ref'] = control_model
    args['--plot'] = os.path.join(cfg['plot_dir'], cfg['plot_name'] + '.png')
    args['--title'] = cfg['plot_title']
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
