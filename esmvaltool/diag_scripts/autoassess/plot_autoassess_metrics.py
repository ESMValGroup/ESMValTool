"""Standard MO metrics plotter."""
import os
import logging
import sys

import iris
import yaml
from esmvaltool.diag_scripts.autoassess._plot_mo_metrics import (
    read_model_metrics, read_obs_metrics, plot_nac)

logger = logging.getLogger(__name__)

# Diagnostic that takes two datasets (control_model and exp_model
# and observational data (ERA-Interim and MERRA);
# plotting OBS is not yet supported; it will be, hold your horses


def get_cfg():
    """Read diagnostic script configuration from settings.yml."""
    settings_file = sys.argv[1]
    with open(settings_file) as file:
        cfg = yaml.safe_load(file)
    return cfg


def main():
    """Call the plotting script via command line."""
    cfg = get_cfg()
    logger.setLevel(cfg['log_level'].upper())

    control_model = cfg['control_model']
    exp_model = cfg['exp_model']

    vsloc = exp_model + '_vs_' + control_model
    file_exp = os.path.join(
        os.path.dirname(os.path.dirname(cfg['plot_dir'])), cfg['diag_tag'],
        cfg['diag_name'], vsloc, cfg['area'], exp_model, 'metrics.csv')
    file_ref = os.path.join(
        os.path.dirname(os.path.dirname(cfg['plot_dir'])), cfg['diag_tag'],
        cfg['diag_name'], vsloc, cfg['area'], control_model, 'metrics.csv')

    plot_title = ' '.join([cfg['area'], control_model, 'vs', exp_model])
    # Read metrics files
    # metrics = read_order_metrics(args.file_ord)
    ref = read_model_metrics(file_ref)
    tests = [read_model_metrics(file_exp)]
    # var = read_model_metrics(args.file_var)
    obs, acc = None, None
    if 'additional_metrics' in cfg:
        # choose the obs file to get the metrics from
        file_obs = os.path.join(
            os.path.dirname(os.path.dirname(cfg['plot_dir'])), cfg['diag_tag'],
            cfg['diag_name'], vsloc, cfg['area'], cfg['error_metric'],
            'metrics.csv')
        (obs, acc) = read_obs_metrics(file_obs)

    # Produce plot
    plot_nac(
        control_model, [exp_model],
        ref,
        tests,
        metrics=None,
        var=None,
        obs=obs,
        acc=acc,
        extend_y=False,
        title=plot_title,
        ofile=os.path.join(cfg['plot_dir'], cfg['plot_name'] + '.png'))


if __name__ == '__main__':
    iris.FUTURE.netcdf_promote = True
    logging.basicConfig(format="%(asctime)s [%(process)d] %(levelname)-8s "
                        "%(name)s,%(lineno)s\t%(message)s")
    main()
