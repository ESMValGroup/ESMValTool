"""A diagnostic that calculates consecutive dry days."""
import logging
import os
from copy import deepcopy

import cmocean.cm
import iris
import numpy as np

from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    save_data,
    save_figure,
)
from esmvaltool.diag_scripts.shared.plot import global_pcolormesh

logger = logging.getLogger(os.path.basename(__file__))


def save_results(cfg, cube, basename, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    basename = basename + '_' + cube.var_name
    provenance = {
        'caption': cube.long_name.replace('\n', ' '),
        'statistics': ['other'],
        'domains': ['global'],
        'authors': ['berg_peter'],
        'references': ['acknow_project'],
        'ancestors': ancestor_files,
    }
    save_data(basename, provenance, cfg, cube)
    kwargs = dict(cfg.get('plot', {}))
    cmap_name = kwargs.get('cmap', 'rain')
    if cmap_name in cmocean.cm.cmap_d:
        kwargs['cmap'] = cmocean.cm.cmap_d[cmap_name]
    global_pcolormesh(cube, **kwargs)
    save_figure(basename, provenance, cfg)


def main(cfg):
    """Calculate drought indices."""
    for filename, attributes in cfg['input_data'].items():
        logger.info("Processing variable %s from dataset %s",
                    attributes['standard_name'], attributes['dataset'])
        logger.debug("Loading %s", filename)
        cube = iris.load_cube(filename)
        drymaxcube, fqthcube = droughtindex(cube, cfg)
        basename = os.path.splitext(os.path.basename(filename))[0]
        save_results(cfg, drymaxcube, basename, ancestor_files=[filename])
        save_results(cfg, fqthcube, basename, ancestor_files=[filename])


def droughtindex(cube, cfg):
    """Calculate drought stats."""
    if cfg['dryindex'] == 'cdd':
        plim = float(cfg['plim']) / 86400.  # units of kg m-2 s-1
        frlim = float(cfg['frlim'])
        precip = deepcopy(cube.data)
        precip[cube.data < plim] = 1
        precip[cube.data >= plim] = 0
        cube.data[0, :, :] = precip[0, :, :]
        for ttt in range(1, cube.data.shape[0]):
            cube.data[ttt, :, :] = (
                (precip[ttt, :, :] + cube.data[ttt - 1, :, :]) *
                precip[ttt, :, :])
        dif = cube.data[0:-1, :, :] - cube.data[1:cube.data.shape[0], :, :]
        whh = np.where(dif != cube.data[0:-1])
        cube.data[whh] = 0
        # Longest consecutive period
        drymaxcube = cube.collapsed('time', iris.analysis.MAX)
        drymaxcube.long_name = (
            'The greatest number of consecutive days per time period\n'
            'with daily precipitation amount below {plim} mm.').format(**cfg)
        drymaxcube.var_name = 'drymax'
        drymaxcube.standard_name = None
        drymaxcube.units = 'days'

        whth = np.where(cube.data > frlim)
        cube.data = cube.data * 0
        cube.data[whth] = 1
        fqthcube = cube.collapsed('time', iris.analysis.SUM)
        fqthcube.long_name = (
            'The number of consecutive dry day periods of at least {frlim} '
            'days\nwith precipitation below {plim} mm each day.').format(**cfg)
        fqthcube.var_name = 'dryfreq'
        fqthcube.standard_name = None
        fqthcube.units = None

    return drymaxcube, fqthcube


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
