"""A diagnostic that calculates consecutive dry days."""
import logging
import os
from copy import deepcopy

import iris
import numpy as np

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared.plot import quickplot
from esmvaltool.diag_scripts.shared._base import (ProvenanceLogger,
                                                  get_diagnostic_filename)

logger = logging.getLogger(os.path.basename(__file__))


def write_provenance_record(cfg, diagnostic_file, caption, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'statistics': ['other'],
        'domains': ['global'],
        'authors': ['berg_pe'],
        'references': ['acknow_project'],
        'ancestors': ancestor_files,
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, record)


def main(cfg):
    """Calculate drought indices."""
    for filename, attributes in cfg['input_data'].items():
        logger.info("Processing variable %s from dataset %s",
                    attributes['standard_name'], attributes['dataset'])
        logger.debug("Loading %s", filename)
        cube = iris.load_cube(filename)
        drymaxcube, dmcap, fqthcube, fqcap = droughtindex(cube, cfg)
        name = os.path.splitext(os.path.basename(filename))[0]
        # Write to file
        if cfg['write_netcdf']:
            drymax_file = get_diagnostic_filename(name + '_drymax', cfg)
            iris.save(drymaxcube, target=drymax_file)
            write_provenance_record(
                cfg, drymax_file, dmcap, ancestor_files=[filename])

            fqth_file = get_diagnostic_filename(name + '_dryfreq', cfg)
            iris.save(fqthcube, target=fqth_file)
            write_provenance_record(
                cfg, fqth_file, fqcap, ancestor_files=[filename])
        if cfg['write_plots'] and cfg.get('quickplot'):
            path = os.path.join(
                cfg['plot_dir'],
                name + '.' + cfg['output_file_type'],
            )
            logger.debug("Plotting analysis results to %s", path)
            quickplot(drymaxcube, filename=path, **cfg['quickplot'])


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
        dmlong_name = ('The greatest number of consecutive days ' +
                       'per time period with daily precipitation ' +
                       'amount below ' + str(cfg['plim']) + ' mm.')
        drymaxcube.long_name = dmlong_name
        whth = np.where(cube.data > frlim)
        cube.data = cube.data * 0
        cube.data[whth] = 1
        fqthcube = cube.collapsed('time', iris.analysis.SUM)
        fqthlong_name = (
            'The number of consecutive dry day periods ' + 'of at least ' +
            str(cfg['frlim']) + ' days ' + 'with precipitation below ' + str(
                cfg['plim']) + ' mm each day.')
        fqthcube.long_name = fqthlong_name
    return drymaxcube, dmlong_name, fqthcube, fqthlong_name


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
