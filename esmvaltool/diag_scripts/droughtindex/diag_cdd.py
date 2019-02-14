"""Diagnostic to select grid points within a shapefile."""
import logging
import os
from copy import deepcopy

import iris
import numpy as np

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared.plot import quickplot
from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger, get_diagnostic_filename)

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(cfg, basename, caption, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'statistics': ['other'],
        'domains': ['global'],
        'authors': ['berg_pe'],
        'references': ['acknow_project'],
        'ancestors': ancestor_files,
    }
    diagnostic_file = get_diagnostic_filename(basename, cfg)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, record)
    return record


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
            path = os.path.join(
                cfg['work_dir'],
                name + '_drymax.nc',
            )
            output_basename = name + '_drymax'
            iris.save(drymaxcube, target=path)
            provenance_record = get_provenance_record(
                cfg, output_basename, dmcap, ancestor_files=[path])
            path = os.path.join(
                cfg['work_dir'],
                name + '_dryfreq.nc',
            )
            output_basename = name + '_dryfreq'
            iris.save(fqthcube, target=path)
            provenance_record = get_provenance_record(
                cfg, output_basename, fqcap, ancestor_files=[path])
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
        for t in range(1, cube.data.shape[0]):
            cube.data[t, :, :] = (
                (precip[t, :, :] + cube.data[t - 1, :, :]) * precip[t, :, :])
        dif = cube.data[0:-1, :, :] - cube.data[1:cube.data.shape[0], :, :]
        wh = np.where(dif != cube.data[0:-1])
        cube.data[wh] = 0
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
        fqthlong_name = ('The number of consecutive dry day periods ' +
                         'of at least ' + str(cfg['frlim']) + ' days ' +
                         'with precipitation below ' + str(cfg['plim']) +
                         ' mm each day.')
        fqthcube.long_name = fqthlong_name
    return drymaxcube, dmlong_name, fqthcube, fqthlong_name


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
