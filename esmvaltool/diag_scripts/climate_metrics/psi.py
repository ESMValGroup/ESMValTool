#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to temperature variability metric psi (Cox et al., 2018).

Description
-----------
Calculate global temperature variability metric psi following Cox et al.
(2018).

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recips
-------------------------------
window_length : int, optional (default: 55)
    Number of years used for the moving window average.
lag : int, optional(default: 1)
    Lag (in years) for the autocorrelation function.

"""

import logging
import os

import cf_units
import iris
import numpy as np
from scipy import stats

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            save_iris_cube, save_scalar_data,
                                            variables_available)

logger = logging.getLogger(os.path.basename(__file__))


def calculate_psi(cube, cfg):
    """Calculate temperature variability metric psi for a given cube."""
    window_length = cfg.get('window_length', 55)
    lag = cfg.get('lag', 1)
    psi_years = []
    psis = []

    # Moving average
    for yr_idx in range(cube.shape[0] - window_length):
        slc = slice(yr_idx, yr_idx + window_length)
        years = cube.coord('year').points[slc]
        tas = np.copy(cube.data[slc])

        # De-trend data
        reg = stats.linregress(years, tas)
        tas -= reg.slope * years + reg.intercept

        # Autocorrelation
        norm = np.sum(np.square(tas))
        [autocorr] = np.correlate(tas[:-lag], tas[lag:], mode='valid') / norm

        # Psi
        psi_years.append(years[-1])
        psis.append(np.std(tas) / np.sqrt(-np.log(autocorr)))

    # Return new cube
    year_coord = iris.coords.DimCoord(
        np.array(psi_years),
        var_name='year',
        long_name='Year',
        units=cf_units.Unit('year'))
    psi_cube = iris.cube.Cube(
        np.array(psis),
        var_name='psi',
        long_name='Temperature variability metric',
        units=cf_units.Unit('1'),
        dim_coords_and_dims=[(year_coord, 0)],
        attributes={
            'window_length': window_length,
            'lag': lag
        })
    return psi_cube


def main(cfg):
    """Run the diagnostic."""
    input_data = cfg['input_data'].values()

    # Check if tas is available
    if not variables_available(cfg, ['tas']):
        raise ValueError("This diagnostics needs 'tas' variable")

    # Calculate psi for every dataset
    grouped_data = group_metadata(input_data, 'dataset')
    psis = {}
    for (dataset, [data]) in grouped_data.items():
        cube = iris.load_cube(data['filename'])
        cube = cube.aggregated_by('year', iris.analysis.MEAN)
        psi_cube = calculate_psi(cube, cfg)
        psi_cube.attributes.update({
            'project': data['project'],
            'dataset': dataset,
        })

        # Save psi for every dataset
        path = os.path.join(cfg['work_dir'], 'psi_{}.nc'.format(dataset))
        save_iris_cube(psi_cube, path, cfg)

        # Save averaged psi
        psis[dataset] = np.mean(psi_cube.data)

    # Save averaged psis for every dataset in one file
    path = os.path.join(cfg['work_dir'], 'psi.nc')
    cube_atts = {
        'var_name': 'psi',
        'long_name': 'Temperature variability metric',
        'units': cf_units.Unit('1'),
    }
    save_scalar_data(psis, path, cfg, **cube_atts)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
