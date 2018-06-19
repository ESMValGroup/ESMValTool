"""Calculate interannual variability diagnostic."""
import logging
import os

import iris

import numpy as np

from esmvaltool.diag_scripts.shared import run_diagnostic

logger = logging.getLogger(os.path.basename(__file__))

def main(cfg):
    logger.setLevel(cfg['log_level'].upper())
    if not os.path.exists(cfg['plot_dir']):
        os.makedirs(cfg['plot_dir'])
    if not os.path.exists(cfg['work_dir']):
        os.makedirs(cfg['work_dir'])
    
    """Compute IAV for each input model."""
    for filename, attributes in cfg['input_data'].items():
        logger.info("Processing variable %s from model %s",
                    attributes['standard_name'], attributes['model'])

        logger.debug("Loading %s", filename)
        cube = iris.load_cube(filename)

        logger.debug("Cut first 100 years and check min length > 500 years")

        
        logger.debug("Calculate annual or seasonal means")
        time_avg = cfg['time_avg']
        time = cube.time
        new_start_yr = cfg['start_year'] + 100
        end_yr = cfg['end_year']
        if cfg['periodlength']:
            length_of_period = cfg['periodlength']
            nr_periods = np.floor(time.shape - 12 * 100) / length_of_period * 12
            yr_possible = (time.shape - 12 * 100) / (length_of_period * 12)
            rest = yr_possible - nr_periods
            start_years = range(new_start_year, end_yr, length_of_period)
        if time_avg == "annualclim":
            # calculate annual means over periods
            new_cube = cube.aggregated_by('year', iris.analysis.MEAN)
        elif time_avg == "seasonalclim":
            # remove seasonal cycle

            # calculate seasonal means over periods
            iris.coord_categorisation.add_season(cube, 'time',
                                                 name = 'clim_season')
            iris.coord_categorisation.add_season_year(cube, 'time',
                                                      name = 'season_year')
            new_cube = cube.aggregated_by(['clim_season', 'season_year'],
                                          iris.analysis.MEAN)

        
        name = os.path.splitext(os.path.basename(filename))[0] + '_IAV'
        if cfg['write_netcdf']:
            path = os.path.join(
                cfg['work_dir'],
                name + '.nc',
            )
            logger.debug("Saving analysis results to %s", path)
            iris.save(cube, target = path)

if __name__ == '__main__':

    with run_diagnostic() as config:
	main(config)

