import logging
from pathlib import Path
from pprint import pformat
import iris

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
    select_metadata,
    sorted_metadata,
)
from esmvaltool.diag_scripts.shared.plot import quickplot
from scipy import stats
import numpy as np
import iris.coord_categorisation as cc  # we could avoid using this


logger = logging.getLogger(Path(__file__).stem)


def extract_cube(data, variable_group):
    """
    Selects records from data based on the variable_group,
    asserts that there is only one matching record,
    then returns an iris cube for that record
    """

    # Load any data in the variable_group
    selection = select_metadata(data, variable_group=variable_group)

    # Ensure there is only one file in the list
    assert len(selection) == 1, f'Too many matching files found for {variable_group}'

    # Load the cube, [0] is because selection returns a list
    cube = iris.load_cube(selection[0]['filename'])

    return cube


# Stolen from Ed Blockley
def calc_trend(times, timeseries, slope_only=True):
    """
    Calculate linear trend for a timeseries
    Returns a numpy array containing the linear fit for that trend, which can
    be subtracted from another run to de-trend it.
    """
    # Use SciPy stats to calculate the slope
    slope, intercept, r, p, stderr = stats.linregress(times, timeseries)

    # Return either the slope or the array [0, 1*slope, 2*slope, 3*slope, ...]
    if slope_only:
        return slope
    else:
        return slope * np.arange(len(times))


def calculate_sensitivity(cfg):
    """
    Calculates a first order best fit gradient for change in sea ice area
    over change in global mean temperature
    """
    # Get data from the configuration object
    input_data = cfg['input_data'].values()

    # Load the preprocessed cubes
    si_cube = extract_cube(input_data, 'arctic_sea_ice')
    tas_cube = extract_cube(input_data, 'avg_ann_global_temp')

    # Add years to si cube (tas already has years)
    cc.add_year(si_cube, 'time', name='year')
    # Check that the years match
    assert (si_cube.coord('year').points == tas_cube.coord('year').points).all()
    # Return the years for later use
    years = tas_cube.coord('year').points

    # Calculate trends and sensitivity, both via time as in Ed's code
    si_trend = calc_trend(years, si_cube.data)
    tas_trend = calc_trend(years, tas_cube.data)
    sensitivity = si_trend / tas_trend

    # Still printing rather than returning for now  TODO: return
    print(f'{sensitivity:.5e} {si_cube.units} per {tas_cube.units}')


def main(cfg):
    calculate_sensitivity(cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
