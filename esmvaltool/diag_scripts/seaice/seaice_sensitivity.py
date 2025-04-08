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
from numpy.polynomial.polynomial import polyfit  # TODO: is this OK?


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


def calculate_sensitivity(cfg):
    """
    Calculates a first order best fit gradient for change in sea ice area
    over change in global mean temperature
    """
    # Get data from the configuration object
    input_data = cfg['input_data'].values()

    # Load the preprocessed cubes
    tas_cube = extract_cube(input_data, 'avg_ann_global_temp')
    si_cube = extract_cube(input_data, 'arctic_sea_ice')

    # Use Numpy polyfit to calculate the gradient
    x = tas_cube.data
    y = si_cube.data
    # c isn't used but it felt wrong to use _
    c, m = polyfit(x, y, 1)

    # Print the gradient to log for now
    print('-' * 40)
    print(f'{m:.2e} {si_cube.units} per {tas_cube.units}')


def main(cfg):
    calculate_sensitivity(cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
