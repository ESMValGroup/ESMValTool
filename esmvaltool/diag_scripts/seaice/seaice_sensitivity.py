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

logger = logging.getLogger(Path(__file__).stem)


def extract_value(data, preprocessor):
    """
    Selects records from data based on the preprocessor,
    asserts that there is only one matching record,
    then returns that record's value
    """
    # TODO: assert that there is only one value in the cube

    # loading any data that uses the preprocessor
    selection = select_metadata(data, preprocessor=preprocessor)

    # ensure there was only one file
    assert len(selection) == 1, f'Too many matching files found for {preprocessor}'

    # loading the cube
    cube = iris.load_cube(selection[0]['filename'])

    # accessing the value
    value = cube.data

    return value


def calculate_sensitivity(cfg):
    """
    Calculates the siconc difference divided by the tas difference
    """
    # TODO: build in a group by dataset step

    input_data = cfg['input_data'].values()

    # Sea ice for the reference period
    si_ref_value = extract_value(input_data, 'pp_arctic_sea_ice_ref')

    # Sea ice for the test period
    si_test_value = extract_value(input_data, 'pp_arctic_sea_ice_test')

    # Mean temperature for the reference period
    tas_ref_value = extract_value(input_data, 'pp_avg_global_temp_ref')

    # Mean temperature for the test period
    tas_test_value = extract_value(input_data, 'pp_avg_global_temp_test')

    # calculating temperature difference
    tas_diff = tas_test_value - tas_ref_value

    # calculating sea ice difference
    si_diff = si_test_value - si_ref_value

    # calculating the ratio
    sensitivity = si_diff / tas_diff

    # Printing the value for now instead of plotting
    print(f'Loss per degree of warming: {sensitivity}')
    print(f'written in standard index form: {sensitivity:.5e}')


def main(cfg):
    calculate_sensitivity(cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
