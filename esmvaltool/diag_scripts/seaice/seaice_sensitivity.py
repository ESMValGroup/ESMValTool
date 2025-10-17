"""Diagnostic that shows the sensitivity of sea ice area to global warming."""

import logging
from pathlib import Path

import iris
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import Normalize
from scipy.stats import linregress

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    run_diagnostic,
)

logger = logging.getLogger(Path(__file__).stem)


# This is stolen from the example in AutoAssess _plot_mo_metrics.py
def get_provenance_record(cfg, caption):
    """Create a provenance record describing the diagnostic data and plot."""
    filenames = [item["filename"] for item in cfg["input_data"].values()]

    record = {
        "caption": caption,
        "plot_type": "other",
        "authors": [
            "sellar_alistair",
            "parsons_naomi",
        ],
        "ancestors": filenames,
    }

    return record


def create_category_dict(cfg):
    """Create a structured dictionary for adding values to later on."""
    # Create blank dictionary of correct structure
    category_dict = {
        'models': {},
        'tasa_obs': {},
        'siconc_obs': {},
    }

    # Initialize the models as a set to avoid duplication
    models_set = set()

    # Read the data from the config object
    input_data = cfg["input_data"].values()

    # Iterate over the datasets to add to the dictionary
    for input in input_data:

        # Check for tasa observations
        if input['variable_group'] == 'tasa_obs':
            category_dict['tasa_obs'][input['dataset']] = {}

        # Check for siconc observations
        elif input['variable_group'] == 'siconc_obs':
            category_dict['siconc_obs'][input['dataset']] = {}

        # Everything else should be a model
        else:
            models_set.add(input['dataset'])

    # Add the models to the dictionary
    for model in models_set:
        category_dict['models'][model] = {}

    return category_dict


def fetch_cube(dataset, variable, cfg):
    """Fetch a data cube for a dataset and variable using info from the config."""
    # Read the data from the config object
    input_data = cfg["input_data"].values()

    # Find the correct filepath for the dataset
    for input in input_data:
        if input['dataset'] == dataset:
            # Also check the variable matches as models have two entries
            # Only matching the first three letters to avoid issues with sic vs siconc
            if input['short_name'][:3] == variable[:3]:
                filepath = input['filename']
                break

    # Load the cube using iris
    cube = iris.load_cube(filepath, variable)
    return cube


def calculate_annual_trend(cube):
    """Calculate the linear trend of a cube over time using scipy.stats.linregres."""
    # Depending on preprocessor, coord may be 'year' or 'time'
    if 'year' in cube.coords():
        years = cube.coord("year").points
    else:
        years = cube.coord("time").points

    # slope, intercept, rvalue, pvalue, stderr = linregress(independent, dependent)
    trend = linregress(years, cube.data)

    # Only the slope is needed in this code
    return trend.slope


def calculate_direct_stats(dataset, cfg):
    """Calculate the direct sensitivity of siconc to tas for a given dataset."""
    # Fetch the required cubes
    siconc_cube = fetch_cube(dataset, 'siconc', cfg)
    tas_cube = fetch_cube(dataset, 'tas', cfg)

    # Calculate direct regression (tas as independent)
    direct_sensitivity = linregress(tas_cube.data, siconc_cube.data)

    # direct_sensitivity = slope, intercept, rvalue, pvalue, stderr
    return direct_sensitivity


def write_values_to_dict(data_dict, cfg):
    """Calculate and write values to the structured dictionary."""
    # Calculate all the values for the models
    for model_dataset in data_dict['models'].keys():
        # Calculate annual tas trend
        tas_cube = fetch_cube(model_dataset, 'tas', cfg)
        ann_tas_trend = calculate_annual_trend(tas_cube)
        data_dict['models'][model_dataset]['annual_tas_trend'] = ann_tas_trend

        # Calculate annual siconc trend
        siconc_cube = fetch_cube(model_dataset, 'siconc', cfg)
        ann_siconc_trend = calculate_annual_trend(siconc_cube)
        data_dict['models'][model_dataset]['annual_siconc_trend'] = ann_siconc_trend

        # Calculate direct sensitivity of siconc to tas
        direct_sensitivity = calculate_direct_stats(model_dataset, cfg)
        data_dict['models'][model_dataset]['direct_sensitivity'] = direct_sensitivity.slope
        data_dict['models'][model_dataset]['direct_r_value'] = direct_sensitivity.rvalue
        data_dict['models'][model_dataset]['direct_p_value'] = direct_sensitivity.pvalue

    # Calculate just the tasa trend for the tasa observations
    for obs_dataset in data_dict['tasa_obs'].keys():
        # Calculate annual tas trend
        tasa_cube = fetch_cube(obs_dataset, 'tasa', cfg)
        ann_tasa_trend = calculate_annual_trend(tasa_cube)
        # Still labelling in final dictionary as tas for consistency with models
        data_dict['tasa_obs'][obs_dataset]['annual_tas_trend'] = ann_tasa_trend

    # Calculate just the siconc trend for the siconc observations
    for obs_dataset in data_dict['siconc_obs'].keys():
        # Calculate annual siconc trend
        siconc_cube = fetch_cube(obs_dataset, 'siconc', cfg)
        ann_siconc_trend = calculate_annual_trend(siconc_cube)
        data_dict['siconc_obs'][obs_dataset]['annual_siconc_trend'] = ann_siconc_trend

    return data_dict


def write_dictionary_to_csv(sub_dict, filename, cfg):
    """
    Output a section of data dictionary to a csv file using Pandas.

    Only sections of the dictionary should be written at a time as otherwise
    the structure is too complex to easily convert to a DataFrame.
    """
    # Create the csv filepath using info from the config
    csv_filepath = f"{cfg['work_dir']}/{filename}.csv"

    # Write the data to a csv file (via a Pandas DataFrame)
    dataframe = pd.DataFrame.from_dict(sub_dict, orient="index")
    dataframe.to_csv(csv_filepath)
    logger.info("Wrote data to %s", csv_filepath)

    # Create a provenance record for the csv file
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(
            csv_filepath,
            get_provenance_record(cfg, "Annual (not decadal) figures"),
        )


def main(cfg):
    # Create the structured dictionary
    data_dict = create_category_dict(cfg)

    # Calculate and write values to the dictionary
    data_dict = write_values_to_dict(data_dict, cfg)

    # Write the model and obs dictionaries to csv files
    write_dictionary_to_csv(data_dict['models'], 'models_values', cfg)
    write_dictionary_to_csv(data_dict['tasa_obs'], 'tasa_obs_values', cfg)
    write_dictionary_to_csv(data_dict['siconc_obs'], 'siconc_obs_values', cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
