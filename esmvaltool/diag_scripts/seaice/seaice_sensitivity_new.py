"""Diagnostic that shows the sensitivity of sea ice area to global warming."""

import logging
from pathlib import Path

import iris
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from scipy.stats import linregress

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    run_diagnostic,
    save_figure,
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


def create_dataset_dict(cfg):
    # Initialise a dictionary
    dataset_dict = {}

    tasa_obs = []
    siconc_obs = []

    # Iterate over the data in the cfg
    for section in cfg["input_data"].values():

        # Read relevant facets
        dataset = section["dataset"]
        group = section["variable_group"]

        # Create the dictionary's section for the dataset if it doesn't exist
        dataset_dict.setdefault(dataset, {})

        # Add the models, they are in tas and siconc
        if group == "tas":
            dataset_dict[dataset]["type"] = "model"

            # Add whether or not to label the models
            if section.get("label_dataset"):
                dataset_dict[dataset]["label"] = "to_label"
            else:
                dataset_dict[dataset]["label"] = "unlabelled"

        # List the GMST obs
        elif group == "tasa_obs":
            tasa_obs.append(dataset)

        # List the SIA obs
        elif group == "siconc_obs":
            siconc_obs.append(dataset)

    # Add the pairs of obs
    pairs = [f"{t}_v_{s}" for t in tasa_obs for s in siconc_obs]
    for pair in pairs:
        dataset_dict[pair]["type"] = "multi-obs"

        # Don't label the obs now, but change here later if needed
        dataset_dict[pair]["label"] = "unlabelled"

    return dataset_dict


def create_df_columns(cfg):
    # This feeds through from the recipe in both diagnostics
    data_start = cfg['observations']['data_period']['start_year']
    data_end = cfg['observations']['data_period']['end_year']

    # This is only present in the arctic diagnostic
    if 'observation_period' in cfg['observations']:
        obs_start = cfg['observations']['observation_period']['start_year']
        obs_end = cfg['observations']['observation_period']['end_year']
        if not obs_start == data_start or not obs_end == data_end:
            periods = [f"{obs_start}-{obs_end}", f"{data_start}-{data_end}"]
    else:
        periods = [f"{data_start}-{data_end}"]

    # In future the historical period will extend beyond 2014
    num_periods = len(periods)

    # Sections for types of regression
    regressions = ["gmst_over_time", "sia_over_time", "sia_over_gmst"]
    num_regressions = len(regressions)

    # Things that will be caluclated
    statistics = ["slope", "r_value", "p_value", "std_err_slope"]
    num_stats = len(statistics)

    # Produces 1979-2014 * 12 then 1979-2021 * 12
    period = [p for p in periods for _ in range(num_regressions * num_stats)]

    # Produces gmst_over_time * 4 then sia_over_time * 4 then sia_over_gmst * 4, two times
    regression = [r for _ in range(num_periods) for r in regressions for _ in range(num_stats)]

    # Produces slope then r_value then p_value then std_err_slope, six times
    statistic = [s for _ in range(num_periods * num_regressions) for s in statistics]

    # Create main column headers
    data_columns = pd.MultiIndex.from_arrays([period, regression, statistic], names=["period", "regression", "statistic"])

    # Add single-level columns for dataset info
    first_columns = pd.Index(["label", "type"])

    # Concatenate columns
    columns = first_columns.append(data_columns)

    return columns


def create_blank_dataframe(cfg, dataset_dict, columns):
    # Create DataFrame from dataset_dict
    df = pd.DataFrame.from_dict(dataset_dict, orient="index")

    # Reindex to ensure all required columns are present
    df = df.reindex(columns=columns, fill_value=np.nan)

    # Fill missing "label" and "type" with empty string
    df[["label", "type"]] = df[["label", "type"]].fillna("")
    return df


if __name__ == "__main__":
    with run_diagnostic() as config:
        print('-------------')
        print(config)
        print('-------------')
        datasets = create_dataset_dict(config)
        print(datasets)
        print('-------------')
        columns = create_df_columns(config)
        print(columns)
        print('--------------')
        print(create_blank_dataframe(config, datasets, columns))
