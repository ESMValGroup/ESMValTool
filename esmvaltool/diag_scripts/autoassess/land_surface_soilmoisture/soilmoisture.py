"""Run module for soil moisture metrics."""

import csv
import logging
import os
from collections.abc import Iterable

import iris
import numpy as np
from esmvalcore.preprocessor import regrid

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
)
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger

# Order of seasons must agree with preprocessor definition in recipe
SEASONS = ("djf", "mam", "jja", "son")

logger = logging.getLogger(__name__)


def get_provenance_record(caption, ancestor_filenames):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        "caption": caption,
        "statistics": ["mean"],
        "domains": ["global"],
        "plot_type": "metrics",
        "authors": [
            "rumbold_heather",
            "sellar_alistair",
        ],
        "references": [
            "esacci-soilmoisture",
            "dorigo17rse",
            "gruber19essd",
        ],
        "ancestors": ancestor_filenames,
    }

    return record


def write_metrics(output_dir, metrics, config, ancestors):
    """Write metrics to CSV file.

    The CSV file will have the name ``metrics.csv`` and can be
    used for the normalised metric assessment plot.

    Parameters
    ----------
    output_dir : string
        The full path to the directory in which the CSV file will be written.
    metrics : dictionary of metric,value pairs
        The seasonal data to write.
    config : dictionary
        ESMValTool configuration object
    ancestors : list
        Filenames of input files for provenance
    """
    os.makedirs(output_dir, exist_ok=True)

    file_name = "metrics.csv"
    file_path = os.path.join(output_dir, file_name)

    with open(file_path, "w", newline="", encoding="utf-8") as csvfile:
        csv_writer = csv.writer(csvfile)
        for line in metrics.items():
            csv_writer.writerow(line)

    record_provenance(file_path, config, ancestors)


def volumetric_soil_moisture(model_file, constr_season):
    """
    Read moisture mass content and convert to volumetric soil moisture.

    Parameters
    ----------
    model_file : string
        Path to model file
    constr_season : iris constraint
        Constraint on season to load

    Returns
    -------
    vol_sm1_run : cube
        Volumetric soil moisture
    """
    # Constant: density of water
    rhow = 1000.0

    # m01s08i223
    # CMOR name: mrsos (soil moisture in top model layer kg/m2)
    mrsos = iris.load_cube(
        model_file,
        "mass_content_of_water_in_soil_layer" & constr_season,
    )

    # Set soil moisture to missing data where no soil (moisture=0)
    np.ma.masked_where(mrsos.data == 0, mrsos.data, copy=False)

    # first soil layer depth
    dz1 = mrsos.coord("depth").bounds[0, 1] - mrsos.coord("depth").bounds[0, 0]

    # Calculate the volumetric soil moisture in m3/m3
    # volumetric soil moisture = volume of water / volume of soil layer
    # = depth equivalent of water / thickness of soil layer
    # = (soil moisture content (kg m-2) / water density (kg m-3) )  /
    #      soil layer thickness (m)
    # = mosrs / (rhow * dz1)
    vol_sm1_run = mrsos / (rhow * dz1)
    vol_sm1_run.units = "m3 m-3"
    vol_sm1_run.long_name = "Top layer Soil Moisture"

    return vol_sm1_run


def flatten(list_of_lists):
    """
    Convert list of lists into a flat list, allowing some items to be non-list.

    Parameters
    ----------
    list_of_lists : list
        List containing iterables to flatten, plus optionally non-list items

    Returns
    -------
    flattened : list
        Flattened list with one level of nesting removed
    """
    flattened = []
    for item in list_of_lists:
        if isinstance(item, Iterable) and not isinstance(item, str | bytes):
            flattened.extend(item)
        else:
            flattened.append(item)

    return flattened


def land_sm_top(clim_file, model_file, model_dataset, config, ancestors):
    """
    Calculate median absolute errors for soil mosture against CCI data.

    Parameters
    ----------
    clim_file : string
        Path to observation climatology file
    model_file : list
        Paths to model files
    model_dataset : string
        Name of model dataset
    config : dict
        ESMValTool configuration object
    ancestors : list
        Filenames of input files for provenance

    Returns
    -------
    metrics: dict
        a dictionary of metrics names and values
    """
    # Work through each season
    metrics = {}
    for index, season in enumerate(SEASONS):
        constr_season = iris.Constraint(season_number=index)
        ecv_clim = iris.load_cube(clim_file, constr_season)

        vol_sm1_run = volumetric_soil_moisture(model_file, constr_season)

        # update the coordinate system ECV data with a WGS84 coord system
        # unify coord systems for regridder
        vol_sm1_run.coord(
            "longitude",
        ).coord_system = iris.coord_systems.GeogCS(
            semi_major_axis=6378137.0,
            inverse_flattening=298.257223563,
        )
        vol_sm1_run.coord("latitude").coord_system = iris.coord_systems.GeogCS(
            semi_major_axis=6378137.0,
            inverse_flattening=298.257223563,
        )
        ecv_clim.coord("longitude").coord_system = iris.coord_systems.GeogCS(
            semi_major_axis=6378137.0,
            inverse_flattening=298.257223563,
        )
        ecv_clim.coord("latitude").coord_system = iris.coord_systems.GeogCS(
            semi_major_axis=6378137.0,
            inverse_flattening=298.257223563,
        )

        # Interpolate to the grid of the climatology and form the difference
        vol_sm1_run = regrid(vol_sm1_run, ecv_clim, "linear")

        # mask invalids
        vol_sm1_run.data = np.ma.masked_invalid(vol_sm1_run.data)
        ecv_clim.data = np.ma.masked_invalid(ecv_clim.data)

        # diff the cubes
        dff = vol_sm1_run - ecv_clim

        # save output and populate metric
        caption = f"{model_dataset} minus CCI soil moisture clim for {season}"
        provenance_record = get_provenance_record(caption, ancestors)
        save_data(
            f"soilmoist_diff_{model_dataset}_{season}",
            provenance_record,
            config,
            dff,
        )

        name = f"soilmoisture MedAbsErr {season}"
        metrics[name] = float(np.ma.median(np.ma.abs(dff.data)))

    return metrics


def record_provenance(diagnostic_file, config, ancestors):
    """Record provenance."""
    caption = f"Autoassess soilmoisture MedAbsErr for {SEASONS}"
    provenance_record = get_provenance_record(caption, ancestors)
    with ProvenanceLogger(config) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)


def main(config):
    """
    Top-level function for soil moisture metrics.

    Parameters
    ----------
    config : dict
        The ESMValTool configuration.
    """
    input_data = config["input_data"]

    # Separate OBS from model datasets
    # (and check there is only one obs dataset)
    obs = [v for v in input_data.values() if v["project"] == "OBS"]
    if len(obs) != 1:
        msg = f"Expected exactly 1 OBS dataset: found {len(obs)}"
        raise RuntimeError(msg)
    clim_file = obs[0]["filename"]

    models = group_metadata(
        [v for v in input_data.values() if v["project"] != "OBS"],
        "dataset",
    )

    for model_dataset, group in models.items():
        # 'model_dataset' is the name of the model dataset.
        # 'group' is a list of dictionaries containing metadata.
        logger.info("Processing data for %s", model_dataset)
        model_file = [item["filename"] for item in group]

        # Input filenames for provenance
        ancestors = flatten([model_file, clim_file])

        # Calculate metrics
        metrics = land_sm_top(
            clim_file,
            model_file,
            model_dataset,
            config,
            ancestors,
        )

        # Write metrics
        metrics_dir = os.path.join(
            config["plot_dir"],
            f"{config['exp_model']}_vs_{config['control_model']}",
            config["area"],
            model_dataset,
        )

        write_metrics(metrics_dir, metrics, config, ancestors)


if __name__ == "__main__":
    with run_diagnostic() as CONFIG:
        main(CONFIG)
