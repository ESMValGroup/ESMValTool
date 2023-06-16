"""Run module for soil moisture metrics."""

import os
import logging
import csv
import numpy as np
import iris
from esmvalcore.preprocessor import regrid
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic

# Order of seasons must agree with preprocessor definition in recipe
SEASONS = ("djf", "mam", "jja", "son")

logger = logging.getLogger(__name__)


def get_provenance_record(caption, filenames):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_type': 'metrics',
        'authors': [
            'rumbold_heather',
            'sellar_alistair',
        ],
        'references': [
            'esacci-soilmoisture',
            'dorigo17rse',
            'gruber19essd',
        ],
        "ancestors": filenames,
    }

    return record


def write_metrics(output_dir, metrics, config):
    """Write metrics to CSV file.

    The CSV file will have the name ``metrics.csv`` and can be
    used for the normalised metric assessment plot.

    Parameters
    ----------
    output_dir : string
        The full path to the directory in which the CSV file will be written.
    metrics : dictionary of metric,value pairs
        The seasonal data to write.
    """
    os.makedirs(output_dir, exist_ok=True)

    #TODO use esmvaltool.diag_scripts.shared.get_diagnostic_filename("metrics", config, extension=csv)
    file_name = "metrics.csv"
    file_path = os.path.join(output_dir, file_name)

    with open(file_path, "w", newline="", encoding="utf-8") as csvfile:
        csv_writer = csv.writer(csvfile)
        for line in metrics.items():
            csv_writer.writerow(line)

    record_provenance(file_path, config)


def land_sm_top(clim_file, model_file, work_dir):
    """
    Calculate median absolute errors for soil mosture against CCI data.

    Parameters
    ----------
    clim_file : string
        Path to observation climatology file
    model_file : string
        Path to model file
    work_dir : string
        Directory for intermediate working files

    Returns
    -------
    metrics: dict
        a dictionary of metrics names and values
    """
    # Constant: density of water
    rhow = 1000.

    # Work through each season
    metrics = {}
    for index, season in enumerate(SEASONS):

        constr_season = iris.Constraint(season_number=index)
        ecv_clim = iris.load_cube(clim_file, constr_season)

        # m01s08i223
        # CMOR name: mrsos (soil moisture in top model layer kg/m2)
        mrsos = iris.load_cube(
            model_file,
            "mass_content_of_water_in_soil_layer" & constr_season
        )

        # Set soil moisture to missing data on ice points (i.e. no soil)
        np.ma.masked_where(mrsos.data == 0, mrsos.data, copy=False)

        # first soil layer depth
        dz1 = mrsos.coord('depth').bounds[0, 1] - \
            mrsos.coord('depth').bounds[0, 0]

        # Calculate the volumetric soil moisture in m3/m3
        # volumetric soil moisture = volume of water / volume of soil layer
        # = depth equivalent of water / thickness of soil layer
        # = (soil moisture content (kg m-2) / water density (kg m-3) )  /
        #      soil layer thickness (m)
        # = mosrs / (rhow * dz1)
        vol_sm1_run = mrsos / (rhow * dz1)
        vol_sm1_run.units = "m3 m-3"
        vol_sm1_run.long_name = "Top layer Soil Moisture"

        # update the coordinate system ECV data with a WGS84 coord system
        # unify coord systems for regridder
        vol_sm1_run.coord('longitude').coord_system = \
            iris.coord_systems.GeogCS(semi_major_axis=6378137.0,
                                      inverse_flattening=298.257223563)
        vol_sm1_run.coord('latitude').coord_system = \
            iris.coord_systems.GeogCS(semi_major_axis=6378137.0,
                                      inverse_flattening=298.257223563)
        ecv_clim.coord('longitude').coord_system = \
            iris.coord_systems.GeogCS(semi_major_axis=6378137.0,
                                      inverse_flattening=298.257223563)
        ecv_clim.coord('latitude').coord_system = \
            iris.coord_systems.GeogCS(semi_major_axis=6378137.0,
                                      inverse_flattening=298.257223563)

        # Interpolate to the grid of the climatology and form the difference
        vol_sm1_run = regrid(vol_sm1_run, ecv_clim, 'linear')

        # mask invalids
        vol_sm1_run.data = np.ma.masked_invalid(vol_sm1_run.data)
        ecv_clim.data = np.ma.masked_invalid(ecv_clim.data)

        # diff the cubes
        dff = vol_sm1_run - ecv_clim

        # save output and populate metric
        iris.save(dff, os.path.join(work_dir,
                                    f"soilmoist_diff_{season}.nc"))
        name = f"soilmoisture MedAbsErr {season}"
        dffs = dff.data
        dffs = np.ma.abs(dffs)
        metrics[name] = float(np.ma.median(dffs))

    return metrics


def record_provenance(diagnostic_file, config):
    """Record provenance."""
    caption = f"Autoassess soilmoisture MedAbsErr for {SEASONS}"
    filenames = [item["filename"] for item in config["input_data"].values()]
    provenance_record = get_provenance_record(caption, filenames)
    cfg = {}
    cfg["run_dir"] = config["run_dir"]  #TODO needed?
    with ProvenanceLogger(cfg) as provenance_logger:
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
        "dataset")

    for model_dataset, group in models.items():
        # 'model_dataset' is the name of the model dataset.
        # 'group' is a list of dictionaries containing metadata.
        logger.info("Processing data for %s", model_dataset)
        model_file = [item["filename"] for item in group]
        metrics = land_sm_top(clim_file, model_file, config["work_dir"])

        # Write metrics
        metrics_dir = os.path.join(
            config["plot_dir"],
            f"{config['exp_model']}_vs_{config['control_model']}",
            config["area"],
            model_dataset,
        )

        write_metrics(metrics_dir, metrics, config)


if __name__ == "__main__":
    with run_diagnostic() as CONFIG:
        main(CONFIG)
