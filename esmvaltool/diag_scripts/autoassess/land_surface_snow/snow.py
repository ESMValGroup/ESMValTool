"""Module for snow metrics."""

import os

import iris
import numpy as np
from esmvalcore.preprocessor import regrid

from esmvaltool.diag_scripts.shared._supermeans import get_supermean


def land_swe_top(run):
    """
    Compute median-absolute difference of SWE against GlobSnow.

    Arguments:
        run - dictionary containing model run metadata
              (see auto_assess/model_run.py for description)

    Returns
    -------
        metrics - dictionary of metrics names and values

    """
    supermean_data_dir = os.path.join(
        run["data_root"],
        run["runid"],
        run["_area"] + "_supermeans",
    )

    snow_seasons = ["son", "djf", "mam"]

    # Calculate rms errors for seasons with snow.
    metrics = dict()
    for season in snow_seasons:
        clim_file = os.path.join(run["climfiles_root"], f"SWE_clm_{season}.pp")
        swe_clim = iris.load_cube(clim_file)
        swe_clim.data = np.ma.masked_array(
            swe_clim.data,
            mask=(swe_clim.data == -1e20),
        )

        # snowfall
        swe_run = get_supermean(
            "surface_snow_amount",
            season,
            supermean_data_dir,
        )

        # Force same coord_system
        swe_run.coord("longitude").coord_system = swe_clim.coord(
            "longitude",
        ).coord_system
        swe_run.coord("latitude").coord_system = swe_clim.coord(
            "latitude",
        ).coord_system

        # Force the units for SWE to match the model
        swe_clim.units = swe_run.units

        # form the difference
        # active regridding here
        swe_run = regrid(swe_run, swe_clim, "linear")
        dff = swe_run - swe_clim
        iris.save(
            dff,
            os.path.join(run["dump_output"], f"snow_diff_{season}.nc"),
        )

        #  Calculate median absolute error of the difference
        name = f"snow MedAbsErr {season}"
        metrics[name] = float(np.ma.median(np.ma.abs(dff.data)))

    return metrics
