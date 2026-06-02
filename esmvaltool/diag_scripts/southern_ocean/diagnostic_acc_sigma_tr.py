"""
Look at this module for guidance how to write your own.

Read the README_PERSONAL_DIAGNOSTIC file associated with this example;

Module for personal diagnostics (example).
Internal imports from exmvaltool work e.g.:

from esmvalcore.preprocessor import regrid
from esmvaltool.diag_scripts.shared.supermeans import get_supermean

Pipe output through logger;

Please consult the documentation for help with esmvaltool's functionalities
and best coding practices.
"""
# place your module imports here:
import gsw
import xarray as xr
import numpy as np

# import cmocean

import logging

# operating system manipulations (e.g. path constructions)
import os
import sys
from pathlib import Path

# to manipulate iris cubes
# import iris
import matplotlib.pyplot as plt
from esmvalcore.preprocessor import area_statistics

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic, ProvenanceLogger, save_figure

# reuse tools already developed for ocean diagnostics
from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools

# This part sends debug statements to stdout
# logger = logging.getLogger(os.path.basename(__file__))
# logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

# This outputs to the run/transect/script1/log.txt file
logger = logging.getLogger(Path(__file__).stem)

# module checks
logger.info("GSW version: %s", gsw.__version__)


def _get_data(cfg):

    '''
        Get the data from the input files specified in the recipe.
        Open the files and return a dictionary of xarray datasets.
    '''

    #TODO needs adjusting for multi dataset recipes

    input_data = cfg["input_data"]

    ds = {}
    for filename, attributes in input_data.items():
        logger.info("Loading %s", filename)
        variable_name = attributes["short_name"]
        ds[variable_name] = xr.open_dataset(filename)

    return ds

def _compute_sigma(ds):
    ''' 
        Compute the sigma2 variable from the input dataset.
        Uses gsw.density.sigma2 from the gsw toolbox 
            https://teos-10.github.io/GSW-Python/_modules/gsw/_wrapped_ufuncs.html#sigma2
    '''

    # logger.info("Numpy so is: %s", ds["so"]["so"].to_numpy())
    # logger.info("Numpy thetao is: %s", ds["thetao"]["thetao"].to_numpy())

    # pull the data and coords as numpy arrays
    so = ds["so"]["so"].to_numpy()              # (lev, lat)
    thetao = ds["thetao"]["thetao"].to_numpy()  # (lev, lat)
    depth = ds["so"]["lev"].to_numpy()          # (lev,)
    lat = ds["so"]["lat"].to_numpy()            # (lat,)
    lon = float(ds["so"]["lon"])                # scalar transect longitude

    # broadcast depth and lat onto the full (lev, lat) transect grid so all
    # operands match the 2D salinity/temperature fields (gsw broadcasts
    # every argument together element-wise).
    depth2d, lat2d = np.meshgrid(depth, lat, indexing="ij")  # both (lev, lat)

    # calculate pressure from depth (latitude-dependent)
    p = gsw.p_from_z(-depth2d, lat2d)
    logger.info("Pressure shape: %s", p.shape)

    # TEOS-10 chain: SA (absolute salinity) -> CT (conservative temp) -> sigma2.
    # calculate absolute salinity from practical salinity
    SA = gsw.SA_from_SP(so, p, lon, lat2d)
    logger.info("Absolute salinity shape: %s", SA.shape)

    # calculate conservative temperature from potential temperature
    CT = gsw.CT_from_pt(SA, thetao)
    logger.info("Conservative temperature shape: %s", CT.shape)

    # compute sigma2
    sigma2 = gsw.sigma2(SA, CT)

    # create a new xarray dataset with the same coordinates as the input datasets
    sigma2_ds = xr.Dataset(
        {
            "sigma2": (["lev", "lat"], sigma2)
        },
        coords={
            "time": ds["so"]["time"],
            "lev": ds["so"]["lev"],
            "lat": ds["so"]["lat"],
            "lon": ds["so"]["lon"]
        },
        attrs={
            "short_name": "sigma2",
            "long_name": "Potential density anomaly referenced to 2000 dbar",
            "units": "kg/m^3",
        }
    )

    return sigma2_ds

def _plot_transect(ds, sigma, cfg):
    '''
        Plot the transect of zonal velocity with sigma2 contours overlaid.

        Pass in the ds input like ds["uo"]
    '''

    # extract data arrays
    sigma2 = sigma["sigma2"]
    uo = ds["uo"]

    # create the figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))

    uo.plot.contourf(ax=ax, x="lat", y="lev", 
                     yincrease=False, 
                     cmap="RdBu_r")
    
    sigma2.plot.contour(ax=ax, x="lat", y="lev", 
                        yincrease=False,
                        colors="k")

    ax.set_title("Drake Passage Zonal Velocity Transect with Sigma2 Contours")

    provenance = {
        "caption": "Zonal velocity transect across Drake Passage with sigma2 contours.",
        "statistics": ["mean"],
        "domains": ["drake_passage_transect"],
        "plot_types": ["contourf", "contour"],
        "authors": ["Thomas Wilder"],
        "ancestors": list(cfg["input_data"].keys()),
    }

    # save the figure
    save_figure("drakepassage_sigma2_transect", provenance, cfg)


def main(cfg):
    # get the data
    datasets = _get_data(cfg)

    logger.info("Datasets: %s", datasets)

    sigma2_ds = _compute_sigma(datasets)

    logger.info("Sigma2 dataset: %s", sigma2_ds)
    

    _plot_transect(datasets["uo"], sigma2_ds, cfg)



if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
