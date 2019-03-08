import os

import numpy as np
import xarray as xr


def load_betas_intercept(stat, inpath):
    if not type(stat) == list:
        stat = [stat]

    betas = {}
    intercept = {}

    filename_int = 'intercepts_all_stations.nc'
    filename_bs = 'betas_all_stations.nc'
    fin = os.path.join(inpath, filename_bs)
    betas_in = xr.open_dataarray(fin)

    fin = os.path.join(inpath, filename_int)
    interc = xr.open_dataarray(fin)

    for s in stat:
        betas[s] = betas_in.sel(station=s).values
        intercept[s] = interc.sel(station=s).values

    return betas, intercept
