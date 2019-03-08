import os

import numpy as np
import xarray as xr


def load_monmean_srgclim(stat, data_dir):
    filename = 'monmean_ERAintWAQUA_surge_1979-2016_speed.nc'
    ncpath = os.path.join(data_dir, filename)
    nc_srg_anom = xr.open_dataset(ncpath)
    srg_anom = nc_srg_anom.WAQUA_surge

    monmean_srgclim = srg_anom.sel(stations=np.bytes_(stat.ljust(8)))

    #ncpath = os.path.join(data_dir, filename + str(ifname) + '.nc')
    #nc = Dataset(ncpath, 'r')
    #nc.variables['WAQUA_surge'][:]
    #nc.close()
    #
    return monmean_srgclim
