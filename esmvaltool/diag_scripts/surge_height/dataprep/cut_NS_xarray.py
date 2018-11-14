import numpy as np
import xarray as xr


def cut_NS(psl, ua, va):
    ymin, ymax = [47.5, 65.6]
    xmin, xmax = [166.5, 193.5] #[-13.5 + 180, 13.5 + 180]

    ## ------------------------------------
    ## I. Interpolate to ERA-Interim grid
    ## ------------------------------------
    # remap manually
    lons = np.arange(0,360,0.75)
    lats = np.arange(-89.75,90,0.75)
    psl = psl.interp(lon=lons,lat=lats,method='linear', kwargs={'fill_value': None})
    ua = ua.interp(lon=lons,lat=lats,method='linear', kwargs={'fill_value': None})
    va = va.interp(lon=lons,lat=lats,method='linear', kwargs={'fill_value': None})

    # -------------------
    # III. Cut NS region
    # -------------------
    pslNS = psl.sel(lon=slice(xmin, xmax), lat=slice(ymin, ymax))
    uaNS = ua.sel(lon=slice(xmin, xmax), lat=slice(ymin, ymax))
    vaNS = va.sel(lon=slice(xmin, xmax), lat=slice(ymin, ymax))

    return pslNS, uaNS, vaNS
