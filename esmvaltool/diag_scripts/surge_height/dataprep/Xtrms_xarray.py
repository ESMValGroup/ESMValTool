import numpy as np
import xarray as xr

def Xtrms(psl, ua, va):

    dmin_psl = psl.resample(time='D').min('time')
    dmax_ua = ua.resample(time='D').max('time')
    dmax_va = va.resample(time='D').max('time')

    for t in range(dmin_psl.shape[0]):
        if (np.isnan(dmin_psl[t])).any(): #while?
            rem_idx = t
            dmin_psl = dmin_psl.drop(dmin_psl.time.values[t], dim='time')
            dmax_ua = dmin_ua.drop(dmin_ua.time.values[t], dim='time')
            dmax_va = dmin_va.drop(dmin_va.time.values[t], dim='time')

    return dmin_psl, dmax_ua, dmax_va
