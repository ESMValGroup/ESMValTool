import numpy as np
import xarray as xr

def calc_monanom(psl, ua, va):
    # ------------------------------
    # Calculate monthly climatology
    # ------------------------------
    #climatology = ds.groupby('time.month').mean('time')
    # PROBLEM: This only caluclates the climatology for the input file not the full model run!!!
    climatology_va = va.groupby('time.month').mean('time')
    climatology_ua = ua.groupby('time.month').mean('time')
    climatology_psl = psl.groupby('time.month').mean('time')

    # ----------------------------
    # Determine monthly anomalies
    # ----------------------------
    #anomalies = ds.groupby('time.month') - climatology
    apsl_array = xr.apply_ufunc(lambda x, m: x - m, psl.groupby('time.month'),
                                climatology_psl)
    aua_array = xr.apply_ufunc(lambda x, m: x - m, ua.groupby('time.month'),
                               climatology_ua)
    ava_array = xr.apply_ufunc(lambda x, m: x - m, va.groupby('time.month'),
                               climatology_va)

    # ----------------------------
    # Convert into list object
    # ----------------------------
    ##apsl = apsl_array.values
    ##aua  = aua_array.values
    ##ava  = ava_array.values
    #apsl = np.stack([data_array for data_array in apsl_array.values])
    #aua = np.stack([data_array for data_array in aua_array.values])
    #ava = np.stack([data_array for data_array in ava_array.values])

    return apsl_array, aua_array, ava_array
