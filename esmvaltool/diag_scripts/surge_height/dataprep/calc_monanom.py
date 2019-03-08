import numpy as np
import xarray as xr


def calc_monanom(psl, uas, vas):
    # ------------------------------
    # Calculate monthly climatology
    # ------------------------------
    climatology_vas = vas.groupby('time.month').mean('time')
    climatology_uas = uas.groupby('time.month').mean('time')
    climatology_psl = psl.groupby('time.month').mean('time')

    # ----------------------------
    # Determine monthly anomalies
    # ----------------------------
    apsl_array = xr.apply_ufunc(lambda x, m: x - m, psl.groupby('time.month'),
                                climatology_psl)
    auas_array = xr.apply_ufunc(lambda x, m: x - m, uas.groupby('time.month'),
                                climatology_uas)
    avas_array = xr.apply_ufunc(lambda x, m: x - m, vas.groupby('time.month'),
                                climatology_vas)

    return apsl_array, auas_array, avas_array
