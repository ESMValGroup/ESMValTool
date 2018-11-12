import numpy as np
from eofs.standard import Eof 
import sys
import os
import pickle
import ConfigParser
from datetime import datetime
from netCDF4 import Dataset
import xarray as xr

import surge_estimator
from surge_estimator import surge_estimator_main

# ==================
# I. Load test data
# ==================
print 'Loading test data...'
ECPATH = '/usr/people/ridder/Documents/0_models/ESMValTool_backup/diag_scripts/lib/python/surge_height/test_data/'

SLPf  = 'daymin_monthly_anom_psl_6h_ECEarth_PD_s01r15_2035'
uf    = 'daymax_monthly_anom_uas_6h_ECEarth_PD_s01r15_2035'
vf    = 'daymax_monthly_anom_vas_6h_ECEarth_PD_s01r15_2035'

nc_mslp = xr.open_dataset(ECPATH + SLPf + '.nc')
#if not 'psl' in va.data_vars.keys():
nc_mslp = nc_mslp.rename({'var151':'psl'})
psl     = nc_mslp.psl

#
nc_uas  = xr.open_dataset(ECPATH + uf + '.nc')
nc_uas  = nc_uas.rename({'var151':'ua'})
uas     = nc_uas.ua

#
nc_vas  = xr.open_dataset(ECPATH + vf + '.nc')
nc_vas  =nc_vas.rename({'var151':'va'})
vas     = nc_vas.va
#exit()

# ======================================
# II. Call surge estimator main script
# ======================================
print 'Calling surge estimator...'
surge_estimator_main(psl, uas, vas) #nc_psl,nc_uas,nc_vas)

