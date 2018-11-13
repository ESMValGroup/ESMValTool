#;;############################################################################# 
#;; Surge estimator 
#;; Author: Nina Nadine Ridder (KNMI, The Netherlands) 
#;; C3 MAGIC project 
#;;############################################################################# 
#;; Description 
#;;    Diagnostic to estimate surge heights from sea level pressure (msp) 
#;;    and u- and v- component of 10m wind
#;;    Additional description of the diagnostic 
#;;    Method to estimate surge described in  et al. () 
#;; 
#;; Required diag_script_info attributes (diagnostics specific) 
#;;    att1: short description 
#;;          keep the indentation if more lines are needed 
#;;    att2: short description 
#;; 
#;; Optional diag_script_info attributes (diagnostic specific) 
#;;    att1: short description 
#;;    att2: short description 
#;; 
#;; Required variable_info attributes (variable specific) 
#;;    att1: short description 
#;;    att2: short description 
#;; 
#;; Optional variable_info attributes (variable specific) 
#;;    att1: short description 
#;;    att2: short description 
#;; 
#;; Caveats 
##;;    assumes psl, ua & va are provided on ERA-Interim grid 
##;;    surge data generated using WAQUA/DCSMv5, which underestimates extreme
##;;    Features to-be-implemented shall also be mentioned here 
##;; 
##;; Modification history 
##;;    YYYYMMDD-A_X4Y4: extended... 
##;;    20180927-A_RN: bug-fixed... 
##;;    20180522-A_RN: adapted to include module structure 
##;;    20180104-A_RN: written. 
##;; 
##;; #############################################################################
   
#load ... 
#load ...   


#begin      

import numpy as np
from eofs.standard import Eof 
import sys
import os
import pickle
import ConfigParser
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
from load.load_config import load_config
from load.load_EOFs import load_EOFs
from load.load_betas_intercept import load_betas_intercept
from dataprep.grad_psl import grad_psl
from estimate.build_predictX import build_predictX 
from estimate.estimate_srg import estimate_srg
from output.save_netCDF import save_netCDF
from output.plot_map import plot_map
from output.plot_map_cartopy import plot_map_cartopy
from output.plot_tseries import plot_tseries
from dataprep.cut_NS_xarray import cut_NS
from dataprep.Xtrms_xarray import Xtrms
from dataprep.calc_monanom import calc_monanom
import load.load_EOFs   as llE
import load.load_betas_intercept as llbi
import dataprep.grad_psl as dpgrd
import estimate.build_predictX as ebX
import estimate.estimate_srg as ees
import output.save_netCDF as osn
import output.plot_map as opm
import output.plot_tseries as opt


def surge_estimator_main(psl_in, ua_in, va_in, cfg, dataset):
	# ---------------
	# I. Definitions
	# ---------------
	# I.1 Read config-file
	#ld.load_config.load_config()
	#load_config()
	# I.2 Stations
	# I.2a Names
	allstats = ["aberdeen", "delfzijl", "europlat", "holyhead", "kornwerd", "northcor", "roompotb",
	  "stornowa", "westkape", "aukalpha", "denhelde", "f3", "huibertg", "lauwerso", "northshi", 
	  "scarboro", "terschel", "westters", "bg2", "denoever", "felixsto", "husum", "leith", "offharwi", 
	  "scheveni", "texelnoo", "weymouth", "borkums", "devonpor", "goeree", "ijmuiden", "lerwick", 
	  "oostende", "scillyis", "torsmind", "wick", "bremerha", "dover", "harlinge", "ilfracom", 
	  "lowestof", "os11", "sheernes", "tregde", "zeebrugg", "cadzand", "duinkerk", "helgeroa", 
	  "immingha", "meetpost", "os15", "southend", "vidaa", "cromer", "ekofisk", "helgolan", 
	  "innerdow", "newhaven", "oscarsbo", "stavange", "vlaktevd", "cuxhaven", "esbjerg", "hoekvanh", 
	  "k13a", "newlyn", "portsmou", "stmarys", "vlissing"]
	#
	if cfg["coastal_map"]:	
		stat       = allstats
		dates_map  = [datetime(cfg['t0'].year,cfg['t0'].month,cfg['t0'].day,0,0)]
	#
	if cfg["plt_tseries"]:
		if cfg["SOIname"] in allstats:
			stat  = [cfg["SOIname"]]
		else:
			#logger.info(
			print('Station ' + str(cfg["plt_tseries"]) + 
					' is not available -> timeseries plot cannot be generated.')
			cfg["plt_tseries"] = False
	#
	# ---------------------------------------------------------------------------
	# II. Determine time series range and dates as if input was daily timeseries
	# ---------------------------------------------------------------------------
        dates = []
	ns = 1e-9 	# number of seconds in a nanosecond
	for t in psl_in.time.values:
		tmp = datetime.utcfromtimestamp(t.astype(int) * ns)
		dates.append(datetime(tmp.year,tmp.month,tmp.day,0,0))

        dates = list(set(dates))
	
	# -------------------------------------------
	# III. Load solver & regression coefficients
	# -------------------------------------------
        #logger.debug("Loading EOFs and regression coefficients")
	load_EOFs()
	load_betas_intercept(stat)
	
	# ----------------------------------------------------
	# IV. Cut NS box & calculate daily maxes & anomalies
	# ----------------------------------------------------
        #logger.debug
	print("Preprocessing input data")

        # remap manually
	lons = np.arange(0,360,0.75) 
	lats = np.arange(-89.75,90,0.75) 
	psl_in = psl_in.interp(lon=lons,lat=lats,method='linear', kwargs={'fill_value': None})
	ua_in = ua_in.interp(lon=lons,lat=lats,method='linear', kwargs={'fill_value': None})
	va_in = va_in.interp(lon=lons,lat=lats,method='linear', kwargs={'fill_value': None})
	#print('after regrid')
	#print(np.isnan(psl_in.values.flatten()).any())
	#exit()

	pslNS, uaNS, vaNS    = cut_NS(psl_in, ua_in, va_in)
	print('after NSbox')
	print(np.isnan(pslNS.values.flatten()).any())

	xpslNS, xuaNS, xvaNS = Xtrms(pslNS, uaNS, vaNS)
	print('after xtrms')
	print(np.isnan(xpslNS.values.flatten()).any())

	psl, ua, va = calc_monanom(xpslNS, xuaNS, xvaNS)
	print('after preprocessing')
	print(np.isnan(psl.values.flatten()).any())
	
	# -----------------------------
	# V. Calculate SLP gradients 
	# -----------------------------
        #logger.debug("Calculating gradient of psl")
	grad_psl(psl)

	# -----------------------------------------
	# VI. Project fields onto ERA-Interim EOFs 
	# -----------------------------------------
        #logger.debug("Generating PCs")
	pseudo_pcs_SLP        = llE.SLPsolver.projectField(psl.values)
	pseudo_pcs_gradlatSLP = llE.gradlatsolver.projectField(dpgrd.gradlatSLP)
	pseudo_pcs_gradlonSLP = llE.gradlonsolver.projectField(dpgrd.gradlonSLP)
	pseudo_pcs_u          = llE.usolver.projectField(ua.values)
	pseudo_pcs_v          = llE.vsolver.projectField(va.values)

	# -------------------------------
	# VII. Generate predictor array
	# -------------------------------
        #logger.debug("Generating predictor array")
	X = {}
	for s in stat:
		X[s] = build_predictX(dates,pseudo_pcs_SLP, pseudo_pcs_gradlatSLP, #...
					pseudo_pcs_gradlonSLP, pseudo_pcs_u, pseudo_pcs_v)

	# ----------------------------
	# VIII. Apply regression model
	# ----------------------------
        #logger.debug("Estimating surge heights")
	estimate_srg(X,dates,stat)

	# ------------------------
	# IX. Save surge to file
	# ------------------------
        if not "fout_name" in cfg.keys():
	     cfg["fout_name"] = 'surge'
        elif cfg["fout_name"] == "":
             cfg["fout_name"] = 'surge'
	if cfg["write_netcdf"]:
		#logger.debug('saving data to netCDF file ')
                save_netCDF(dates,stat,ees.srg_est_full, cfg, dataset)

	# -----------
	# X. Plot
	# -----------
	# if write_plots:
	if cfg["coastal_map"]: # generate geographical map with surge levels on day specified in config file
        	#logger.info("Plotting and saving geographical map with surge heights")	
		plot_map_cartopy(dates_map[0], ees.srg_est_full,dates.index(dates_map[0]),cfg, dataset)
	#
	if cfg["plt_tseries"]: # generate timeseries plot
        	#logger.info("Plotting and saving surge timeseries")
		for s in ees.srg_est_full.keys():
			plot_tseries(dates, ees.srg_est_full[s], stat, cfg, dataset)
	#
	#plt.show()


if __name__ == '__main__':
	print 'Is main'
	import sys
	surge_estimator_main(sys.argv[1],sys.argv[2],sys.argv[3]) 

