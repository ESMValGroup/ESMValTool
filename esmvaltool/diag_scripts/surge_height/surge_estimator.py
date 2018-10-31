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

config = ConfigParser.ConfigParser()
config.read('/usr/people/ridder/Documents/0_models/ESMValTool/nml/cfg_srg_estim/cfg_srg_estim.conf')


def surge_estimator_main(psl_in, ua_in, va_in):#???
	from load.load_config import load_config
	from load.load_EOFs import load_EOFs
	from load.load_betas_intercept import load_betas_intercept
	from dataprep.grad_psl import grad_psl
	from estimate.build_predictX import build_predictX 
	from estimate.estimate_srg import estimate_srg
	from output.save_netCDF import save_netCDF
	from output.plot_map import plot_map
	from output.plot_tseries import plot_tseries
	from dataprep.cut_NS_xarray import cut_NS
	from dataprep.Xtrms_xarray import Xtrms
	from dataprep.calc_monanom import calc_monanom
	import load.load_config as llc
	import load.load_EOFs   as llE
	import load.load_betas_intercept as llbi
	import dataprep.grad_psl as dpgrd
	import estimate.build_predictX as ebX
	import estimate.estimate_srg as ees
	import output.save_netCDF as osn
	import output.plot_map as opm
	import output.plot_tseries as opt
	# ---------------
	# I. Definitions
	# ---------------
	# I.1 Read config-file
	#ld.load_config.load_config()
	load_config()
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
	# WAQUA lon/lat points
	lons = np.arange(-12.00,-12.00+(0.125*201),0.125).tolist()
	lats = np.arange(48.00,48.00+(0.08333*173),0.08333).tolist()
	#
	if llc.coastal_map:	
		stat       = allstats
		dates_map  = [llc.t0]
	#
	if llc.plt_tseries:
		if llc.SOIname in allstats:
			stat  = [llc.SOIname]
		else:
			print('Station not available -> timeseries plot cannot be generated.')
			llc.plt_tseries = False
	#
	tlen  = (llc.tend - llc.tstart).total_seconds()/60./60./24.
	dates = pd.date_range(llc.tstart,periods = tlen+1).tolist()
	dates = list(map(pd.Timestamp.to_pydatetime,dates))

	# ----------------------------------------------------------------
	# II. Test if requested dates for plotting are provided dataset?
	# ----------------------------------------------------------------
	# probably unnessecary in ESMValTool context as time range is requested using definition in nml/config file
	dates_in = []
	ns = 1e-9 	# number of seconds in a nanosecond
	for t in psl_in.time.values:
		tmp = datetime.utcfromtimestamp(t.astype(int) * ns)
		dates_in.append(datetime(tmp.year,tmp.month,tmp.day,0,0))

	for x in dates: 
		if not x in dates_in:
			print ('WARNING: Selected time period not in provided data! ' +
				'Using full time range provided in dataset instead.')
			dates = dates_in
	
	# -------------------------------------------
	# III. Load solver & regression coefficients
	# -------------------------------------------
	load_EOFs()
	load_betas_intercept(stat)
	
	# ----------------------------------------------------
	# IV. Cut NS box & calculate daily maxes & anomalies
	# ----------------------------------------------------
	#ymin, ymax = [47.5, 65.6]
	#xmin, xmax = [-13.5,13.5]
	#pslNS = esmvaltool.preprocessor.extract_region(psl_in, xmin, xmax, ymin, ymax)
	#...
	pslNS, uaNS, vaNS    = cut_NS(psl_in, ua_in, va_in)

	xpslNS, xuaNS, xvaNS = Xtrms(pslNS, uaNS, vaNS)

	psl, ua, va = calc_monanom(xpslNS, xuaNS, xvaNS)
	
	# -----------------------------
	# V. Calculate SLP gradients 
	# -----------------------------
	grad_psl(psl)

	# -----------------------------------------
	# VI. Project fields onto ERA-Interim EOFs 
	# -----------------------------------------
	pseudo_pcs_SLP        = llE.SLPsolver.projectField(psl.values)
	pseudo_pcs_gradlatSLP = llE.gradlatsolver.projectField(dpgrd.gradlatSLP)
	pseudo_pcs_gradlonSLP = llE.gradlonsolver.projectField(dpgrd.gradlonSLP)
	pseudo_pcs_u          = llE.usolver.projectField(ua.values)
	pseudo_pcs_v          = llE.vsolver.projectField(va.values)

	# -------------------------------
	# VII. Generate predictor array
	# -------------------------------
	X = {}
	for s in stat:
		X[s] = build_predictX(dates,pseudo_pcs_SLP, pseudo_pcs_gradlatSLP, #...
					pseudo_pcs_gradlonSLP, pseudo_pcs_u, pseudo_pcs_v)

	# ----------------------------
	# VIII. Apply regression model
	# ----------------------------
	estimate_srg(X,dates,stat)

	# ------------------------
	# IX. Save surge to file
	# ------------------------
	#if  write_netcdf:
	save_netCDF(dates,stat,ees.srg_est_full)

	# -----------
	# X. Plot
	# -----------
	# if write_plots:
	if llc.coastal_map: # generate geographical map with surge levels on day specified in config file	
		plot_map(dates_map, ees.srg_est_full,dates.index(dates_map[0]))
	#
	if llc.plt_tseries: # generate timeseries plot
		for s in ees.srg_est_full.keys():
			plot_tseries(dates, ees.srg_est_full[s], stat)
	#
	#plt.show()


if __name__ == '__main__':
	print 'Is main'
	import sys
	surge_estimator_main(sys.argv[1],sys.argv[2],sys.argv[3]) 

