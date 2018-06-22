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
##;;    YYYYMMDD-A_X3Y3: bug-fixed... 
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

config = ConfigParser.ConfigParser()
config.read('/usr/people/ridder/Documents/0_models/ESMValTool/nml/cfg_srg_estim/cfg_srg_estim.conf')


def surge_estimator_main(psl, ua, va):#???
	from load.load_config import load_config
	from load.load_EOFs import load_EOFs
	from load.load_betas_intercept import load_betas_intercept
	from dataprep.grad_psl import grad_psl
	from estimate.build_predictX import build_predictX 
	from estimate.estimate_srg import estimate_srg
	from output.save_netCDF import save_netCDF
	from output.plot_map import plot_map
	from output.plot_tseries import plot_tseries
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
	# I.2b coordinates of stations
	coords = {'aberdeen': [81, 111], 'aukalpha': [114, 102], 'bg2': [126, 46], 'borkums': [150, 68],  
		   'bremerha': [165, 68],'cadzand': [123, 42], 'cromer': [108, 60], 'cuxhaven': [167, 72],  
		   'delfzijl': [153, 65], 'denhelde': [138, 65],'denoever': [137, 61], 'devonpor': [64, 29],  
		   'dover': [108, 39], 'duinkerk': [116, 38], 'ekofisk': [123, 104], 'esbjerg': [163, 91],  
		   'europlat': [123, 49], 'f3': [134, 83], 'felixsto': [108, 48], 'goeree': [126, 48], 
		   'harlinge': [140, 64], 'helgeroa': [175, 133], 'helgolan': [160, 75], 'hoekvanh': [130, 49], 
		   'holyhead': [60, 66], 'huibertg': [148, 68], 'husum': [168, 78], 'ijmuiden': [133, 54], 
		   'ilfracom': [64, 40], 'immingha': [98, 69], 'innerdow': [101, 65], 'k13a': [123, 64],  
		   'kornwerd': [138, 61], 'lauwerso': [147, 66], 'leith': [72, 97], 'lerwick': [89, 151],  
		   'lowestof': [111, 55], 'meetpost': [131, 53], 'newhaven': [96, 34], 'newlyn': [53, 26],  
		   'northcor': [106, 160], 'northshi': [86, 85], 'offharwi': [109, 47], 'oostende': [120, 40], 
		   'os11': [125, 45], 'os15': [125, 44], 'oscarsbo': [182, 139], 'portsmou': [88, 35],  
		   'roompotb': [126, 45], 'scarboro': [94, 77], 'scheveni': [131, 51], 'scillyis': [45, 23],  
		   'sheernes': [104, 42], 'southend': [103, 43],'stavange': [141, 133], 'stmarys': [47, 24],  
		   'stornowa': [51, 122], 'terschel': [139, 66], 'texelnoo': [135, 63], 'torsmind': [161, 101], 
		   'tregde': [158, 121], 'vidaa': [166, 85], 'vlaktevd': [122, 43], 'vlissing': [125, 42], 
		   'westkape': [124, 43],'westters': [138, 65], 'weymouth': [77, 32], 'wick': [73, 126], 
		   'zeebrugg': [122, 41]}
	# WAQUA lon/lat points
	lons = np.arange(-12.00,-12.00+(0.125*201),0.125).tolist()
	lats = np.arange(48.00,48.00+(0.08333*173),0.08333).tolist()

	#
	if llc.coastal_map:			# defined in config file
		stat   = allstats
		dates  = [llc.t0]
	else:
		if llc.SOIname in allstats:
			stat  = [llc.SOIname] # defined in config file
			tlen  = (llc.tend - llc.tstart).total_seconds()/60./60./24.
			dates = pd.date_range(llc.tstart,periods = tlen+1).tolist()
			dates = list(map(pd.Timestamp.to_pydatetime,dates))
		else:
			exit('Station not available. For a list of available stations refer to config-file.')
	
	# ------------------------------------------
	# II. Load solver & regression coefficients
	# ------------------------------------------
	load_EOFs()
	load_betas_intercept(stat)

	# -----------------------------
	# III. Calculate SLP gradients 
	# -----------------------------
	grad_psl(psl)

	# -----------------------------------------
	# IV. Project fields onto ERA-Interim EOFs 
	# -----------------------------------------
	pseudo_pcs_SLP        = llE.SLPsolver.projectField(psl)
	pseudo_pcs_gradlatSLP = llE.gradlatsolver.projectField(dpgrd.gradlatSLP)
	pseudo_pcs_gradlonSLP = llE.gradlonsolver.projectField(dpgrd.gradlonSLP)
	pseudo_pcs_u          = llE.usolver.projectField(ua)
	pseudo_pcs_v          = llE.vsolver.projectField(va)

	# -----------------------------
	# V. Generate predictor array
	# -----------------------------
	X = {}
	for s in stat:
		X[s] = build_predictX(dates,pseudo_pcs_SLP, pseudo_pcs_gradlatSLP, #...
					pseudo_pcs_gradlonSLP, pseudo_pcs_u, pseudo_pcs_v)

	# ---------------------------
	# VI. Apply regression model
	# ---------------------------
	estimate_srg(X,dates,stat)

	# ------------------------
	# VII. Save surge to file
	# ------------------------
	#save_netCDF(dates,stat,coords,lons,lats)

	# -----------
	# VIII. Plot
	# -----------
	if llc.coastal_map:	
		plot_map(dates, ees.srg_est_full)
	else:
		for s in ees.srg_est_full.keys():
			plot_tseries(dates, ees.srg_est_full[s], s)

	plt.show()


if __name__ == '__main__':
	print 'Is main'
	import sys
	surge_estimator_main(sys.argv[1],sys.argv[2],sys.argv[3]) 

