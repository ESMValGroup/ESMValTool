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

import ConfigParser
import os
import pickle
import sys
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from eofs.standard import Eof

#import .dataprep.grad_psl as dpgrd
#import .estimate.build_predictX as ebX
#import .estimate.estimate_srg as ees
#import .load.load_betas_intercept as llbi
#import .load.load_EOFs
#import .output.plot_map as opm
#import .output.plot_tseries as opt
#import .output.save_netCDF as osn
from dataprep.calc_monanom import calc_monanom
from dataprep.cut_NS_xarray import cut_NS
from dataprep.grad_psl import grad_psl
from dataprep.Xtrms_xarray import Xtrms
from estimate.build_predictX import build_predictX
from estimate.estimate_srg import estimate_srg
from load.check_solver import check_solver
from load.load_betas_intercept import load_betas_intercept
from load.load_config import load_config
from load.load_EOFs import load_eofs
from output.plot_map import plot_map
from output.plot_map_cartopy import plot_map_cartopy
from output.plot_tseries import plot_tseries
from output.save_netCDF import save_netCDF


def surge_estimator_main(psl_in, ua_in, va_in, cfg, dataset):
    # ---------------
    # I. Definitions
    # ---------------
    # I.1 Read config-file
    #ld.load_config.load_config()
    #load_config()
    # I.2 Stations
    # I.2a Names
    allstats = [
        "aberdeen", "delfzijl", "europlat", "holyhead", "kornwerd", "northcor",
        "roompotb", "stornowa", "westkape", "aukalpha", "denhelde", "f3",
        "huibertg", "lauwerso", "northshi", "scarboro", "terschel", "westters",
        "bg2", "denoever", "felixsto", "husum", "leith", "offharwi",
        "scheveni", "texelnoo", "weymouth", "borkums", "devonpor", "goeree",
        "ijmuiden", "lerwick", "oostende", "scillyis", "torsmind", "wick",
        "bremerha", "dover", "harlinge", "ilfracom", "lowestof", "os11",
        "sheernes", "tregde", "zeebrugg", "cadzand", "duinkerk", "helgeroa",
        "immingha", "meetpost", "os15", "southend", "vidaa", "cromer",
        "ekofisk", "helgolan", "innerdow", "newhaven", "oscarsbo", "stavange",
        "vlaktevd", "cuxhaven", "esbjerg", "hoekvanh", "k13a", "newlyn",
        "portsmou", "stmarys", "vlissing"
    ]
    #
    if cfg["coastal_map"]:
        stat = allstats
        dates_map = [
            datetime(cfg['t0'].year, cfg['t0'].month, cfg['t0'].day, 0, 0)
        ]
    #
    if cfg["plt_tseries"]:
        if cfg["SOIname"] in allstats:
            stat = [cfg["SOIname"]]
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
    ns = 1e-9  # number of seconds in a nanosecond
    for t in psl_in.time.values:
        tmp = datetime.utcfromtimestamp(t.astype(int) * ns)
        dates.append(datetime(tmp.year, tmp.month, tmp.day, 0, 0))

    dates = list(set(dates))

    # ----------------------------------------------------
    # III. Cut NS box & calculate daily maxes & anomalies
    # ----------------------------------------------------
    #logger.debug
    print("Preprocessing input data")

    pslNS, uaNS, vaNS = cut_NS(psl_in, ua_in, va_in)

    xpslNS, xuaNS, xvaNS = Xtrms(pslNS, uaNS, vaNS)
    print('after xtrms')
    print(np.isnan(xpslNS.values.flatten()).any())

    psl, ua, va = calc_monanom(xpslNS, xuaNS, xvaNS)
    print('after preprocessing')
    print(np.isnan(psl.values.flatten()).any())

    # -----------------------------
    # IV. Calculate SLP gradients
    # -----------------------------
    #logger.debug("Calculating gradient of psl")
    gradlatpsl, gradlonpsl = grad_psl(psl)

    # -----------------------------------------------------------
    # V. Load or determine solver & load regression coefficients
    # -----------------------------------------------------------
    #logger.debug("Loading EOFs and regression coefficients")
    dirname = os.path.dirname(__file__)
    data_dir = os.path.join(dirname, 'data')

    solvers_exist = check_solver(data_dir)

    if cfg['perf_regres'] or not solvers_exist:
        [psl_solver, gradlon_solver, gradlat_solver, 
          ua_solver, va_solver] = calc_eofs(cfg['path2eraint'], data_dir)
    else:
        psl_solver, gradlon_solver, gradlat_solver, ua_solver, va_solver = load_eofs(data_dir)

    betas, intercept = load_betas_intercept(stat, data_dir)

    # -----------------------------------------
    # VI. Project fields onto ERA-Interim EOFs
    # -----------------------------------------
    #logger.debug("Generating PCs")
    pseudo_pcs_slp = psl_solver.projectField(psl.values)
    pseudo_pcs_gradlatpsl = gradlat_solver.projectField(gradlatpsl)
    pseudo_pcs_gradlonpsl = gradlon_solver.projectField(gradlonpsl)
    pseudo_pcs_u = ua_solver.projectField(ua.values)
    pseudo_pcs_v = va_solver.projectField(va.values)

    # -------------------------------
    # VII. Generate predictor array
    # -------------------------------
    #logger.debug("Generating predictor array")
    X = {}
    for s in stat:
        X[s] = build_predictX(
            dates,
            pseudo_pcs_slp,
            pseudo_pcs_gradlatpsl,  #...
            pseudo_pcs_gradlonpsl,
            pseudo_pcs_u,
            pseudo_pcs_v)

    # ----------------------------
    # VIII. Apply regression model
    # ----------------------------
    #logger.debug("Estimating surge heights")
    srg_estim = estimate_srg(X, dates, stat, betas, intercept)

    # ------------------------
    # IX. Save surge to file
    # ------------------------
    if not "fout_name" in cfg.keys():
        cfg["fout_name"] = 'surge'
    elif cfg["fout_name"] == "":
        cfg["fout_name"] = 'surge'
    if cfg["write_netcdf"]:
        #logger.debug('saving data to netCDF file ')
        save_netCDF(dates, stat, srg_estim, cfg, dataset)

    # -----------
    # X. Plot
    # -----------
    # if write_plots:
    if cfg["coastal_map"]:  # generate geographical map with surge levels on day specified in config file
        #logger.info("Plotting and saving geographical map with surge heights")
        plot_map_cartopy(dates_map[0], srg_estim,
                         dates.index(dates_map[0]), cfg, dataset)
    #
    if cfg["plt_tseries"]:  # generate timeseries plot
        #logger.info("Plotting and saving surge timeseries")
        for s in srg_estim.keys():
            plot_tseries(dates, srg_estim[s], stat, cfg, dataset)
    #
    #plt.show()


if __name__ == '__main__':
    print 'Is main'
    import sys
    surge_estimator_main(sys.argv[1], sys.argv[2], sys.argv[3])
