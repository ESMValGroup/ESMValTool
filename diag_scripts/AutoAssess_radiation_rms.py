"""
;;#############################################################################
;; AutoAssess_radiation_rms.py
;; Author: Yoko Tsushima (Met Office, UK)
;; CMUG project
;;#############################################################################
;; Description
;;    This script is the RMS error metric script of 
;;    AutoAssess radiation
;;
;; Required diag_script_info attributes (diagnostics specific)
;;
;; [general]
;; plot_clouds:           Switch to plot cloud diags
;; plot_fluxes:           Switch to plot fluxes diags
;; plot_radiation:        Switch to plot radiation diags
;; plot_scatter:          Switch to plot scatter diags
;; plot_background_grid:  Switch to plot background grid
;;
;; plot_total_cover:     Sub map switch for plotting
;; plot_liquid_path:     Sub map switch for plotting
;; plot_ice_path:        Sub map switch for plotting
;; plot_optical_depth:   Sub map switch for plotting
;; plot_flux_maps:       Sub map switch for plotting
;; plot_radiation_maps:  Sub map switch for plotting
;;
;; plot_lat_averages:        General switch valid for all plots
;; plot_lon_averages:        General switch valid for all plots
;; plot_monthly_averages:    General switch valid for all plots
;; plot_sub_areas:           General switch valid for all plots
;;
;; mask_unwanted_values: Switch for masking values
;; mask_limit_low:       Mask values below this value
;; mask_limit_high:      Mask values above this value
;;
;; [SouthernHemisphere]
;; areas:           Control what plots (regions/seasons) are generated
;; sub_areas:       Control what plots (regions/seasons) are generated
;; scatter_areas:   Control what plots (regions/seasons) are generated
;; seasons:         Control what plots (regions/seasons) are generated
;;
;; [SouthernHemisphere_default]
;; # Latitudes [-90, 90 degrees]: 10S = -10, 10N = 10; longtitudes [0, 360 degrees]
;; lat_min: Default area
;; lat_max: Default area
;; lon_min: Default area
;; lon_max: Default area
;;
;; stride:      Color difference interval
;; maxshades:   Color difference interval
;;
;; contour_limits_clt:   (min, max, diff, [dev_min, dev_max]
;; contour_limits_clivi: (min, max, diff, [dev_min, dev_max]
;; contour_limits_clwvi: (min, max, diff, [dev_min, dev_max]
;; contour_limits_hfls:  (min, max, diff, [dev_min, dev_max]
;; contour_limits_hfss:  (min, max, diff, [dev_min, dev_max]
;; contour_limits_rlut:  (min, max, diff, [dev_min, dev_max]
;; contour_limits_rsut:  (min, max, diff, [dev_min, dev_max]
;; contour_limits_rlds:  (min, max, diff, [dev_min, dev_max]
;; contour_limits_rsds:  (min, max, diff, [dev_min, dev_max]
;;
;; colourmap_clouds: Python matplotlib colormaps, '_r' indicates a reversed colormap
;; colourmap_model:  Python matplotlib colormaps, '_r' indicates a reversed colormap
;; colourmap_diff:   Python matplotlib colormaps, '_r' indicates a reversed colormap
;; colourmap_dev:    Python matplotlib colormaps, '_r' indicates a reversed colormap
;;
;;
;; # Define area specifications for sub_areas
;; [SouthernHemisphere_northern]
;; lat_min: Define areas for northern parts
;; lat_max: Define areas for northern parts
;; lon_min: Define areas for northern parts
;; lon_max: Define areas for northern parts
;;
;; [SouthernHemisphere_southern]
;; lat_min: Define areas for southern parts
;; lat_max: Define areas for southern parts
;; lon_min: Define areas for southern parts
;; lon_max: Define areas for southern parts
;;
;; # Months to use for each season - 1 is January and so forth.
;; [SouthernHemisphere_season_DJF]
;; season_months: Months to use, 1 i January and so forth
;;
;; [SouthernHemisphere_season_MAM]
;; season_months: Months to use, 1 i January and so forth
;;
;; [SouthernHemisphere_season_JJA]
;; season_months: Months to use, 1 i January and so forth
;;
;; [SouthernHemisphere_season_SON]
;; season_months: Months to use, 1 i January and so forth
;;
;; # Define configuration for cloud vs radiation scatter plots
;; [SouthernHemisphere_scatter_default]
;; lat_min: Configuration for cloud vs radiation scatter plots
;; lat_max: Configuration for cloud vs radiation scatter plots
;; lon_min: Configuration for cloud vs radiation scatter plots
;; lon_max: Configuration for cloud vs radiation scatter plots
;; points:  Configuration for cloud vs radiation scatter plots
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;
;; Required variable_info attributes (variable specific)
;;    long_name: Name displayed in plot
;;    units:     Units
;;
;; Optional variable_info attributes (variable specific)
;;
;; Caveats
;;   landsea.nc from NCLA for land/sea mask and resolution. This does not have coord system info.
;;
;; Modification history
;;    20170323-_AutoAssess_radiation_rms: Test finished.
;;    20160819-_test_AutoAssess_radiation_rms: written based on calc_rms code.
;;
;; #############################################################################
"""

# Basic Python packages
import ConfigParser
import os
import pdb
import sys
import numpy as np
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# ESMValTool defined Python packages
sys.path.append("./interface_scripts")
from esmval_lib import ESMValProject
from auxiliary import info

# Force matplotlib to not use any Xwindows backend.
import matplotlib
matplotlib.use('Agg')

# Common Python packages
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Python netcdf handler
import netCDF4 as nc

# Deduce code directory (where is this file)
code_dir=os.path.dirname(os.path.realpath(__file__))
#print 'code_dir', code_dir
# Before continuing we need to be able to import python code in the general/cma_python directory.
# To do this we need to add the full path to your system path (removing any other instances of cma_python in the process).
lib_dir = os.path.join(code_dir, 'aux/AutoAssess')
print 'lib_dir', lib_dir
sys.path.append(lib_dir)

# Import extra routines that we couldn't import at the start
import utility.valmod_radiation as vm
import utility.globalvar as globalvar
import utility.rms as rms
globalvar.debug=True

def main(project_info):
    #print '>>>entering radiation_rms.py <<<<<<<<<<<<'
    # create instance of a wrapper that allows easy access to data
    E = ESMValProject(project_info)
    config_file = E.get_configfile()
    datakeys = E.get_currVars()
    plot_dir = E.get_plot_dir()
    verbosity = E.get_verbosity()

    diag_script = E.get_diag_script_name()
    print 'diag_script', diag_script

# A_laue_ax+
    res = E.write_references(diag_script,              # diag script name
                             ["A_tsus_yo"],            # authors
                             ["A_read_si"],            # contributors
                             [""],                     # diag_references
                             [""],                     # obs_references
                             ["P_cmug"],               # proj_references
                             project_info,
                             verbosity,
                             False)
# A_laue_ax-

    # Load some extra data needed for the grids and land fractions used in the RMS calculations
    print 'Loading supplementary data'
    extras_dict=vm.read_info_file(
        os.path.join(
            lib_dir,
            'extras_file.dat'
        )
    )
    vm.load_extra_data(extras_dict)

    # Initialise data_dict
    data_dict={}

    for datakey in datakeys:
        print 'datakey', datakey
        # Decide on an output directory for your CSV files
        #work_dir = os.path.join(code_dir,'../work/')
        work_dir = E.get_work_dir()
        if not os.path.exists(work_dir): os.mkdir(work_dir)
        print 'work_dir', work_dir
        summary_dir = os.path.join(work_dir,'AutoAssess_radiation_rms_summary')
        if not os.path.exists(summary_dir): os.mkdir(summary_dir)
        print 'summary_dir', summary_dir
        csv_dir = os.path.join(summary_dir,datakey)
        if not os.path.exists(csv_dir): os.mkdir(csv_dir)
        print 'csv_dir', csv_dir

        # Observations location
        obsmodel, obs_loc, models = E.get_clim_model_and_obs_filenames(datakey)
        print 'obsmodel, obs_loc, models', obsmodel, obs_loc, models
        #Output dirs based on reference data
        #csv_dir = os.path.join(code_dir,'../work/'+datakey+'/'+obsmodel)
        #if not os.path.exists(csv_dir): os.mkdir(csv_dir)

        cube_mon_obs=iris.load_cube(obs_loc)
        # A-laue_ax+
        E.add_to_filelist(obs_loc)
        # A-laue_ax-

        # supermean
        obs_variable = cube_mon_obs.collapsed(['time'], iris.analysis.MEAN)
        # We want to put this into a dictionary for passing to perform_equation
        data_dict['obs_variable']=obs_variable

        for model in models:
            # Start an rms list to hold the rms data
            print 'model', model

            model_id = model
            rms_list = rms.start(model_id,model_id)

            cube_mon_model=iris.load_cube(models[model])
            # A-laue_ax+
            E.add_to_filelist(models[model])
            # A-laue_ax-

            # supermean
            exper_variable = cube_mon_model.collapsed(['time'], iris.analysis.MEAN)

            # We want to put this into a dictionary for passing to perform_equation
            data_dict['exper_variable']=exper_variable

            # Run perform equation
            print 'Performing calculations'
            model_vs_obs_variable, rms_exper_variable = vm.perform_equation('exper_variable - obs_variable',data_dict,datakey,'a',rms_list=rms_list,calc_rms=True)


            # ************** End calculations and print out final RMS values ******************
            # Display one of the plots to show that we have reasonable figures
            # Close the RMS list and output to CSV files
            #qplt.contourf(model_vs_obs_variable, extend='both')
            #plt.gca().coastlines()
            #plt.title(datakey+'  '+obsmodel+'  '+model)
            #plt.savefig('/home/users/ytsushima/plot/'+datakey+'_'+obsmodel+'_'+model+'.png')                                               
            #plt.show()
            #plt.close()
            rms.end(rms_list,csv_dir)
            print 'Finished outputting csv files into {0}.'.format(csv_dir)

# End of main
