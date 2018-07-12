"""
;; Original Description from Version 1 Diagnostic
;; Ported to Version 2 with implementation of v2-specific changes
;; Valeriu Predoi, UREAD, July 2018
;;###########################################################################
;; AutoAssess_radiation_rms.py
;; Author: Yoko Tsushima (Met Office, UK)
;; CMUG project
;;###########################################################################
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
;; # Latitudes [-90, 90 degrees]: 10S = -10, 10N = 10;
;; longtitudes [0, 360 degrees]
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
;; colourmap_clouds: Python matplotlib colormaps, '_r' reversed colormap
;; colourmap_model:  Python matplotlib colormaps, '_r' reversed colormap
;; colourmap_diff:   Python matplotlib colormaps, '_r' reversed colormap
;; colourmap_dev:    Python matplotlib colormaps, '_r' reversed colormap
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
;;   landsea.nc from NCLA for land/sea mask and resolution.
;;   This does not have coord system info.
;;
;; Modification history
;;    20180712- autoassess_radiation_rms: porting to v2
;;    20170323-_AutoAssess_radiation_rms: Test finished.
;;    20160819-_test_AutoAssess_radiation_rms: written based on calc_rms code.
;;
;; ###########################################################################
"""

import os
import logging
import sys
import numpy as np
import iris
from esmvaltool.diag_scripts.autoassess.autoassess_source import rms
from esmvaltool.diag_scripts.autoassess.autoassess_source import valmod_radiation as vm
from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)

logger = logging.getLogger(os.path.basename(__file__))

# dict of model var-to-obs name
VAR_DICT = {'pr': 'MPI-ESM-MR'}


def main(cfg):
    """Execute the radiation rms diag"""
    logger.setLevel(cfg['log_level'].upper())
    if not os.path.exists(cfg['plot_dir']):
        os.makedirs(cfg['plot_dir'])
    if not os.path.exists(cfg['work_dir']):
        os.makedirs(cfg['work_dir'])


    # Load some extra data needed for the grids and land fractions used in the RMS calculations
    #print 'Loading supplementary data'
    #extras_dict=vm.read_info_file(
    #    os.path.join(
    #        lib_dir,
    #        'extras_file.dat'
    #    )
    #)
    #vm.load_extra_data(extras_dict)

    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(
        input_data, 'short_name', sort='dataset')

    # select variables and their corresponding
    # obs files
    for short_name in grouped_input_data:
        dataset_selection = select_metadata(input_data,
                                            short_name=short_name,
                                            project='CMIP5')
        obs_selection = select_metadata(input_data,
                                        short_name=short_name,
                                        dataset=VAR_DICT[short_name])
        logger.info("Processing variable %s", short_name)
        # this is a length=1 list
        obs = obs_selection[0]
        for model in dataset_selection:
            logger.info("Processing dataset %s", model['dataset'])
            logger.info("Processing obs data %s", obs['dataset'])

            # load files
            model_file = model['filename']
            obs_file = obs['filename']
            logger.info("Loading %s", model_file)
            logger.info("Loading %s", obs_file)

            # load cubes
            model_cube = iris.load_cube(model_file)
            obs_cube = iris.load_cube(obs_file)
            model_cube = model_cube.collapsed('time', iris.analysis.MEAN)
            logger.debug("Time-averaged %s", model_cube)
            obs_cube = obs_cube.collapsed('time', iris.analysis.MEAN)
            logger.debug("Time-averaged %s", obs_cube)

            # apply rms
            rms_list = rms.start(model['dataset'], obs['dataset'])
            data_dict['exper_variable'] = model_cube
            data_dict['obs_variable'] = obs_cube
            model_vs_obs_variable, rms_exper_variable = vm.perform_equation('exper_variable - obs_variable',data_dict,datakey,'a',rms_list=rms_list,calc_rms=True)
            rms.end(rms_list, cfg['work_dir'])

if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
