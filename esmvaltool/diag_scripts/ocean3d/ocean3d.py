"""
*********************************************************************
APPLICATE/TRR Ocean Diagnostics
*********************************************************************
"""
import logging
import os
import joblib
from collections import OrderedDict
import iris

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))

import inspect
from netCDF4 import Dataset
import numpy as np
import os
import matplotlib as mpl
mpl.use('agg')  #noqa
import matplotlib.pylab as plt
import math
from matplotlib import cm
from netCDF4 import num2date
#import seawater as sw
from collections import OrderedDict
from cdo import Cdo
import cmocean.cm as cmo
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import addcyclic
import pandas as pd
import pyresample
from scipy.interpolate import interp1d
import ESMF
import pyproj
#import seaborn as sns
import palettable
import itertools
#plt.style.context('seaborn-talk')
#sns.set_context("paper")

from esmvaltool.diag_scripts.ocean3d.utils import get_clim_model_filenames, get_fx_filenames, find_observations_name, shiftedColorMap, genfilename, timmean 
from esmvaltool.diag_scripts.ocean3d.regions import hofm_regions, transect_points
from esmvaltool.diag_scripts.ocean3d.getdata import hofm_data, meanminmax_data, transect_data, aw_core, tsplot_data
from esmvaltool.diag_scripts.ocean3d.plotting import hofm_plot, hofm_plot2, tsplot_plot, plot_profile, meanminmax_plot, plot2d_original_grid, plot2d_bias, plot2d_speed, plot2d_original_grid_AWdepth, plot_aw_core_stat, transect_plot, plot2d_bias2


def main(cfg):
    """Compute the time average for each input model."""
    # print(cfg['testing'])
    joblib.dump(cfg, 'cfg_NK.joblib')
    plotdir = cfg['plot_dir']
    log_level = cfg['log_level']
    plot_type = cfg['output_file_type']
    diagworkdir = cfg['work_dir']
    # diagworkdir = '/mnt/lustre01/work/ab0995/a270088/ESMV2/DUMP'
    diagplotdir = cfg['plot_dir']

    logger.info("Starting APPLICATE Arctic Ocean diagnostics")

    model_filenames_thetao = get_clim_model_filenames(cfg, 'thetao')
    model_filenames_thetao = OrderedDict(sorted(model_filenames_thetao.items(), key=lambda t: t[0]))
    print(model_filenames_thetao)

    model_filenames_so = get_clim_model_filenames(cfg, 'so')
    model_filenames_so = OrderedDict(sorted(model_filenames_so.items(), key=lambda t: t[0]))
    print(model_filenames_so)

    observations = find_observations_name(cfg, 'thetao')
    
    gg = genfilename(diagworkdir, basis='arctic_ocean', variable=None,
                mmodel=None, region=None, data_type=None, extension=None)
    print(gg)
    gg = genfilename(diagworkdir, basis='arctic_ocean', variable='thetao',
                mmodel=None, region=None, data_type='hofm', extension='.npy')
    print(gg)

    ### to be replaces with the function that get fx file information
    areacello_fx = get_fx_filenames(cfg, 'thetao', 'areacello')
    print(areacello_fx)

    # Extract data from Hovmoeller diagram
    if cfg['hofm_data']:
        for hofm_var in cfg['hofm_vars']:
            model_filenames = get_clim_model_filenames(cfg, hofm_var)
            model_filenames = OrderedDict(sorted(model_filenames.items(),
                                                     key=lambda t: t[0]))
            for mmodel, region in itertools.product(model_filenames,
                                                    cfg['hofm_regions']):
                hofm_data(model_filenames,
                          mmodel, 
                          hofm_var,
                          areacello_fx,
                          cfg['hofm_depth'],
                          region,
                          diagworkdir)

                    

    
        ############# Hofm data extraction #############################
        # for mmodel in model_filenames_thetao:
        #     hofm_data(model_filenames_thetao, mmodel, 'thetao', areacello_fx,
        #             5000, 'EB', diagworkdir)
        # for mmodel in model_filenames_thetao:
        #     hofm_data(model_filenames_thetao, mmodel, 'thetao', areacello_fx,
        #             5000, 'AB', diagworkdir)
        # for mmodel in model_filenames_so:
        #     hofm_data(model_filenames_so, mmodel, 'so', areacello_fx,
        #             5000, 'EB', diagworkdir)
        # for mmodel in model_filenames_so:
        #     hofm_data(model_filenames_so, mmodel, 'so', areacello_fx,
        #             5000, 'AB', diagworkdir)
        ############# END Hofm data extraction #########################

        ############## Create time mean #################################
    # for mmodel in model_filenames_thetao:
    #     timmean(model_filenames_thetao, mmodel, 'thetao', diagworkdir,
    #             observations=observations)
    # for mmodel in model_filenames_so:
    #     timmean(model_filenames_so, mmodel, 'so', diagworkdir, 
    #     observations=observations)
    ############## END Create time mean #############################

    ############# Hofm plot ########################################
    # hofm_plot(model_filenames_thetao, 'thetao', 5000,
    #           'EB', diagworkdir, diagplotdir,
    #           levels=np.round(np.linspace(-2, 2.3, 41), 1),
    #           ncols=3, cmap=cm.Spectral_r, observations=observations)

    # hofm_plot(model_filenames_thetao, 'thetao', 5000,
    #           'AB', diagworkdir, diagplotdir,
    #           levels=np.round(np.linspace(-2, 2.3, 41), 1),
    #           ncols=3, cmap=cm.Spectral_r, observations=observations)

    # custom_cmap = shiftedColorMap(palettable.cubehelix.cubehelix3_16.mpl_colormap, start=0, midpoint=0.89, stop=0.9, name='shiftedcmap')

    # hofm_plot(model_filenames_so, 'so', 5000,
    #           'EB', diagworkdir, diagplotdir,
    #           levels=np.round(np.linspace(30.5, 35.1, 47), 2),
    #           ncols=3, cmap=custom_cmap, observations=observations)

    # hofm_plot(model_filenames_so, 'so', 5000,
    #           'AB', diagworkdir, diagplotdir,
    #           levels=np.round(np.linspace(29, 36.5, 41), 1),
    #           ncols=3, cmap=cm.Spectral_r, observations=observations)
    ############# END Hofm plot ########################################

    ################# Plot Hofm anomalies ####################
    # hofm_plot2(model_filenames_thetao, 'thetao', 1500,
    #           'EB', diagworkdir, diagplotdir,
    #           levels=np.round(np.linspace(-0.5, 0.5, 41), 2),
    #           ncols=3, cmap=cmo.balance, observations=observations)

    # hofm_plot2(model_filenames_thetao, 'thetao', 1500,
    #           'AB', diagworkdir, diagplotdir,
    #           levels=np.round(np.linspace(-0.5, 0.5, 41), 2),
    #           ncols=3, cmap=cmo.balance, observations=observations)

    # hofm_plot2(model_filenames_so, 'so', 1500,
    #           'EB', diagworkdir, diagplotdir,
    #           levels=np.round(np.linspace(-0.1, 0.1, 41), 3),
    #           ncols=3, cmap=cmo.balance, observations=observations)

    # hofm_plot2(model_filenames_so, 'so', 1500,
    #           'AB', diagworkdir, diagplotdir,
    #           levels=np.round(np.linspace(-0.1, 0.1, 41), 3),
    #           ncols=3, cmap=cmo.balance, observations=observations)
    ################# END Plot Hofm anomalies ####################


    ############# plot vertical profiles ##########################################
    # plot_profile(model_filenames_thetao, 'thetao',
    #              5000, 'EB', diagworkdir, diagplotdir,
    #              cmap=cm.Set2, dpi=100, observations=observations)
    # plot_profile(model_filenames_thetao, 'thetao',
    #              5000, 'AB', diagworkdir, diagplotdir,
    #              cmap=cm.Set2, dpi=100, observations=observations)

    # plot_profile(model_filenames_so, 'so',
    #              5000, 'EB', diagworkdir, diagplotdir,
    #              cmap=cm.Set2, dpi=100, observations=observations)
    # plot_profile(model_filenames_so, 'so',
    #              5000, 'AB', diagworkdir, diagplotdir,
    #              cmap=cm.Set2, dpi=100, observations=observations)
############# END vertical plot profiles ########################################## 

#    ############# plot 2d original grid ###########################
    # for ddepth in [10, 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000 ]:
    #     plot2d_original_grid(model_filenames_thetao, 'thetao', 
    #                      ddepth, 'AO', diagworkdir, diagplotdir, dpi=100)
#     ############# END plot 2d original grid #######################
    # for ddepth in [10, 100, 200, 300, 400, 500 ]:
    # for ddepth in [10]:
        # plot2d_bias(model_filenames_thetao, 'thetao', 
        #         ddepth, 'AO', diagworkdir, diagplotdir, contours=[-3, 3, 21], dpi=100, observations = observations)
    # for ddepth in [300]:
    # for ddepth in [10]:
        # plot2d_bias2(model_filenames_thetao, 'thetao', 
        #         ddepth, 'AO', diagworkdir, diagplotdir, contours=[-3, 3, 21], dpi=100, observations = observations)

#     for mmodel in model_filenames_thetao:
#         transect_data(mmodel, 'thetao', 1500,
#                   'AWpath', diagworkdir, observations=observations)
    
#     transect_plot(model_filenames_thetao, 'thetao',
#               1500, 'AWpath', diagworkdir, diagplotdir,
#               levels=np.round(np.linspace(0, 5, 21), 1), ncols=3, cmap=cm.Spectral_r)
    
    ############# Calculate AW parameters and plot #########################
    # aw_core_parameters = aw_core(model_filenames_thetao, diagworkdir, 'EB', 'thetao')
    # plot_aw_core_stat(aw_core_parameters, diagplotdir)
    # plot2d_original_grid_AWdepth(model_filenames_thetao, 'thetao', 
    #                      aw_core_parameters, 'AO', diagworkdir, diagplotdir,
    #                      dpi=100, observations=observations)
    ############# END Calculate AW parameters and plot #########################
    # tsplot_data_clim(1500, 'EB', diagworkdir)
    # for mmodel in model_filenames_thetao:
    #     tsplot_data(mmodel, 1500, 'EB', diagworkdir, observations=observations)
    # tsplot_plot(model_filenames_thetao, 1500, 'EB', diagworkdir, diagplotdir,
    #             ncols=3, cmap=cm.Set1, observations = observations)
# def main(project_info):

    
    # tsplot_data_clim(1500, 'EB', diagworkdir)
    # for mmodel in model_filenames:
    #     tsplot_data(mmodel, 1500, 'EB', diagworkdir)
    # tsplot_plot(model_filenames, 1500, 'EB', diagworkdir, diagplotdir,
    #             ncols=3, cmap=cm.Set1)

    # tsplot_data_clim(1500, 'AB', diagworkdir)
    # for mmodel in model_filenames:
    #     tsplot_data(mmodel, 1500, 'AB', diagworkdir)
    # tsplot_plot(model_filenames, 1500, 'AB', diagworkdir, diagplotdir,
    #             ncols=3, cmap=cm.Set1)
       
    # for mmodel in model_filenames:
    #     meanminmax_data(model_filenames, mmodel, 'thetao', areacello_fx,
    #             0, 100, 'Barents_sea', diagworkdir)

    # meanminmax_plot(model_filenames,'thetao',
    #           -1.5, 19, 7, 15, 'Barents_sea' , diagworkdir, diagplotdir,
    #           (-3, 20), ncols=3)

    # for mmodel in model_filenames:
    #     meanminmax_data(model_filenames, mmodel, 'thetao', areacello_fx,
    #             0, 50, 'North_sea', diagworkdir)

    # meanminmax_plot(model_filenames,'thetao',
    #           -1.5, 19, 7, 15, 'North_sea' , diagworkdir, diagplotdir,
    #           (-3, 23), ncols=3)
    # phc_dict = OrderedDict([('PHC3','./')])
    # phc_dict.update(model_filenames)

    # plot_hist(phc_dict, 'thetao', 0, 100, 'Barents_sea', 
    #              range(-2, 15), diagworkdir, diagplotdir, ncols = 4, dpi=100)

    # plot_hist(phc_dict, 'thetao', 0, 100, 'North_sea', 
    #              range(5, 15), diagworkdir, diagplotdir, ncols = 4, dpi=100)


    ############# plot 2d speed ##################################  
    # plot2d_speed(model_filenames_u, 
    #              10, 0.05,  diagworkdir, diagplotdir, dpi=100)
    # plot2d_speed(model_filenames_u, 
    #              100, 0.05,  diagworkdir, diagplotdir, dpi=100)
    # plot2d_speed(model_filenames_u, 
    #              300, 0.03,  diagworkdir, diagplotdir, dpi=100)
    # plot2d_speed(model_filenames_u, 
    #              500, 0.01,  diagworkdir, diagplotdir, dpi=100)
    # ############# END plot 2d speed ##############################

    ############# Calculate AW parameters and plot #########################
    # aw_core_parameters = aw_core(model_filenames, diagworkdir, 'EB', 'thetao')
    # print(aw_core_parameters)
    # plot_aw_core_stat(aw_core_parameters, diagplotdir)
    # plot2d_original_grid_AWdepth(model_filenames, 'thetao', 
    #                      aw_core_parameters, 'AO', diagworkdir, diagplotdir, dpi=100)
    ############# END Calculate AW parameters and plot #########################


  
 
if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)