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

from esmvaltool.diag_scripts.ocean3d.utils import get_clim_model_filenames, get_fx_filenames, find_observations_name, get_cmap, genfilename, timmean 
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

    # model_filenames_thetao = get_clim_model_filenames(cfg, 'thetao')
    # model_filenames_thetao = OrderedDict(sorted(model_filenames_thetao.items(), key=lambda t: t[0]))
    # print(model_filenames_thetao)

    # model_filenames_so = get_clim_model_filenames(cfg, 'so')
    # model_filenames_so = OrderedDict(sorted(model_filenames_so.items(), key=lambda t: t[0]))
    # print(model_filenames_so)

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

    # Extract data for Hovmoeller diagrams
    if cfg['hofm_data']:
        for hofm_var in cfg['hofm_vars']:
            print(cfg['hofm_vars'])
            print(hofm_var)
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

    # Plot Hovmoeller diagrams for each variable
    if cfg['hofm_plot']:
        for var_number, hofm_var in enumerate(cfg['hofm_vars']):
            model_filenames = get_clim_model_filenames(cfg, hofm_var)
            model_filenames = OrderedDict(sorted(model_filenames.items(),
                                                     key=lambda t: t[0]))
            if cfg['hofm_cmap']:
                cmap = get_cmap(cfg['hofm_cmap'][var_number])
            else:
                cmap = get_cmap('Spectral_r')
            
            if cfg['hofm_ncol']:
                ncols = cfg['hofm_ncol']
            else:
                ncols = 3

            vmin, vmax, sstep, roundlimit = cfg['hofm_limits'][var_number]

            for mmodel, region in itertools.product(model_filenames,
                                                    cfg['hofm_regions']):   
                print(mmodel)
                print(region)
                print(hofm_var)                
                hofm_plot(model_filenames,
                          hofm_var, 
                          cfg['hofm_depth'],
                          region,
                          diagworkdir,
                          diagplotdir,
                          levels=np.round(np.linspace(vmin, vmax, sstep),
                                          roundlimit),
                          ncols=ncols, 
                          cmap=cmap, 
                          observations=observations)

    # Create timemean (to be replaced by ESMValTool function)
    if cfg['mean']:
        for hofm_var in cfg['hofm_vars']:
            model_filenames = get_clim_model_filenames(cfg, hofm_var)
            model_filenames = OrderedDict(sorted(model_filenames.items(),
                                                     key=lambda t: t[0]))
            for model in model_filenames:
                timmean(model_filenames, model, hofm_var, diagworkdir,
                        observations=observations)
    

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

    if cfg['profiles']:
        for var_number, hofm_var in enumerate(cfg['hofm_vars']):
            model_filenames = get_clim_model_filenames(cfg, hofm_var)
            model_filenames = OrderedDict(sorted(model_filenames.items(),
                                                     key=lambda t: t[0]))
            for region in cfg['hofm_regions']:
                plot_profile(model_filenames,
                             hofm_var,
                             cfg['hofm_depth'],
                             region,
                             diagworkdir,
                             diagplotdir,
                             cmap=cm.Set2,
                             dpi=100,
                             observations=observations)

    if cfg['plot2d']:
        for var_number, plot2d_var in enumerate(cfg['plot2d_vars']):
            model_filenames = get_clim_model_filenames(cfg, plot2d_var)
            model_filenames = OrderedDict(sorted(model_filenames.items(),
                                                     key=lambda t: t[0]))
            if cfg['plot2d_cmap']:
                cmap = get_cmap(cfg['plot2d_cmap'][var_number])
            else:
                cmap = get_cmap('Spectral_r')
            
            if cfg['plot2d_ncol']:
                ncols = cfg['plot2d_ncol']
            else:
                ncols = 3

            vmin, vmax, sstep, roundlimit = cfg['plot2d_limits'][var_number]
            for depth in cfg['plot2d_depths']:
                plot2d_original_grid(model_filenames,
                                     plot2d_var, 
                                     depth,
                                     levels = np.round(np.linspace(vmin,
                                                                   vmax,
                                                                   sstep),
                                                       roundlimit),
                                     region = 'AO',
                                     diagworkdir = diagworkdir,
                                     diagplotdir = diagplotdir,
                                     dpi=100)

    if cfg['plot2d_bias']:
        for var_number, plot2d_bias_var in enumerate(cfg['plot2d_bias_vars']):
            model_filenames = get_clim_model_filenames(cfg, plot2d_bias_var)
            model_filenames = OrderedDict(sorted(model_filenames.items(),
                                                     key=lambda t: t[0]))
            if cfg['plot2d_bias_cmap']:
                cmap = get_cmap(cfg['plot2d_bias_cmap'][var_number])
            else:
                cmap = get_cmap('Spectral_r')
            
            if cfg['plot2d_bias_ncol']:
                ncols = cfg['plot2d_bias_ncol']
            else:
                ncols = 3

            vmin, vmax, sstep, roundlimit = cfg['plot2d_bias_limits'][var_number]
            for depth in cfg['plot2d_bias_depths']:
                plot2d_bias(model_filenames,
                            plot2d_bias_var, 
                            depth,
                            'AO',
                            diagworkdir,
                            diagplotdir,
                            contours=np.round(np.linspace(vmin,
                                              vmax,
                                              sstep),
                                              roundlimit),
                            dpi=100,
                            observations=observations)

    if cfg['transects']:
        for var_number, trans_var in enumerate(cfg['transects_vars']):
            model_filenames = get_clim_model_filenames(cfg, hofm_var)
            model_filenames = OrderedDict(sorted(model_filenames.items(),
                                                 key=lambda t: t[0]))
            for mmodel, region in itertools.product(model_filenames,
                                                    cfg['transects_regions']):
                transect_data(mmodel,
                              trans_var,
                              cfg['transects_depth'],
                              region,
                              diagworkdir, 
                              observations=observations)

            if cfg['transects_cmap']:
                cmap = get_cmap(cfg['transects_cmap'][var_number])
            else:
                cmap = get_cmap('Spectral_r')
            
            if cfg['transects_ncol']:
                ncols = cfg['transects_ncol']
            else:
                ncols = 3

            vmin, vmax, sstep, roundlimit = cfg['transects_limits'][var_number]
            for region in cfg['transects_regions']:
                transect_plot(model_filenames, 
                              trans_var,
                              cfg['transects_depth'],
                              region, diagworkdir, diagplotdir,
                              levels=np.round(np.linspace(vmin,
                                              vmax,
                                              sstep),
                                              roundlimit),
                              ncols=ncols, 
                              cmap=cmap)

    # Will change to more general definition of the water mass core
    if cfg['AW_core']:
        model_filenames = get_clim_model_filenames(cfg, 'thetao')
        model_filenames = OrderedDict(sorted(model_filenames.items(),
                                                 key=lambda t: t[0]))
        aw_core_parameters = aw_core(model_filenames, diagworkdir, 'EB', 'thetao')
        plot_aw_core_stat(aw_core_parameters, diagplotdir)

    if cfg['AW_core_2d']:
        model_filenames = get_clim_model_filenames(cfg, 'thetao')
        model_filenames = OrderedDict(sorted(model_filenames.items(),
                                                 key=lambda t: t[0]))
        plot2d_original_grid_AWdepth(model_filenames, 'thetao', 
                         aw_core_parameters, 'AO', diagworkdir, diagplotdir,
                         dpi=100, observations=observations)
    
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