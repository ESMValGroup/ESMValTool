# -*- coding: utf-8 -*-
"""Script to calculate Arctic Ocean diagnostics.

Description
-----------
The main focus of this diagnostics is evaluation of ocean components
 of climate models in the Arctic Ocean, however most of the diagnostics
 are implemented in a way that can be easily expanded to other parts
 of the World Ocean. Most of the diagnostics aim at model comparison
 to climatological data (PHC3), so we target historical CMIP simulations.
 However scenario runs also can be analysed to have an impression
 of how Arcti Ocean hydrography will chnage in the future.

Author
------
Nikolay Koldunov (MARUM/AWI, Germany)

Project
-------
TRR181/APPLICATE

Configuration options in recipe
-------------------------------
See documentation

"""
import cartopy.crs as ccrs
import itertools
import logging
import os
from collections import OrderedDict

import cartopy.crs as ccrs
import matplotlib as mpl
from matplotlib import cm
import numpy as np
import pickle

from esmvaltool.diag_scripts.arctic_ocean.getdata import (aw_core, hofm_data,
                                                          transect_data,
                                                          tsplot_data)
from esmvaltool.diag_scripts.arctic_ocean.plotting import (
    hofm_plot, plot2d_bias, plot2d_original_grid, plot_aw_core_stat,
    plot_profile, transect_map, transect_plot, tsplot_plot)
from esmvaltool.diag_scripts.arctic_ocean.utils import (
    find_observations_name, get_clim_model_filenames, get_cmap,
    get_fx_filenames, timmean)
from esmvaltool.diag_scripts.shared import run_diagnostic

mpl.use('agg')  # noqa

logger = logging.getLogger(os.path.basename(__file__))


def run_hofm_data(cfg):
    """Extract data for Hovmoeller diagrams.

    Parameters
    ----------
    cfg: dict
        configuration dictionary ESMValTool format.
    areacello_fx : dict
        configuration dictionary with names of the areacello_fx
        files associated to dataset names.
    diagworkdir: str
        path to the diagnostic work directory.
    """
    logger.info("The `hofm_data` is True, going \
                 to extract monthly values for `hofm_regions`")

    logger.info("`hofm_vars` are: %s", cfg['hofm_vars'])
    # doing the loop for every variable
    for hofm_var in cfg['hofm_vars']:
        logger.info("Processing %s", hofm_var)
        # get dictionary with model names as key and path to the
        # preprocessed file as a value
        model_filenames = get_clim_model_filenames(cfg, hofm_var)
        model_filenames = OrderedDict(
            sorted(model_filenames.items(), key=lambda t: t[0]))
        # loop over regions and models
        for mmodel, region in itertools.product(model_filenames,
                                                cfg['hofm_regions']):
            # actual extraction of the data for specific model and region
            hofm_data(cfg, model_filenames, mmodel, hofm_var, region)


def hofm_plot_params(cfg, hofm_var, var_number, observations):
    model_filenames = get_clim_model_filenames(cfg, hofm_var)
    model_filenames = OrderedDict(
        sorted(model_filenames.items(), key=lambda t: t[0]))
    # remove "model" that contain observations,
    # since there will be no monthly data
    model_filenames = model_filenames.copy()
    if observations:
        del model_filenames[observations]
    # set the color map if not default
    if cfg['hofm_cmap']:
        cmap = get_cmap(cfg['hofm_cmap'][var_number])
    else:
        cmap = get_cmap('Spectral_r')
    # set the number of columns in the output figure
    # if defined
    if cfg['hofm_ncol']:
        ncols = cfg['hofm_ncol']
    else:
        ncols = 3
    # get the levels for plots of this variable
    vmin, vmax, sstep, roundlimit = cfg['hofm_limits'][var_number]
    plot_params = {}
    plot_params['variable'] = hofm_var
    plot_params['model_filenames'] = model_filenames
    plot_params['cmap'] = cmap
    plot_params['ncols'] = ncols
    plot_params['levels'] = np.round(np.linspace(vmin, vmax, sstep),
                                     roundlimit)
    plot_params['observations'] = observations

    return plot_params


def run_hofm_plot(cfg, observations):
    """Plot Hovmoeller diagrams for each variable.

    Parameters
    ----------
    cfg: dict
        configuration dictionary ESMValTool format.
    observations: str
        name of the observation data set
    """
    # loop over variables
    for var_number, hofm_var in enumerate(cfg['hofm_vars']):
        plot_params = hofm_plot_params(cfg, hofm_var, var_number, observations)

        # loop over models and regions
        for region in cfg['hofm_regions']:
            logger.info("Plotting Hovmoeller: for Region: %s, Variable %s",
                        region, hofm_var)
            plot_params['region'] = region
            hofm_plot(cfg, plot_params)


def run_mean(cfg, diagworkdir, observations):
    """Create time mean.

    Parameters
    ----------
    cfg: dict
        configuration dictionary ESMValTool format.
    diagworkdir: str
        path to the diagnostic work directory.
    observations: str
        name of the observation data set
    """
    # loop over variables
    for hofm_var in cfg['hofm_vars']:
        model_filenames = get_clim_model_filenames(
            cfg,
            hofm_var,
        )
        model_filenames = OrderedDict(
            sorted(model_filenames.items(), key=lambda t: t[0]))
        # loop over models
        for model in model_filenames:
            timmean(cfg,
                    model_filenames,
                    model,
                    hofm_var,
                    observations=observations)


def plot_profile_params(cfg, hofm_var, observations):
    model_filenames = get_clim_model_filenames(cfg, hofm_var)
    model_filenames = OrderedDict(
        sorted(model_filenames.items(), key=lambda t: t[0]))

    plot_params = {}

    plot_params['variable'] = hofm_var
    plot_params['model_filenames'] = model_filenames
    plot_params['cmap'] = cm.Set2
    plot_params['dpi'] = 100
    plot_params['observations'] = observations

    return plot_params


def run_profiles(cfg, diagworkdir, diagplotdir, observations):
    """Plot average vertical profiles for regions.

    Parameters
    ----------
    cfg: dict
        configuration dictionary ESMValTool format.
    diagworkdir: str
        path to the diagnostic work directory.
    diagplotdir: str
        path to the diagnostic plot directory.
    observations: str
        name of the observation data set
    """
    # loop over variables
    for hofm_var in cfg['hofm_vars']:
        plot_params = plot_profile_params(cfg, hofm_var, observations)
        # loop over regions
        for region in cfg['hofm_regions']:
            plot_params['region'] = region
            plot_profile(cfg, plot_params)


def plot2d_params(cfg, plot2d_var, var_number):
    model_filenames = get_clim_model_filenames(cfg, plot2d_var)
    model_filenames = OrderedDict(
        sorted(model_filenames.items(), key=lambda t: t[0]))
    # set the color map
    if cfg['plot2d_cmap']:
        cmap = get_cmap(cfg['plot2d_cmap'][var_number])
    else:
        cmap = get_cmap('Spectral_r')
    # set the number of columns in the plot
    if cfg['plot2d_ncol']:
        ncols = cfg['plot2d_ncol']
    else:
        ncols = 4
    # create color limits for the plot
    vmin, vmax, sstep, roundlimit = cfg['plot2d_limits'][var_number]
    # loop over depths
    plot_params = {}
    plot_params['variable'] = plot2d_var
    plot_params['model_filenames'] = model_filenames
    plot_params['cmap'] = cmap
    plot_params['ncols'] = ncols
    plot_params['levels'] = np.round(np.linspace(vmin, vmax, sstep),
                                     roundlimit)
    plot_params['dpi'] = 100
    plot_params['explicit_depths'] = None
    plot_params['projection'] = ccrs.NorthPolarStereo()
    plot_params['bbox'] = (-180, 180, 60, 90)

    return plot_params


def run_plot2d(cfg, diagworkdir, diagplotdir):
    """Plot 2d maps on original grid.

    Parameters
    ----------
    cfg: dict
        configuration dictionary ESMValTool format.
    diagworkdir: str
        path to the diagnostic work directory.
    diagplotdir: str
        path to the diagnostic plot directory.
    """
    # loop over variables
    for var_number, plot2d_var in enumerate(cfg['plot2d_vars']):

        plot_params = plot2d_params(cfg, plot2d_var, var_number)

        for depth in cfg['plot2d_depths']:
            plot_params['depth'] = depth
            plot2d_original_grid(cfg, plot_params)



def run_plot2d_bias(cfg, diagworkdir, diagplotdir, observations):
    """Plot model biases over depth.

    Parameters
    ----------
    cfg: dict
        configuration dictionary ESMValTool format.
    diagworkdir: str
        path to the diagnostic work directory.
    diagplotdir: str
        path to the diagnostic plot directory.
    observations: str
        name of the observation data set
    """
    # loop over variables
    for var_number, plot2d_bias_var in enumerate(cfg['plot2d_bias_vars']):

        model_filenames = get_clim_model_filenames(cfg, plot2d_bias_var)
        model_filenames = OrderedDict(
            sorted(model_filenames.items(), key=lambda t: t[0]))
        # setup the color map
        if cfg['plot2d_bias_cmap']:
            cmap = get_cmap(cfg['plot2d_bias_cmap'][var_number])
        else:
            cmap = get_cmap('Spectral_r')
        # setup the number of columns
        if cfg['plot2d_bias_ncol']:
            ncols = cfg['plot2d_bias_ncol']
        else:
            ncols = 3
        # setup color limits
        vmin, vmax, sstep, roundlimit = cfg['plot2d_bias_limits'][var_number]
        # loop over depths
        for depth in cfg['plot2d_bias_depths']:
            plot2d_bias(model_filenames,
                        plot2d_bias_var,
                        depth,
                        diagworkdir,
                        diagplotdir,
                        cmap,
                        levels=np.round(np.linspace(vmin, vmax, sstep),
                                        roundlimit),
                        dpi=100,
                        observations=observations,
                        ncols=ncols)


def run_transects(cfg, diagworkdir, diagplotdir):
    """Plot transects.

    Parameters
    ----------
    cfg: dict
        configuration dictionary ESMValTool format.
    diagworkdir: str
        path to the diagnostic work directory.
    diagplotdir: str
        path to the diagnostic plot directory.
    """
    # First plot the map woth transect points for each "region"
    for region in cfg['transects_regions']:
        transect_map(region,
                     diagplotdir,
                     projection=ccrs.NorthPolarStereo(),
                     bbox=[-180, 180, 60, 90],
                     mult=2)
    # loop over variables
    for var_number, trans_var in enumerate(cfg['transects_vars']):
        model_filenames = get_clim_model_filenames(cfg, trans_var)
        model_filenames = OrderedDict(
            sorted(model_filenames.items(), key=lambda t: t[0]))
        # loop over regions
        for mmodel, region in itertools.product(model_filenames,
                                                cfg['transects_regions']):
            # ploting a transect
            transect_data(cfg, mmodel, trans_var, region)
        # setup a color map
        if cfg['transects_cmap']:
            cmap = get_cmap(cfg['transects_cmap'][var_number])
        else:
            cmap = get_cmap('Spectral_r')
        # setup number of columns
        if cfg['transects_ncol']:
            ncols = cfg['transects_ncol']
        else:
            ncols = 3
        # setup color limits
        vmin, vmax, sstep, roundlimit = cfg['transects_limits'][var_number]
        # loop over regions
        for region in cfg['transects_regions']:
            transect_plot(model_filenames,
                          trans_var,
                          cfg['transects_depth'],
                          region,
                          diagworkdir,
                          diagplotdir,
                          levels=np.round(np.linspace(vmin, vmax, sstep),
                                          roundlimit),
                          ncols=ncols,
                          cmap=cmap)


def run_aw_core(cfg, diagworkdir, diagplotdir):
    """Calculate depth and temperature of the Atlantic Water core.

    Parameters
    ----------
    cfg: dict
        configuration dictionary ESMValTool format.
    diagworkdir: str
        path to the diagnostic work directory.
    diagplotdir: str
        path to the diagnostic plot directory.
    """
    model_filenames = get_clim_model_filenames(cfg, 'thetao')
    model_filenames = OrderedDict(
        sorted(model_filenames.items(), key=lambda t: t[0]))
    aw_core_parameters = aw_core(model_filenames, diagworkdir, 'EB', 'thetao')
    plot_aw_core_stat(aw_core_parameters, diagplotdir)
    return aw_core_parameters


def run_aw_core_2d(cfg, diagworkdir, diagplotdir, aw_core_parameters):
    """Plot temperature spatial distribution at AW core depth.

    Parameters
    ----------
    cfg: dict
        configuration dictionary ESMValTool format.
    diagworkdir: str
        path to the diagnostic work directory.
    diagplotdir: str
        path to the diagnostic plot directory.
    aw_core_parameters: dict
        dictionary that contain AW core parameters generated
        by run_aw_core function.
    """
    model_filenames = get_clim_model_filenames(cfg, 'thetao')
    model_filenames = OrderedDict(
        sorted(model_filenames.items(), key=lambda t: t[0]))
    aw_core_parameters = aw_core(model_filenames, diagworkdir, 'EB', 'thetao')
    # this is now just using plot2d_original_grid with
    # additional `explicit_depths` parameter
    plot_params = {}
    plot_params['variable'] = 'thetao'
    plot_params['model_filenames'] = model_filenames
    plot_params['depth'] = 0
    plot_params['cmap'] = cm.Spectral_r
    plot_params['ncols'] = 4
    plot_params['levels'] = np.round(np.linspace(-2, 2.3, 41), 1)
    plot_params['dpi'] = 100
    plot_params['explicit_depths'] = aw_core_parameters
    plot_params['projection'] = ccrs.NorthPolarStereo()
    plot_params['bbox'] = (-180, 180, 60, 90)

    plot2d_original_grid(cfg, plot_params)



def run_tsdiag(cfg, diagworkdir, diagplotdir, observations):
    """Plot TS diagrams.

    Parameters
    ----------
    cfg: dict
        configuration dictionary ESMValTool format.
    diagworkdir: str
        path to the diagnostic work directory.
    diagplotdir: str
        path to the diagnostic plot directory.
    observations: str
        name of the observation data set
    """
    # get the dictionary with model file names
    model_filenames = get_clim_model_filenames(cfg, 'thetao')
    model_filenames = OrderedDict(
        sorted(model_filenames.items(), key=lambda t: t[0]))
    # loop over models and regions
    for mmodel, region in itertools.product(model_filenames,
                                            cfg['tsdiag_regions']):
        # this function will generate files with T and S points
        # selected from the region untill `tsdiag_depth` for
        # every model.
        tsplot_data(cfg, mmodel, region, observations=observations)
    # setting the number of columns for the plot
    if cfg['tsdiag_ncol']:
        ncols = cfg['tsdiag_ncol']
    else:
        ncols = 3
    # actually plot TS diagrams
    for region in cfg['tsdiag_regions']:
        tsplot_plot(
            model_filenames,
            cfg['tsdiag_depth'],
            region,
            diagworkdir,
            diagplotdir,
            ncols=ncols,
            cmap=cm.Set1,
        )


def main(cfg):
    """Compute the time average for each input model."""
    # for debuging save the configuration in a pickle file
    with open('cfg_NK.joblib', 'wb') as handle:
        pickle.dump(cfg, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # sets the directories
    diagworkdir = cfg['work_dir']
    diagplotdir = cfg['plot_dir']

    logger.info("Starting APPLICATE/TRR Arctic Ocean diagnostics")

    # find the name of the observational dataset
    observations = find_observations_name(cfg)
    logger.info("Name of the Observations: %s{}", observations)

    # get the names of fx filenames (for now are the same for
    # all variables (this is why "thetao" is hardcoded))
    areacello_fx = get_fx_filenames(cfg, 'areacello')
    logger.info("areacello_fx files: %s", areacello_fx)

    # Extract data for Hovmoeller diagrams
    # run_hofm_data(cfg)

    # # Plot Hovmoeller diagrams for each variable
    # run_hofm_plot(cfg, observations)

    # # # Create timemean
    # run_mean(cfg, diagworkdir, observations)

    # # Plot average vertical profiles for regions
    # run_profiles(cfg, diagworkdir, diagplotdir, observations)

    # Plot 2d maps on original grid
    run_plot2d(cfg, diagworkdir, diagplotdir)

    # # Plot model biases over depth
    # run_plot2d_bias(cfg, diagworkdir, diagplotdir, observations)

    # # Plot transects
    # run_transects(cfg, diagworkdir, diagplotdir)

    # # Calculate depth and temperature of the Atlantic Water core
    # # and make plots.
    aw_core_parameters = run_aw_core(cfg, diagworkdir, diagplotdir)

    # # Plot temperature spatial distribution at the depth of the
    # # atlantic water core in different models
    run_aw_core_2d(cfg, diagworkdir, diagplotdir, aw_core_parameters)

    # # Plot TS diagrams
    # run_tsdiag(cfg, diagworkdir, diagplotdir, observations)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
