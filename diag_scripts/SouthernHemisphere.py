"""
;;#############################################################################
;; SouthernHemisphere.py
;; Author: Jarmo Makela (FMI, Finland)
;; EMBRACE project
;;#############################################################################
;; Description
;;    This script is the diagnostics and plotting script for
;;    Southern Hemisphere radiation and fluxes maps and graphs
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
;;    See TODO in 'esmval_lib' about hardcoded CMIP5 project class
;;    assumptions.
;;
;; Modification history
;;    20151029-A_laue_ax: added output of acknowledgements + processed files
;;                        to log-file
;;    20150415-A_maek_ja: written
;;
;; #############################################################################
"""

# Basic Python packages
import ConfigParser
import os
import pdb
import sys
import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import interp1d

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

def main(project_info):
    """Diagnostics and plotting script for Southern Hemisphere radiation."""

    # ESMValProject provides some easy methods to access information in
    # project_info but you can also access it directly (as in main.py)
    # First we get some general configurations (ESMValProject)
    # and then we read in the model configuration file

    E = ESMValProject(project_info)
    config_file = E.get_configfile()
    datakeys = E.get_currVars()
    plot_dir = E.get_plot_dir()
    verbosity = E.get_verbosity()

    # A-laue_ax+
    diag_script = E.get_diag_script_name()
    res = E.write_references(diag_script,              # diag script name
                             ["A_maek_ja"],            # authors
                             ["A_eval_ma", "A_jone_co"], # contributors
                             [""],                     # diag_references
                             [""],                     # obs_references
                             ["P_embrace"],            # proj_references
                             project_info,
                             verbosity,
                             False)
    # A-laue_ax-

    modelconfig = ConfigParser.ConfigParser()
    modelconfig.read(config_file)
    E.ensure_directory(plot_dir)

    # Check at which stage of the program we are
    clouds = False
    fluxes = False
    radiation = False
    if ('clt' in datakeys or 'clivi' in datakeys or 'clwvi' in datakeys):
        clouds = True
    elif ('hfls' in datakeys or 'hfss' in datakeys):
        fluxes = True
    else:
        radiation = True

    # Check which parts of the code to run
    if (modelconfig.getboolean('general', 'plot_clouds') and clouds is True):
        info("Starting to plot clouds", verbosity, 2)
        process_clouds(E, modelconfig)

    if (modelconfig.getboolean('general', 'plot_fluxes') and fluxes is True):
        info("Starting to plot turbulent fluxes", verbosity, 2)
        process_fluxes(E, modelconfig)

    if (modelconfig.getboolean('general', 'plot_radiation') and radiation is True):
        info("Starting to plot radiation graphs", verbosity, 2)
        process_radiation(E, modelconfig)


# Further diagnostic separation and which parts of the program to run

def process_clouds(E, modelconfig):
    """Wrapper for cloud graphs"""
    experiment = 'SouthernHemisphere'
    specifier = 'cloud'
    datakeys = E.get_currVars()

    # Check settings what is to be plotted
    plot_sub_areas = modelconfig.getboolean('general', 'plot_sub_areas')
    plot_mon_avg = modelconfig.getboolean('general', 'plot_monthly_averages')
    plot_lat_avg = modelconfig.getboolean('general', 'plot_lat_averages')
    plot_lon_avg = modelconfig.getboolean('general', 'plot_lon_averages')
    plot_total_cover = modelconfig.getboolean('general', 'plot_total_cover')
    plot_liq_path = modelconfig.getboolean('general', 'plot_liquid_path')
    plot_ice_path = modelconfig.getboolean('general', 'plot_ice_path')

    # Call plots based on variable - setup is only one key at a time
    # so we should only have one variable at a time
    for datakey in datakeys:

        # First call common averaging graphs
        if (plot_lat_avg):
            process_mean_plots(E, modelconfig, datakey, 'lats', specifier)

        if (plot_lon_avg):
            process_mean_plots(E, modelconfig, datakey, 'lons', specifier)

        if (plot_mon_avg or plot_sub_areas):
            process_mean_plots(E, modelconfig, datakey, 'monthly', specifier)

        # And next call contour map plots
        if (plot_total_cover and datakey == 'clt'):
            process_simple_maps(E, modelconfig, datakey, specifier + '-map')

        if (plot_liq_path and datakey == 'clwvi'):
            process_simple_maps(E, modelconfig, datakey, specifier + '-map')

        if (plot_liq_path and datakey == 'clivi'):
            process_simple_maps(E, modelconfig, datakey, specifier + '-map')

        # FIXME: when optical depth becomes available, you can remove thei
        # quotes below and add the correct datakey to enable optical depth
        # diagnostics
        """
        if (modelconfig.getboolean('general', 'plot_optical_depth') and
            datakey == 'plot_optical_depth'):
            process_simple_maps(E, modelconfig, datakey, specifier + '_map')
        """


def process_fluxes(E, modelconfig):
    """Wrapper for turbulent fluxes"""
    experiment = 'SouthernHemisphere'
    specifier = 'flux'
    datakeys = E.get_currVars()

    # Check settings what is to be plotted
    plot_sub_areas = modelconfig.getboolean('general', 'plot_sub_areas')
    plot_mon_avg = modelconfig.getboolean('general', 'plot_monthly_averages')
    plot_lat_avg = modelconfig.getboolean('general', 'plot_lat_averages')
    plot_lon_avg = modelconfig.getboolean('general', 'plot_lon_averages')
    plot_flux_maps = modelconfig.getboolean('general', 'plot_flux_maps')

    for datakey in datakeys:
        # Call averaging plots
        if (plot_lat_avg):
            process_mean_plots(E, modelconfig, datakey, 'lats', specifier)

        if (plot_lon_avg):
            process_mean_plots(E, modelconfig, datakey, 'lons', specifier)

        if (plot_mon_avg or plot_sub_areas):
            process_mean_plots(E, modelconfig, datakey, 'monthly', specifier)

        # Call turbulent flux seasonal maps
        if (plot_flux_maps):
            process_simple_maps(E, modelconfig, datakey, specifier + '-map')


def process_radiation(E, modelconfig):
    """Wrapper for radation maps and graphs"""
    experiment = 'SouthernHemisphere'
    specifier = 'radiation'
    datakeys = E.get_currVars()
    # Separate all datakeys to clear sky variants and others
    keyscs, keys = separate_list(datakeys, 'cs')

    # Check settings what is to be plotted
    plot_sub_areas = modelconfig.getboolean('general', 'plot_sub_areas')
    plot_mon_avg = modelconfig.getboolean('general', 'plot_monthly_averages')
    plot_lat_avg = modelconfig.getboolean('general', 'plot_lat_averages')
    plot_lon_avg = modelconfig.getboolean('general', 'plot_lon_averages')
    plot_rad_maps = modelconfig.getboolean('general', 'plot_radiation_maps')

    for datakey in datakeys:
        # Call averaging plots
        if (plot_lat_avg):
            process_mean_plots(E, modelconfig, datakey, 'lats', specifier)

        if (plot_lon_avg):
            process_mean_plots(E, modelconfig, datakey, 'lons', specifier)

        if (plot_mon_avg or plot_sub_areas):
            process_mean_plots(E, modelconfig, datakey, 'monthly', specifier)

        # Call radiation seasonal maps for those that have datakey + 'cs'
        if (plot_rad_maps):
            if (datakey + 'cs' in keyscs):
                process_radiation_maps(E, modelconfig, datakey, specifier + '-map')


# All callable functions below are in alphabetical order
# E. style functions are in the general python file esmval_lib.py

def extract_seasonal_mean_values(modelconfig, data, experiment, season):
    """Returns the season specific mean values  for each lat, lon from the data
    We assume the usual indexing of time, lat, lon"""
    season_key = experiment + '_season_' + season
    data_shape = data.shape

    if (season == 'annual'):
        # For annual season we merely copy the data
        masked_values = data
    else:
        # For a specific season we mask the undesired values
        season_months = modelconfig.get(season_key, 'season_months').split()
        mask = np.ones((data_shape))

        for month in season_months:
            month_loc = int(month) - 1
            mask[month_loc::12, :] = 0
        masked_values = np.ma.masked_array(data, mask)

    mean_values = np.zeros((data_shape[1], data_shape[2]))
    # Seasonal mean values
    for lat in xrange(data_shape[1]):
        for lon in xrange(data_shape[2]):
            mean_values[lat, lon] = masked_values[:, lat, lon].mean()
    return mean_values


def get_contour_config(modelconfig, config_file, area_key, datakey):
    """Returns area specific contour limits. """
    contour_key = 'contour_limits_' + datakey
    if   (datakey in ['clt', 'clivi', 'clwvi']):
        required = 3
        cm_model = modelconfig.get(area_key, 'colourmap_clouds')
    elif (datakey in ['hfls', 'hfss']):
        required = 3
        cm_model = modelconfig.get(area_key, 'colourmap_model')
    elif (datakey in ['rlut', 'rsut', 'rlds', 'rsds']):
        required = 5
        cm_model = modelconfig.get(area_key, 'colourmap_model')

    limits = []
    # Check the existance of needed limits - otherwise we let python decide
    if modelconfig.has_option(area_key, contour_key):
        # Check the limits and thei number
        str_limits = modelconfig.get(area_key, contour_key).split()
        if (len(str_limits) == required):
            limits = np.zeros(required)
            for limit in xrange(required):
                limits[limit] = int(str_limits[limit])
        # Wrong number of limits present
        else:
            print("PY  ERROR: Misdefined contour limits for : '" + contour_key + "'")
            print("PY  ERROR: You should use " + required + " integers for the limits")
            print("PY  ERROR: Check your configuration file: " + config_file)
            print("PY  ERROR: Reverting to automated limits")
            limits = []
    # The entire section is missing
    else:
        print("PY  ERROR: Misdefined contour limits in : '[" + area_key + "]'")
        print("PY  ERROR: Couldn't find key: '" + contour_key + "'")
        print("PY  ERROR: Check your configuration file: " + config_file)
        print("PY  ERROR: Reverting to automated limits")
        limits = []

    stride = modelconfig.getint(area_key, 'stride')
    cm_diff = modelconfig.get(area_key, 'colourmap_diff')

    # if stride = 0 we define strides based on maxshades
    if (stride == 0):
        shades = modelconfig.getint(area_key, 'maxshades')
        lev = np.linspace(limits[0], limits[1], shades - 1)
        diff = np.linspace(-limits[2], limits[2], shades - 1)
    else:
        lev = np.arange(limits[0], limits[1] + 1, stride)
        diff = np.arange(-limits[2], limits[2] + 1, stride)

    # Last to check is that radiation uses more limits
    # This affects the return values
    if (required == 3):
        return lev, diff, cm_model, cm_diff
    else:
        if (limits[3] < 0):
            cm_dev = modelconfig.get(area_key, 'colourmap_dev').split()[0]
        else:
            cm_dev = modelconfig.get(area_key, 'colourmap_dev').split()[1]
        if (stride == 0):
            diff2 = np.linspace(limits[3], limits[4], shades - 1)
        else:
            diff2 = np.arange(limits[3], limits[4] + 1, stride)
        return lev, diff, diff2, cm_model, cm_diff, cm_dev


def interpolate_data_grid(data, lats, lons, target_lats, target_lons):
    """Interpolates the data values to a specific lat/lon grid.
    This function should only be used for 2D arrays (no time indeces etc.)"""

    # First check if the coordinates are the same, otherwise interpolate
    if (np.array_equal(lats, target_lats) and np.array_equal(lons, target_lons)):
        return data

    else:
        grid_lats, grid_lons = np.meshgrid(target_lats, target_lons)
        datapoints = np.zeros((len(lats) * len(lons), 2))
        datavalues = data.reshape(-1)
        i = 0
        for lat in lats:
            datapoints[i:i + len(lons), 0] = lat
            datapoints[i:i + len(lons), 1] = lons
            i += len(lons)

        # Calculate the values on the new grid
        gridvalues = griddata(datapoints,
                              datavalues,
                              (grid_lats, grid_lons),
                              method='cubic')
        nearest = griddata(datapoints,
                           datavalues,
                           (grid_lats, grid_lons),
                           method='nearest')

        # In some cases griddata cannot export values on the border of the grid
        for i in xrange(gridvalues.shape[0]):
            for j in xrange(gridvalues.shape[1]):
                if (np.isnan(gridvalues[i, j])):
                    gridvalues[i, j] = nearest[i, j]
        return gridvalues.T


def process_mean_plots(E, modelconfig, datakey, orientation, mean_name):
    """Generates latitudal / lontitudal / monthly mean plots.
    Seasonal values are extracted for lat/lon plots. """
    experiment = 'SouthernHemisphere'
    plot_dir = E.get_plot_dir()
    verbosity = E.get_verbosity()
    areas = []
    plot_grid = modelconfig.getboolean('general', 'plot_background_grid')

    # Extract required keys based on datakey
    if (orientation == 'monthly'):
        if (modelconfig.getboolean('general', 'plot_monthly_averages')):
            areas += modelconfig.get(experiment, 'areas').split()
        if (modelconfig.getboolean('general', 'plot_sub_areas')):
            areas += modelconfig.get(experiment, 'sub_areas').split()
        seasons = ['monthly']
    else:
        areas += modelconfig.get(experiment, 'areas').split()
        seasons = modelconfig.get(experiment, 'seasons').split()

    # Observations location
    obsmodel, obs_loc, models = E.get_clim_model_and_obs_filenames(datakey)
    obsfile = nc.Dataset(obs_loc, 'r')
    # A-laue_ax+
    E.add_to_filelist(obs_loc)
    # A-laue_ax-

    scale_cloud = False
    # Ugly fix for scaling cloud ice/liquid water path values
    if (datakey == 'clivi' or datakey == 'clwvi'):
        scale_cloud = True
        scale = 1E3
        scale_str = ' * 1E3'

    for area in areas:
        # Value configurations
        area_key = experiment + '_' + area
        lat_min = modelconfig.getint(area_key, 'lat_min')
        lat_max = modelconfig.getint(area_key, 'lat_max')
        lon_min = modelconfig.getint(area_key, 'lon_min')
        lon_max = modelconfig.getint(area_key, 'lon_max')
        lat_area = E.get_ticks_labels(np.array([lat_min, lat_max]), 'lats')

        # First extract the observations and then start looping over seasons
        obslats, obslons, obsdata = E.get_model_data(modelconfig,
                                                     experiment,
                                                     area,
                                                     datakey,
                                                     obsfile)
        for season in seasons:
            # Plot layout configuration
            plt.clf()
            fig, axs = plt.subplots(2, 1, figsize=(15, 10))
            fig.subplots_adjust(top=0.9)
            fig.subplots_adjust(right=0.67)
            fig.subplots_adjust(hspace=0.35)

            if   (orientation == 'lats'):
                odata = E.extract_seasonal_mean_values(modelconfig,
                                                       obsdata,
                                                       experiment,
                                                       season,
                                                       monthly=True)
                odata = E.average_data(odata, 1)
                xobs = obslats
            elif (orientation == 'lons'):
                odata = E.extract_seasonal_mean_values(modelconfig,
                                                       obsdata,
                                                       experiment,
                                                       season,
                                                       monthly=True)
                odata = E.average_data(odata, 2)
                xobs = obslons
            elif (orientation == 'monthly'):
                odata = E.average_data(obsdata, 'monthly')
                xobs = np.arange(0, 12, 1)

            # Bad cloud values fix
            if (scale_cloud):
                odata = odata * scale

            # Plot model values to first graph
            ocolor, odashes, owidth = E.get_model_plot_style(obsmodel)
            if (len(odashes) == 0):
                axs[0].plot(xobs,
                            odata,
                            color=ocolor,
                            linewidth=owidth,
                            label=obsmodel + ' (obs)')
            else:
                line1, = axs[0].plot(xobs,
                                     odata,
                                     '--',
                                     color=ocolor,
                                     linewidth=owidth,
                                     label=obsmodel + ' (obs)')
                line1.set_dashes(odashes)

            multimodelmean_initialized = False
            for model in models:
                # Read in model specific data
                datafile = nc.Dataset(models[model], 'r')
                # A-laue_ax+
                E.add_to_filelist(models[model])
                # A-laue_ax-
                data_units = datafile.variables[datakey].units
                if (orientation == 'monthly'):
                    lats, lons, data = E.get_model_data(modelconfig,
                                                        experiment,
                                                        area,
                                                        datakey,
                                                        datafile)
                else:
                    lats, lons, data = E.get_model_data(modelconfig,
                                                        experiment,
                                                        area,
                                                        datakey,
                                                        datafile,
                                                        extend=orientation)
                datafile.close()

                # Bad cloud values fix
                if (scale_cloud):
                    data = data * scale
                    if (data_units == 'kg m-2'):
                        data_units = r'$\mu$m'
                    else:
                        data_units += scale_str

                # Process depending on orientation
                if   (orientation == 'lats'):
                    data = E.extract_seasonal_mean_values(modelconfig,
                                                          data,
                                                          experiment,
                                                          season,
                                                          monthly=True)
                    data = E.average_data(data, 1)
                    xdata = lats
                    interpolate = interp1d(xdata, data, kind='cubic', bounds_error=False)
                    idata = interpolate(xobs)
                elif (orientation == 'lons'):
                    data = E.extract_seasonal_mean_values(modelconfig,
                                                          data,
                                                          experiment,
                                                          season,
                                                          monthly=True)
                    data = E.average_data(data, 2)
                    xdata = lons
                    interpolate = interp1d(xdata, data, kind='cubic', bounds_error=False)
                    idata = interpolate(xobs)
                elif (orientation == 'monthly'):
                    data = E.average_data(data, 'monthly')
                    xdata = xobs
                    idata = data

                # plotting custom dashes requires some extra effort (else)
                # with empty dashes the format is default
                rmse = round(np.mean((idata - odata) ** 2) ** 0.5, 1)
                mlabel = model + " (RMSE: " + str(rmse) + ")"
                color, dashes, width = E.get_model_plot_style(model)
                if (len(dashes) == 0):
                    axs[0].plot(xdata,
                                data, color=color,
                                linewidth=width,
                                label=mlabel)
                    axs[1].plot(xobs,
                                idata - odata,
                                color=color,
                                linewidth=width,
                                label=mlabel)
                else:
                    line1, = axs[0].plot(xdata,
                                         data,
                                         '--',
                                         color=color,
                                         linewidth=width,
                                         label=mlabel)
                    line1.set_dashes(dashes)
                    line2, = axs[1].plot(xobs,
                                         idata - odata,
                                         color=color,
                                         linewidth=width,
                                         label=mlabel)
                    line2.set_dashes(dashes)

                # Next the multimodel mean. We can just sum it up and divide
                if not (multimodelmean_initialized):
                    multimodelmean_initialized = True
                    mmmean = idata
                else:
                    mmmean += idata

            # Plot the multimodel mean out if required
            if (len(models) > 1):
                mmmean = mmmean / len(models)
                rmse = round(np.mean((mmmean - odata) ** 2) ** 0.5, 1)
                mlabel = "Multimodel mean (RMSE: " + str(rmse) + ")"

                color, dashes, width = E.get_model_plot_style('model_mean')
                if (len(dashes) == 0):
                    axs[0].plot(xobs,
                                mmmean,
                                color=color,
                                linewidth=width,
                                label=mlabel)
                    axs[1].plot(xobs,
                                mmmean - odata,
                                color=color,
                                linewidth=width,
                                label=mlabel)
                else:
                    line1, = axs[0].plot(xobs,
                                         mmmean,
                                         '--',
                                         color=color,
                                         linewidth=width,
                                         label=mlabel)
                    line1.set_dashes(dashes)
                    line2, = axs[1].plot(xobs,
                                         mmmean - odata,
                                         color=color,
                                         linewidth=width,
                                         label=mlabel)
                    line2.set_dashes(dashes)

            # Now to make the plot pretty i.e. headers, ticks, title etc.
            suptitle = E.get_title_basename(datakey)
            specifier = mean_name
            if   (orientation == 'lats'):
                title = "Latitudinal mean"
                specifier += '-latitudinal-mean'
                x_min = lat_min
                x_max = lat_max
                xticks = E.get_ticks(8, np.array([x_min, x_max]))
                lon_info = E.get_ticks_labels(np.array([lon_min, lon_max]),
                                              'lons')
                if (int(lon_max) == 360):
                    suptitle += " (lons [" + lon_info[0] + ":360" + "])"
                else:
                    suptitle += " (lons [" + lon_info[0] + ":" + lon_info[1] + "])"

            elif (orientation == 'lons'):
                title = "Longitudinal mean"
                specifier += '-longitudinal-mean'
                x_min = lon_min
                x_max = lon_max
                xticks = E.get_ticks(8, np.array([x_min, x_max]))
                lat_info = E.get_ticks_labels(np.array([lat_min, lat_max]),
                                              'lats')
                suptitle += " (lats [" + lat_info[0] + ":" + lat_info[1] + "])"

            elif (orientation == 'monthly'):
                title = "Monthly mean"
                specifier += '-yearly-cycle-' + area
                x_min = 0
                x_max = 11
                xticks = xobs
                lat_info = E.get_ticks_labels(np.array([lat_min, lat_max]),
                                              'lats')
                lon_info = E.get_ticks_labels(np.array([lon_min, lon_max]),
                                              'lons')
                suptitle += " (lats [" + lat_info[0] + ":" + lat_info[1] + "], "
                if (lon_max == 360):
                    suptitle += " lons [" + lon_info[0] + ":360" + "])"
                else:
                    suptitle += " lons [" + lon_info[0] + ":" + lon_info[1] + "])"

            if (season == 'monthly'):
                pass
            else:
                suptitle += " for " + season
                specifier = specifier + "-" + season
            plt.suptitle(suptitle, fontsize=20)
            labels = E.get_ticks_labels(xticks, orientation)
            axs[0].grid(plot_grid)
            axs[0].set_xlim(x_min, x_max)
            axs[0].xaxis.set_ticks(xticks)
            axs[0].set_xticklabels(labels)
            axs[0].locator_params(axis='y', nbins=5)
            axs[0].set_title(title + " [" + data_units + "]")
            axs[1].grid(plot_grid)
            axs[1].set_xlim(x_min, x_max)
            axs[1].xaxis.set_ticks(xticks)
            axs[1].set_xticklabels(labels)
            axs[1].locator_params(axis='y', nbins=5)
            axs[1].set_title(title + " (model - obs) [" + data_units + "]")

            # Saving the plot and letting the user know what just happened
            diag_name = E.get_diag_script_name()
            output_file = E.get_plot_output_filename(diag_name=diag_name,
                                                     variable=datakey,
                                                     specifier=specifier)
            output_dir = os.path.join(plot_dir, diag_name)
            E.ensure_directory(output_dir)
            handles, labels = axs[0].get_legend_handles_labels()
            plt.legend(handles,
                       labels,
                       loc='center left',
                       bbox_to_anchor=(0.7, 0.5),
                       bbox_transform=plt.gcf().transFigure)
            plt.savefig(os.path.join(output_dir, output_file))
            info("", verbosity, 1)
            info("Created image: ", verbosity, 1)
            info(output_file, verbosity, 1)
    obsfile.close()


def process_radiation_maps(E, modelconfig, datakey, specifier):
    """
        Main script for gathering radiation seasonal
        map values and plotting them.
    """
    config_file = E.get_configfile()
    experiment = 'SouthernHemisphere'
    areas = modelconfig.get(experiment, 'areas').split()
    seasons = modelconfig.get(experiment, 'seasons').split()
    plot_dir = E.get_plot_dir()
    verbosity = E.get_verbosity()
    work_dir = E.get_work_dir()
    datakeycs = datakey + 'cs'

    # Check that we have similar model sets (obs should be the same)
    obsmodel, obs_loc, models = E.get_clim_model_and_obs_filenames(datakey)
    obsmodel, obscs_loc, modelscs = E.get_clim_model_and_obs_filenames(datakeycs)
    E.check_model_instances(models, modelscs)

    # Loop over areas and seasons - generate season specific maps
    obsfile = nc.Dataset(obs_loc, 'r')
    obscsfile = nc.Dataset(obscs_loc, 'r')
    # A-laue_ax+
    E.add_to_filelist(obs_loc)
    E.add_to_filelist(obscs_loc)
    # A-laue_ax-

    # Basic structure of the figure is same for all plots,
    # we clear it at intervals
    for area in areas:
        # Get area specific configuration and observations
        coords = E.get_area_coordinates(modelconfig, experiment, area)
        xticks = E.get_ticks(5, coords[2:])
        yticks = E.get_ticks(3, coords[:2])
        xlabels = E.get_ticks_labels(xticks, 'lons')
        ylabels = E.get_ticks_labels(yticks, 'lats')
        area_key = experiment + '_' + area
        obslats, obslons, obsdata = E.get_model_data(modelconfig,
                                                     experiment,
                                                     area,
                                                     datakey,
                                                     obsfile)
        obscslats, obscslons, obscsdata = E.get_model_data(modelconfig,
                                                           experiment,
                                                           area,
                                                           datakeycs,
                                                           obscsfile)

        for model in models:
            # Get model and area specific values
            datafile = nc.Dataset(models[model], 'r')
            datafilecs = nc.Dataset(modelscs[model], 'r')
            # A-laue_ax+
            E.add_to_filelist(models[model])
            E.add_to_filelist(modelscs[model])
            # A-laue_ax-
            mlats, mlons, mdata = E.get_model_data(modelconfig,
                                                   experiment,
                                                   area,
                                                   datakey,
                                                   datafile)
            mcslats, mcslons, mcsdata = E.get_model_data(modelconfig,
                                                         experiment,
                                                         area,
                                                         datakeycs,
                                                         datafilecs)
            units = datafile.variables[datakey].units

            for season in seasons:
                # Extract seasonal specific values from area values
                data = E.extract_seasonal_mean_values(modelconfig,
                                                      mdata,
                                                      experiment,
                                                      season)
                datacs = E.extract_seasonal_mean_values(modelconfig,
                                                        mcsdata,
                                                        experiment,
                                                        season)
                obs = E.extract_seasonal_mean_values(modelconfig,
                                                     obsdata,
                                                     experiment,
                                                     season)
                obscs = E.extract_seasonal_mean_values(modelconfig,
                                                       obscsdata,
                                                       experiment,
                                                       season)

                # Interpolate to obs grid
                data = interpolate_data_grid(data,
                                             mlats,
                                             mlons,
                                             obslats,
                                             obslons)
                datacs = interpolate_data_grid(datacs,
                                               mcslats,
                                               mcslons,
                                               obslats,
                                               obslons)
                obscs = interpolate_data_grid(obscs,
                                              obscslats,
                                              obscslons,
                                              obslats,
                                              obslons)

                # Get contour specific plotting limits
                a, b, c, d, e, f = get_contour_config(modelconfig,
                                                      config_file,
                                                      area_key,
                                                      datakey)
                lev, diff, diff2, cm_model, cm_diff, cm_dev = a, b, c, d, e, f

                # Plot the figures
                fig, axs = plt.subplots(3, 2, figsize=(21, 15))
                fig.subplots_adjust(hspace=1.0)

                # Same Basemap (coastlines) for all subplots
                # This is what usually takes time to process
                for ax in axs.flat:
                    map_ax = Basemap(ax=ax,
                                     fix_aspect=False,
                                     llcrnrlat=coords[0], urcrnrlat=coords[1],
                                     llcrnrlon=coords[2], urcrnrlon=coords[3])
                    map_ax.drawcoastlines()

                # Plot the contours and colourmaps
                dash = '  -  '  # Could also use u"\u2014" for long dash

                c00 = axs[0, 0].contourf(obslons,
                                         obslats,
                                         data,
                                         levels=lev,
                                         cmap=cm_model,
                                         extend='both',
                                         latlon=True)
                cax00 = fig.add_axes([0.15, 0.69, 0.3, 0.01])
                cbar00 = plt.colorbar(c00, cax=cax00, orientation='horizontal')
                cbar00.set_label(units)
                axs[0, 0].set_title(model)

                c01 = axs[0, 1].contourf(obslons,
                                         obslats,
                                         datacs,
                                         levels=lev,
                                         cmap=cm_model,
                                         extend='both',
                                         latlon=True)
                cax01 = fig.add_axes([0.575, 0.69, 0.3, 0.01])
                cbar01 = plt.colorbar(c01, cax=cax01, orientation='horizontal')
                cbar01.set_label(units)
                axs[0, 1].set_title(model + ' (clear sky)')

                c10 = axs[1, 0].contourf(obslons,
                                         obslats,
                                         data - obs,
                                         levels=diff,
                                         cmap=cm_diff,
                                         extend='both',
                                         latlon=True)
                cax10 = fig.add_axes([0.15, 0.37, 0.3, 0.01])
                cbar10 = plt.colorbar(c10, cax=cax10, orientation='horizontal')
                cbar10.set_label(units)
                axs[1, 0].set_title(model + dash + obsmodel)

                c11 = axs[1, 1].contourf(obslons,
                                         obslats,
                                         datacs - obscs,
                                         levels=diff,
                                         cmap=cm_diff,
                                         extend='both',
                                         latlon=True)
                cax11 = fig.add_axes([0.575, 0.37, 0.3, 0.01])
                cbar11 = plt.colorbar(c11, cax=cax11, orientation='horizontal')
                cbar11.set_label(units)
                axs[1, 1].set_title(model + ' (cs)' + dash + obsmodel + ' (cs)')

                c20 = axs[2, 0].contourf(obslons,
                                         obslats,
                                         data - datacs,
                                         levels=diff2,
                                         cmap=cm_dev,
                                         extend='both',
                                         latlon=True)
                cax20 = fig.add_axes([0.15, 0.05, 0.3, 0.01])
                cbar20 = plt.colorbar(c20, cax=cax20, orientation='horizontal')
                cbar20.set_label(units)
                axs[2, 0].set_title(model + dash + model + ' (cs)')

                c21 = axs[2, 1].contourf(obslons,
                                         obslats,
                                         obs - obscs,
                                         levels=diff2,
                                         cmap=cm_dev,
                                         extend='both',
                                         latlon=True)
                cax21 = fig.add_axes([0.575, 0.05, 0.3, 0.01])
                cbar21 = plt.colorbar(c21, cax=cax21, orientation='horizontal')
                cbar21.set_label(units)
                axs[2, 1].set_title(obsmodel + dash + obsmodel + ' (cs)')

                # Labels
                for row in xrange(3):
                    for col in xrange(2):
                        axs[row, col].xaxis.set_ticks(xticks)
                        axs[row, col].set_xticklabels(xlabels)
                        axs[row, col].yaxis.set_ticks(yticks)
                        axs[row, col].set_yticklabels(ylabels)

                # Title and filename based on variables
                suptitle = E.get_title_basename(datakey) + " for months: " + season
                plt.suptitle(suptitle, fontsize=20)
                variable = datakey
                specifier = specifier + "-" + season
                diag_name = E.get_diag_script_name()
                output_file = E.get_plot_output_filename(variable=variable,
                                                         specifier=specifier,
                                                         model=model)
                output_dir = os.path.join(plot_dir, diag_name)
                E.ensure_directory(output_dir)
                plt.savefig(os.path.join(output_dir, output_file))
                plt.clf()
                info("", verbosity, 1)
                info("Created image: ", verbosity, 1)
                info(output_file, verbosity, 1)
            datafile.close()
            datafilecs.close()
    obsfile.close()
    obscsfile.close()


def separate_list(inlist, rule):
    """
    Separates a list in two based on the rule.
    The rule is a string that is mached for last characters of list elements.
    """
    lrule = len(rule)
    list1 = []
    list2 = []
    for item in inlist:
        if (item[-lrule:] == rule):
            list1.append(item)
        else:
            list2.append(item)
    return list1, list2


def process_simple_maps(E, modelconfig, datakey, specifier):
    """ Main script for gathering seasonal map values and plotting them. """
    config_file = E.get_configfile()
    experiment = 'SouthernHemisphere'
    areas = modelconfig.get(experiment, 'areas').split()
    seasons = modelconfig.get(experiment, 'seasons').split()
    plot_dir = E.get_plot_dir()
    verbosity = E.get_verbosity()
    work_dir = E.get_work_dir()

    # Get the observations
    obsmodel, obs_loc, models = E.get_clim_model_and_obs_filenames(datakey)
    obsfile = nc.Dataset(obs_loc, 'r')
    # A-laue_ax+
    E.add_to_filelist(obs_loc)
    # A-laue_ax-

    scale_cloud = False
    # Ugly fix for scaling cloud ice/liquid water path values
    if (datakey == 'clivi' or datakey == 'clwvi'):
        scale_cloud = True
        scale = 1E3
        scale_str = ' * 1E3'

    # Basic structure of the figure is same for all plots,
    # we clear it at intervals
    for area in areas:
        # Get area specific configuration and observations
        coords = E.get_area_coordinates(modelconfig, experiment, area)
        xticks = E.get_ticks(5, coords[2:])
        yticks = E.get_ticks(3, coords[:2])
        xlabels = E.get_ticks_labels(xticks, 'lons')
        ylabels = E.get_ticks_labels(yticks, 'lats')
        area_key = experiment + '_' + area
        obslats, obslons, obsdata = E.get_model_data(modelconfig,
                                                     experiment,
                                                     area,
                                                     datakey,
                                                     obsfile)
        if (scale_cloud):
            obsdata = obsdata * scale

        for model in models:
            # Get model and area specific values
            datafile = nc.Dataset(models[model], 'r')
            # A-laue_ax+
            E.add_to_filelist(models[model])
            # A-laue_ax-
            mlats, mlons, mdata = E.get_model_data(modelconfig,
                                                   experiment,
                                                   area,
                                                   datakey,
                                                   datafile)
            units = datafile.variables[datakey].units
            if (scale_cloud):
                mdata = mdata * scale
                if (units == 'kg m-2'):
                    units = r'$\mu$m'
                else:
                    units += scale_str

            for season in seasons:
                # Extract seasonal specific values from area values
                data = E.extract_seasonal_mean_values(modelconfig,
                                                      mdata,
                                                      experiment,
                                                      season)
                obs = E.extract_seasonal_mean_values(modelconfig,
                                                     obsdata,
                                                     experiment,
                                                     season)
                # Interpolate to obs grid
                data = interpolate_data_grid(data,
                                             mlats,
                                             mlons,
                                             obslats,
                                             obslons)
                # Get contour specific plotting limits
                lev, diff, cm_model, cm_diff = get_contour_config(modelconfig,
                                                                  config_file,
                                                                  area_key,
                                                                  datakey)

                # Plot the figures
                fig, axs = plt.subplots(2, 1, figsize=(15, 10))
                fig.subplots_adjust(top=0.9)
                fig.subplots_adjust(hspace=0.5)

                # Same Basemap (coastlines) for all subplots
                # This is what usually takes time to process
                for ax in axs.flat:
                    map_ax = Basemap(ax=ax,
                                     fix_aspect=False,
                                     llcrnrlat=coords[0], urcrnrlat=coords[1],
                                     llcrnrlon=coords[2], urcrnrlon=coords[3])
                    map_ax.drawcoastlines()

                # Plot the contours and colourmaps
                dash = '  -  '  # Could also use u"\u2014" for long dash

                c0 = axs[0].contourf(obslons, obslats, data,
                                     levels=lev, cmap=cm_model,
                                     extend='both', latlon=True)
                cax0 = fig.add_axes([0.36, 0.525, 0.3, 0.01])
                cbar0 = plt.colorbar(c0, cax=cax0, orientation='horizontal')
                cbar0.set_label(units)
                axs[0].set_title(model)

                c1 = axs[1].contourf(obslons, obslats, data - obs,
                                     levels=diff, cmap=cm_diff,
                                     extend='both', latlon=True)
                cax1 = fig.add_axes([0.36, 0.05, 0.3, 0.01])
                cbar1 = plt.colorbar(c1, cax=cax1, orientation='horizontal')
                cbar1.set_label(units)
                axs[1].set_title(model + dash + obsmodel)

                # Labels
                for row in xrange(2):
                    axs[row].xaxis.set_ticks(xticks)
                    axs[row].set_xticklabels(xlabels)
                    axs[row].yaxis.set_ticks(yticks)
                    axs[row].set_yticklabels(ylabels)

                # Title and filename based on variables
                suptitle = E.get_title_basename(datakey) + " for months: " + season
                plt.suptitle(suptitle, fontsize=20)
                variable = datakey
                specifier = specifier + "-" + season
                diag_name = E.get_diag_script_name()
                output_file = E.get_plot_output_filename(diag_name=diag_name,
                                                         variable=variable,
                                                         specifier=specifier,
                                                         model=model)
                # Save the image and let the user know what happened
                output_dir = os.path.join(plot_dir, diag_name)
                E.ensure_directory(output_dir)
                plt.savefig(os.path.join(output_dir, output_file))
                plt.clf()
                info("", verbosity, 1)
                info("Created image: ", verbosity, 1)
                info(output_file, verbosity, 1)
            datafile.close()
    obsfile.close()
