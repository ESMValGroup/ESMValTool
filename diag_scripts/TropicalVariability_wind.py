"""
;;#############################################################################
;; TropicalVariablity_wind.py
;; Author: Jarmo Makela (FMI, Finland)
;; EMBRACE project
;;#############################################################################
;; Description
;;    This script is the preprocessing script for Tropical Variability wind
;;    divergence plots (part of equatorial plots)
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    plot_equatorial:      Switch for equatorial plots
;;    plot_scatter:         Switch for scatter plots
;;    plot_zonal_means:     Switch for zonal means plots
;;    mask_unwanted_values: Mask values outside given range
;;    mask_limit_low:       lower mask limit
;;    mask_limit_high:      uppper mask limit
;;    plot_grid:            provides a background grid for relavant plots
;;
;;    [equatorial]
;;    areas:                One of "Atlantic" "Indian" "Pacific"
;;
;;    [equatorial_Atlantic]
;;    lat_min:  Spatial extent
;;    lat_max:  Spatial extent
;;    lon_min:  Spatial extent
;;    lon_max:  Spatial extent
;;    prec_min: Range of values
;;    prec_max: Range of values
;;    temp_min: Range of values
;;    temp_max: Range of values
;;    wind_min: Range of values
;;    wind_max: Range of values
;;    div_min:  Range of values
;;    div_max:  Range of values
;;
;;    [equatorial_Indian]
;;    lat_min:  Spatial extent
;;    lat_max:  Spatial extent
;;    lon_min:  Spatial extent
;;    lon_max:  Spatial extent
;;    prec_min: Range of values
;;    prec_max: Range of values
;;    temp_min: Range of values
;;    temp_max: Range of values
;;    wind_min: Range of values
;;    wind_max: Range of values
;;    div_min:  Range of values
;;    div_max:  Range of values
;;
;;    [equatorial_Pacific]
;;    lat_min:  Spatial extent
;;    lat_max:  Spatial extent
;;    lon_min:  Spatial extent
;;    lon_max:  Spatial extent
;;    prec_min: Range of values
;;    prec_max: Range of values
;;    temp_min: Range of values
;;    temp_max: Range of values
;;    wind_min: Range of values
;;    wind_max: Range of values
;;    div_min:  Range of values
;;    div_max:  Range of values
;;
;;    [scatter]
;;    areas:           One of "West-Pacific" "Central-Pacific" "East-Pacific"
;;    seasons:         One of "annual" "DJF" "MAM" "JJA" "SON"
;;    seasonal_limits: True/False if you want to use your own limits
;;
;;    [scatter_West-Pacific]
;;    lat_min:              Spatial extent for West-Pacific
;;    lat_max:              Spatial extent for West-Pacific
;;    lon_min:              Spatial extent for West-Pacific
;;    lon_max:              Spatial extent for West-Pacific
;;    season_limits_annual: Spatial extent for annual for West-Pacific
;;    season_limits_DJF:    Spatial extent for DJF for West-Pacific
;;    season_limits_MAM:    Spatial extent for MAM for West-Pacific
;;    season_limits_JJA:    Spatial extent for JJA for West-Pacific
;;    season_limits_SON:    Spatial extent for SON for West-Pacific
;;
;;    [scatter_Central-Pacific]
;;    lat_min               Spatial extent for Central-Pacific
;;    lat_max               Spatial extent for Central-Pacific
;;    lon_min               Spatial extent for Central-Pacific
;;    lon_max               Spatial extent for Central-Pacific
;;    season_limits_annual  Spatial extent for annual for Central-Pacific
;;    season_limits_DJF     Spatial extent for DJF for Central-Pacific
;;    season_limits_MAM     Spatial extent for MAM for Central-Pacific
;;    season_limits_JJA     Spatial extent for JJA for Central-Pacific
;;    season_limits_SON     Spatial extent for SON for Central-Pacific
;;
;;    [scatter_East-Pacific]
;;    lat_min                       Spatial extent for East-Pacific
;;    lat_max                       Spatial extent for East-Pacific
;;    lon_min                       Spatial extent for East-Pacific
;;    lon_max                       Spatial extent for East-Pacific
;;    season_limits_annual          Spatial extent for annual for East-Pacific
;;    season_limits_DJF             Spatial extent for DJF for East-Pacific
;;    season_limits_MAM             Spatial extent for MAM for East-Pacific
;;    season_limits_JJA             Spatial extent for JJA for East-Pacific
;;    season_limits_SON             Spatial extent for SON for East-Pacific
;;
;;    [scatter_season_DJF]
;;    season_months: DJF months
;;
;;    [scatter_season_MAM]
;;    season_months: MAM months
;;
;;    [scatter_season_JJA]
;;    season_months: JJA months
;;
;;    [scatter_season_SON]
;;    season_months: SON months
;;
;;    [zonal_means]
;;    areas: One of "Pacific" "Atlantic" "Indian"
;;
;;    [zonal_means_Atlantic]
;;    lat_min: Spatial extent for zonal mean in the Atlantic Ocean
;;    lat_max: Spatial extent for zonal mean in the Atlantic Ocean
;;    lon_min: Spatial extent for zonal mean in the Atlantic Ocean
;;    lon_max: Spatial extent for zonal mean in the Atlantic Ocean
;;
;;    [zonal_means_Indian]
;;    lat_min: Spatial extent for zonal mean in the Indian Ocean
;;    lat_max: Spatial extent for zonal mean in the Indian Ocean
;;    lon_min: Spatial extent for zonal mean in the Indian Ocean
;;    lon_max: Spatial extent for zonal mean in the Indian Ocean
;;
;;    [zonal_means_Pacific]
;;    lat_min: Spatial extent for zonal mean in the Pacific Ocean
;;    lat_max: Spatial extent for zonal mean in the Pacific Ocean
;;    lon_min: Spatial extent for zonal mean in the Pacific Ocean
;;    lon_max: Spatial extent for zonal mean in the Pacific Ocean
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
;;    Dependency - must be run after "diag_script/TropicalVariability.py"
;;
;; Modification history
;;    20151029-A_laue_ax: added output of acknowledgements + processed files
;;                        to log-file
;;    20150515-A_maek_ja: written
;;
;; #############################################################################
"""

# Basic Python packages
import ConfigParser
import os
import pdb
import sys
import numpy as np

# Force matplotlib to not use any Xwindows backend.
import matplotlib
matplotlib.use('Agg')

# Common Python packages
import matplotlib.pyplot as plt

from scipy.interpolate import griddata

# ESMValTool defined Python packages
sys.path.append("./interface_scripts")
from esmval_lib import ESMValProject
from auxiliary import info

from mpl_toolkits.basemap import Basemap

# Python extensions and netcdf handler
import netCDF4 as nc

def main(project_info):
    """Preprocessing script for Tropical Variability equatorial divergence
    winds. """

    # ESMValProject provides some easy methods to access information in
    # project_info but you can also access it directly (as in main.py)
    # First we get some general configurations (ESMValProject)
    # and then we read in the model configuration file

    E = ESMValProject(project_info)
    config_file = E.get_configfile()
    plot_dir = E.get_plot_dir()
    verbosity = E.get_verbosity()

    # A-laue_ax+
    diag_script = E.get_diag_script_name()
    res = E.write_references(diag_script,              # diag script name
                             ["A_maek_ja"],            # authors
                             ["A_eval_ma", "A_jone_co"], # contributors
                             ["D_li14jclim"],          # diag_references
                             [""],                     # obs_references
                             ["P_embrace"],            # proj_references
                             project_info,
                             verbosity,
                             False)
    # A-laue_ax-

    modelconfig = ConfigParser.ConfigParser()
    modelconfig.read(config_file)
    E.ensure_directory(plot_dir)

    # Here we check and process only desired parts of the diagnostics
    if (modelconfig.getboolean('general', 'plot_equatorial')):
        info("Starting to gather values for equatorial divergence plots",
             verbosity, 2)
        process_divergence(E, modelconfig)


# All callable functions below are in alphabetical order
# E. style functions are in the general python file esmval_lib.py

def calculate_derivatives(E, lats, lons, data, key):
    """Calculates pointwise derivatives."""
    # First we average over time
    means = np.zeros((data.shape[1], data.shape[2]))
    for lat in xrange(data.shape[1]):
        for lon in xrange(data.shape[2]):
            means[lat, lon] = np.mean(data[:, lat, lon])
    # For ua the derivate is longitudinal
    if (key == 'ua'):
        output = np.zeros((data.shape[1], data.shape[2] - 2))
        for lat in xrange(len(lats)):
            for lon in xrange(len(lons) - 2):
                # Delta is lat/lon specific (we don't assume regular grid)
                delta = convert_latlon_meters(lats[lat], lats[lat],
                                              lons[lon], lons[lon + 2])
                output[lat, lon] = (means[lat, lon + 2] - means[lat, lon]) / delta
        # return mean values and the derivatives
        return means[:, 1:-1], output
    # For va the derivate is latitudinal
    elif (key == 'va'):
        output = np.zeros((data.shape[1] - 2, data.shape[2]))
        for lat in xrange(len(lats) - 2):
            for lon in xrange(len(lons)):
                # Delta is lat/lon specific (we don't assume regular grid)
                delta = convert_latlon_meters(lats[lat], lats[lat + 2],
                                              lons[lon], lons[lon])
                output[lat, lon] = (means[lat + 2, lon] - means[lat, lon]) / delta
        # return mean values and the derivatives
        return means[1:-1, :], output


def convert_latlon_meters(lat1, lat2, lon1, lon2):
    """Converts the distance between two coordinates to meters"""
    return np.arccos(np.sin(lat1 * np.pi / 180) * np.sin(lat2 * np.pi / 180)
                     + np.cos(lat1 * np.pi / 180) * np.cos(lat2 * np.pi / 180)
                     * np.cos(lon2 * np.pi / 180 - lon1 * np.pi / 180)) * 6371000


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
        gridvalues = griddata(datapoints, datavalues, (grid_lats, grid_lons),
                              method='cubic')
        nearest = griddata(datapoints, datavalues, (grid_lats, grid_lons),
                           method='nearest')

        # In some cases griddata cannot export values on the border of the grid
        for i in xrange(gridvalues.shape[0]):
            for j in xrange(gridvalues.shape[1]):
                if (np.isnan(gridvalues[i, j])):
                    gridvalues[i, j] = nearest[i, j]
        return gridvalues.T


def process_divergence(E, modelconfig):
    """This script preprocessess equatorial divergence values.
    The data is extracted and saved to a npy-file. Plotting is done in
    TropicalVariability_EQ.py. """
    experiment = 'equatorial'
    datakeys = E.get_currVars()
    verbosity = E.get_verbosity()
    plot_dir = E.get_plot_dir()
    work_dir = E.get_work_dir()
    areas = modelconfig.get(experiment, 'areas').split()

    if ('ua' in datakeys[0]):
        ua_key = datakeys[0]
        va_key = datakeys[1]
    elif ('ua' in datakeys[1]):
        ua_key = datakeys[1]
        va_key = datakeys[0]
    else:
        print("PY  ERROR: Unexpected variables for divergence.")
        print("PY  ERROR: I was expecting only 'ua' and 'va' (or similar).")
        print("PY  ERROR: What I got was: " + datakeys[0] + ", " + datakeys[1])
        print ("PY ERROR: Stopping the script and exiting")
        sys.exit()

    # Starting to extract required areas
    for area in areas:
        area_key = experiment + '_' + area
        lat_min = modelconfig.getint(area_key, 'lat_min')
        lat_max = modelconfig.getint(area_key, 'lat_max')
        lon_min = modelconfig.getint(area_key, 'lon_min')
        lon_max = modelconfig.getint(area_key, 'lon_max')

        # create a python array file for each area, model  and datakey
        # overwrite if exists
        base_name = work_dir\
                    + 'TropicalVariability_'\
                    + area_key\
                    + '_divergence_model_'
        ua_models = E.get_clim_model_filenames(variable=ua_key)
        va_models = E.get_clim_model_filenames(variable=va_key)
        E.check_model_instances(ua_models, va_models)

        # Get obs model name and extract these first
        for model in ua_models:
            if (E.get_model_id(model) == 'obs'):
                obs = model
                out_file = base_name + model + '.npy'
                ua_datafile = nc.Dataset(ua_models[model], 'r')
                va_datafile = nc.Dataset(va_models[model], 'r')
                # A-laue_ax+
                E.add_to_filelist(ua_models[model])
                E.add_to_filelist(va_models[model])
                # A-laue_ax-
                # Add ghost layers for ua longtitudes
                # va doesn't need them. Lats and lons are same for both
                ulats, ulons, ua_data = E.get_model_data(modelconfig, experiment,
                                                         area, ua_key, ua_datafile,
                                                         extend='lons')
                vlats, vlons, va_data = E.get_model_data(modelconfig, experiment,
                                                         area, va_key, va_datafile,
                                                         extend='lats')
                ulons = E.ensure_looping(ulons)
                vlons = E.ensure_looping(vlons)
                ua_datafile.close()
                va_datafile.close()
                ua_obs, ua_derv = calculate_derivatives(E, ulats, ulons, ua_data, 'ua')
                va_obs, va_derv = calculate_derivatives(E, vlats, vlons, va_data, 'va')
                diverg = np.zeros(len(vlons))
                for lon in xrange(diverg.shape[0]):
                    diverg[lon] = np.mean(ua_derv[:, lon] + va_derv[:, lon])
                # Now we have everything so only output needed
                outdata = vlons, diverg
                np.save(out_file, outdata)
                obscolor, obsdashes, obswidth = E.get_model_plot_style(model)
                obslats, obslons = ulats, vlons
        del ua_models[obs]

        # Now we do the same for all models and plot the wind vector field
        for model in ua_models:
            out_file = base_name + model + '.npy'
            ua_datafile = nc.Dataset(ua_models[model], 'r')
            va_datafile = nc.Dataset(va_models[model], 'r')
            # A-laue_ax+
            E.add_to_filelist(ua_models[model])
            E.add_to_filelist(va_models[model])
            # A-laue_ax-
            # Add ghost layers for ua longtitudes
            # va doesn't need them. Lats and lons are same for both
            ulats, ulons, ua_data = E.get_model_data(modelconfig, experiment,
                                                     area, ua_key, ua_datafile,
                                                     extend='lons')
            vlats, vlons, va_data = E.get_model_data(modelconfig, experiment,
                                                     area, va_key, va_datafile,
                                                     extend='lats')
            ulons = E.ensure_looping(ulons)
            vlons = E.ensure_looping(vlons)
            ua_datafile.close()
            va_datafile.close()
            ua_mean, ua_derv = calculate_derivatives(E, ulats, ulons, ua_data, 'ua')
            va_mean, va_derv = calculate_derivatives(E, vlats, vlons, va_data, 'va')
            diverg = np.zeros(len(vlons))
            for lon in xrange(diverg.shape[0]):
                diverg[lon] = np.mean(ua_derv[:, lon] + va_derv[:, lon])
            # Now we have everything so only output needed
            outdata = vlons, diverg
            np.save(out_file, outdata)

            # Interpolate model data to obs grid
            ua_mean = interpolate_data_grid(ua_mean, ulats, vlons, obslats, obslons)
            va_mean = interpolate_data_grid(va_mean, ulats, vlons, obslats, obslons)
            color, dashes, width = E.get_model_plot_style(model)

            # Ensure proper handling for lons
            if (lon_max < lon_min):
                lon_max += 360
            # Start the plots
            plt.clf()
            fig, axs = plt.subplots(2, 1, figsize=(5 + (lon_max - lon_min) / 10, 10))
            fig.subplots_adjust(hspace=0.3)

            for ax in axs.flat:
                map_ax = Basemap(ax=ax, fix_aspect=False,
                                 llcrnrlat=obslats[0], urcrnrlat=obslats[-1],
                                 llcrnrlon=obslons[0], urcrnrlon=obslons[-1])
                map_ax.drawcoastlines()
            Q0 = axs[0].quiver(obslons, obslats, va_mean, ua_mean, color=color,
                               pivot='middle', units='xy', scale=5.0, width=0.1,
                               headwidth=4, headlength=3, headaxislength=4)
            Q1 = axs[1].quiver(obslons, obslats, va_obs, ua_obs, color=obscolor,
                               pivot='middle', units='xy', scale=5.0, width=0.1,
                               headwidth=4, headlength=3, headaxislength=4)
            qk0 = plt.quiverkey(Q0, 1.05, 0.5, 5, r'$5 \frac{m}{s}$')
            qk1 = plt.quiverkey(Q1, 1.05, 0.5, 5, r'$5 \frac{m}{s}$')

            # Ticks, labels, header etc.
            xticks = E.get_ticks(10, np.array([lon_min, lon_max]))
            yticks = E.get_ticks(5, np.array([lat_min, lat_max]))
            xlabels = E.get_ticks_labels(xticks, 'lons')
            ylabels = E.get_ticks_labels(yticks, 'lats')
            axs[0].set_title(model)
            axs[0].set_xlim(lon_min, lon_max)
            axs[0].set_ylim(lat_min, lat_max)
            axs[0].xaxis.set_ticks(xticks)
            axs[0].yaxis.set_ticks(yticks)
            axs[0].set_xticklabels(xlabels)
            axs[0].set_yticklabels(ylabels)
            axs[1].set_title(obs)
            axs[1].set_xlim(lon_min, lon_max)
            axs[1].set_ylim(lat_min, lat_max)
            axs[1].xaxis.set_ticks(xticks)
            axs[1].yaxis.set_ticks(yticks)
            axs[1].set_xticklabels(xlabels)
            axs[1].set_yticklabels(ylabels)

            # Get output filename and save the figure
            diag_name = E.get_diag_script_name()
            output_dir = os.path.join(plot_dir, diag_name)
            E.ensure_directory(output_dir)
            plt.suptitle(area + " ocean mean wind direction and strength", fontsize=16)
            variable_string = ''.join([ua_key, va_key])
            output_file = E.get_plot_output_filename(diag_name=diag_name,
                                                     model=model,
                                                     variable=variable_string,
                                                     specifier=area)
            plt.savefig(os.path.join(output_dir, output_file))

            info("", verbosity, 1)
            info("Created image: ", verbosity, 1)
            info(output_file, verbosity, 1)
