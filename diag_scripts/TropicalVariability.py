"""
;;#############################################################################
;; TropicalVariablity_wind.py
;; Author: Jarmo Makela (FMI, Finland)
;; EMBRACE project
;;#############################################################################
;; Description
;;    This script is the diagnostics and plotting script for
;;    Tropical Variability zonal means and scatterplots
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
;; Modification history
;;    20151029-A_laue_ax: added output of acknowledgements + processed files
;;                        to log-file
;;    20150115-A_maek_ja: written
;;
;; #############################################################################
"""
# Basic Python packages
import ConfigParser
import sys
import numpy as np
import os
import pdb

# ESMValTool defined Python packages
sys.path.append("./interface_scripts")
from esmval_lib import ESMValProject
from auxiliary import info

# Force matplotlib to not use any Xwindows backend.
import matplotlib
matplotlib.use('Agg')

# Common Python packages
import matplotlib.pyplot as plt

# Python extensions and netcdf handler
import netCDF4 as nc

def main(project_info):
    """Diagnostics and plotting script for Tropical Variability.
    We use ts as a proxy for Sea Surface Temperature. """

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
    if (modelconfig.getboolean('general', 'plot_zonal_means')):
        info("Starting to plot zonal means", verbosity, 2)
        process_zonal_means(E, modelconfig)

    if (modelconfig.getboolean('general', 'plot_scatter')):
        info("Starting scatterplot of temperature and precipitation", verbosity, 2)
        process_scatterplot(E, modelconfig)

    if (modelconfig.getboolean('general', 'plot_equatorial')):
        info("Starting to gather values for equatorial means", verbosity, 2)
        process_equatorial_means(E, modelconfig)


# All callable functions below are in alphabetical order
# E. style functions are in the general python file esmval_lib.py

def get_scatterplot_limits(modelconfig, config_file, area, season):
    """Returns area specific seasonal limits for scatterplotting. """
    area_key = 'scatter_' + area
    season_key = 'season_limits_' + season
    limits = []

    # Check the existance of needed limits - otherwise we let python decide
    if modelconfig.has_option(area_key, season_key):
        # Check the limits and thei number
        str_limits = modelconfig.get(area_key, season_key).split()
        if (len(str_limits) == 4):
            limits = np.zeros(4)
            for limit in xrange(4):
                limits[limit] = int(str_limits[limit])
        # Wrong number of limits present
        else:
            print("PY  ERROR: Misdefined seasonal limits for : '" + season_key + "'")
            print("PY  ERROR: You should use four (4) integers for the limits")
            print("PY  ERROR: Check your configuration file: " + config_file)
            print("PY  ERROR: Reverting to automated limits")
            limits = []
    # The entire section is missing
    else:
        print("PY  ERROR: Misdefined seasonal limits in : '[" + area_key + "]'")
        print("PY  ERROR: Couldn't find key: '" + season_key + "'")
        print("PY  ERROR: Check your configuration file: " + config_file)
        print("PY  ERROR: Reverting to automated limits")
        limits = []
    return limits


def get_scatterplot_values(modelconfig, season, data):
    """Extracts the values to plot from all model specific values. """
    season_key = 'scatter_season_' + season

    # We transform monthly values to yearly means or just take the monthly ones
    if (season == 'annual'):
        values = np.zeros(len(data) / 12)
        for year in xrange(len(values)):
            values[year] = np.mean(data[12 * year: 12 * (year + 1), :, :])
    else:
        season_months = modelconfig.get(season_key, 'season_months').split()
        months = len(season_months)
        values = np.zeros(months * len(data) / 12)
        for year in xrange(len(data) / 12):
            cm = 0
            for month in season_months:
                month_loc = int(month) - 1
                values[year * months + cm]\
                    = np.mean(data[year * 12 + month_loc, :, :])
                cm += 1
    return values


def process_equatorial_means(E, modelconfig):
    """This script preprocessess equatorial means for precipitation and
    temperature - the data is extracted and saved to a temporary npy-file.
    Plotting and the wind processing are in TropicalVariability_EQ.py. """
    experiment = 'equatorial'
    datakeys = E.get_currVars()
    work_dir = E.get_work_dir()
    areas = modelconfig.get(experiment, 'areas').split()

    # This should be redundant (no more than one equatorial area) but we'll do
    # the looping similarly to other parts of the code
    for area in areas:
        area_key = experiment + '_' + area
        for datakey in datakeys:
            # create a python array file for each area, model  and datakey
            # overwrite if exists
            base_name = work_dir\
                        + 'TropicalVariability_'\
                        + area_key + '_'\
                        + datakey\
                        + '_model_'
            models = E.get_clim_model_filenames(variable=datakey)
            for model in models:
                out_file = base_name + model + '.npy'
                datafile = nc.Dataset(models[model], 'r')
                # A-laue_ax+
                E.add_to_filelist(models[model])
                # A-laue_ax-
                lats, lons, model_data = E.get_model_data(modelconfig,
                                                          experiment,
                                                          area,
                                                          datakey,
                                                          datafile)
                datafile.close()
                model_data = E.average_data(model_data, 2)
                # Now we have everything so only output neede
                outdata = lons, model_data
                np.save(out_file, outdata)


def process_scatterplot(E, modelconfig):
    """Main script for scatterplots. Outputs model specific scatterplots with
    specified observations. Also prints out comparable statistics. """
    experiment = 'scatter'
    config_file = E.get_configfile()
    datakeys = E.get_currVars()
    plot_dir = E.get_plot_dir()
    verbosity = E.get_verbosity()
    plot_grid = modelconfig.getboolean('general', 'plot_grid')
    areas = modelconfig.get(experiment, 'areas').split()
    seasons = modelconfig.get(experiment, 'seasons').split()

    # Scatterplots are only for temperature/precipitation
    if 'ua' in datakeys:
        datakeys = datakeys.remove('ua')

    # We extract the explicit precipitation datakey since precip can be in
    # mmday or kg/m2s - this affects the datakey (pr or pr-mmday)
    # For consistency we do this also for ts
    for datakey in datakeys:
        if   (datakey == 'ts'):
            ts_key = datakey
        elif (datakey == 'pr' or datakey == 'pr-mmday'):
            pr_key = datakey

    # Get the required filenames and check consistency of model filenames
    # The model keys should be same - the key values (paths) are different
    pr_model, pr_obs_loc, pr_models = E.get_clim_model_and_obs_filenames(pr_key)
    ts_model, ts_obs_loc, ts_models = E.get_clim_model_and_obs_filenames(ts_key)
    E.check_model_instances(pr_models, ts_models)

    # Initialize observations
    pr_obs_file = nc.Dataset(pr_obs_loc, 'r')
    ts_obs_file = nc.Dataset(ts_obs_loc, 'r')
    # A-laue_ax+
    E.add_to_filelist(pr_obs_loc)
    E.add_to_filelist(ts_obs_loc)
    # A-laue_ax-

    # We are outputting model specific scatterplots with areas and seasons
    # So some looping to do. Remember that keys in pr_models and ts_models
    # are the same (we checked this) - the paths differ.
    for model in pr_models:
        # Initializing model specific files
        pr_file = nc.Dataset(pr_models[model], 'r')
        ts_file = nc.Dataset(ts_models[model], 'r')

        # A-laue_ax+
        E.add_to_filelist(pr_models[model])
        E.add_to_filelist(ts_models[model])
        # A-laue_ax-

        pr_units = pr_file.variables[pr_key].units
        ts_units = ts_file.variables[ts_key].units

        # Extract model specific plotting style (we only need the colour)
        modelcolor, dashes, width = E.get_model_plot_style(model)

        plt.clf()
        fig, axs = plt.subplots(len(areas), len(seasons),
                                figsize=(3 * len(seasons), 1 + 3 * (len(areas))))
        fig.subplots_adjust(top=0.83)
        fig.subplots_adjust(right=0.85)
        fig.subplots_adjust(hspace=0.4)
        fig.subplots_adjust(wspace=0.3)

        # Next we extract the required values for each area
        for area in areas:
            pr_data = E.get_model_data(modelconfig,
                                       experiment,
                                       area,
                                       pr_key,
                                       pr_file)

            ts_data = E.get_model_data(modelconfig,
                                       experiment,
                                       area,
                                       ts_key,
                                       ts_file)

            pr_obs = E.get_model_data(modelconfig,
                                      experiment,
                                      area,
                                      pr_key,
                                      pr_obs_file)

            ts_obs = E.get_model_data(modelconfig,
                                      experiment,
                                      area,
                                      ts_key,
                                      ts_obs_file)

            row = areas.index(area)
            area_key = experiment + '_' + area
            lat_min = modelconfig.getint(area_key, 'lat_min')
            lat_max = modelconfig.getint(area_key, 'lat_max')
            lon_min = modelconfig.getint(area_key, 'lon_min')
            lon_max = modelconfig.getint(area_key, 'lon_max')
            lat_area = E.get_ticks_labels(np.array([lat_min, lat_max]), 'lats')
            lon_area = E.get_ticks_labels(np.array([lon_min, lon_max]), 'lons')
            c_area = "[" + lat_area[0] + ":" + lat_area[1] + ", "\
                     + lon_area[0] + ":" + lon_area[1] + "]"

            # And one more loop for different seasons
            for season in seasons:
                col = seasons.index(season)
                pr_season = get_scatterplot_values(modelconfig, season, pr_data)
                ts_season = get_scatterplot_values(modelconfig, season, ts_data)
                pr_obs_se = get_scatterplot_values(modelconfig, season, pr_obs)
                ts_obs_se = get_scatterplot_values(modelconfig, season, ts_obs)

                # plot
                axs[row, col].scatter(ts_season, pr_season,
                                      facecolors=modelcolor,
                                      edgecolors=modelcolor,
                                      label='modelled')

                axs[row, col].scatter(ts_obs_se, pr_obs_se,
                                      facecolors='none',
                                      edgecolors='black',
                                      label='observed')

                # statitics: errors in precip/temp and both
                mse_pr = np.mean((pr_season - pr_obs_se) ** 2)
                mse_ts = np.mean((ts_season - ts_obs_se) ** 2)
                mse = (mse_pr + mse_ts) ** 0.5
                rmse = np.round(mse, 2)
                if (row == 0):
                    axs[row, col].set_title(season + "\n rmse: " + str(rmse))
                else:
                    axs[row, col].set_title("rmse: " + str(rmse))

                # Get specific plotting limits
                if modelconfig.getboolean(experiment, 'seasonal_limits'):
                    # If we use user defined limits
                    limits = get_scatterplot_limits(modelconfig, config_file,
                                                    area, season)
                    if (len(limits) == 4):
                        axs[row, col].set_xlim(limits[0], limits[1])
                        axs[row, col].set_ylim(limits[2], limits[3])
                        xticks = E.get_ticks(5, limits[:2])
                        yticks = E.get_ticks(5, limits[2:])
                    else:
                        xticks = E.get_ticks(5, ts_season, ts_obs_se)
                        yticks = E.get_ticks(5, pr_season, pr_obs_se)
                    # Setting the ticks (where to place the numbers
                    axs[row, col].xaxis.set_ticks(xticks)
                    axs[row, col].yaxis.set_ticks(yticks)
                else:
                    # Let the Gods decide
                    ts_min = int(np.min(ts_obs_se - 0.7))
                    ts_max = int(np.max(ts_obs_se + 1.7))
                    pr_min = int(np.min(pr_obs_se - 1.7))
                    pr_max = int(np.max(pr_obs_se + 2.7))
                    axs[row, col].set_xlim(ts_min, ts_max)
                    axs[row, col].set_ylim(pr_min, pr_max)
                    if (ts_max - ts_min > 4):
                        axs[row, col].locator_params(axis='x', nbins=5)
                    else:
                        axs[row, col].locator_params(axis='x',
                                                     nbins=ts_max - ts_min)
                    if (pr_max - pr_min > 6):
                        axs[row, col].locator_params(axis='y', nbins=7)
                    else:
                        axs[row, col].locator_params(axis='y',
                                                     nbins=pr_max - pr_min)

                axs[row, -1].yaxis.set_label_position("right")
                axs[row, -1].set_ylabel(area + "\n" + c_area,
                                        rotation=0, horizontalalignment='left')
                axs[row, col].grid(plot_grid)

        # Closing the model datafiles
        pr_file.close()
        ts_file.close()

        axs[-1, 0].set_ylabel("Precipitation [" + pr_units + "]")
        axs[-1, 0].set_xlabel("SST [" + ts_units + "]")
        plt.suptitle("Seasonal (annual/monthly) mean values of "
                     + "SST and precipitation for \n"
                     + model
                     + " (coloured) vs "
                     + ts_model
                     + " (sst obs) and "
                     + pr_model
                     + " (pr obs)", fontsize=20)
        # Saving the plot and letting the user know what just happened
        variable = ''.join([ts_key,
                             pr_key])

        diag_name = E.get_diag_script_name()
        output_file = E.get_plot_output_filename(diag_name=diag_name,
                                                 variable=variable,
                                                 model=model,
                                                 specifier=experiment)
        output_dir = os.path.join(plot_dir, diag_name)
        E.ensure_directory(output_dir)
        plt.savefig(os.path.join(output_dir, output_file))
        info("", verbosity, 1)
        info("Created image: ", verbosity, 1)
        info(output_file, verbosity, 1)

    # After all the loops remember to close the observation files aswell
    pr_obs_file.close()
    ts_obs_file.close()


def process_zonal_means(E, modelconfig):
    """Main script for plotting zonal means. Outputs variable specific plots
    containing all models for that variable (including observations). """
    experiment = 'zonal_means'
    datakeys = E.get_currVars()
    plot_dir = E.get_plot_dir()
    verbosity = E.get_verbosity()
    areas = modelconfig.get(experiment, 'areas').split()
    plot_grid = modelconfig.getboolean('general', 'plot_grid')

    # We don't do the zonal means for winds (equatorial elsewhere)
    if 'ua' in datakeys:
        datakeys.remove('ua')

    for area in areas:
        plt.clf()
        fig, axs = plt.subplots(2, 2, figsize=(12, 8))
        fig.subplots_adjust(top=0.87)
        fig.subplots_adjust(right=0.67)
        fig.subplots_adjust(hspace=0.2)
        fig.subplots_adjust(wspace=0.4)

        for datakey in datakeys:
            col = datakeys.index(datakey)
            # We output a file for each variable
            model_filenames = E.get_clim_model_filenames(variable=datakey)

            # Start looping and generating subplots
            ymin, ymax = False, False
            lon_min = modelconfig.getint(experiment + '_' + area, 'lon_min')
            lon_max = modelconfig.getint(experiment + '_' + area, 'lon_max')
            lon_area = E.get_ticks_labels(np.array([lon_min, lon_max]), 'lons')
            ymin = -9999
            ymax = -9999
            sd = np.zeros(len(model_filenames))
            mindex = 0
            for model in model_filenames:
                # read in specified data and concatenate
                datafile = nc.Dataset(model_filenames[model], 'r')
                # A-laue_ax+
                E.add_to_filelist(model_filenames[model])
                # A-laue_ax-
                data_units = datafile.variables[datakey].units
                lats, data = E.get_model_data(modelconfig,
                                              experiment,
                                              area,
                                              datakey,
                                              datafile)
                data = E.average_data(data, 1)
                datafile.close()
                norm = data.mean()

                # get y min and max iteratively

                if (ymin == -9999):
                    ymin = data.min()
                    ymax = data.max()
                else:
                    ymin = np.min([data.min(), ymin])
                    ymax = np.max([data.max(), ymax])
                sd[mindex] = (data / norm).std()
                mindex += 1

                # read model specific plotting style
                color, dashes, width = E.get_model_plot_style(model)

                # plotting custom dashes requires some extra effort (else)
                # with empty dashes the format is default
                if (len(dashes) == 0):
                    axs[0, col].plot(lats, data, color=color,
                                     linewidth=width, label=model)
                    axs[1, col].plot(lats, data / norm, color=color,
                                     linewidth=width, label=model)
                else:
                    line1, = axs[0, col].plot(lats, data, '--', color=color,
                                              linewidth=width, label=model)
                    line1.set_dashes(dashes)
                    line2, = axs[1, col].plot(lats, data / norm, color=color,
                                              linewidth=width, label=model)
                    line2.set_dashes(dashes)

            # Now to make the plot pretty i.e. headers, ticks etc.
            # axs[-1, col].yaxis.set_label_position("right")
            axs[0, col].grid(plot_grid)
            axs[1, col].grid(plot_grid)

            area_key = experiment + '_' + area
            xmin = modelconfig.getint(area_key, 'lat_min')
            xmax = modelconfig.getint(area_key, 'lat_max')
            axs[0, col].set_xlim(xmin, xmax)
            axs[1, col].set_xlim(xmin, xmax)

            xticks = E.get_ticks(5, np.array([xmin, xmax]))
            labels = E.get_ticks_labels(xticks, 'lats')

            axs[0, col].xaxis.set_ticks(xticks)
            axs[0, col].set_xticklabels(labels)
            axs[1, col].xaxis.set_ticks(xticks)
            axs[1, col].set_xticklabels(labels)

            yticks = E.get_ticks(8, np.array([ymin, ymax]))
            axs[0, col].yaxis.set_ticks(yticks)
            axs[1, col].yaxis.set_ticks([1.0])
            axs[1, col].set_yticklabels(['$\mu$'])
            #axs[1, col].yaxis.set_ticks([1.0-sd.mean(), 1.0, 1.0+sd.mean()])
            #axs[1, col].set_yticklabels(['$\mu$-$\sigma$', '$\mu$', '$\mu$+$\sigma$'])

            if   (datakey == 'pr' or datakey == 'pr-mmday'):
                axs[0, col].set_title("Precipitation")
                axs[0, col].set_ylabel("[" + data_units + "]")
                axs[1, col].set_ylabel("normalized")

            elif (datakey == 'ts'):
                axs[0, col].set_title("Temperature")
                axs[0, 0].set_ylabel("[" + data_units + "]")
                axs[1, 0].set_ylabel("normalized")

        plt.suptitle(area
                     + " ocean ["
                     + lon_area[0]
                     + ":"
                     + lon_area[1]
                     + "] seasonal mean", fontsize=16)
        # Adding the legend
        handles0, labels0 = axs[0, 0].get_legend_handles_labels()
        handles1, labels1 = axs[0, 1].get_legend_handles_labels()
        handles_all = handles0 + handles1
        labels_all = labels0 + labels1
        labels_sorted = zip(*sorted(zip(labels_all, handles_all)))[0]
        handles_sorted = zip(*sorted(zip(labels_all, handles_all)))[1]
        handles = []
        labels = []
        for item in xrange(len(labels_sorted)):
            if (labels_sorted[item] not in labels):
                labels.append(labels_sorted[item])
                handles.append(handles_sorted[item])
        plt.legend(handles, labels, loc='center left',
                   bbox_to_anchor=(0.7, 0.5),
                   bbox_transform=plt.gcf().transFigure)
        # And save the plot
        diag_name = E.get_diag_script_name()
        output_dir = os.path.join(plot_dir, diag_name)
        E.ensure_directory(output_dir)
        specifier = area + '-seasonal-mean'
        variable = "".join(E.get_currVars())
        output_file = E.get_plot_output_filename(diag_name=diag_name,
                                                 variable=variable,
                                                 specifier=specifier)
        plt.savefig(os.path.join(output_dir, output_file))
        info("", verbosity, 1)
        info("Created image: ", verbosity, 1)
        info(output_file, verbosity, 1)
