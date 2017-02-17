"""
;;#############################################################################
;; TropicalVariablity_wind.py
;; Author: Jarmo Makela (FMI, Finland)
;; EMBRACE project
;;#############################################################################
;; Description
;;    This script is the diagnostics and plotting script for
;;    Tropical Variability Equatorial means of
;;    precipitation, sea surface temperature and surface winds
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
;;    lat_min:               Spatial extent for Central-Pacific
;;    lat_max:               Spatial extent for Central-Pacific
;;    lon_min:               Spatial extent for Central-Pacific
;;    lon_max:               Spatial extent for Central-Pacific
;;    season_limits_annual:  Spatial extent for annual for Central-Pacific
;;    season_limits_DJF:     Spatial extent for DJF for Central-Pacific
;;    season_limits_MAM:     Spatial extent for MAM for Central-Pacific
;;    season_limits_JJA:     Spatial extent for JJA for Central-Pacific
;;    season_limits_SON:     Spatial extent for SON for Central-Pacific
;;
;;    [scatter_East-Pacific]
;;    lat_min:               Spatial extent for East-Pacific
;;    lat_max:               Spatial extent for East-Pacific
;;    lon_min:               Spatial extent for East-Pacific
;;    lon_max:               Spatial extent for East-Pacific
;;    season_limits_annual:  Spatial extent for annual for East-Pacific
;;    season_limits_DJF:     Spatial extent for DJF for East-Pacific
;;    season_limits_MAM:     Spatial extent for MAM for East-Pacific
;;    season_limits_JJA:     Spatial extent for JJA for East-Pacific
;;    season_limits_SON:     Spatial extent for SON for East-Pacific
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

;;    See TODO in 'process_equatorial_means' about flexibilty wrt obs
;;
;;    Dependency - must be run after "diag_script/TropicalVariability_wind.py"
;;                 which in turn depends on "TropicalVariability_wind.py".
;;
;; Modification history
;;    20151029-A_laue_ax: added output of acknowledgements + processed files
;;                        to log-file
;;    20150515-A_maek_ja: written
;;
;; #############################################################################
"""
"""
*********************************************************************
 TropicalVariablity_EQ.py
*********************************************************************
 PYTHON script
 Diagnostics: TropicalVariability
 Date:        January 2015
 Author:      Jarmo Makela (FMI, Finland, jarmo.makela@fmi.fi)
*******************************************************************************
 This script is the diagnostics and plotting script for
 Tropical Variability Equatorial means of
 precipitation, sea surface temperature and surface winds
*******************************************************************************
"""
# Basic Python packages
import ConfigParser
import os
import pdb
import sys
import numpy as np
from glob import glob

# Force matplotlib to not use any Xwindows backend.
import matplotlib
matplotlib.use('Agg')

# Common Python packages
import matplotlib.pyplot as plt

# Python netcdf handler
import netCDF4 as nc

# ESMValTool defined Python packages
sys.path.append("./interface_scripts")
from esmval_lib import ESMValProject
from auxiliary import info


def main(project_info):
    """Diagnostics and plotting script for Tropical Variability Equatorial.
    We use ts as a proxy for Sea Surface Temperature. """

    # ESMValProject provides some easy methods to access information in
    # project_info but you can also access it directly (as in main.py)
    # We need to ensure that needed files for precipitation and temperature
    # are present before we can proceed

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
        info("Starting to plot equatorial means", verbosity, 2)
        process_equatorial_means(E, modelconfig)
        post_process(E)


# All callable functions below are in alphabetical order
# E. style functions are in the general python file esmval_lib.py

def check_data_existance(E, modelconfig):
    """Checks that all needed datafiles exist for temp / precip. """
    experiment = 'equatorial'
    areas = modelconfig.get(experiment, 'areas').split()
    work_dir = E.get_work_dir()
    base_name = 'TropicalVariability_' + experiment + '_' + areas[0] + '_'

    additional_keys = []
    # First test temperature
    ts_test = glob(work_dir + base_name + 'ts*')
    if (len(ts_test) > 0):
        additional_keys.append('ts')
    # Next precipitation
    pr_test = glob(work_dir + base_name + 'pr*')
    if (len(pr_test) > 0):
        if ('pr-mmday' in pr_test[0]):
            additional_keys.append('pr-mmday')
        else:
            additional_keys.append('pr')
    # Finally divergence
    divergence_test = glob(work_dir + base_name + 'divergence*')
    if (len(divergence_test) > 0):
        additional_keys.append('divergence')

    if (len(ts_test) > 0 and len(pr_test) > 0 and len(divergence_test) > 0):
        return additional_keys
    else:
        print("PY  ERROR: I'm missing some of the needed files for equatorial plots")
        print("PY  ERROR: These should've been automatically generated by")
        print("PY  ERROR: TropicalVariability.py")
        print("PY  ERROR: There's no point for me to proceed, exiting")
        sys.exit()


def post_process(E):
    """This script is for post processing of equatorial means.
    We erase temporary npy files from memory."""
    base_names = E.get_work_dir()
    base_names += 'TropicalVariability_equatorial_*'
    erase_files = glob(base_names)
    for erase in erase_files:
        os.remove(erase)


def process_equatorial_means(E, modelconfig):
    """Main script for plotting equatorial means. Outputs one image with
    four subplots (precipitation, temperature, equatorial winds and divergence).
    TODO: This script is not too flexible - since we're using ncl for preprocessing
    we cant exclude multiple obs files - so we can use only two different
    observation models (so they can exclude one another). """
    experiment = 'equatorial'
    datakeys = E.get_currVars()
    plot_dir = E.get_plot_dir()
    verbosity = E.get_verbosity()
    work_dir = E.get_work_dir()
    areas = modelconfig.get(experiment, 'areas').split()
    plot_grid = modelconfig.getboolean('general', 'plot_grid')
    ua_key = datakeys[0]

    # Get the model paths etc. including the obs file for winds
    models = E.get_clim_model_filenames(ua_key)
    # Get other variables and check data existance
    datakeys.extend(check_data_existance(E, modelconfig))

    # Looping over areas etc.
    for area in areas:
        area_key = experiment + '_' + area
        lat_min = modelconfig.getint(area_key, 'lat_min')
        lat_max = modelconfig.getint(area_key, 'lat_max')
        lon_min = modelconfig.getint(area_key, 'lon_min')
        lon_max = modelconfig.getint(area_key, 'lon_max')
        lat_area = E.get_ticks_labels(np.array([lat_min, lat_max]), 'lats')
        # Initial figure config
        fig, axs = plt.subplots(4, 1, figsize=(12, 11))
        fig.subplots_adjust(top=0.88)
        fig.subplots_adjust(right=0.67)
        fig.subplots_adjust(hspace=0.4)
        for datakey in datakeys:
            # Output precipitation / temperature and winds
            # Adjust titles accordingly
            suptitle = area + " ocean equatorial ["
            suptitle += lat_area[0] + ":" + lat_area[1] + "] mean"
            plt.suptitle(suptitle, fontsize=16)
            if (datakey == 'pr' or datakey == 'pr-mmday'):
                row = 0
                title = "Precipitation"
            elif (datakey == 'ts'):
                row = 1
                title = "Temperature"
            elif (datakey == 'ua' or datakey == 'ua-1000'):
                row = 2
                title = "Eastward wind"
            elif (datakey == 'divergence'):
                row = 3
                title = "Wind divergence"

            # The data is handled differently for winds and the rest
            # For winds we have to extract the data for each model
            # For the rest we only have one file / datakey
            scaling = ''
            if (datakey == 'ua' or datakey == 'ua-1000'):
                for model in models:
                    # extract the model specific data
                    datafile = nc.Dataset(models[model], 'r')
                    # A-laue_ax+
                    E.add_to_filelist(models[model])
                    # A-laue_ax-
                    data_units = datafile.variables[datakey].units
                    lats, lons, data = E.get_model_data(modelconfig, experiment,
                                                        area, datakey, datafile)
                    datafile.close()
                    lons = E.ensure_looping(lons)

                    # Take the means and get the plot style
                    data = E.average_data(data, 2)
                    color, dashes, width = E.get_model_plot_style(model)

                    if (len(dashes) == 0):
                        axs[row].plot(lons, data, color=color,
                                      linewidth=width, label=model)
                    else:
                        line, = axs[row].plot(lons, data, '--', color=color,
                                              linewidth=width, label=model)
                        line.set_dashes(dashes)

                ymin = modelconfig.getint(area_key, 'wind_min')
                ymax = modelconfig.getint(area_key, 'wind_max')
                yticks = E.get_ticks(8, np.array([ymin, ymax]))

            # Next the datakeys with all data in one file (per datakey)
            # These will be read from external files previously written by
            # TropicalVariability.py
            elif (datakey == 'pr' or datakey == 'pr-mmday' or datakey == 'ts'):
                base_name = work_dir + 'TropicalVariability_'\
                            + area_key + '_' + datakey + '_model_'
                datafiles = glob(base_name + '*')
                key_models = []
                for dfile in datafiles:
                    start = dfile.find('_model_') + 7
                    key_models.append(str(dfile)[start:-4])

                # now we can do the loop over the models
                for model in key_models:
                    # extract model specific data
                    alldata = np.load(base_name + model + '.npy')
                    lons = alldata[0]
                    data = alldata[1]
                    color, dashes, width = E.get_model_plot_style(model)
                    lons = E.ensure_looping(lons)

                    # get the style and plot it out
                    if (len(dashes) == 0):
                        axs[row].plot(lons, data, color=color,
                                      linewidth=width, label=model)
                    else:
                        line, = axs[row].plot(lons, data, '--', color=color,
                                              linewidth=width, label=model)
                        line.set_dashes(dashes)

                # Some additional figure config
                if   (datakey == 'pr' or datakey == 'pr-mmday'):
                    ymin = modelconfig.getint(area_key, 'prec_min')
                    ymax = modelconfig.getint(area_key, 'prec_max')
                    yticks = E.get_ticks(6, np.array([ymin, ymax]))
                    if (datakey == 'pr'):
                        data_units = 'kg/m2s'
                    else:
                        data_units = 'mm/day'
                elif (datakey == 'ts'):
                    ymin = modelconfig.getint(area_key, 'temp_min')
                    ymax = modelconfig.getint(area_key, 'temp_max')
                    yticks = E.get_ticks(6, np.array([ymin, ymax]))
                    data_units = 'K'

            # Next the wind divergence that is also read from external files
            # previously written by TropicalVariability_wind.py
            elif (datakey == 'divergence'):
                base_name = work_dir + 'TropicalVariability_'\
                            + area_key + '_' + datakey + '_model_'
                datafiles = glob(base_name + '*')
                key_models = []
                for dfile in datafiles:
                    start = dfile.find('_model_') + 7
                    key_models.append(str(dfile)[start:-4])

                # now we can do the loop over the models
                for model in key_models:
                    # extract model specific data
                    alldata = np.load(base_name + model + '.npy')
                    lons = alldata[0]
                    data = 1E6 * alldata[1]
                    color, dashes, width = E.get_model_plot_style(model)
                    lons = E.ensure_looping(lons)

                    # get the style and plot it out
                    if (len(dashes) == 0):
                        axs[row].plot(lons, data, color=color,
                                      linewidth=width, label=model)
                    else:
                        line, = axs[row].plot(lons, data, '--', color=color,
                                              linewidth=width, label=model)
                        line.set_dashes(dashes)

                ymin = modelconfig.getint(area_key, 'div_min')
                ymax = modelconfig.getint(area_key, 'div_max')
                yticks = E.get_ticks(6, np.array([ymin, ymax]))
                data_units = '1/s'
                scaling = '1E6'

            # This part is the same for all subplots
            axs[row].set_ylim(ymin, ymax)
            axs[row].yaxis.set_ticks(yticks)
            xticks = E.get_ticks(9, np.array([lon_min, lon_max]))
            labels = E.get_ticks_labels(xticks, 'lons')
            axs[row].xaxis.set_ticks(xticks)
            axs[row].set_xticklabels(labels)
            if (lon_min < lon_max):
                axs[row].set_xlim(lon_min, lon_max)
            else:
                axs[row].set_xlim(lon_min, 360 + lon_max)
            if (len(scaling) > 0):
                axs[row].set_title(title + ' ' + scaling + ' * [' + data_units + ']')
            else:
                axs[row].set_title(title + ' [' + data_units + ']')
            axs[row].grid(plot_grid)

        # Adding the legend
        handles0, labels0 = axs[0].get_legend_handles_labels()
        handles1, labels1 = axs[1].get_legend_handles_labels()
        handles2, labels2 = axs[2].get_legend_handles_labels()
        handles3, labels3 = axs[3].get_legend_handles_labels()
        handles_all = handles0 + handles1 + handles2 + handles3
        labels_all = labels0 + labels1 + labels2 + labels3
        # Sort both lists based on labels
        labels_sorted = zip(*sorted(zip(labels_all, handles_all)))[0]
        handles_sorted = zip(*sorted(zip(labels_all, handles_all)))[1]
        # Now we remove duplicates
        handles = []
        labels = []
        for item in xrange(len(labels_sorted)):
            if (labels_sorted[item] not in labels):
                labels.append(labels_sorted[item])
                handles.append(handles_sorted[item])
        # Plotting the legend at its proper place
        plt.legend(handles, labels, loc='center left',
                   bbox_to_anchor=(0.7, 0.5),
                   bbox_transform=plt.gcf().transFigure)
        # Saving the plot and printing the image location
        diag_name = E.get_diag_script_name()
        output_dir = os.path.join(plot_dir, diag_name)
        E.ensure_directory(output_dir)
        plt.suptitle(area + " ocean mean wind direction and strength", fontsize=16)
        specifier = '-'.join([area,
                              experiment,
                              'mean'])
        variable = "".join(E.get_currVars())
        output_file = E.get_plot_output_filename(diag_name=diag_name,
                                                 variable=variable,
                                                 specifier=specifier)
        plt.savefig(os.path.join(output_dir, output_file))
        info("", verbosity, 1)
        info("Created image: ", verbosity, 1)
        info(output_file, verbosity, 1)
