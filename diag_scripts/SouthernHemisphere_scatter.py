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
;; colourmap_clouds: Python matplotlib colormaps, suffix '_r' = reversed cmap
;; colourmap_model:  Python matplotlib colormaps, suffix '_r' = reversed cmap
;; colourmap_diff:   Python matplotlib colormaps, suffix '_r' = reversed cmap
;; colourmap_dev:    Python matplotlib colormaps, suffix '_r' = reversed cmap
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
import sys
import os
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

    # Check which parts of the code to run
    if (modelconfig.getboolean('general', 'plot_scatter')):
        info("Starting scatter plot", verbosity, 2)
        process_scatter(E, modelconfig)


# All callable functions below are in alphabetical order
# E. style functions are in the general python file esmval_lib.py

def calculate_scatterplot_values(modelconfig,
                                 area,
                                 cl_key,
                                 cloud,
                                 radiation,
                                 bins=''):
    """Calculate cloud vs radiation scatterplot values from input."""
    # Check if limits are already defined
    if (len(bins) > 0):
        nbins = len(bins) - 1
    else:
        nbins = modelconfig.getint('SouthernHemisphere_scatter_' + area,
                                   'points')
        if (cl_key == 'clt'):
            cl_min = 0
            cl_max = 100
        else:
            cl_min = cloud.min()
            cl_max = cloud.max()
        if not (cl_min > cl_max):
            cl_min = 0
            cl_max = 100
        nbins = 20
        bins = np.linspace(cl_min, cl_max, nbins + 1)
    # Mask unwanted values (open ended at end limits)
    cl_fraction = np.zeros(nbins)
    cl_values = np.zeros(nbins)
    rd_values = np.zeros(nbins)
    for i in xrange(nbins):
        # Open ended for end points
        if (i == 0):
            cl_masked = np.ma.masked_greater(cloud, bins[i + 1])
        elif (i == nbins - 1):
            cl_masked = np.ma.masked_less(cloud, bins[i])
        else:
            cl_masked = np.ma.masked_outside(cloud, bins[i], bins[i + 1])
        rd_masked = np.ma.masked_array(radiation, cl_masked.mask)
        cl_values[i] = cl_masked.mean()
        rd_values[i] = rd_masked.mean()
        cl_fraction[i] = np.ma.count(cl_masked)
    cl_sum = cl_fraction.sum()
    cl_fraction[:] = cl_fraction[:] / cl_sum
    return bins, cl_fraction, cl_values, rd_values


def interpolate_3d(data, lats, lons, target_lats, target_lons):
    """Interpolate time/lat/lon grid to specific lat/lon coordinates."""
    new_data = np.zeros((data.shape[0], len(target_lats), len(target_lons)))
    for time in xrange(data.shape[0]):
        new_data[time] = interpolate_data_grid(data[time, :, :],
                                               lats,
                                               lons,
                                               target_lats,
                                               target_lons)
    return new_data


def interpolate_data_grid(data, lats, lons, target_lats, target_lons):
    """Interpolates the data values to a specific lat/lon grid.
    This function should only be used for 2D arrays (no time indeces etc.)"""

    # First check if the coordinates are the same, otherwise interpolate
    if (np.array_equal(lats, target_lats)
        and np.array_equal(lons, target_lons)):
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


def process_scatter(E, modelconfig):
    """
    Main script for gathering cloud vs radiation values and plotting them.
    """
    config_file = E.get_configfile()
    experiment = 'SouthernHemisphere'
    areas = modelconfig.get(experiment, 'scatter_areas').split()
    seasons = modelconfig.get(experiment, 'seasons').split()
    plot_grid = modelconfig.getboolean('general', 'plot_background_grid')
    plot_dir = E.get_plot_dir()
    verbosity = E.get_verbosity()
    work_dir = E.get_work_dir()
    datakeys = E.get_currVars()

    # Check which cloud diagnostic we are working with
    if   ('clt' in datakeys):
        cl_key = 'clt'
    elif ('clivi' in datakeys):
        cl_key = 'clivi'
    elif ('clwvi' in datakeys):
        cl_key = 'clwvi'

    # Extract cloud observations
    datakeys.remove(cl_key)
    cl_obs, cl_loc, cl_models = E.get_clim_model_and_obs_filenames(cl_key)
    cl_obsfile = nc.Dataset(cl_loc, 'r')
    # A-laue_ax+
    E.add_to_filelist(cl_loc)
    # A-laue_ax-

    # Basic structure of the figure is same for all
    # plots - we clear it at intervals
    for area in areas:
        # Get area specific configuration and observations
        coords = E.get_area_coordinates(modelconfig, experiment, area)
        area_key = experiment + '_scatter_' + area
        cl_lats, cl_lons, cl_data = E.get_model_data(modelconfig, experiment,
                                                     area, cl_key, cl_obsfile)
        cl_data = E.average_data(cl_data, 'annual')
        for datakey in datakeys:
            # Extract radiation observations
            rd_obs, rd_loc, rd_models = E.get_clim_model_and_obs_filenames(datakey)
            rd_obsfile = nc.Dataset(rd_loc, 'r')
            # A-laue_ax+
            E.add_to_filelist(rd_loc)
            # A-laue_ax-
            rd_lats, rd_lons, rd_data = E.get_model_data(modelconfig,
                                                         experiment,
                                                         area,
                                                         datakey,
                                                         rd_obsfile,
                                                         extend='lons')
            # Check that we have similar model sets
            E.check_model_instances(cl_models, rd_models)

            # Average and interpolate to cloud obsevation grid
            rd_data = E.average_data(rd_data, 'annual')
            rd_data = interpolate_3d(rd_data,
                                     rd_lats,
                                     rd_lons,
                                     cl_lats,
                                     cl_lons)
            for model in rd_models:
                # Get model and area specific values
                cl_datafile = nc.Dataset(cl_models[model], 'r')
                rd_datafile = nc.Dataset(rd_models[model], 'r')
                # A-laue_ax+
                E.add_to_filelist(cl_models[model])
                E.add_to_filelist(rd_models[model])
                # A-laue_ax-
                mcl_lats, mcl_lons, mcl_data = E.get_model_data(modelconfig,
                                                                experiment,
                                                                area,
                                                                cl_key,
                                                                cl_datafile)
                mrd_lats, mrd_lons, mrd_data = E.get_model_data(modelconfig,
                                                                experiment,
                                                                area,
                                                                datakey,
                                                                rd_datafile,
                                                                extend='lons')
                cl_datafile.close()
                rd_datafile.close()
                # Average and interpolate model data to cloud data grid
                mcl_data = E.average_data(mcl_data, 'annual')
                mrd_data = E.average_data(mrd_data, 'annual')
                mrd_data = interpolate_3d(mrd_data,
                                          mrd_lats,
                                          mrd_lons,
                                          mcl_lats,
                                          mcl_lons)

                # One plot for all seasons
                plt.clf()
                fig, axs = plt.subplots(2, (len(seasons)),
                                        figsize=(4 * len(seasons), 8))
                fig.subplots_adjust(top=0.8)
                fig.subplots_adjust(bottom=0.15)
                fig.subplots_adjust(right=0.75)
                fig.subplots_adjust(hspace=0.4)
                fig.subplots_adjust(wspace=0.4)
                for season in seasons:
                    col = seasons.index(season)
                    # Mask unwanted seasonal values
                    obs_cl = E.extract_seasonal_mean_values(modelconfig,
                                                            cl_data,
                                                            experiment,
                                                            season,
                                                            monthly=True)
                    obs_rd = E.extract_seasonal_mean_values(modelconfig,
                                                            rd_data,
                                                            experiment,
                                                            season,
                                                            monthly=True)
                    data_cl = E.extract_seasonal_mean_values(modelconfig,
                                                             mcl_data,
                                                             experiment,
                                                             season,
                                                             monthly=True)
                    data_rd = E.extract_seasonal_mean_values(modelconfig,
                                                             mrd_data,
                                                             experiment,
                                                             season,
                                                             monthly=True)
                    # Calculate scatterplot values
                    obs_out = calculate_scatterplot_values(modelconfig,
                                                           area,
                                                           cl_key,
                                                           obs_cl,
                                                           obs_rd)
                    model_out = calculate_scatterplot_values(modelconfig,
                                                             area,
                                                             cl_key,
                                                             data_cl,
                                                             data_rd,
                                                             bins=obs_out[0])
                    centers = np.zeros(len(obs_out[0]) - 1)
                    for i in xrange(len(centers)):
                        centers[i] = (obs_out[0][i + 1] + obs_out[0][i]) / 2.

                    # Plot and make it pretty
                    modelcolor, dashes, width = E.get_model_plot_style(model)
                    obs_label = cl_obs + ', \n' + rd_obs
                    axs[0, col].scatter(model_out[2],
                                        model_out[3],
                                        facecolors=modelcolor,
                                        edgecolors=modelcolor,
                                        label=model)
                    axs[0, col].scatter(obs_out[2],
                                        obs_out[3],
                                        facecolors='none',
                                        edgecolors='black',
                                        label=obs_label)

                    if (len(dashes) == 0):
                        axs[1, col].plot(centers,
                                         obs_out[1],
                                         color=modelcolor,
                                         linewidth=width,
                                         label=model)
                    else:
                        line1, = axs[1, col].plot(centers,
                                                  model_out[1],
                                                  '--',
                                                  color=modelcolor,
                                                  linewidth=width,
                                                  label=model)
                        line1.set_dashes(dashes)
                    axs[1, col].plot(centers,
                                     obs_out[1],
                                     color='black',
                                     label=obs_label)

                    # If we get nan as values we convert it to zero
                    obs_x = np.ma.masked_array(obs_out[2])
                    obs_y = np.ma.masked_array(obs_out[3])
                    model_x = np.ma.masked_array(model_out[2])
                    model_y = np.ma.masked_array(model_out[3])
                    rmse = (np.mean((obs_x - model_x) ** 2)
                            + np.mean((obs_y - model_y) ** 2)) ** 0.5
                    rmse = round(rmse, 1)
                    xmin = int(obs_out[0][0])
                    xmax = int(np.ceil(obs_out[0][-1]))

                    # Make it prettier
                    axs[0, col].set_title(season
                                          + " (RMSE: "
                                          + str(rmse) + ")",
                                          fontsize=12)
                    axs[0, col].set_xlim(xmin, xmax)
                    axs[0, col].locator_params(axis='x', nbins=5)
                    axs[0, col].locator_params(axis='y', nbins=5)
                    axs[0, col].grid(plot_grid)
                    axs[1, col].set_xlim(xmin, xmax)
                    axs[1, col].locator_params(axis='x', nbins=5)
                    axs[1, col].locator_params(axis='y', nbins=5)
                    axs[1, col].grid(plot_grid)

                xlabel = E.get_title_basename(cl_key)
                ylabel = E.get_title_basename(datakey)
                axs[0, 0].set_xlabel(xlabel)
                axs[0, 0].set_ylabel('Radiation')

                axs[1, 0].set_xlabel('Cloud cover percentage')
                axs[1, 0].set_ylabel('Proportion of values')

                # Title and legend
                suptitle = ylabel + " sensitivity to " + xlabel
                plt.suptitle(suptitle, fontsize=20)
                handles0, labels0 = axs[0, 0].get_legend_handles_labels()
                handles1, labels1 = axs[1, 0].get_legend_handles_labels()
                handles = handles0 + handles1
                labels = labels0 + labels1
                plt.legend(handles, labels, loc='center left',
                           bbox_to_anchor=(0.78, 0.5),
                           bbox_transform=plt.gcf().transFigure)

                # Define ouput name and save
                variable = cl_key + '-' + datakey
                diag_name = E.get_diag_script_name()
                output_file = E.get_plot_output_filename(diag_name=diag_name,
                                                         variable=variable,
                                                         model=model)
                output_dir = os.path.join(plot_dir, diag_name)
                E.ensure_directory(output_dir)
                plt.savefig(os.path.join(output_dir, output_file))
                plt.clf()
                info("", verbosity, 1)
                info("Created image: ", verbosity, 1)
                info(output_file, verbosity, 1)
    cl_obsfile.close()
    rd_obsfile.close()
