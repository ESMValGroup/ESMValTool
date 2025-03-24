"""
Model 1 vs Model 2 vs Observations diagnostics.
===============================================

Diagnostic to produce an image showing four maps, based on a comparison of two
different models results against an observational dataset. This process is
often used to compare a new iteration of a model under development against
a previous version of the same model. The four map plots are:

* Top left: model 1
* Top right: model 1 minus model 2
* Bottom left: model 2 minus obs
* Bottom right: model 1 minus obs

All four plots show latitude vs longitude and the cube value is used as the
colour scale.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, a small number of depth layers,
and a latitude and longitude coordinates.

This diagnostic also requires the ``exper_model``, ``control_model`` and
``observational_dataset`` keys in the recipe.

This tool is part of the ocean diagnostic tools package in the ESMValTool,
and was based on the plots produced by the Ocean Assess/Marine Assess toolkit.

Original script:
Author: Lee de Mora (PML)
        ledm@pml.ac.uk

Refactored script: (Alphabetical order)
Author: Sophie Hall (Met Office)
        sophie.hall@metoffice.gov.uk
Author: Emma Hogan (Met Office)
        emma.hogan@metoffice.gov.uk
Author: Dave Storkey (Met Office)
        dave.storkey@metoffice.gov.uk
"""

import logging
import os
import sys

import cartopy.crs as ccrs
import iris
import iris.analysis.cartography
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic, save_figure

# Create a logger object.
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def main(config):
    """Loading the configuration file and running through the order of
    operation.

    Parameters
    ----------
    config : dictionary
        configuration dictionary that contains all the necessary information
        for the function to run. It includes details about the models,
        observational datasets, file paths, and other settings.
    """
    # Call the load_data function
    experiment, control, observation = load_data(config)

    # Call create_plotting_data
    (exp, exp_minus_ctr, ctr_minus_obs,
     exp_minus_obs) = create_plotting_data(control, experiment, observation)

    # If 'Sea Surface' set only 1 level to plot.
    # Otherwise set array of depth levels to loop through.
    if experiment.long_name in [
            'Sea Surface Temperature', 'Sea Surface Salinity'
    ]:
        levels = [1]
    else:
        # Levels are the chosen depths we wish to plot.
        levels = [2, 50, 100, 300]

    # Loop over levels to extract_single_level.
    for level in levels:
        exp_single_level = extract_global_single_level(exp, level)
        exp_minus_ctr_single_level = extract_global_single_level(
            exp_minus_ctr, level)
        ctr_minus_obs_single_level = extract_global_single_level(
            ctr_minus_obs, level)
        exp_minus_obs_single_level = (extract_global_single_level(
            exp_minus_obs, level))

        # Call create_quadmap, which contains plot_global_single_level
        create_quadmap(exp_single_level, exp_minus_ctr_single_level,
                       ctr_minus_obs_single_level, exp_minus_obs_single_level)

    # After successfully generating plots, function logs a success message.
    logger.info('Success')


def load_data(config):
    """Loads in all necessary data to output experiment, control, observation.
    Fixing variable names, fixing units then passed into create_plotting_data()

    Parameters
    ----------
    config : dictionary
        configuration dictionary that contains all the necessary information
        for the function to run. It includes details about the
        models, observational datasets, file paths, and other settings.

    Returns
    -------
    experiment : iris cube
        Data as defined as exper_model from the recipe.
    control : iris cube
        Data as defined as control_model from the recipe.
    observation : iris cube
        Data as defined as observational_dataset from the recipe.
    """
    # The datasets defined by:
    # exper_model, control_model, observational_datasets in the recipe
    # determine what is loaded in this function.
    exp_label_from_recipe = 'exper_model'
    ctl_label_from_recipe = 'control_model'
    obs_label_from_recipe = 'observational_dataset'

    # Getting all input files from configuration.
    input_files = diagtools.get_input_files(config)

    # Get experiment input files from config
    exp_filename = diagtools.match_model_to_key(exp_label_from_recipe,
                                                config[exp_label_from_recipe],
                                                input_files)
    # Get control input files from config
    ctl_filename = diagtools.match_model_to_key(ctl_label_from_recipe,
                                                config[ctl_label_from_recipe],
                                                input_files)
    # Get observation input files from config
    obs_filename = diagtools.match_model_to_key(obs_label_from_recipe,
                                                config[obs_label_from_recipe],
                                                input_files)

    # Set variable names to filename cubes above.
    experiment = iris.load_cube(exp_filename)
    control = iris.load_cube(ctl_filename)
    observation = iris.load_cube(obs_filename)

    # Fixing all units
    experiment = diagtools.bgc_units(experiment,
                                     input_files[exp_filename]['short_name'])
    control = diagtools.bgc_units(control,
                                  input_files[ctl_filename]['short_name'])
    observation = diagtools.bgc_units(observation,
                                      input_files[obs_filename]['short_name'])

    return experiment, control, observation


def fix_cube(cube_to_fix, cube_with_good_values):
    """!!! ONLY REQUIRED WHEN WORKING WITH NEMO DATASETS !!! The latitude and
    longitude coordinates of the cubes are adjusted to ensure the grids align
    for performing operations later. Only required with working with NEMO
    datasets due to different grids.

    Parameters
    ----------
    cube_to_fix : iris cube
        The cubes experiment and control with bad coord values.
    cube_with_good_values : iris cube
        The cube observation with the good coord values.
    """
    # This is necessary for subtracting one cube from another.
    cube_to_fix.coord("latitude").points[:] = cube_with_good_values.coord(
        "latitude").points
    cube_to_fix.coord("longitude").points[:] = cube_with_good_values.coord(
        "longitude").points


def remove_extra_time_axis(cube):
    """!!! ONLY REQUIRED WHEN WORKING WITH NEMO DATASETS !!! Remove the extra
    time axis from the input provided by the cube parameter.

    Parameters
    ----------
    cube : iris cube
        Contains the climate data to be plotted. Including information like:
        temperature values, latitude, longitude, and depth.
    """
    # Counting the number of time coordinates in the cube.
    time_coords = [
        coord for coord in cube.coords()
        if iris.util.guess_coord_axis(coord) == 'T'
    ]
    # If the cube has multiple time coordinates with the same standard name,
    # the function removes the auxiliary time coordinates.
    time_axes_names = [coord.standard_name for coord in time_coords]

    if len(time_coords) >= 2 and len(set(time_axes_names)) == 1:
        for aux_coord in time_coords:
            if not aux_coord.is_dim_coord():
                cube.remove_coord(aux_coord)
    else:
        time_counter_coord = next(
            (coord
             for coord in cube.coords() if coord.var_name == 'time_counter'),
            None)
        # If the cube has a time_counter coordinate, the function removes it.
        if time_counter_coord:
            cube.remove_coord(time_counter_coord)


def create_plotting_data(control, experiment, observation):
    """Calculate the difference between the control, experiment and observation
    datasets to prepare data for plotting. Set as new variables.

    Parameters
    ----------
    control : iris cube
        Data as defined as control_model from the recipe.
    experiment : iris cube
        Data as defined as exper_model from the recipe.
    observation : iris cube
        Data as defined as observational_dataset from the recipe.

    Returns
    -------
    exp : iris cube
        Untouched experimental input.
    exp_minus_ctr : iris cube
        Experiment model minus control model.
    ctr_minus_obs : iris cube
        Control model minus observational dataset.
    exp_minus_obs : iris cube
        Experimental model minus observational dataset.
    """
    # The data for models and the obs dataset is loaded into Iris cubes.
    # These cubes contain the climate data that will be plotted.
    exp = experiment
    exp_minus_ctr = experiment - control
    ctr_minus_obs = control - observation
    exp_minus_obs = experiment - observation

    # Fixing the long_name of the cubes
    exp_minus_ctr.long_name = experiment.long_name
    ctr_minus_obs.long_name = experiment.long_name
    exp_minus_obs.long_name = experiment.long_name

    # Fixing exp_minus_ctr model_id for the plot title used later
    exp_minus_ctr.attributes['model_id'] = (
        f"{experiment.attributes['model_id']} minus "
        f"{control.attributes['model_id']}")

    # Fixing ctr_minus_obs model_id for the plot title used later
    ctr_minus_obs.attributes['model_id'] = (
        f"{control.attributes['model_id']} minus "
        f"{observation.attributes['source_id']}")

    # Fixing exp_minus_obs model_id for the plot title used later
    exp_minus_obs.attributes['model_id'] = (
        f"{experiment.attributes['model_id']} minus "
        f"{observation.attributes['source_id']}")

    return (exp, exp_minus_ctr, ctr_minus_obs, exp_minus_obs)


def extract_global_single_level(cube, level):
    """Extract a single level from the cube. Making all cubes 2D for plotting.

    Parameters
    ----------
    cube : iris cube
        The input data cube.
    level : float
        The depth level to extract.

    Returns
    -------
    single_level: iris cube
        The extracted single level cube.
    """
    # If cube.shape is already 2D make no change
    if len(cube.shape) == 2:
        single_level = cube
    # If the cube.shape is (1, 1, 1207, 1442) or (1, 1207, 1442): squeeze.
    elif len(cube.shape) > 2 and cube.shape[0] == 1 and (cube.shape[1] == 1 or
                                                         len(cube.shape) == 3):
        single_level = iris.util.squeeze(cube)
    # 3D cube - Takes cubes with depth levels and extracts one (called in main)
    else:
        constraint = iris.Constraint(depth=level)
        single_level = cube.extract(constraint)
        single_level = iris.util.squeeze(single_level)

    return single_level


def plot_global_single_level(axis, cube, contour_levels, title):
    """Create each individual plot before being added to create_quadmap.

    Parameters
    ----------
    axis: matplotlib 'ax'
        Represents one (sub-)plot in a figure.
        It contains the plotted data, axis ticks, labels, title, legend, etc.
    cube : iris cube
        This is a data structure that contains the climate data to be plotted.
        Including information like temperature values, latitude, longitude,
        and depth.
    contour_levels : numpy.array
        Used to set the ticks on the colour bar and used to define
        levels for the contour plot.
    title : str
        This is a string that will be used as the title of the subplot.
    """
    # Setting the colour of axis 1 always to viridis.
    if title == "HadGEM2-ES":
        cmap = 'viridis'
    # Setting the colour of all other plots to bwr.
    else:
        cmap = 'bwr'

    # This step transforms the data so it can be displayed as 2D
    new_cube, extent = iris.analysis.cartography.project(cube,
                                                         ccrs.PlateCarree(),
                                                         nx=400,
                                                         ny=200)

    # Sets the current Axes instance to the specified axis
    plt.sca(axis)

    # Creates a filled contour plot of the projected data.
    contour_result = iplt.contourf(new_cube,
                                   levels=contour_levels,
                                   linewidth=0,
                                   cmap=plt.cm.get_cmap(cmap))

    # Converts contour_levels to a numpy array.
    contour_levels = np.array(contour_levels)

    # Checks if the contour plot was created successfully
    if contour_result is None:
        raise ValueError(
            "Failed to create contour plot. The plt object is None.")

    # A color bar is added to the plot
    colorbar = plt.colorbar(contour_result, orientation='horizontal')
    # Colour scale dependent on range of data.
    colorbar.set_ticks([
        contour_levels.min(),
        (contour_levels.max() + contour_levels.min()) / 2.,
        contour_levels.max()
    ])
    # Adding latitude & longitude axis on each plot.
    fontsize = 7
    ax = plt.gca()
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
    gl.xlines = False
    gl.ylines = False
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': fontsize, 'color': 'gray'}
    gl.ylabel_style = {'size': fontsize, 'color': 'gray'}
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Coastlines are added to the map to provide geographical context.
    plt.gca().coastlines()

    # Sets the title of the plot
    plt.title(title)

    # Display the plots
    iplt.show()


def create_quadmap(exp_single_level, exp_minus_ctr_single_level,
                   ctr_minus_obs_single_level, exp_minus_obs_single_level):
    """Add all subplots to a main plot with positions of pre-set subplots.

    Parameters
    ----------
    exp_single_level : iris cube
        Extracted single level of experiment cube.
    exp_minus_ctr_single_level : iris cube
        Extracted single level of exp_minus_ctr cube.
    ctr_minus_obs_single_level : iris cube
        Extracted single level of ctr_minus_obs cube.
    exp_minus_obs_single_level : iris cube
        Extracted single level of exp_minus_obs cube.
    level : str
        Set depth that we wish to plot.

    Returns
    -------
    quadmap :
        Make the four pane model vs model vs obs comparison plot
    """
    # Setting zrange1 that is used for axis1
    zrange1 = diagtools.get_cube_range([exp_single_level])

    # Setting zrange dependent on the plot produced.
    if exp_single_level.long_name in [
            'Sea Water Salinity', 'Sea Surface Salinity'
    ]:
        # zrange1 set for average salinity range
        # zrange2 set for salinity at -2,2 due to a smaller range of values.
        zrange1 = [20.0, 40.0]
        zrange2 = [-2.0, 2.0]
    else:
        # zrange1 set for average temperature range
        # zrange2 for all other plots. Set for consistency
        zrange1 = [-2.0, 32]
        zrange2 = [-8.86, 8.86]

    # Generates 12 evenly spaced values between zrange1[0] and zrange1[1]
    linspace1 = np.linspace(zrange1[0], zrange1[1], 12, endpoint=True)
    # Generates 12 evenly spaced values between zrange2[0] and zrange2[1]
    linspace2 = np.linspace(zrange2[0], zrange2[1], 12, endpoint=True)

    # Prepare image and figure with a 2x2 grid of subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9, 6))

    # Set the figure title for 'Sea Surface' plots
    if exp_single_level.long_name in [
            'Sea Surface Temperature', 'Sea Surface Salinity'
    ]:
        fig.suptitle("Annual Mean:" + exp_single_level.long_name)
        formatted_depth = 0
    # Set the figure for depth plots to include the depth value.
    else:
        depth = exp_single_level.coords("depth")[0].points[0]
        formatted_depth = str(f"{depth:.3f}")
        fig.suptitle("Annual Mean:" + exp_single_level.long_name + ' at ' +
                     formatted_depth + 'm',
                     fontsize=14)

    # Calling the plot_global_single_level plot with set parameters
    plot_global_single_level(ax1, exp_single_level, linspace1,
                             exp_single_level.attributes['model_id'])
    plot_global_single_level(ax2, exp_minus_ctr_single_level, linspace2,
                             exp_minus_ctr_single_level.attributes['model_id'])
    plot_global_single_level(ax3, ctr_minus_obs_single_level, linspace2,
                             ctr_minus_obs_single_level.attributes['model_id'])
    plot_global_single_level(ax4, exp_minus_obs_single_level, linspace2,
                             exp_minus_obs_single_level.attributes['model_id'])

    # Prepare to save the figure
    fn_list = [exp_single_level.long_name, str(formatted_depth)]
    input_files = diagtools.get_input_files(config)
    image_extention = diagtools.get_image_format(config)

    # Construct the file path for saving the plot
    path = diagtools.folder(
        config['plot_dir']) + '_'.join(fn_list) + str(formatted_depth)
    path = path.replace(' ', '') + image_extention
    logger.info('Saving plots to %s', path)

    # Prepare provenance record for the plot
    provenance_record = diagtools.prepare_provenance_record(
        config,
        caption=f"Quadmap models comparison against observation level="
        f"{formatted_depth})",
        statistics=[
            'mean',
            'diff',
        ],
        domain=['global'],
        plot_type=['map'],
        ancestors=list(input_files.keys()),
    )

    # Save the figure and close
    save_figure('_'.join(fn_list), provenance_record, config, fig, close=True)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
