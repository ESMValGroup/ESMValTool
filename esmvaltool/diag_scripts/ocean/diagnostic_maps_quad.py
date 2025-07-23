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
Author: Lee de Mora (PML): ledm@pml.ac.uk

Refactored script: (Alphabetical order)
Author: Sophie Hall (Met Office): sophie.hall@metoffice.gov.uk
Author: Emma Hogan (Met Office): emma.hogan@metoffice.gov.uk
Author: Dave Storkey (Met Office): dave.storkey@metoffice.gov.uk
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
from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    save_data,
    save_figure,
)

# Create a logger object.
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def main(config):
    """Load the config file and running through the order of operation.

    Parameters
    ----------
    config : dictionary
        configuration dictionary that contains all the necessary information
        for the function to run. It includes details about the models,
        observational datasets, file paths, and other settings.
    """
    # Call the load_data function
    experiment, control, observation = load_data(config)

    # Call the create_plotting_data function
    (exp, exp_minus_ctr, ctr_minus_obs, exp_minus_obs) = create_plotting_data(
        control, experiment, observation
    )

    # If true, set the lists to contain only one level
    if experiment.long_name in [
        "Sea Surface Temperature",
        "Sea Surface Salinity",
    ]:
        exp_list = [exp]
        exp_minus_ctr_list = [exp_minus_ctr]
        ctr_minus_obs_list = [ctr_minus_obs]
        exp_minus_obs_list = [exp_minus_obs]
    else:
        # Set the lists to contain slices of data over different depth levels
        exp_list = exp.slices_over("depth")
        exp_minus_ctr_list = exp_minus_ctr.slices_over("depth")
        ctr_minus_obs_list = ctr_minus_obs.slices_over("depth")
        exp_minus_obs_list = exp_minus_obs.slices_over("depth")

    # Combine the lists into a single list of tuples for easier iteration
    cube_list = list(
        zip(
            exp_list,
            exp_minus_ctr_list,
            ctr_minus_obs_list,
            exp_minus_obs_list,
            strict=True,
        )
    )

    # Iterate through each depth level
    for (
        exp_single_level,
        exp_minus_ctr_single_level,
        ctr_minus_obs_single_level,
        exp_minus_obs_single_level,
    ) in cube_list:
        # Call create_quadmap, which contains plot_global_single_level
        create_quadmap(
            exp_single_level,
            exp_minus_ctr_single_level,
            ctr_minus_obs_single_level,
            exp_minus_obs_single_level,
            config,
        )

    # After successfully generating plots, function logs a success message.
    logger.info("Success")


def load_data(config):
    """Load all necessary data to output experiment, control, observation.

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
    exp_label_from_recipe = "exper_model"
    ctl_label_from_recipe = "control_model"
    obs_label_from_recipe = "observational_dataset"

    # Getting all input files from configuration.
    input_files = diagtools.get_input_files(config)

    # Get experiment input files from config
    exp_filename = diagtools.match_model_to_key(
        exp_label_from_recipe, config[exp_label_from_recipe], input_files
    )
    # Get control input files from config
    ctl_filename = diagtools.match_model_to_key(
        ctl_label_from_recipe, config[ctl_label_from_recipe], input_files
    )
    # Get observation input files from config
    obs_filename = diagtools.match_model_to_key(
        obs_label_from_recipe, config[obs_label_from_recipe], input_files
    )

    # Set variable names to filename cubes above.
    experiment = iris.load_cube(exp_filename)
    control = iris.load_cube(ctl_filename)
    observation = iris.load_cube(obs_filename)

    # Fixing all units
    experiment = diagtools.bgc_units(
        experiment, input_files[exp_filename]["short_name"]
    )
    control = diagtools.bgc_units(
        control, input_files[ctl_filename]["short_name"]
    )
    observation = diagtools.bgc_units(
        observation, input_files[obs_filename]["short_name"]
    )

    return experiment, control, observation


def create_plotting_data(control, experiment, observation):
    """Calculate the diff between the ctr, exp & obs to prepare data for plots.

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

    # Fixing exp_minus_ctr source_id for the plot title used later
    exp_minus_ctr.attributes["source_id"] = (
        f"{experiment.attributes['source_id']} minus "
        f"{control.attributes['source_id']}"
    )

    # Fixing ctr_minus_obs source_id for the plot title used later
    ctr_minus_obs.attributes["source_id"] = (
        f"{control.attributes['source_id']} minus "
        f"{observation.attributes['short_name']}"
    )

    # Fixing exp_minus_obs source_id for the plot title used later
    exp_minus_obs.attributes["source_id"] = (
        f"{experiment.attributes['source_id']} minus "
        f"{observation.attributes['short_name']}"
    )

    return (exp, exp_minus_ctr, ctr_minus_obs, exp_minus_obs)


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
    # Setting the colour of axis1 always to viridis.
    if title == "UKESM1-0-LL":
        cmap = "viridis"
    # Setting the colour of all other plots to bwr.
    else:
        cmap = "bwr"

    # This step transforms the data so it can be displayed as 2D

    new_cube, _ = iris.analysis.cartography.project(
        cube, ccrs.PlateCarree(), nx=400, ny=200
    )

    # Sets the current Axes instance to the specified axis
    plt.sca(axis)

    # Converts contour_levels to a numpy array.
    contour_levels = np.array(contour_levels)
    # Creates a filled contour plot of the projected data.
    contour_result = iplt.contourf(
        new_cube,
        levels=contour_levels,
        linewidth=0,
        cmap=plt.cm.get_cmap(cmap),
    )
    # Checks if the contour plot was created successfully
    if contour_result is None:
        raise ValueError("Failed to create contour plot. plt object is None.")

    # A color bar is added to the plot
    colorbar = plt.colorbar(contour_result, orientation="horizontal")

    # Colour scale dependent on range of data.
    colorbar.set_ticks(
        [
            contour_levels.min(),
            (contour_levels.max() + contour_levels.min()) / 2.0,
            contour_levels.max(),
        ]
    )

    # Adding latitude & longitude axis on each plot.
    fontsize = 7
    axis = plt.gca()
    grid_lines = axis.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
    grid_lines.xlines = False
    grid_lines.ylines = False
    grid_lines.top_labels = False
    grid_lines.right_labels = False
    grid_lines.xlabel_style = {"size": fontsize, "color": "gray"}
    grid_lines.ylabel_style = {"size": fontsize, "color": "gray"}
    grid_lines.xformatter = LONGITUDE_FORMATTER
    grid_lines.yformatter = LATITUDE_FORMATTER

    # Coastlines are added to the map to provide geographical context.
    plt.gca().coastlines()

    # Sets the title of the plot
    plt.title(title)

    # Display the plots
    iplt.show()


def save_cube(cube, field_name, config, ancestors):
    """
    Produces a provenance record and saves data for each cube.

    Parameters
    ----------
    cube : iris cube
        This is a data structure that contains the climate data to be plotted.
        Including information like temperature values, latitude, longitude,
        and depth.
    field_name : str
        A string that contains the cube name with the corresponding extracted
        depth level.
    config : dictionary
        configuration dictionary that contains all the necessary information
        for the function to run. It includes details about the models,
        observational datasets, file paths, and other settings.
    ancestors : list
        A list of keys from the input_files dictionary, representing the
        provenance of the data. This list helps track the origin and
        transformation history of the data used in the cube

    """
    # Prepare provenance record for the plot
    provenance_record = diagtools.prepare_provenance_record(
        config,
        caption=field_name,
        statistics=["mean", "diff"],
        domain=["global"],
        plot_type=["map"],
        ancestors=ancestors,
    )
    save_data(
        field_name,
        provenance_record,
        config,
        cube,
    )


def create_quadmap(
    exp_single_level,
    exp_minus_ctr_single_level,
    ctr_minus_obs_single_level,
    exp_minus_obs_single_level,
    config,
):
    """
    Add all subplots to a main plot, positions of subplots are pre-set.

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
    config : dictionary
        configuration dictionary that contains all the necessary information
        for the function to run. It includes details about the models,
        observational datasets, file paths, and other settings.

    Returns
    -------
    quadmap :
        Make the four pane model vs model vs obs comparison plot
    """
    # Setting zrange dependent on the plot produced.
    if exp_single_level.long_name in [
        "Sea Water Salinity",
        "Sea Surface Salinity",
    ]:
        # zrange1 set for average salinity range
        # zrange2 set for salinity at -2,2 due to a smaller range of values.
        zrange1 = [20.0, 40.0]
        zrange2 = [-2.0, 2.0]
    else:
        # zrange1 set for average temperature range
        # zrange2 for all other plots. Set for consistency
        zrange1 = [-2.0, 32]
        zrange2 = [-5.0, 5.0]

    # Generates 12 evenly spaced values between zrange1[0] and zrange1[1]
    linspace1 = np.linspace(zrange1[0], zrange1[1], 12, endpoint=True)
    # Generates 12 evenly spaced values between zrange2[0] and zrange2[1]
    linspace2 = np.linspace(zrange2[0], zrange2[1], 12, endpoint=True)

    # Prepare image and figure with a 2x2 grid of subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9, 6))

    # Set the figure title for 'Sea Surface' plots
    if exp_single_level.long_name in [
        "Sea Surface Temperature",
        "Sea Surface Salinity",
    ]:
        fig.suptitle(f"Annual Mean: {exp_single_level.long_name}")
        formatted_depth = 0
    # Set the figure for depth plots to include the depth value.
    else:
        depth = exp_single_level.coords("depth")[0].points[0]
        # Making the depth a string and 3pd
        formatted_depth = str(f"{int(depth):04d}")
        depth_title = str(f"{depth:.1f}")
        fig.suptitle(
            f"Annual Mean: {exp_single_level.long_name} at {depth_title}m",
            fontsize=14,
        )

    # Calling the plot_global_single_level plot with set parameters
    plot_global_single_level(
        ax1,
        exp_single_level,
        linspace1,
        exp_single_level.attributes["source_id"],
    )
    plot_global_single_level(
        ax2,
        exp_minus_ctr_single_level,
        linspace2,
        exp_minus_ctr_single_level.attributes["source_id"],
    )
    plot_global_single_level(
        ax3,
        ctr_minus_obs_single_level,
        linspace2,
        ctr_minus_obs_single_level.attributes["source_id"],
    )
    plot_global_single_level(
        ax4,
        exp_minus_obs_single_level,
        linspace2,
        exp_minus_obs_single_level.attributes["source_id"],
    )

    input_files = diagtools.get_input_files(config)
    ancestors = list(input_files.keys())
    # Calling save_cube for each cube.
    save_cube(
        exp_single_level, f"experiment_{formatted_depth}", config, ancestors
    )
    save_cube(
        exp_minus_ctr_single_level,
        f"experiment_minus_control_{formatted_depth}",
        config,
        ancestors,
    )
    save_cube(
        ctr_minus_obs_single_level,
        f"control_minus_observation_{formatted_depth}",
        config,
        ancestors,
    )
    save_cube(
        exp_minus_obs_single_level,
        f"experiment_minus_observation_{formatted_depth}",
        config,
        ancestors,
    )

    # Prepare to save the figure
    fn_list = [exp_single_level.long_name, str(formatted_depth)]
    image_extention = diagtools.get_image_format(config)

    # Construct the file path for saving the plot
    path = (
        f"{diagtools.folder(config['plot_dir'])}"
        f"{'_'.join(fn_list)}"
        f"{formatted_depth}"
    )
    path = f"{path.replace(' ', '')}{image_extention}"
    logger.info("Saving plots to %s", path)

    # Prepare provenance record for the plot
    provenance_record = diagtools.prepare_provenance_record(
        config,
        caption=f"Quadmap models comparison against observation level="
        f"{formatted_depth})",
        statistics=["mean", "diff"],
        domain=["global"],
        plot_type=["map"],
        ancestors=ancestors,
    )

    # Save the figure and close
    save_figure("_".join(fn_list), provenance_record, config, fig, close=True)


if __name__ == "__main__":
    with run_diagnostic() as CONFIG:
        main(CONFIG)
