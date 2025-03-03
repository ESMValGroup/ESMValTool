import logging
import os
import sys

import cartopy.crs as ccrs
import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic, save_figure

logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def main(config):
    """Diagnostic to produce an image showing four maps, based on a comparison
    of two different models results against an observational dataset.

    Parameters
    ----------
    config : dictionary
        configuration dictionary that contains
        all the necessary informationfor the function to run.
        It includes details about the
        models, observational datasets, file paths, and other settings.
    """

    # function starts by loading the config files.
    # this contains info about the models, obs, file paths and settings.

    control, experiment, observation = load_data(config)

    fix_cube(cube_to_fix=control, cube_with_good_values=observation)
    fix_cube(cube_to_fix=experiment, cube_with_good_values=observation)

    remove_extra_time_axis(control)
    remove_extra_time_axis(experiment)
    remove_extra_time_axis(observation)

    (experiment, experiment_minus_control, control_minus_observation,
     experiment_minus_observation) = create_plotting_data(
         control, experiment, observation)

    # Pick a value for level.
    levels = [2, 50, 100, 300]
    for level in levels:
        experiment_single_level = extract_global_single_level(
            experiment, level)
        experiment_minus_control_single_level = extract_global_single_level(
            experiment_minus_control, level)
        control_minus_observation_single_level = extract_global_single_level(
            control_minus_observation, level)
        experiment_minus_observation_single_level = (
            extract_global_single_level(experiment_minus_observation, level))

        create_quadmap(experiment_single_level,
                       experiment_minus_control_single_level,
                       control_minus_observation_single_level,
                       experiment_minus_observation_single_level, level)

    # After successfully generating plots, function logs a success message.
    logger.info('Success')


def load_data(config):
    """Loading in all necessary data to output control, experiment, observation
    to be passed into create_plotting_data()

    Parameters
    ----------
    config : dictionary
        configuration dictionary that contains all the necessary information
        for the function to run. It includes details about the
        models, observational datasets, file paths, and other settings.

    Returns
    -------
    control : iris cube
        Data as defined as control_model from the recipe.
    experiment : iris cube
        Data as defined as exper_model from the recipe.
    observation : iris cube
        Data as defined as observational_dataset from the recipe.

    Notes
    -----
    For each file listed in the input_files section of the configuration,
    the function logs the filename. This helps in tracking which files are
    being processed.

    The function retrieves the input files using the
    diagtools.get_input_files function.
    This prepares the necessary data files for further processing.
    """

    # The datasets defined by
    # control_model, exper_model, observational_datasets in the recipe
    # determine what is loaded in this function.
    ctl_label_from_recipe = 'control_model'
    exp_label_from_recipe = 'exper_model'
    obs_label_from_recipe = 'observational_dataset'

    # Getting all input files from configuration.
    input_files = diagtools.get_input_files(config)

    # Get control input files from config
    ctl_filename = diagtools.match_model_to_key(ctl_label_from_recipe,
                                                config[ctl_label_from_recipe],
                                                input_files)
    # Get experiment input files from config
    exp_filename = diagtools.match_model_to_key(exp_label_from_recipe,
                                                config[exp_label_from_recipe],
                                                input_files)
    # Get observation input files from config
    obs_filename = diagtools.match_model_to_key(obs_label_from_recipe,
                                                config[obs_label_from_recipe],
                                                input_files)
    print(f'{ctl_filename=}')
    control = iris.load_cube(ctl_filename)
    print(f'{exp_filename=}')
    experiment = iris.load_cube(exp_filename)
    print(f'{obs_filename=}')
    observation = iris.load_cube(obs_filename)

    print(f'{control=}')
    print(f'{experiment=}')
    print(f'{observation=}')
    print("load_data success")
    return control, experiment, observation


def fix_cube(cube_to_fix, cube_with_good_values):
    """The latitude and longitude coordinates of the cubes are adjusted to
    ensure they match cube with good values.

    Parameters
    ----------
    cube_to_fix :
        The cubes experiment and control with bad coord values.
    cube_with_good_values :
        The cube observation with the good coord values.
    """
    # This is necessary for subtracting one cube from another.
    cube_to_fix.coord("latitude").points[:] = cube_with_good_values.coord(
        "latitude").points
    cube_to_fix.coord("longitude").points[:] = cube_with_good_values.coord(
        "longitude").points
    print("fix_cube success")


def remove_extra_time_axis(cube):
    """Remove the extra time axis from the |input variable| provided by the
    ``cube`` parameter.

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
    print("remove_extra_time success")


def create_plotting_data(control, experiment, observation):
    """Calculating the differences between the control, experiment and
    observation datasets to prepare data for plotting.

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
    experiment : iris cube
        Untouched experimental input.
    experiment_minus_control : iris cube
        Experiment model minus control model.
    control_minus_observation : iris cube
        Control model minus observational dataset.
    experiment_minus_observation : iris cube
        Experimental model minus observational dataset.
    """
    # The data for models and the obs dataset is loaded into Iris cubes.
    # These cubes contain the climate data that will be plotted.
    experiment_minus_control = experiment - control
    control_minus_observation = control - observation
    experiment_minus_observation = experiment - observation

    # Fixing the long_name of the cubes
    experiment_minus_control.long_name = experiment.long_name
    control_minus_observation.long_name = experiment.long_name
    experiment_minus_observation.long_name = experiment.long_name

    return (experiment, experiment_minus_control, control_minus_observation,
            experiment_minus_observation)


def extract_global_single_level(cube, level):
    """Extracts a single level from the cube.

    Parameters
    ----------
    cube : iris cube
        The input data cube.
    level : float
        The depth level to extract.

    Returns
    -------
    iris cube
        The extracted single level cube.
    """
    print(f'{cube.shape=}')
    if len(cube.coord('depth').points) == 1:
        # 2D cube - Sea Surface
        single_level = cube
    else:
        # 3D cube - select relevant level
        constraint = iris.Constraint(depth=level)
        single_level = cube.extract(constraint)

    single_level = iris.util.squeeze(single_level)
    print(f'{single_level.shape=}')
    return single_level


def plot_global_single_level(axis, cube, contour_levels, title):
    """Creating each individual plot before being added to create_quadmap.

    Parameters
    ----------
    cube : iris cube
        This is a data structure that contains the climate data to be plotted.
        Including information like temperature values, latitude, longitude,
        and depth.
    cmap : str
        This is a string that specifies the color map to be used for the plot.
    title : str
        This is a string that will be used as the title of the subplot.
    nspace : numpy.array
        nspace is used to set the ticks on the colour bar and used to define
        levels for the contour plot.

    Returns
    -------
    plot of single depth (called four times)

    Notes
    -----
    The plots are then saved as image files in the specified directory.
    """
    if title == "experiment":
        cmap = 'viridis'
    else:
        cmap = 'bwr'

    # This step transforms the data so it can be displayed as 2D
    new_cube, extent = iris.analysis.cartography.project(cube,
                                                         ccrs.PlateCarree(),
                                                         nx=400,
                                                         ny=200)
    # Set at the top of the function.s
    plt.sca(axis)

    # The function then creates a filled contour plot of the projected data.
    contour_result = iplt.contourf(new_cube,
                                   levels=contour_levels,
                                   linewidth=0,
                                   cmap=plt.cm.get_cmap(cmap))

    contour_levels = np.array(contour_levels)

    if contour_result is None:
        raise ValueError(
            "Failed to create contour plot. The plt object is None.")

    # A color bar is added to the plot to show the range of values.
    colorbar = plt.colorbar(contour_result, orientation='horizontal')
    colorbar.set_ticks([
        contour_levels.min(),
        (contour_levels.max() + contour_levels.min()) / 2.,
        contour_levels.max()
    ])
    # Coastlines are added to the map to provide geographical context.
    plt.gca().coastlines()
    # title plotted
    plt.title(title)
    iplt.show()


#    return qplot


def create_quadmap(experiment_single_level,
                   experiment_minus_control_single_level,
                   control_minus_observation_single_level,
                   experiment_minus_observation_single_level, level):
    """The function starts by specifying the position of each of the plots in
    the quadmap. We want this to be set.

    Parameters
    ----------
    experiment_plot : iris cube
        Untouched experimental input.
    experiment_minus_control_plot : iris cube
        Experiment model minus control model.
    control_minus_observation_plot : iris cube
        Control model minus observational dataset.
    experiment_minus_observation_plot : iris cube
        Experimental model minus observational dataset.

    Returns
    -------
    quadmap :
        Make the four pane model vs model vs obs comparison plot

    Notes
    -----
    Set the number that tells the function where to place the map within the
    larger figure. 224 means the map will be placed in the fourth position
    of a 2x2 grid. We want this to be unchanged.

    The function matches the mo    ProvenanceLogger,
    group_metadata,dels to their respective keys using the
    information in the configuration dictionary.
    """
    # Setting cmap and nspace.

    # fig.title and plt.title

    zrange1 = diagtools.get_cube_range([experiment_single_level])
    if (experiment_single_level.long_name == 'Sea Surface Salinity'
            or experiment_single_level.long_name == 'Sea Water Salinity'):
        zrange2 = [-2.0, 2.0]
    else:
        zrange2 = [-5.0, 5.0]

    linspace1 = np.linspace(zrange1[0], zrange1[1], 12, endpoint=True)
    linspace2 = np.linspace(zrange2[0], zrange2[1], 12, endpoint=True)
    # prepare image and figure
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9, 6))
    level = str(level)
    fig.suptitle(experiment_single_level.long_name + ' at ' + level + 'm',
                 fontsize=14)
    plot_global_single_level(ax1, experiment_single_level, linspace1,
                             "experiment")
    plot_global_single_level(ax2, experiment_minus_control_single_level,
                             linspace2, "experiment minus control")
    plot_global_single_level(ax3, control_minus_observation_single_level,
                             linspace2, "control minus observation")
    plot_global_single_level(ax4, experiment_minus_observation_single_level,
                             linspace2, "experiment minus observation")

    # Saving files:
    fn_list = [experiment_single_level.long_name, str(level)]
    input_files = diagtools.get_input_files(config)
    image_extention = diagtools.get_image_format(config)

    path = diagtools.folder(
        config['plot_dir']) + '_'.join(fn_list) + str(level)
    path = path.replace(' ', '') + image_extention
    logger.info('Saving plots to %s', path)
    provenance_record = diagtools.prepare_provenance_record(
        config,
        caption=f'Quadmap models comparison against observation level={level}',
        statistics=[
            'mean',
            'diff',
        ],
        domain=['global'],
        plot_type=['map'],
        ancestors=list(input_files.keys()),
    )
    save_figure('_'.join(fn_list), provenance_record, config, fig, close=True)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
