import logging
import sys

import iris

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic


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

    Notes
    -----
    After successfully generating the plots,
    the function logs a success message.

    The if statement checks if the script is being run directly
    (as opposed to being imported as a module).
    If it is, the code inside this block will be executed.

    The function starts by loading the configuration dictionary.
    This dictionary contains information about the models,
    observational datasets, file paths, and other settings.
    """

    create_logger()

    control, experiment, observation = load_data(config)

    fix_cube(cube_to_fix=control, cube_with_good_values=observation)
    fix_cube(cube_to_fix=experiment, cube_with_good_values=observation)

    remove_extra_time_axis(control)
    remove_extra_time_axis(experiment)
    remove_extra_time_axis(observation)

    (experiment, experiment_minus_control, control_minus_observation,
     experiment_minus_observation) = create_plotting_data(
         control, experiment, observation)

    #    level = 2
    #    single_level_experiment = extract_global_single_level(
    #    experiment, level)
    #    single_level_control = extract_global_single_level(
    #        experiment_minus_control, level)
    #    single_level_observation = extract_global_single_level(
    #        control_minus_observation, level)
    #    single_level_observation = extract_global_single_level(
    #        experiment_minus_observation, level)

    level = 2
    extract_global_single_level(experiment, level)
    extract_global_single_level(experiment_minus_control, level)
    extract_global_single_level(control_minus_observation, level)
    extract_global_single_level(experiment_minus_observation, level)


#    experiment_plot = plot_global_single_level(
#    experiment)
#    experiment_minus_control_plot = plot_global_single_level(
#    experiment_minus_control)
#    control_minus_observation_plot = plot_global_single_level(
#    control_minus_observation)
#    experiment_minus_observation_plot = plot_global_single_level(
#    experiment_minus_observation)

#    create_quadmap(
#    experiment_plot, experiment_minus_control_plot,
#    control_minus_observation_plot, experiment_minus_observation_plot)


def create_logger():
    """The script configures the logging system to create a logger named after
    the current file and to print log messages to the console.

    This helps in monitoring the script's progress and debugging any
    issues.
    """
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


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
    The function starts by loading the configuration dictionary.
    This dictionary contains information about the
    models, observational datasets, file paths, and other settings.

    The input data is grouped by dataset. This helps in organising the data for
    easier processing and plotting.

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

    input_files = diagtools.get_input_files(config)

    ctl_filename = diagtools.match_model_to_key(ctl_label_from_recipe, config,
                                                input_files)

    exp_filename = diagtools.match_model_to_key(exp_label_from_recipe, config,
                                                input_files)

    obs_filename = diagtools.match_model_to_key(obs_label_from_recipe, config,
                                                input_files)

    control = iris.load_cube(ctl_filename)
    experiment = iris.load_cube(exp_filename)
    observation = iris.load_cube(obs_filename)

    print(f'{control=}')
    print(f'{experiment=}')
    print(f'{observation=}')

    return control, experiment, observation


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

    Notes
    -----
    The data for each model and the observational dataset is loaded into
    Iris cubes. These cubes contain the climate data that will be plotted.
    """

    experiment = experiment

    experiment_minus_control = experiment - control

    control_minus_observation = control - observation

    experiment_minus_observation = experiment - observation

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

    Notes
    -----
    Pick a value for layer and loop over the cubes.
    """
    if len(cube.coord('depth').points) == 1:
        # 2D cube
        print('aaaaaaaaaaasection 1')
        return iris.util.squeeze(cube)
    else:
        # 3D cube - select relevant layer
        slices = [slice(None)] * len(cube.shape)
        coord_dim = cube.coord_dims('depth')[0]
        slices[coord_dim] = level
        print('aaaaaaaaaaaasection 2')
        return iris.util.squeeze(cube[tuple(slices)])

    print('aaaaaaaacube:', cube)
    print('aaaaaaaaaaalevel:', cube)


def plot_global_single_level(cube, cmap, title, contour_levels):
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
    contour_levels : numpy.array
        nspace is used to set the ticks on the colour bar and used to define
        levels for the contour plot.

    Returns
    -------
    plot of single depth (called four times)

    Notes
    -----
    The climate data is projected onto a 2D map using a specific map projection
    (PlateCarree). This step transforms the data so it can be displayed on a
    flat map.

    The function then creates a filled contour plot of the projected data.
    The array is used to set the color scale, and the input determines the
    color scheme.

    A color bar is added to the plot to show the range of values represented by
    different colors. The ticks on the color bar are set to the minimum,
    midpoint, and maximum values of nspace.

    Coastlines are added to the map to provide geographical context.
    The title of the subplot is set using the title input.

    The overall title of the figure is set based on the data being plotted.

    Titles are added to each plot.

    The plots are then saved as image files in the specified directory.

    Plots are then passed as parameters in create_quadmap().
    """


def create_quadmap(experiment_plot, experiment_minus_control_plot,
                   control_minus_observation_plot,
                   experiment_minus_observation_plot):
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
    Top left: Experiment model
    Top right: Experiment model minus control model
    Bottom left: Control model minus observational dataset
    Bottom right: Experiment model minus observational dataset

    Set the number that tells the function where to place the map within the
    larger figure. 224 means the map will be placed in the fourth position
    of a 2x2 grid. We want this to be unchanged.

    The function matches the models to their respective keys using the
    information in the configuration dictionary.
    """


def fix_cube(cube_to_fix, cube_with_good_values):
    """The latitude and longitude coordinates of the cubes are adjusted to
    ensure they match cube with good values.

    This is necessary for subtracting one cube from another to calculate
    differences.
    """
    cube_to_fix.coord("latitude").points[:] = cube_with_good_values.coord(
        "latitude").points
    cube_to_fix.coord("longitude").points[:] = cube_with_good_values.coord(
        "longitude").points


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


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
