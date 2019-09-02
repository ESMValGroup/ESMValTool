"""
Look at this module for guidance how to write your own.

Read the README_PERSONAL_DIAGNOSTIC file associated with this example;

Module for personal diagnostics (example).
Internal imports from exmvaltool work e.g.:

from esmvalcore.preprocessor import regrid
from esmvaltool.diag_scripts.shared.supermeans import get_supermean

Pipe output through logger;

Please consult the documentation for help with esmvaltool's functionalities
and best coding practices.
"""
# place your module imports here:

# operating system manipulations (e.g. path constructions)
import os

# to manipulate iris cubes
import iris
import matplotlib.pyplot as plt

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvalcore.preprocessor import area_statistics


def _plot_time_series(cfg, cube, dataset):
    """
    Example of personal diagnostic plotting function.

    Arguments:
        cfg - nested dictionary of metadata
        cube - the cube to plot
        dataset - name of the dataset to plot

    Returns:
        string; makes some time-series plots

    Note: this function is private; remove the '_'
    so you can make it public.
    """
    # custom local paths for e.g. plots are supported -
    # here is an example
    # root_dir = '/group_workspaces/jasmin2/cmip6_prep/'  # edit as per need
    # out_path = 'esmvaltool_users/valeriu/'   # edit as per need
    # local_path = os.path.join(root_dir, out_path)
    # but one can use the already defined esmvaltool output paths
    local_path = cfg['plot_dir']

    # do the plotting dance
    plt.plot(cube.data, label=dataset)
    plt.xlabel('Time (months)')
    plt.ylabel('Area average')
    plt.title('Time series at (ground level - first level)')
    plt.tight_layout()
    plt.grid()
    plt.legend()
    png_name = 'Time_series_' + dataset + '.png'
    plt.savefig(os.path.join(local_path, png_name))
    plt.close()

    # no need to brag :)
    return 'I made some plots!'


def run_my_diagnostic(cfg):
    """
    Simple example of a diagnostic.

    This is a basic (and rather esotherical) diagnostic that firstly
    loads the needed model data as iris cubes, performs a difference between
    values at ground level and first vertical level, then squares the
    result.

    Before plotting, we grab the squared result (not all operations on cubes)
    and apply an area average on it. This is a useful example of how to use
    standard esmvalcore.preprocessor functionality within a diagnostic, and
    especially after a certain (custom) diagnostic has been run and the user
    needs to perform an operation that is already part of the preprocessor
    standard library of functions.

    The user will implement their own (custom) diagnostics, but this
    example shows that once the preprocessor has finished a whole lot of
    user-specific metrics can be computed as part of the diagnostic,
    and then plotted in various manners.

    Arguments:
        cfg - nested dictionary of metadata

    Returns:
        string; runs the user diagnostic

    """
    # assemble the data dictionary keyed by dataset name
    # this makes use of the handy group_metadata function that
    # orders the data by 'dataset'; the resulting dictionary is
    # keyed on datasets e.g. dict = {'MPI-ESM-LR': [var1, var2...]}
    # where var1, var2 are dicts holding all needed information per variable
    my_files_dict = group_metadata(cfg['input_data'].values(), 'dataset')

    # iterate over key(dataset) and values(list of vars)
    for key, value in my_files_dict.items():
        # load the cube from data files only
        # using a single variable here so just grab the first (and only)
        # list element
        cube = iris.load_cube(value[0]['filename'])

        # the first data analysis bit: simple cube difference:
        # perform a difference between ground and first levels
        diff_cube = cube[:, 0, :, :] - cube[:, 1, :, :]
        # square the difference'd cube just for fun
        squared_cube = diff_cube ** 2.

        # the second data analysis bit (slightly more advanced):
        # compute an area average over the squared cube
        # to apply the area average use a preprocessor function
        # rather than writing your own function
        area_avg_cube = area_statistics(squared_cube, 'mean')

        # finalize your analysis by plotting a time series of the
        # diffed, squared and area averaged cube; call the plot function:
        _plot_time_series(cfg, area_avg_cube, key)

    # that's it, we're done!
    return 'I am done with my first ESMValTool diagnostic!'


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        run_my_diagnostic(config)
