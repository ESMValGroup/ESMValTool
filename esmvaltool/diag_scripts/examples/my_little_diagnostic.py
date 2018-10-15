"""
Look at this module for guidance how to write your own.

Read the README_PERSONAL_DIAGNOSTIC file associated with this example;

Module for personal diagnostics (example).
Internal imports from exmvaltool work e.g.:

from esmvaltool.preprocessor._regrid import regrid
from esmvaltool.diag_scripts.shared._supermeans import get_supermean

Pipe output through logger;

Please consult the documentation for help with esmvaltool's functionalities
and best coding practices.
"""
# place your module imports here:

# operating system manipulations (e.g. path constructions)
import os

# for plotting
import matplotlib
# use this everytime you import matplotlib
# modules; some machines dont have graphical interface (X)
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt

# to manipulate iris cubes
import iris

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.preprocessor._area_pp import area_average


def _get_my_files(cfg):
    """Put files in dicts of datasets and return them."""
    files_dict = {}
    for filename, attributes in cfg['input_data'].items():
        base_file = os.path.basename(filename)
        dataset = base_file.split('_')[1]
        files_dict[dataset] = {}
        files_dict[dataset]['file'] = filename
        if 'fx_files' in attributes:
            for fx_var in attributes['fx_files']:
                files_dict[dataset][fx_var] = attributes['fx_files'][fx_var]

    return files_dict


def _run_diagnostic(cfg):
    """
    Simple example of a diagnostic.

    This is a basic (and rather esotherical) diagnostic that firstly
    loads the needed model data as iris cubes, performs a difference between
    values at ground level and first vertical level, then squares the
    result. The user will implement their own (custom) diagnostics, but this
    example shows that once the preprocessor has finished a whole lot of
    user-specific metrics can be computed as part of the diagnostic.

    """
    # assemble the data dictionary keyed by dataset name
    my_files_dict = _get_my_files(cfg)
    # assemble an empty dict keyed by dataset
    # but that will have cubes as values
    diagnostic_dict = {}
    # iterate over key(dataset) and values(dict of files and fx_files)
    for key, value in my_files_dict.items():
        # load the cube from files only
        cube = iris.load_cube(value['file'])
        # perform a difference between ground and first levels
        diff_cube = cube[:, 0, :, :] - cube[:, 1, :, :]
        # square the difference'd cube and populate the output dictionary
        diagnostic_dict[key] = diff_cube ** 2.

    return diagnostic_dict


def plot_time_series(cfg):
    """
    Example of personal diagnostic plotting function.

    Before plotting, we grab the output of the diagnostic (_run_diagnostic)
    and apply an area average on it. This is a useful example of how to use
    standard esmvaltool-preprocessor functionality within a diagnostic, and
    especially after a certain (custom) diagnostic has been run and the user
    needs to perform an operation that is already part of the preprocessor
    standard library of functions.

    Arguments:
        run - dictionary of data files

    Returns:
        string; makes some time-series plots

    """
    # local path for e.g. plots: user input
    root_dir = '/group_workspaces/jasmin2/cmip6_prep/'  # edit as per need
    out_path = 'esmvaltool_users/valeriu/'   # edit as per need
    local_path = os.path.join(root_dir, out_path)

    # run your diagnostic - get the squared data
    my_data_dict = _run_diagnostic(cfg)

    # iterate through data for plotting
    for key, value in my_data_dict.items():

        # apply an area average (using a preprocessor function
        # rather than writing your own function)
        area_avg_cube = area_average(value, 'latitude', 'longitude')

        # do the plotting dance
        plt.plot(area_avg_cube.data, label=key)
        plt.xlabel('Time (months)')
        plt.ylabel('The squared difference')
        plt.title('Time series at ground level')
        plt.tight_layout()
        plt.grid()
        plt.legend()
        png_name = 'Time_series_' + key + '.png'
        plt.savefig(os.path.join(local_path, png_name))
        plt.close()

    # no need to brag :)
    return 'I made some plots!'


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        plot_time_series(config)
