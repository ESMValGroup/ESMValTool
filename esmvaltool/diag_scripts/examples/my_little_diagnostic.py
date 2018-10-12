"""
Look at this module for guidance how to write your own.

Module for personal diagnostics (example).
Internal imports from exmvaltool work e.g.:

from esmvaltool.preprocessor._regrid import regrid
from esmvaltool.diag_scripts.shared._supermeans import get_supermean

Pipe output through logger;
"""
import os

import matplotlib
# use this everytime you import matplotlib
# modules; some machines dont have graphical interface (X)
matplotlib.use('Agg')  # noqa

import iris
import matplotlib.pyplot as plt

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


def plot_time_series(cfg):
    """
    Example of personal diagnostic function.

    Arguments:
        run - dictionary of data files

    Returns:
        string; makes some time-series plots

    """
    # local path for e.g. plots: user input
    root_dir = '/group_workspaces/jasmin2/cmip6_prep/'  # edit as per need
    out_path = 'esmvaltool_users/valeriu/'   # edit as per need
    local_path = os.path.join(root_dir, out_path)

    # get the files (simple case, one-variable only)
    my_files_dict = _get_my_files(cfg)

    # iterate through preprocessed model data
    for key, value in my_files_dict.items():
        cube = iris.load_cube(value['file'])
        area_avg_cube = area_average(cube, 'latitude', 'longitude')
        plt.plot(area_avg_cube.data[:, 0], label=key)
        plt.xlabel('Time (months)')
        plt.ylabel(cube.standard_name)
        plt.title('Time series at ground level')
        plt.tight_layout()
        plt.grid()
        plt.legend()
        png_name = 'Time_series_' + key + '.png'
        plt.savefig(os.path.join(local_path, png_name))
        plt.close()

    return 'I made some plots!'


if __name__ == '__main__':

    with run_diagnostic() as config:
        plot_time_series(config)
