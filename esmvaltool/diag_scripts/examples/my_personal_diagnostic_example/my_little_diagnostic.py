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
matplotlib.use('Agg')  # noqa

import iris
import matplotlib.pyplot as plt

from esmvaltool.preprocessor._area_pp import area_average


def plot_time_series(my_files_dict):
    """
    Example of personal diagnostic function.

    Arguments:
        run - dictionary of data files

    Returns:
        string; makes some time-series plots

    """
    # local path for e.g. plots
    root_dir = '/group_workspaces/jasmin2/cmip6_prep/'
    out_path = 'esmvaltool_users/valeriu/'
    local_path = os.path.join(root_dir, out_path)

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
