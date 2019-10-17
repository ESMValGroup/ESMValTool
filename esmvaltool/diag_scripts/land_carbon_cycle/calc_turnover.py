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

import numpy as np
# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import get_diagnostic_filename, group_metadata, select_metadata, run_diagnostic, extract_variables
from esmvalcore.preprocessor import area_statistics

# esmvaltool.diag_scripts.shared.extract_variables

def _plot_map(cfg, cube, dataset):
    """
    makes the maps of variables 

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
    plt.imshow(cube.data, label=dataset)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Total Carbon Stock')
    plt.tight_layout()
    plt.grid()
    plt.colorbar()
    png_name = 'map_' + dataset + '.png'
    plt.savefig(os.path.join(local_path, png_name))
    plt.close()

    # no need to brag :)
    return 'I made some plots!'


def load_variable(metadata, var_name):
    candidates = select_metadata(metadata, short_name=var_name)
    assert len(candidates) == 1
    filename = candidates[0]['filename']
    cube = iris.load_cube(filename)
    return cube
            


def get_cTau(cfg):
    """
    A diagnostic function to calculate the total carbon stock from the list of variables
    passed to the diagnostic script.

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
    varNames=list(extract_variables(cfg).keys())

    my_files_dict = group_metadata(cfg['input_data'].values(), 'dataset')
    # import pdb; pdb.set_trace()

    # iterate over key(dataset) and values(list of vars)
    for key, value in my_files_dict.items():
        # load the cube from data files only
        # using a single variable here so just grab the first (and only)
        # list element

        c_total = load_variable(value, 'ctotal')
        gpp = load_variable(value, 'gpp')

        c_tau = (c_total/gpp)
        c_tau.var_name = 'cTau'
        c_tau.standard_name = None 
        c_tau.long_name = 'Ecosystem Turnover time of carbon'
        c_tau.units = 'years'
        ofilename = get_diagnostic_filename(key+'_cTau', cfg)
        iris.save(c_tau, ofilename)
        # print(vars()[_val['short_name']])
        # print(varNames)
        # print('----')
        # cube = iris.load_cube(value[0]['filename'])

        # the first data analysis bit: simple cube difference:
        # perform a difference between ground and first levels
        # diff_cube = cube[:, :, :] #- cube[:, 1, :, :]
        # square the difference'd cube just for fun
        # squared_cube = diff_cube ** 2.

        # the second data analysis bit (slightly more advanced):
        # compute an area average over the squared cube
        # to apply the area average use a preprocessor function
        # rather than writing your own function
        # area_avg_cube = area_statistics(squared_cube, 'mean')

        # finalize your analysis by plotting a time series of the
        # diffed, squared and area averaged cube; call the plot function:
        _plot_map(cfg, c_tau, key)

    # that's it, we're done!
    return 'I am done with my first ESMValTool diagnostic!'


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        get_cTau(config)
