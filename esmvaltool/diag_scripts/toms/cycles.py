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
import logging

# operating system manipulations (e.g. path constructions)
import os

# to manipulate iris cubes
import iris
import matplotlib.pyplot as plt
import iris.quickplot as qplt

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic

logger = logging.getLogger(os.path.basename(__file__))


def _plot_cycle(cfg, data_dict):
    logger.info('Making the plot')
    local_path = cfg['plot_dir']

    # do some plotting
    # loop over the dictionary
    # for each dataset plot the first variable found
    plt.figure()
    for key, value in data_dict.items():
        if value[0].name() == 'precipitation_flux':
            logging.info(f'Converting precip units for {key} to mm day-1')
            value[0].data = value[0].data * 86400
            value[0].units = 'mm day-1'
        qplt.plot(value[0], label=key)

    plt.legend()
    png_name = 'cycle_plot.png'
    plt.savefig(os.path.join(local_path, png_name))
    plt.close()

    return 'Finished'


def process_data(cfg):
    """
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

    # create empty dict to hold the data we will plot
    cycles_dat = {}

    # iterate over key(dataset) and values(list of vars)
    for key, value in my_files_dict.items():
        # load the data that should exist for for each key (dataset)
        # 'value' is a list containing the variables for each key (dataset)
        cycles_dat[key] = []
        # load the data
        for v in value:
            cycles_dat[key].append(iris.load_cube(v['filename']))

    # now make the plot
    _plot_cycle(cfg, cycles_dat)

    # that's it, we're done!
    return 'I am done with my first ESMValTool diagnostic!'


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        process_data(config)
