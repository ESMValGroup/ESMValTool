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

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic

logger = logging.getLogger(os.path.basename(__file__))


def _plot_scatter(cfg, x_data, y_data, labels):
    logger.info('Making the plot')
    local_path = cfg['plot_dir']

    # do some plotting
    plt.scatter(x_data, y_data)

    # now add some labels
    for i, txt in enumerate(labels):
        plt.annotate(txt, (x_data[i], y_data[i]))

    plt.tight_layout()
    png_name = 'scatter_plot.png'
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

    # create empty lists to hold the data we will make a scatter plot of
    sactter_dat = [[], []]
    labels = []

    # get an arbitrary dataset to work out what we're dealing with
    val = list(my_files_dict.values())[0]
    # get the variable names and experiments of each variable
    vars = set([v['short_name'] for v in val])
    exps = set([v['exp'] for v in val])

    logger.info(f'Found the following variables {vars}')
    logger.info(f'Found the following scenarios {exps}')

    # our past experiment should always be 'historical'
    # get the name of the future experiment
    fut_name = exps.copy()
    fut_name.remove('historical')
    fut_name, = fut_name

    # iterate over key(dataset) and values(list of vars)
    for key, value in my_files_dict.items():
        # load the data that should exist for for each key
        # 'value' is a list containing the 4 variables
        # (present and future var1 and var2)
        # for each individual dataset

        # create nested dict to store the info
        data = dict.fromkeys(vars)
        for v in vars:
            data[v] = dict.fromkeys(exps)

        # load the data into our dict
        for v in value:
            data[v['short_name']][v['exp']] = iris.load_cube(v['filename'])

        # compute the future to present day difference
        # the variables should already be single point iris cubes
        for i, v in enumerate(vars):
            delta = data[v][fut_name].data - data[v]['historical'].data
            # convert units if needed
            if v == 'pr':
                delta = delta * 86400
            sactter_dat[i].append(delta)

        labels.append(key)

    # now make the plot
    _plot_scatter(cfg, sactter_dat[1], sactter_dat[0], labels)

    # that's it, we're done!
    return 'I am done with my first ESMValTool diagnostic!'


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        process_data(config)
