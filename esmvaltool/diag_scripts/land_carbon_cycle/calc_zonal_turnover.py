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

def _plot_line(cfg, data,dataset):
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
    plt.plot(data)
    plt.xlabel('Latitude')
    plt.ylabel('tau')
    plt.title('zonal turnover')
    plt.tight_layout()
    plt.ylim(0,1000)
    plt.grid()
    # plt.colorbar()
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
            

def get_zonal_tau(_datgpp,_datcTotal):
    __dat=np.zeros((np.shape(_datgpp)[0]))
    minPoints=20
    bandsize=9
    for li in range(len(__dat)):
        istart=max(0,li-bandsize)
        iend=min(359,li+bandsize+1)
        _datgppZone=_datgpp[istart:iend,:]#* _areaZone
        _datcTotalZone=_datcTotal[istart:iend,:]#* _areaZone
        # _datgppZone=remove_tail_percentiles(_datgppZone,_outlierPerc=10)
        # _datcTotalZone=remove_tail_percentiles(_datcTotalZone,_outlierPerc=10)
        # print(_datZone)
        # print(_datgpp)
        nValids=np.nansum(1-np.ma.getmask(np.ma.masked_invalid(_datgppZone)))
        if nValids > minPoints:
            __dat[li] = (np.nansum(_datcTotalZone)/np.nansum(_datgppZone)) / (86400*365)
        else:
            __dat[li] = np.nan
            print(nValids,'valid points for tau, setting nan')
    return(__dat)
def get_zonal_cTau(cfg):
    """
    A diagnostic function to calculate the total carbon stock from the list of variables
    passed to the diagnostic script.

    Arguments:
        cfg - nested dictionary of metadata

    Returns:
        string; runs the user diagnostic

    """

    my_files_dict = group_metadata(cfg['input_data'].values(), 'dataset')
    # import pdb; pdb.set_trace()

    # iterate over key(dataset) and values(list of vars)
    zonal_tau_models={}
    for key, value in my_files_dict.items():
        # load the cube from data files only
        # using a single variable here so just grab the first (and only)
        # list element

        c_total = load_variable(value, 'ctotal')
        gpp = load_variable(value, 'gpp')
        print(gpp,c_total)

        c_tau = get_zonal_tau(gpp.data, c_total.data)
        _plot_line(cfg, c_tau, key)
        print(key,value)
    # that's it, we're done!
    return 'I am done with my first ESMValTool diagnostic!'


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        get_zonal_cTau(config)
