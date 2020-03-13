import iris
import numpy as np

from esmvaltool.diag_scripts.shared import select_metadata


def _apply_gpp_threshold(gpp_dat, fig_config):
    '''
    returns the input array with values below the threshold of gpp set to nan
    '''
    # converting gC m-2 yr-1 to kgC m-2 s-1
    gpp_thres = fig_config["gpp_threshold"] / (86400.0 * 365 * 1000.)
    gpp_dat = np.ma.masked_less(gpp_dat,
                                gpp_thres).filled(fig_config["fill_value"])
    return gpp_dat


def _load_variable(metadata, var_name):
    '''
    load the data for the variable based on the metadata of the diagnostic
    variable

    Arguments:
        metadata - nested dictionary of metadata

    Returns:
        iris cube of the data
    '''
    candidates = select_metadata(metadata, short_name=var_name)
    assert len(candidates) == 1
    filename = candidates[0]['filename']
    cube = iris.load_cube(filename)
    return cube
