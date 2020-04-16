import os
import iris
import numpy as np

from esmvaltool.diag_scripts.shared import select_metadata


def _apply_common_mask(*args):
    '''
    apply common mask to all arrays passed as arguments

    Arguments:
        arrays

    Returns:
        an array with size of nargs x common size of all input arrays
    '''
    nargs = len(args)
    for ar in range(nargs):
        _dat = args[ar]
        vars()['datMask' + str(ar)] = np.ones((np.shape(_dat)))
        vars()['datMaskInv' + str(ar)] = np.ma.masked_invalid(_dat).mask
        vars()['datMask' + str(ar)][vars()['datMaskInv' + str(ar)]] = 0
    datMask = vars()['datMask0']
    for ar in range(nargs):
        datMask = datMask * vars()['datMask' + str(ar)]
    mask_where = np.ma.getmask(np.ma.masked_less(datMask, 1.))
    odat = []
    for ar in range(nargs):
        _dat = args[ar].astype(np.float)
        _dat[mask_where] = np.nan
        odat = np.append(odat, np.ma.masked_invalid(_dat))
    odat = odat.reshape(nargs, _dat.shape[0], _dat.shape[1])
    return odat


def _apply_gpp_threshold(gpp_dat, fig_config):
    '''
    returns the input array with values below the threshold of gpp set to nan
    '''
    # converting gC m-2 yr-1 to kgC m-2 s-1
    gpp_thres = fig_config["gpp_threshold"] / (86400.0 * 365 * 1000.)
    gpp_dat = np.ma.masked_less(gpp_dat,
                                gpp_thres).filled(fig_config["fill_value"])
    return gpp_dat


def _get_obs_data_zonal(diag_config):
    '''
    Get and handle the observations of turnover time from Carvalhais 2014.

    Arguments:
        diag_config - nested dictionary of metadata

    Returns:
        dictionary with observation data with different variables as keys
    '''
    if not diag_config.get('obs_variable'):
        raise ValueError('The observation variable needs to be specified in '
                         'the recipe (see recipe description for details)')
    else:
        obs_dir = os.path.join(diag_config['auxiliary_data_dir'],
                               diag_config['obs_info']['obs_data_subdir'])

    all_data = {}
    var_list = diag_config.get('obs_variable')

    input_files = []
    for _var in var_list:
        var_list = np.append(var_list, '{var}_{perc:d}'.format(var=_var,
                                                               perc=5))
        var_list = np.append(var_list, '{var}_{perc:d}'.format(var=_var,
                                                               perc=95))
        obs_filename = (f'{_var}_{{frequency}}_{{source_label}}_'
                        f'{{variant_label}}_{{grid_label}}z.nc'.format(
                            **diag_config['obs_info']))
        input_files = np.append(input_files,
                                os.path.join(obs_dir, obs_filename))

    for v_ind in range(len(var_list)):
        var_obs = var_list[v_ind]
        variable_constraint = iris.Constraint(
            cube_func=(lambda c: c.var_name == var_obs))
        cube = iris.load_cube(input_files, constraint=variable_constraint)
        all_data[var_obs] = cube
    for coord in cube.coords():
        all_data[coord.name()] = coord
    return all_data


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


def _remove_invalid(tmp, fill_value=-9999.):
    '''
    removes the invalid non-numeric values from the input array and fills it
    with fill_value. Also removes all large and small values with magnitude
    beyond 1e15
    '''
    tmp = np.ma.masked_outside(tmp, -1e15, 1e15).filled(fill_value)
    whereisNan = np.isnan(tmp)
    tmp[whereisNan] = fill_value
    whereisNan = np.isinf(tmp)
    tmp[whereisNan] = fill_value
    return tmp
