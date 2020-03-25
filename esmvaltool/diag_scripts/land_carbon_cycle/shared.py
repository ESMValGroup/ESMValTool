import os
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

def _get_obs_data_zonal(diag_config):
    '''
    Get and handle the observations of turnover time from Carvalhais 2014.

    Arguments:
        diag_config - nested dictionary of metadata

    Returns:
        dictionary with observation data with different variables as keys
    '''
    if not diag_config.get('obs_variable'):
        raise ValueError((
            'The observation variable needs to be specified in the recipe (see recipe description for details)'
        ))
    else:
        obs_dir = os.path.join(diag_config['auxiliary_data_dir'],
                               diag_config['obs_info']['obs_data_subdir'])

    all_data = {}
    var_list = diag_config.get('obs_variable')

    input_files = []
    for _var in var_list:
        var_list = np.append(var_list, '{var}_{perc:d}'.format(var=_var, perc=5))
        var_list = np.append(var_list, '{var}_{perc:d}'.format(var=_var, perc=95))
        obs_filename = '{variable}_{frequency}_{source_label}_{variant_label}_{grid_label}z.nc'.format(
            variable=_var,
            frequency=diag_config['obs_info']['frequency'],
            source_label=diag_config['obs_info']['source_label'],
            variant_label=diag_config['obs_info']['variant_label'],
            grid_label=diag_config['obs_info']['grid_label'])
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