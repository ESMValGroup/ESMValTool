"""Personal diagnostics for Arctic-midlatitude research."""

# import libraries
import iris
# import numpy as np
# import seaborn as sns

from esmvalcore.preprocessor import (
    anomalies,
    area_statistics,
    meridional_statistics,
    zonal_statistics,
)

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
)
from esmvaltool.diag_scripts.shared._base import get_diagnostic_filename


def calculate_polar_vortex(dict_item):
    """Calculate polar vortex."""
    var = iris.load_cube(dict_item['filename'])
    var = var.collapsed('air_pressure', iris.analysis.MEAN)
    # change the sign of polar vortex so the positive values
    # (negative geopotential height anomalies) stand for
    # the strong polar vortex, similarly to
    # Kretschmer et al., 2016 and Galytska et al., 2023
    var.data *= -1
    var.var_name = 'PV'
    return var


def calculate_arctic_tas(dict_item):
    """Read Arctic temperature data."""
    var = iris.load_cube(dict_item['filename'])
    var.var_name = 'Arctic_temperature'
    return var


def calculate_slp(dict_item):
    """Get surface pressure."""
    var = iris.load_cube(dict_item['filename'])
    # calculate hPa from Pa.
    var.data /= 100
    return var


def finalize_bk_ice(dict_item):
    """Read sea ice data (Barents-Kara seas)."""
    var = iris.load_cube(dict_item['filename'])
    var.var_name = 'BK_sic'
    return var


def finalize_ok_ice(dict_item):
    """Read sea ice data (Sea of Okhotsk)."""
    var = iris.load_cube(dict_item['filename'])
    var.var_name = 'Ok_sic'
    return var


def prepare_heat_flux(dict_item):
    """Prepare variables for the heat flux calculations."""
    var = iris.load_cube(dict_item['filename'])
    var_avg = area_statistics(var, operator='mean')
    var_mermean = meridional_statistics(var, operator='mean')
    deviation = var_mermean - var_avg
    return deviation


def calculate_heat_flux(list_va_ta):
    """Calculate eddy poleward heat flux."""
    heat_flux = list_va_ta[0] * list_va_ta[1]
    hf_anom = anomalies(heat_flux, period='monthly')
    hf_anom_zm = zonal_statistics(hf_anom, operator='mean')
    hf_anom_zm.var_name = 'heat_flux'
    return hf_anom_zm


def variable_cases(var, item):
    if var == 'pv':
        out_var  = calculate_polar_vortex(item)
    if var == 'pre_tas':
        out_var = calculate_arctic_tas(item)
    elif var =='pressure_ural':
        out_var = calculate_slp(item)
        out_var.var_name = 'Psl_Ural'
    elif var == 'pressure_sib':
        out_var = calculate_slp(item)
        out_var.var_name = 'Psl_Sib'
    elif var == 'pressure_aleut':
        out_var = calculate_slp(item)
        out_var.var_name = 'Psl_Aleut'
    elif var == 'bk_ice':
        out_var = finalize_bk_ice(item)
    elif var == 'ok_ice':
        out_var = finalize_ok_ice(item)
    elif var == 'heat_flux':
        out_var = prepare_heat_flux(item)
    else:
        raise NotImplementedError(f"Variable '{var}' not supported")
    return out_var

    

def calculate_variables(input_dict):
    """Calculate all necessary variables."""
    dictionary = {}
    for key, value in input_dict.items():
        dictionary.setdefault(key, {})
        tmp_list = []
        for item in value:
            if item['preprocessor'] == "heat_flux":
                tmp_list.append(variable_cases(item['preprocessor'], item))
            else:
                dictionary[key].setdefault(
                    variable_cases(item['preprocessor'], item).var_name,
                    variable_cases(item['preprocessor'], item)
                )

        if key != "HadISST":
            # calculate heat flux for all data sources except HadISST
            heat_flux = calculate_heat_flux(tmp_list)
            dictionary[key].setdefault(heat_flux.var_name, heat_flux)
    return dictionary

# def plot_selected_timeseries(dictionary):
#     """Plot timeseries of indicated variables."""
#     var_names = ['heat_flux', 'PV']
#     for var in var_names:
#         plt.figure(figsize=(14, 4))
#         sns.set_theme()
#         for key in dictionary:
#             time_orig = dictionary[key][var].coord('time')
#             times = np.asarray(time_orig.units.num2date(time_orig.points))
#             time_pts = [t.strftime('%Y-%m') for t in times]
#             plt.plot(time_pts, dictionary[key][var].data)
#             plt.title(var)
#             plt.ylabel('Anomalies')
#             plt.xticks(rotation=45, ha="right", rotation_mode='anchor')


def main(cfg):
    """Calculate and save final variables into .nc files."""
    my_files_dict = group_metadata(cfg['input_data'].values(), 'dataset')
    all_variables = calculate_variables(my_files_dict)
    for key in my_files_dict:
        diagnostic_file = get_diagnostic_filename(key, cfg)
        var = all_variables[key]
        if key == "ERA5":
            cube_list = iris.cube.CubeList([
                var['PV'], var['Arctic_temperature'], var['Psl_Ural'],
                var['Psl_Sib'], var['Psl_Aleut'], var['heat_flux']])
        elif key == "HadISST":
            cube_list = iris.cube.CubeList([
                var['BK_sic'], var['Ok_sic']])
        else:
            cube_list = iris.cube.CubeList([
                var['PV'], var['Arctic_temperature'], var['Psl_Ural'],
                var['Psl_Sib'], var['Psl_Aleut'], var['heat_flux'],
                var['BK_sic'], var['Ok_sic']])
        iris.save(cube_list, diagnostic_file)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
