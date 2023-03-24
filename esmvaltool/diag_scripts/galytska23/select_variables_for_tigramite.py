"""Personal diagnostics for Arctic-midlatitude research."""

# import libraries
import iris
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import logging

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

logger = logging.getLogger(Path(__file__).stem)

plot_timeseries= False

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


def calculate_slp(dict_item):
    """Get surface pressure."""
    var = iris.load_cube(dict_item['filename'])
    # calculate hPa from Pa.
    var.data /= 100
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


def calculate_variables(input_dict):
    """Calculate all necessary variables."""
    dictionary = {}
    for key, value in input_dict.items():
        dictionary.setdefault(key,{})    
        tmp_list = []
        for item in value:
            if item['preprocessor'] == 'pv':
                polar_vortex = calculate_polar_vortex(item)
            elif item['preprocessor'] == 'pre_tas':
                tas = iris.load_cube(item['filename'])
                tas.var_name = 'Arctic_temperature'
            elif item['preprocessor'] == 'pressure_ural':
                psl_ural = calculate_slp(item)
                psl_ural.var_name = 'Psl_Ural'
            elif item['preprocessor'] == 'pressure_sib':
                psl_sib = calculate_slp(item)
                psl_sib.var_name = 'Psl_Sib'
            elif item['preprocessor'] == 'pressure_aleut':
                psl_aleut = calculate_slp(item)
                psl_aleut.var_name = 'Psl_Aleut'
            elif item['preprocessor'] == 'bk_ice':
                sic_bk = iris.load_cube(item['filename'])
                sic_bk.var_name = 'BK_sic'
            elif item['preprocessor'] == 'ok_ice':
                sic_ok = iris.load_cube(item['filename'])
                sic_ok.var_name = 'Ok_sic'
            elif item['preprocessor'] == 'heat_flux':
                var = prepare_heat_flux(item)
                tmp_list.append(var)
        heat_flux = calculate_heat_flux(tmp_list)
        dictionary[key].setdefault ('PV', polar_vortex)
        dictionary[key].setdefault ('Arctic_temperature', tas)
        dictionary[key].setdefault ('Psl_Ural', psl_ural)
        dictionary[key].setdefault ('Psl_Sib', psl_sib)
        dictionary[key].setdefault ('Psl_Aleut', psl_aleut)
        dictionary[key].setdefault ('BK_sic', sic_bk)
        dictionary[key].setdefault ('Ok_sic', sic_ok)
        dictionary[key].setdefault ('heat_flux', heat_flux)
    return dictionary


def plot_selected_timeseries(dictionary):
    """Plot timeseries of indicated variables."""
    #provide the variables in the list below
    var_names =['heat_flux', 'PV']
    for var in var_names: 
        fig = plt.figure(figsize=(14,4))
        sns.set_theme()
        for key in dictionary: 
            time_orig = dictionary[key][var].coord('time')
            times = np.asarray(time_orig.units.num2date(time_orig.points))
            time_pts=[t.strftime('%Y-%m') for t in times]   
            plt.plot(time_pts, dictionary[key][var].data)
            plt.title (var)
            plt.ylabel ('Anomalies')
            plt.xticks (rotation = 45, ha="right", rotation_mode='anchor')
        #save_figure(?,  cfg) 


def run_my_diagnostic(cfg):
    """Calculate and save final variables into .nc files."""
    my_files_dict = group_metadata(cfg['input_data'].values(), 'dataset')
    print ('my_files_dict', my_files_dict)
    all_variables = calculate_variables (my_files_dict)
    for key in my_files_dict.keys():
        diagnostic_file = get_diagnostic_filename(key, cfg)
        var = all_variables[key]
        if key == "ERA5":
            cube_list = iris.cube.CubeList([
                var['PV'], var['Arctic_temperature'], var['Psl_Ural'], 
                var['Psl_Sib'], var['Psl_Aleut'], var['heat_flux']])
        elif key == "HadISST":
            cube_list = iris.cube.CubeList([
                var['BK_sic'], var ['Ok_sic']])
        else:
            cube_list = iris.cube.CubeList([
                var['PV'], var['Arctic_temperature'], var['Psl_Ural'], 
                var['Psl_Sib'], var['Psl_Aleut'], var['heat_flux'],
                var['BK_sic'], var ['Ok_sic']])
        iris.save(cube_list, diagnostic_file)

#    if plot_timeseries == "True":
#

if __name__ == '__main__':
    with run_diagnostic() as config:
        run_my_diagnostic(config)