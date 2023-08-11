"""
Arctic-midlatitude teleconnections Diagnostics.

Diagnostic calculates timeseries needed for Causal Model Evaluation
of Arctic-midlatitude teleconnections.

  Description:
    This diagnostics calculates timeseries of the variables that represent
    Arctic-midlatitude teleconenctions that are further used for the
    Causal Model Evaluation of CMIP6 (Galytska et al., 2023). The output of
    this diagnostics is a .nc file per data source. Optionally this diagnostics
    plots the timeseries of the evolution of each selected variable. If the
    user kept "plot_timeseries: True" in recipe_galytska23jgr.yml, then
    "variable_to_plot:" expects the name of the variable to be plotted.
    Possible options for "variable_to_plot:" are:
    Arctic_temperature
    Psl_Ural
    Psl_Sib
    Psl_Aleut
    PV
    heat_flux
    BK_sic
    Ok_sic
  Author: Evgenia Galytska, IUP-UB
          egalytska@iup.physik.uni-bremen.de
  Project: USMILE
"""
import logging
from pathlib import Path
import iris
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from esmvalcore.preprocessor import (
    anomalies,
    area_statistics,
    meridional_statistics,
    zonal_statistics,
)

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
)
from esmvaltool.diag_scripts.shared._base import (
    get_plot_filename,
)

logger = logging.getLogger(Path(__file__).stem)


def get_provenance_record(ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'authors': ['galytska_evgenia'],
        'ancestors': ancestor_files,
        'projects': ['usmile'],
        'references': [
            'galytska23jgr',
        ],
    }
    return record


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
    """Match preprocessor name and corresponding calculations."""
    if var == 'pv':
        out_var = calculate_polar_vortex(item)
    elif var == 'pre_tas':
        out_var = calculate_arctic_tas(item)
    elif var == 'pressure_ural':
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
    logger.debug("Variables are calculated for the following datasources:%s",
                 input_dict.keys())
    dictionary = {}
    for key, value in input_dict.items():
        logger.debug("Calculating final variables for %s dataset", key)
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
            # Calculate heat flux for all data sources except HadISST
            heat_flux = calculate_heat_flux(tmp_list)
            dictionary[key].setdefault(heat_flux.var_name, heat_flux)
    return dictionary


def plotting_support(cube, key, **kwargs):
    """Help for the pretty plot."""
    if cube.coords('time', dim_coords=True):
        ih.unify_time_coord(cube)
    iris.quickplot.plot(cube, label=key, **kwargs)
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.ylabel("Anomalies, " + str(cube.units))
    plt.title(f"Time series of monthly mean {cube.var_name.upper()} anomalies")
    plt.xticks(rotation=45, ha="right", rotation_mode="anchor")


def plot_timeseries(dictionary, var, cfg):
    """Timeseries plot."""
    fig = plt.figure(figsize=(10, 4))
    sns.set_style('whitegrid')
    colors = plt.cm.viridis(np.linspace(0, 1, len(dictionary.keys())))
    baseplotname = f"Timeseries_{var}_anomalies"
    filename = get_plot_filename(baseplotname, cfg)
    for i, key in enumerate(dictionary.keys()):
        if var not in ('BK_sic', 'Ok_sic'):
            if key == "HadISST":
                continue
            if key != 'ERA5':
                plotting_support(dictionary[key][var], key,
                                 color=colors[i])
            else:
                plotting_support(dictionary[key][var], key,
                                 color='k', linewidth=2)
        else:
            if key == "ERA5":
                continue
            if key != 'HadISST':
                plotting_support(dictionary[key][var], key, color=colors[i])
            else:
                plotting_support(dictionary[key][var], key, color='blue',
                                 linewidth=2)
    fig.savefig(filename, bbox_inches='tight')


def main(cfg):
    """Calculate and save final variables into .nc files."""
    my_files_dict = group_metadata(cfg['input_data'].values(), 'dataset')
    all_variables = calculate_variables(my_files_dict)
    # Check is timeseries should be plotted
    if cfg['plot_timeseries'] is True:
        plot_timeseries(all_variables, cfg['variable_to_plot'], cfg)
    for key in my_files_dict:
        logger.info("Processing final calculations in dataset %s", key)
        prov_record = get_provenance_record([key])
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
        save_data(key, prov_record, cfg, cube_list)
        logger.info("%s data is saved in .nc", key)
    logger.info("Done.")


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
