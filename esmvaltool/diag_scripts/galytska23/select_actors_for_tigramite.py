"""Personal diagnostics for Arctic-midlatitude research."""

# to manipulate iris cubes
import iris
import numpy as np
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
    save_data,
)
from esmvaltool.diag_scripts.shared._base import get_diagnostic_filename


def polar_vortex(dict_item):
    var = iris.load_cube(dict_item['filename'])
    var = var.collapsed('air_pressure', iris.analysis.MEAN)
    var.data *= -1
    var.var_name = 'PV'
    return var


def slp(dict_item):
    var = iris.load_cube(dict_item['filename'])
    var.data /= 100
    return var


def calculate_heat_flux(list_va_ta):
    for i in range(0, len(list_va_ta)):
        hf = list_va_ta[0] * list_va_ta[1]
        hf_anom = anomalies(hf, period='monthly')
        hf_anom_zm = zonal_statistics(hf_anom, operator='mean')
        hf_anom_zm.var_name = 'heat_flux'
        return hf_anom_zm


def run_my_diagnostic(cfg):
    """Save variables into .nc files.

    Returns:
     string; runs the user diagnostic
    """
    # assemble the data dictionary keyed by dataset name
    # via usage of group_metadata function
    my_files_dict = group_metadata(cfg['input_data'].values(), 'dataset')
    #    print ("DICTIONARY", my_files_dict)

    for key, value in my_files_dict.items():
        diagnostic_file = get_diagnostic_filename(key, cfg)
        #        print (diagnostic_file)
        tmp_list = list()
        for item in value:
            if item['preprocessor'] == 'pv':
                pv = polar_vortex(item)
            elif item['preprocessor'] == 'pre_tas':
                tas = iris.load_cube(item['filename'])
                tas.var_name = 'Arctic_temperature'
            elif item['preprocessor'] == 'pre_tas_baffin':
                tas_baffin = iris.load_cube(item['filename'])
                tas_baffin.var_name = 'Temperature_Baffin'
            elif item['preprocessor'] == 'pre_tas_sib':
                tas_sib = iris.load_cube(item['filename'])
                tas_sib.var_name = 'Temperature_Sib'
            elif item['preprocessor'] == 'pressure_ural':
                psl_Ural = slp(item)
                psl_Ural.var_name = 'Psl_Ural'
            elif item['preprocessor'] == 'pressure_sib':
                psl_Sib = slp(item)
                psl_Sib.var_name = 'Psl_Sib'
            elif item['preprocessor'] == 'pressure_aleut':
                psl_Aleut = slp(item)
                psl_Aleut.var_name = 'Psl_Aleut'
            elif item['preprocessor'] == 'zonal_wind':
                zon_wind = iris.load_cube(item['filename'])
                zon_wind.var_name = 'zonal_wind'
            elif item['preprocessor'] == 'bk_ice':
                sic_BK = iris.load_cube(item['filename'])
                sic_BK.var_name = 'BK_sic'
            elif item['preprocessor'] == 'ok_ice':
                sic_Ok = iris.load_cube(item['filename'])
                sic_Ok.var_name = 'Ok_sic'
            elif item['preprocessor'] == 'heat_flux':
                var = iris.load_cube(item['filename'])
                var_avg = area_statistics(var, operator='mean')
                var_mermean = meridional_statistics(var, operator='mean')
                deviation = var_mermean - var_avg
                tmp_list.append(deviation)
        hf = calculate_heat_flux(tmp_list)
        if key == "ERA5":
            cube_list = iris.cube.CubeList([
                pv, tas, tas_baffin, tas_sib, psl_Ural, psl_Sib, psl_Aleut,
                zon_wind, hf
            ])
            iris.save(cube_list, diagnostic_file)
        elif key == "HadISST":
            cube_list = iris.cube.CubeList([sic_BK, sic_Ok])
            iris.save(cube_list, diagnostic_file)
        else:
            cube_list = iris.cube.CubeList([
                pv, tas, tas_baffin, tas_sib, psl_Ural, psl_Sib, psl_Aleut,
                zon_wind, hf, sic_BK, sic_Ok
            ])
            iris.save(cube_list, diagnostic_file)

    return 'Done with saving .nc files'


if __name__ == '__main__':
    with run_diagnostic() as config:
        run_my_diagnostic(config)
