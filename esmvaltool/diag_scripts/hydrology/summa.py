"""summa diagnostic."""
import logging
from pathlib import Path
import warnings
import numpy as np

# import dask.array as da
import iris

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic)
from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(Path(__file__).name)


def get_provenance_record(dataset):
    """Create a provenance record."""
    record = {
        'caption': "Forcings for the SUMMA hydrological model.",
        'domains': ['global'],
        'authors': [
            # 'kalverla_peter',
            # 'alidoost_sarah',
            'camphuijsen_jaro',
        ],
        'projects': [
            'ewatercycle',
        ],
        'references': [
            'acknow_project',
        ],
        'dataset': [dataset],
    }
    return record


def _rename_var(input_dict, output_name):
    """Change the cmor names for SUMMA"""
    variables = {}
    for key in output_name:
        if key in input_dict:
            cube = input_dict[key]
            cube.var_name = output_name[key]
            variables[key] = cube
    return variables


def compute_windspeed(u_component, v_component):
    """ Compute wind speed magnitude based on vector components """
    wind_speed = (u_component**2 + v_component**2)**.5
    return wind_speed


def logarithmic_profile(windspeed, measurement_height):
    """ Convert wind speed from input height to 2m with logarithmic profile """
    warnings.warn('This version of the logarithmic wind profile is violating\
                  general principles of transparancy. Better replace it with\
                  a more general formula incorporating friction velocity or\
                  roughness length')
    cube = windspeed * 4.87/np.log(67.8*measurement_height-5.42)
    return cube

def compute_specific_humidity(dewpoint_temperature, surface_pressure):
    # see e.g. https://github.com/Unidata/MetPy/issues/791
    # TODO confirm this after teleco
    d2m = dewpoint_temperature
    airpres = surface_pressure
    kelvin = 273.15
    vapour_pressure_act = 611 * np.exp((17.76*(d2m-kelvin))/(d2m-29.65))
    vapor_pressure = vapour_pressure_act /1000.
    air_pressure = airpres /1000.
    eps = 0.62196351
    mix_rat = (eps * vapor_pressure) / (air_pressure - vapor_pressure)
    spechum = mix_rat/ (1 + mix_rat)
    return spechum


def _fix_cube(input_dict):
    """Fixing attributes for new cube"""
    for key in input_dict:
        cube = input_dict[key]
        if key in {'windspd'}:
            # remove the height_0 coordinate (10m)
            cube.remove_coord(cube.coord('height'))
            # add the height coordinate (2m)
            utils.add_scalar_height_coord(cube, 2.)
            cube.var_name = 'windspd'
            cube.standard_name = 'wind_speed'
        if key in {'spechum'}:
            cube.var_name = 'spechum'
            cube.standard_name = 'specific_humidity'
            cube.unit = 'g/g'
    return cube


def main(cfg):
    """Process data for use as input to the summa hydrological model """
    input_data = cfg['input_data'].values()
    logger.info(input_data)
    grouped_input_data = group_metadata(input_data,
                                        'standard_name',
                                        sort='dataset')
    variables = {}
    # output variable's name in SUMMA
    output_var_name = {
        'tas':'airtemp',
        'rsds':'SWRadAtm',
        'windspd':'windspd',
        'strd':'LWRadAtm',
        'pr':'pptrate',
        'spechum':'spechum',
        'ps':'airpres'
    }
    for standard_name in grouped_input_data:
        # get the dataset name to use in save function later
        # TODO add support multiple dataset in one diagnostic
        dataset = grouped_input_data[standard_name][0]['alias']
        logger.info("Processing variable %s", standard_name)
        cube_list_all_years = iris.cube.CubeList()
        for attributes in grouped_input_data[standard_name]:
            logger.info("Processing dataset %s", attributes['dataset'])
            input_file = attributes['filename']
            cube = iris.load_cube(input_file)
            cube_list_all_years.append(cube)
        cube_all_years = cube_list_all_years.concatenate_cube()
        key = grouped_input_data[standard_name][0]['short_name']
        variables[key] = cube_all_years

        # Do stuff
        # The data need to be aggregated for each HRU (subcatchment)
        # Inti's `decomposed` function in extract_shape should add
        # this as a dimension to the cubes, so it's just a matter of
        # aggregating latitude and longitude. The resulting cubes
        # will have dimensions 'time' and 'hru'.
        #
        # Lorenz workshop prepared output used metsim to compute spechum
        # and a weird logarithmic wind profile expression... see notebook at:
        # ssh userX@jupyter.ewatercycle.org
        # cd /mnt/data/lorentz-models/SUMMA/summa_era5_scripts/
        #
        # Unit conversion:
        # - precip: kg m-2 s-1
        # - radiation: w m-2
        # - temperature: K
        # - wind speed: m s-1
        # - pressure: Pa
        # - specific humidity: g g-1
        #
        # example output file can also be found on jupyter server.

    # Extract wind component variables from the cube list
    u_component = variables['uas']
    v_component = variables['vas']

    # compute wind speed and convert to 2m
    wind_speed = compute_windspeed(u_component, v_component)
    wind_speed_2m = logarithmic_profile(wind_speed, 10)
    # Add wind speed to cube dict
    variables['windspd'] = wind_speed_2m

    # Fix cube attributes
    variables['windspd'] = _fix_cube(variables)

    # TODO: add specific humidity calculation
    dewpoint_temperature = variables['tdps']
    surface_pressure = variables['ps']
    variables['spechum'] = compute_specific_humidity(dewpoint_temperature,
                                                     surface_pressure)
    # Fix cube attributes
    # TODO check the attributes, avoid code duplication
    variables['spechum'] = _fix_cube(variables)

    # Select and rename the desired variables
    variables = _rename_var(variables, output_var_name)

    # Make a list from all cubes in dictionary
    cube_list_all_vars = iris.cube.CubeList()
    for key in variables:
        cube_list_all_vars.append(variables[key])
    # Save data
    basename = dataset + '_summa'
    output_file = get_diagnostic_filename(basename, cfg)
    iris.save(cube_list_all_vars, output_file, fill_value=1.e20)

    # Store provenance
    provenance_record = get_provenance_record(dataset)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(output_file, provenance_record)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
