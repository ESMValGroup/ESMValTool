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

logger = logging.getLogger(Path(__file__).name)


def get_provenance_record(ancestor_file):
    """Create a provenance record."""
    record = {
        'caption': "Forcings for the summa hydrological model.",
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
        'ancestors': [ancestor_file],
    }
    return record


def main(cfg):
    """Process data for use as input to the summa hydrological model """
    input_data = cfg['input_data'].values()
    logger.info(input_data)
    grouped_input_data = group_metadata(input_data,
                                        'standard_name',
                                        sort='dataset')
    variables = {}
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

    # Add wind speed to cube list and remove old vars
    variables[wind_speed_2m.var_name] = wind_speed_2m
    # for key in variables.items():
    #     if key in {'uas', 'vas'}:
    #         variables.pop(key)

    # TODO: add specific humidity calculation
    # specific_humidity = compute_specific_humidity(dewpoint_temperature, surface_pressure)

    # Select the desired variables
    print('#######################')
    cube_list_all_vars = iris.cube.CubeList()
    for key in variables:
        cube_list_all_vars.append(variables[key])
    # Save data # check the dataset!
    basename = dataset + '_summa'
    output_file = get_diagnostic_filename(basename, cfg)
    iris.save(cube_list_all_vars, output_file, fill_value=1.e20)

    # Store provenance # check this!
    provenance_record = get_provenance_record(input_file)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(output_file, provenance_record)


def compute_windspeed(u_component, v_component):
    """ Compute wind speed magnitude based on vector components """
    wind_speed = (u_component**2 + v_component**2)**.5
    #TODO: add variable names, attributes etc.
    wind_speed.var_name = 'windspd'
    return wind_speed


def logarithmic_profile(windspeed, measurement_height):
    """ Convert wind speed from input height to 2m with logarithmic profile """
    warnings.warn('This version of the logarithmic wind profile is violating\
                  general principles of transparancy. Better replace it with\
                  a more general formula incorporating friction velocity or\
                  roughness length')
    return windspeed * 4.87/np.log(67.8*measurement_height-5.42)


def compute_specific_humidity(dewpoint_temperature, surface_pressure):
    # see e.g. https://github.com/Unidata/MetPy/issues/791
    return

if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
