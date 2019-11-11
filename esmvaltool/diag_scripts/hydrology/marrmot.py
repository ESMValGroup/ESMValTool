"""Marrmot diagnostic."""
import logging
from pathlib import Path

import iris
import scipy.io as sio

from esmvalcore import preprocessor as preproc
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic)

logger = logging.getLogger(Path(__file__).name)

def create_provenance_record():
    """Create a provenance record."""
    record = {
        'caption': "Forcings for the Marrmot hydrological model.",
        'domains': ['global'],
        'authors': [
            'kalverla_peter',
            'camphuijsen_jaro',
            'alidoost_sarah',
        ],
        'projects': [
            'ewatercycle',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': [],
    }
    return record


def tetens_derivative(temp):
    """ Derivative of Teten's formula for saturated vapor pressure.

    Tetens formula (https://en.wikipedia.org/wiki/Tetens_equation) :=
    es(T) = e0 * exp(a * T / (T + b))

    Derivate (checked with Wolfram alpha)
    des / dT = a * b * e0 * exp(a * T / (b + T)) / (b + T)^2
    """
    # Assert temperature is in degC
    temp = preproc.convert_units(temp, 'degC')

    # Saturated vapour pressure at 273 Kelvin
    e0_const = iris.coords.AuxCoord(6.112,
                                    long_name='Saturated vapour pressure',
                                    units='hPa')
    emp_a = 17.67 # empirical constant a

    # Empirical constant b in Tetens formula
    emp_b = iris.coords.AuxCoord(243.5,
                                 long_name='Empirical constant b',
                                 units='degC')
    exponent = iris.analysis.maths.exp(emp_a * temp / (emp_b + temp))
    # return emp_a * emp_b * e0 * exponent / (emp_b + temp)**2
    # iris.exceptions.NotYetImplementedError: coord * coord (emp_b * e0)
    # workaround:
    tmp1 = emp_a * emp_b
    tmp2 = e0_const * exponent / (emp_b + temp)**2
    return tmp1 * tmp2


def get_constants(psl):
    """
    The Definition of rv and rd constants is provided in
    Wallace and Hobbs (2006), 2.6 equation 3.14.
    The Definition of lambda and cp is provided in Wallace and Hobbs 2006.
    The Definition of beta and cs is provided in De Bruin (2016), section 4a.
    """
    # Definition of constants
    rv_const = iris.coords.AuxCoord(461.51,
                                    long_name='Gas constant water vapour',
                                    units='J K-1 kg-1')

    rd_const = iris.coords.AuxCoord(287.0,
                                    long_name='Gas constant dry air',
                                    units='J K-1 kg-1')

    lambda_ = iris.coords.AuxCoord(2.5e6,
                                   long_name='Latent heat of vaporization',
                                   units='J kg-1')

    # Specific heat of dry air constant pressure
    cp_const = iris.coords.AuxCoord(1004,
                                    long_name='Specific heat of dry air',
                                    units='J K-1 kg-1')

    beta = iris.coords.AuxCoord(20,
                                long_name='Correction Constant',
                                units='W m-2')

    cs_const = iris.coords.AuxCoord(110,
                                    long_name='Empirical constant',
                                    units='W m-2')

    # gamma = rv/rd * cp*msl/lambda_
    # iris.exceptions.NotYetImplementedError: coord / coord
    gamma = rv_const.points[0] / rd_const.points[0] * cp_const * psl / lambda_
    return gamma, cs_const, beta, lambda_

def debruin_pet(var_dict):
    """ Determine De Bruin (2016) reference evaporation
    Implement equation 6 from De Bruin (10.1175/JHM-D-15-0006.1)
    """
    # Unit checks:
    psl = preproc.convert_units(var_dict['psl'], 'hPa')
    tas = preproc.convert_units(var_dict['tas'], 'degC')

    # Variable derivation
    delta_svp = tetens_derivative(tas)
    gamma, cs_const, beta, lambda_ = get_constants(psl)

    # the definition of the radiation components according to the paper:
    kdown = var_dict['rsds']
    kdown_ext = var_dict['rsdt']
    # Equation 6
    rad_term = (1-0.23)*kdown - cs_const*kdown/kdown_ext
    ref_evap = delta_svp / (delta_svp + gamma) * rad_term + beta

    pet = ref_evap/lambda_
    pet.var_name = 'potential_evapotranspiration'
    return pet


def get_input_cubes(cfg):
    """ Return a dict with all (preprocessed) input files """
    provenance = create_provenance_record()
    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(input_data,
                                        'short_name',
                                        sort='dataset')
    all_vars = {}
    for short_name in grouped_input_data:
        logger.info("Loading variable %s", short_name)
        input_files = [
            attr['filename'] for attr in grouped_input_data[short_name]
            ]
        allyears = iris.load_cubes(input_files).concatenate_cube()
        all_vars[short_name] = allyears
        provenance['ancestors'].append(input_files)
    return all_vars, provenance


def _get_dataset_name(cfg):
    """ Get the dataset name """
    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(input_data, 'dataset')
    dataset_name = list(grouped_input_data.keys())
    if len(dataset_name) == 1:
        dataset_name = dataset_name[0]
    basename = dataset_name + '_marrmot'
    return basename


def _get_extra_info(cube):
    """ Get the start and end times as an array with length 6
    and get latitude and longitude as an array with length 2
    """
    coord = cube.coord('time')
    time_start_end = []
    for index in 0, -1:
        time_val = coord.units.num2date(coord.points[index])
        time_val = time_val.strftime("%Y %m %d %H %M %S").split()
        time_val = [float(time) for time in time_val]
        time_start_end.append(time_val)

    # Add data_origin
    lat_lon = [cube.coord('latitude').points[0]]
    lat_lon.append(cube.coord('longitude').points[0])
    return time_start_end, lat_lon


def main(cfg):
    """Process data for use as input to the marrmot hydrological model """
    all_vars, provenance = get_input_cubes(cfg)
    # These keys are now available in all_vars:
    # > tas (air_temperature)
    # > pr (precipitation_flux)
    # > psl (air_pressure_at_mean_sea_level)
    # > rsds (surface_downwelling_shortwave_flux_in_air)
    # > rsdt (toa_incoming_shortwave_flux)

    ## Processing temperature
    logger.info("Processing variable tas")
    temp = preproc.area_statistics(all_vars['tas'], operator='mean')
    # convert kelvin to celcius
    temp.convert_units('celsius')

    ## Processing Precipitation (pr)
    logger.info("Processing variable pr")
    precip = preproc.area_statistics(all_vars['pr'], operator='mean')
    # convert kg/m2/s to kg/m2/day (or mm/day)
    precip.convert_units('kg m-2 day-1')

    ## Processing Reference EvapoTranspiration (PET)
    logger.info("Processing variable PET")
    all_vars['pet'] = debruin_pet(all_vars)
    pet = preproc.area_statistics(all_vars['pet'], operator='mean')
    # convert kg/m2/s to kg/m2/day (or mm/day)
    pet.convert_units('kg m-2 day-1')

    # Get the dataset name
    basename = _get_dataset_name(cfg)

    ## Save to matlab structure
    # Get the start and end times and latitude longitude
    time_start_end, lat_lon = _get_extra_info(temp)

    # make data structure
    forcing_dict = {
        'precip': precip.data,
        'temp': temp.data,
        'pet': pet.data,
        'delta_t': 1,  # this could also be extracted from the cube
        'time_unit': 'day'
        }
    output_data = {
        'forcing': forcing_dict,
        'time_start': time_start_end[0],
        'time_end': time_start_end[1],
        'data_origin': lat_lon
        }
    output_name = get_diagnostic_filename(basename, cfg, extension='mat')
    sio.savemat(output_name, output_data)

    # Store provenance
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(output_name, provenance)

    # Do stuff
    # >> A vector with initial values for each of the model stores,
    # >> of size 1x[number of stores].


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
