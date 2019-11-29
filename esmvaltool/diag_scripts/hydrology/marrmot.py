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
    """Compute the derivative of Teten's formula for saturated vapor pressure.

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
    emp_a = 17.67  # empirical constant a

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
    """Define constants to compute De Bruin (2016) reference evaporation.

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

    # Latent heat of vaporization in J kg-1 (or J m-2 day-1)
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
    """Compute De Bruin (2016) reference evaporation.

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
    rad_term = (1 - 0.23) * kdown - cs_const * kdown / kdown_ext
    # the unit is W m-2
    ref_evap = delta_svp / (delta_svp + gamma) * rad_term + beta

    pet = ref_evap / lambda_
    pet.var_name = 'potential_evapotranspiration'
    pet.convert_units('kg m-2 day-1')  # equivalent to mm/day
    return pet


def get_input_cubes(metadata):
    """Return a dict with all (preprocessed) input files."""
    provenance = create_provenance_record()
    all_vars = {}
    for attributes in metadata:
        short_name = attributes['short_name']
        if short_name in all_vars:
            raise ValueError(
                f"Multiple input files found for variable '{short_name}'.")
        filename = attributes['filename']
        logger.info("Loading variable %s", short_name)
        all_vars[short_name] = iris.load_cube(filename)
        provenance['ancestors'].append(filename)
    return all_vars, provenance


def _get_extra_info(cube):
    """Get start/end time and origin of cube.

    Get the start and end times as an array with length 6
    and get latitude and longitude as an array with length 2
    """
    coord = cube.coord('time')
    time_start_end = []
    for index in 0, -1:
        time_val = coord.cell(index).point
        time_val = [
            time_val.year,
            time_val.month,
            time_val.day,
            time_val.hour,
            time_val.minute,
            time_val.second,
        ]
        time_val = [float(time) for time in time_val]
        time_start_end.append(time_val)

    # Add data_origin
    lat_lon = [
        cube.coord(name).points[0] for name in ('latitude', 'longitude')
    ]
    return time_start_end, lat_lon


def main(cfg):
    """Process data for use as input to the marrmot hydrological model.

    These variables are needed in all_vars:
    tas (air_temperature)
    pr (precipitation_flux)
    psl (air_pressure_at_mean_sea_level)
    rsds (surface_downwelling_shortwave_flux_in_air)
    rsdt (toa_incoming_shortwave_flux)
    """
    input_metadata = cfg['input_data'].values()
    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        all_vars, provenance = get_input_cubes(metadata)

        # Processing variables and unit conversion
        # Unit of the fluxes in marrmot should be in kg m-2 day-1 (or mm/day)
        logger.info("Processing variable tas")
        temp = preproc.area_statistics(all_vars['tas'], operator='mean')
        temp.convert_units('celsius')

        logger.info("Processing variable pr")
        precip = preproc.area_statistics(all_vars['pr'], operator='mean')
        precip.convert_units('kg m-2 day-1')  # equivalent to mm/day

        logger.info("Processing variable PET")
        all_vars['pet'] = debruin_pet(all_vars)
        pet = preproc.area_statistics(all_vars['pet'], operator='mean')

        # Get the start and end times and latitude longitude
        time_start_end, lat_lon = _get_extra_info(temp)

        # make data structure
        # delta_t_days could also be extracted from the cube
        output_data = {
            'forcing': {
                'precip': precip.data,
                'temp': temp.data,
                'pet': pet.data,
                'delta_t_days': float(1),
                'time_unit': 'day',
            },
            'time_start': time_start_end[0],
            'time_end': time_start_end[1],
            'data_origin': lat_lon,
        }

        # Save to matlab structure
        basename = '_'.join([
            'marrmot',
            dataset,
            cfg['basin'],
            str(int(output_data['time_start'][0])),
            str(int(output_data['time_end'][0])),
        ])
        output_name = get_diagnostic_filename(basename, cfg, extension='mat')
        sio.savemat(output_name, output_data)

        # Store provenance
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(output_name, provenance)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
