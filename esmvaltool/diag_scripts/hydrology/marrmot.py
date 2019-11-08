"""marrmot diagnostic."""
import logging
from pathlib import Path

# import dask.array as da
import iris
import scipy.io as sio

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic)

logger = logging.getLogger(Path(__file__).name)
from esmvalcore import preprocessor as preproc

def create_provenance_record():
    """Create a provenance record."""
    record = {
        'caption': "Forcings for the wflow hydrological model.",
        'domains': ['global'],
        'authors': [
            'kalverla_peter',
            'camphuijsen_jaro',
            # 'alidoost_sarah',
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

def debruin_pet(tas, psl, rsds, rsdt, **kwargs):
    """ Determine De Bruin (2016) reference evaporation

    Implement equation 6 from De Bruin (10.1175/JHM-D-15-0006.1)
    Implementation using iris.
    """
    # Definition of constants
    rv = iris.coords.AuxCoord(461.51,
        long_name='Gas constant water vapour',
        # source='Wallace and Hobbs (2006), 2.6 equation 3.14',
        units='J K-1 kg-1')

    rd = iris.coords.AuxCoord(287.0,
        long_name='Gas constant dry air',
        # source='Wallace and Hobbs (2006), 2.6 equation 3.14',
        units='J K-1 kg-1')

    lambda_ = iris.coords.AuxCoord(2.5e6,
        long_name='Latent heat of vaporization',
        # source='Wallace and Hobbs 2006',
        units='J kg-1')

    cp = iris.coords.AuxCoord(1004,
        long_name='Specific heat of dry air constant pressure',
        # source='Wallace and Hobbs 2006',
        units='J K-1 kg-1')

    beta = iris.coords.AuxCoord(20,
        long_name='Correction Constant',
        # source='De Bruin (2016), section 4a',
        units='W m-2')

    cs = iris.coords.AuxCoord(110,
        long_name = 'Empirical constant',
        # source = 'De Bruin (2016), section 4a',
        units = 'W m-2')

    def tetens_derivative(temp):
        """ Derivative of Teten's formula for saturated vapor pressure.

        Tetens formula (https://en.wikipedia.org/wiki/Tetens_equation) :=
        es(T) = e0 * exp(a * T / (T + b))

        Derivate (checked with Wolfram alpha)
        des / dT = a * b * e0 * exp(a * T / (b + T)) / (b + T)^2
        """
        # Assert temperature is in degC
        temp = preproc.convert_units(temp,'degC')

        e0 = iris.coords.AuxCoord(6.112,
            long_name='Saturated vapour pressure at 273 Kelvin',
            units='hPa')
        emp_a = 17.67 # empirical constant a
        emp_b = iris.coords.AuxCoord(243.5,
            long_name='Empirical constant b in Tetens formula',
            units='degC')
        exponent = iris.analysis.maths.exp(emp_a * temp / (emp_b + temp))
        # return emp_a * emp_b * e0 * exponent / (emp_b + temp)**2
        # iris.exceptions.NotYetImplementedError: coord * coord (emp_b * e0)
        # workaround:
        tmp1 = emp_a * emp_b
        tmp2 = e0 * exponent / (emp_b + temp)**2
        return tmp1 * tmp2

    # Unit checks:
    psl = preproc.convert_units(psl,'hPa')
    tas = preproc.convert_units(tas,'degC')

    # Variable derivation
    delta_svp = tetens_derivative(tas)
    # gamma = rv/rd * cp*msl/lambda_
    # iris.exceptions.NotYetImplementedError: coord / coord
    gamma = rv.points[0]/rd.points[0] * cp*psl/lambda_

    # Renaming for consistency with paper
    kdown = rsds
    kdown_ext = rsdt

    # Equation 6
    rad_term = (1-0.23)*kdown - cs*kdown/kdown_ext
    ref_evap = delta_svp / (delta_svp + gamma) * rad_term + beta

    return ref_evap/lambda_

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
        input_files = [attr['filename'] for attr in grouped_input_data[short_name]]
        allyears = iris.load_cubes(input_files).concatenate_cube()
        all_vars[short_name] = allyears
        provenance['ancestors'].append(input_files)
    return all_vars, provenance

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
    temp = all_vars['tas']
    temp = preproc.area_statistics(temp, operator='mean')
    # convert kelvin to celcius
    temp.convert_units('celsius')

    ## Processing Precipitation (pr)
    logger.info("Processing variable pr")
    precip = all_vars['pr']
    precip = preproc.area_statistics(precip, operator='mean')
    # convert kg/m2/s to kg/m2/day (or mm/day)
    precip.convert_units('kg m-2 day-1')

    ## Processing Reference EvapoTranspiration (PET)
    logger.info("Processing variable PET")
    pet = debruin_pet(**all_vars)
    pet.var_name = 'potential_evapotranspiration'
    pet = preproc.area_statistics(pet, operator='mean')
    # convert kg/m2/s to kg/m2/day (or mm/day)
    pet.convert_units('kg m-2 day-1')

    # # Save output
    # cubelist = iris.cube.CubeList([pr_accumulated, tas_accumulated, pet_accumulated])
    # # add temp to matlab structure
    # output_file = get_diagnostic_filename('marrmot_input', cfg)
    # iris.save(cubelist, output_file, fill_value=1.e20)

    # Get the dataset name
    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(input_data, 'dataset')
    dataset_name = list(grouped_input_data.keys())
    if len(dataset_name) == 1:
        dataset_name = dataset_name[0]
    basename = dataset_name + '_marrmot'

    ## Save to matlab structure
    # TODO add data_orogin
    # Get the start and end times as an array with lenght 6
    coord = temp.coord(axis='T')
    time_start = coord.units.num2date(coord.points[0])
    time_start = time_start.strftime("%Y %m %d %H %M %S").split()
    time_start = [int(time) for time in time_start]
    time_end = coord.units.num2date(coord.points[-1])
    time_end = time_end.strftime("%Y %m %d %H %M %S").split()
    time_end = [int(time) for time in time_end]

    # make data structure
    mdict = {
        'precip': precip.data,
        'temp': temp.data,
        'pet': pet.data,
        'delta_t': 1,  # this could also be extracted from the cube
        'time_unit': 'day'
        }
    output_data = {
        'forcing': mdict,
        'time_start': time_start,
        'time_end': time_end
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
