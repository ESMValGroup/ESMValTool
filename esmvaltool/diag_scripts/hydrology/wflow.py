"""wflow diagnostic."""
import logging
from pathlib import Path
import os

import iris
from esmvalcore import preprocessor as preproc

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic,
                                            select_metadata)

logger = logging.getLogger(Path(__file__).name)

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

def get_input_cubes(cfg):
    """ Return a dict with all (preprocessed) input files """
    #TODO:
    # - Add concat as option instead of default?
    # - Use short_name instead of standard_name
    # - Load multiple years in one iris statement, not a loop

    provenance = create_provenance_record()
    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(input_data,
                                        'standard_name',
                                        sort='dataset')
    all_vars = {}
    for standard_name in grouped_input_data:
        logger.info("Loading variable %s", standard_name)

        # Combine multiple years into a singe iris cube
        all_years = []
        for attributes in grouped_input_data[standard_name]:
            input_file = attributes['filename']
            cube = iris.load_cube(input_file)
            all_years.append(cube)
            provenance['ancestors'].append(input_file)
        allyears = iris.cube.CubeList(all_years).concatenate_cube()

        all_vars[standard_name] = allyears
    return all_vars, provenance

def geopotential_to_height(geopotential):
    """ Convert geopotential to geopotential height """
    gravity = iris.coords.AuxCoord(9.80665,
        long_name='Acceleration due to gravity',
        units='m s-2')
    return geopotential/gravity

def lapse_rate_correction(height):
    """ Temperature correction over a given height interval """
    gamma = iris.coords.AuxCoord(0.0065,
        long_name='Environmental lapse rate',
        units='K m-1')
    return height*gamma

def regrid_temperature(src_temp, src_height, target_height):
    """ Convert temperature to target grid with lapse rate correction """
    #TODO: Fix iris issue to get rid of workaround

    src_dtemp = lapse_rate_correction(src_height)

    # src_slt = src_temp + dtemp
    # ValueError: This operation cannot be performed as there are differing
    # coordinates (time) remaining which cannot be ignored.
    # Related (?): https://github.com/SciTools/iris/issues/2765

    # Workaround: overwrite data in compatible cube
    src_dtemp_compat = src_temp.collapsed('time', iris.analysis.MEAN)
    src_dtemp_compat.data = src_height.data.copy()

    # Convert 2m temperature to sea-level temperature (slt)
    src_slt = src_temp + src_dtemp_compat

    # Interpolate sea-level temperature to target grid
    target_slt = preproc.regrid(src_slt, target_grid=target_height, scheme='linear')

    # Convert sea-level temperature to new target elevation
    target_dtemp = lapse_rate_correction(target_height)
    target_dtemp_compat = target_slt.collapsed('time', iris.analysis.MEAN)
    target_dtemp_compat.data = target_dtemp.data.copy()

    target_temp = target_slt - target_dtemp_compat
    target_temp.var_name = src_temp.var_name
    return target_temp

def debruin_PET(t2m, msl, ssrd, tisr):
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
    msl = preproc.convert_units(msl,'hPa')
    t2m = preproc.convert_units(t2m,'degC')

    # Variable derivation
    delta_svp = tetens_derivative(t2m)
    # gamma = rv/rd * cp*msl/lambda_
    # iris.exceptions.NotYetImplementedError: coord / coord
    gamma = rv.points[0]/rd.points[0] * cp*msl/lambda_

    # Renaming for consistency with paper
    kdown = ssrd
    kdown_ext = tisr

    # Equation 6
    rad_term = (1-0.23)*kdown - cs*kdown/kdown_ext
    ref_evap = delta_svp / (delta_svp + gamma) * rad_term + beta

    return ref_evap/lambda_

def main(cfg):
    """Process data for use as input to the wflow hydrological model """
    all_vars, provenance = get_input_cubes(cfg)
    # These keys are now available in all_vars:
    # print('############')
    # print(all_vars)
    # > air_temperature
    # > precipitation_flux
    # > air_pressure_at_mean_sea_level
    # > dew_point_temperature
    # > surface_downwelling_shortwave_flux_in_air
    # > toa_incoming_shortwave_flux

    # Interpolating precipitation to the target grid
    # Read the target cube, which contains target grid and target elevation
    dem_path = os.path.join(cfg['auxiliary_data_dir'], cfg['dem_file'])
    dem = iris.load_cube(dem_path)

    # Read source orography (add era5 later) and try to make it cmor compatible
    era_orography_path = os.path.join(cfg['auxiliary_data_dir'], cfg['source_orography'])
    oro = iris.load_cube(era_orography_path)

    ## Processing precipitation
    logger.info("Processing variable precipitation_flux")
    pr = all_vars['precipitation_flux']
    pr_dem = preproc.regrid(pr, target_grid=dem, scheme='linear')

    ## Processing temperature
    logger.info("Processing variable temperature")
    tas = all_vars['air_temperature']
    tas_dem = regrid_temperature(tas, oro, dem)

    ## Processing Reference EvapoTranspiration (PET)
    logger.info("Processing variable PET")
    t2m = all_vars['air_temperature']
    msl = all_vars['air_pressure_at_mean_sea_level']
    ssrd = all_vars['surface_downwelling_shortwave_flux_in_air']
    tisr = all_vars['toa_incoming_shortwave_flux']
    pet = debruin_PET(t2m, msl, ssrd, tisr)
    pet.var_name = 'potential_evapotranspiration'
    pet_dem = preproc.regrid(pet, target_grid=dem, scheme='linear')

    # Save output
    # Output format: "wflow_local_forcing_ERA5_Meuse_1990_2018.nc"
    dataset = iris.cube.CubeList([pr_dem, tas_dem, pet_dem])
    print('####################')
    print(cfg)
    output_file = get_diagnostic_filename('wflow_local_forcing_'
        # + '_' + cfg['dataset'] + '_' + cfg['basin'] + '_' + cfg['start_year'] + '_' + cfg['end_year'], cfg)
        + 'ERA-Interim_Meuse_2000_2001', cfg)
    iris.save(dataset, output_file, fill_value=1.e20)

    # Store provenance
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(output_file, provenance)

    # TODO
    # - Use CMORized orography as recipe variable rather than aux file
    # - Check whether the correct units are used
    # - See whether we can work with wflow pcraster .map files directly
    #   (currently, we use .nc dem files that Jerom converted externally)
    # - Compare output to prepared input during workshop
    # - Save the output files with correct variable- and file-names


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
