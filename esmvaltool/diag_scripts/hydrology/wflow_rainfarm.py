"""wflow diagnostic."""
import logging
from pathlib import Path
import os

import iris
import dask.array as da
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

def get_input_cubes(cfg):
    """ Return a dict with all (preprocessed) input files """
    provenance = create_provenance_record()
    pr_dir = Path(cfg['input_files'][-1])
    print(pr_dir)
    pr_file = next(pr_dir.glob('*pr*_downscaled*.nc'))
    print(pr_file)
    pr_cube = iris.load_cube(str(pr_file))
    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(input_data,
                                        'short_name',
                                        sort='dataset')
    all_vars = {'pr': pr_cube}
    for short_name in grouped_input_data:
        logger.info("Loading variable %s", short_name)
        input_files = [attr['filename'] for attr in grouped_input_data[short_name]]
        allyears = iris.load_cubes(input_files).concatenate_cube()
        all_vars[short_name] = allyears
        provenance['ancestors'].append(input_files)
    return all_vars, provenance

def lapse_rate_correction(height):
    """ Temperature correction over a given height interval """
    gamma = iris.coords.AuxCoord(0.0065,
        long_name='Environmental lapse rate',
        units='K m-1')
    return height*gamma

def to_volume(mass, density=1000.):
    rho = iris.coords.AuxCoord(density,
        long_name='Density',
        units='kg m-3')
    return mass/rho

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
    src_dtemp_compat.data = src_dtemp.data.copy()

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

def debruin_PET(tas, psl, rsds, rsdt, **kwargs):
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
        long_name='Gas  dry air',
        # source='Wallace and Hobbs (2006), 2.6 equation 3.14',
        units='J K-1 kg-1')

    lambda_ = iris.coords.AuxCoord(2.45e6,
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

def main(cfg):
    """Process data for use as input to the wflow hydrological model """
    all_vars, provenance = get_input_cubes(cfg)
    # These keys are now available in all_vars:
    # > tas (air_temperature)
    # > pr (precipitation_flux)
    # > psl (air_pressure_at_mean_sea_level)
    # > rsds (surface_downwelling_shortwave_flux_in_air)
    # > rsdt (toa_incoming_shortwave_flux)
    # > orog (surface_altitude)

    # Interpolating precipitation to the target grid
    # Read the target cube, which contains target grid and target elevation
    dem_path = os.path.join(cfg['auxiliary_data_dir'], cfg['dem_file'])
    dem = iris.load_cube(dem_path)

    ## Processing orography
    logger.info("Processing variable orography")
    orog = all_vars['orog']

    ## Processing precipitation
    logger.info("Processing variable precipitation_flux")
    pr = all_vars['pr']
    pr_dem = preproc.regrid(pr, target_grid=dem, scheme='linear')

    ## Processing temperature
    logger.info("Processing variable temperature")
    tas = all_vars['tas']
    tas_dem = regrid_temperature(tas, orog, dem)
    
    ## Processing air pressure at mean sea level
    logger.info("Processing variable air pressure at mean sea level")
    psl = all_vars['psl']
    psl_dem = preproc.regrid(psl, target_grid=dem, scheme='linear')

    ## Processing air pressure at mean sea level
    logger.info("Processing variable surface downwelling shortwave flux in air")
    rsds = all_vars['rsds']
    rsds_dem = preproc.regrid(rsds, target_grid=dem, scheme='linear')

    ## Processing air pressure at mean sea level
    logger.info("Processing variable surface downwelling shortwave flux in air")
    rsdt = all_vars['rsdt']
    rsdt_dem = preproc.regrid(rsdt, target_grid=dem, scheme='linear')

    ## Processing Reference EvapoTranspiration (PET)
    logger.info("Processing variable PET")
    pet = debruin_PET(tas=tas_dem, psl=psl_dem, rsds=rsds_dem, rsdt=rsdt_dem)
    pet.var_name = 'pet'
    pet_dem = preproc.regrid(pet, target_grid=dem, scheme='linear')

    ## Convert units   
    pet_dem.units = pet_dem.units / 'kg m-3'
    pet_dem.data = pet_dem.core_data() / 1000.
    pet_dem.convert_units('mm day-1')

    pr_dem.units = pr_dem.units / 'kg m-3'
    pr_dem.data = pr_dem.core_data() / 1000.
    pr_dem.convert_units('mm day-1')
    
    tas_dem.convert_units('degC')


    # Round times to integer number of days
    time_coord = pr_dem.coord('time')
    time_coord.points = da.floor(time_coord.core_points())
    time_coord.bounds = None

    time_coord = tas_dem.coord('time')
    time_coord.points = da.floor(time_coord.core_points())
    time_coord.bounds = None

    time_coord = pet_dem.coord('time')
    time_coord.points = da.floor(time_coord.core_points())
    time_coord.bounds = None

    # Save output
    # Output format: "wflow_local_forcing_ERA5_Meuse_1990_2018.nc"
    cubelist = iris.cube.CubeList([pr_dem, tas_dem, pet_dem])
    basin = cfg['basin_name']
    dataset = cfg['dataset']
    startyear = cfg['startyear']
    endyear = cfg['endyear']
    description = '_'.join(['wflow_local_forcing', dataset, basin, startyear, endyear])
    output_file = get_diagnostic_filename(description, cfg)
    iris.save(cubelist, output_file, fill_value=1.e20)

    # Store provenance
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(output_file, provenance)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)

