"""wflow diagnostic."""
import logging
from pathlib import Path

import iris

from esmvalcore.preprocessor import extract_region, regrid
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic)

logger = logging.getLogger(Path(__file__).name)


def create_provenance_record():
    """Create a provenance record."""
    record = {
        'caption':
        "Forcings for the wflow hydrological model.",
        'domains': ['global'],
        'authors': [
            'kalverla_peter',
            'camphuijsen_jaro',
            'alidoost_sarah',
            'aerts_jerom',
            'andela_bouwe',
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


def get_input_cubes(metadata):
    """Create a dict with all (preprocessed) input files."""
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


def save(cubes, dataset, provenance, cfg):
    """Save cubes to file.

    Output format: "wflow_local_forcing_ERA5_Meuse_1990_2018.nc"
    """
    time_coord = cubes[0].coord('time')
    start_year = time_coord.cell(0).point.year
    end_year = time_coord.cell(-1).point.year
    basename = '_'.join([
        'wflow_local_forcing',
        dataset,
        cfg['basin_name'],
        str(start_year),
        str(end_year),
    ])
    output_file = get_diagnostic_filename(basename, cfg)
    logger.info("Saving cubes to file %s", output_file)
    iris.save(cubes, output_file, fill_value=1.e20)

    # Store provenance
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(output_file, provenance)


def lapse_rate_correction(height):
    """Temperature correction over a given height interval."""
    gamma = iris.coords.AuxCoord(0.0065,
                                 long_name='Environmental lapse rate',
                                 units='K m-1')
    return height * gamma


def regrid_temperature(src_temp, src_height, target_height):
    """Convert temperature to target grid with lapse rate correction."""
    # Convert 2m temperature to sea-level temperature (slt)
    src_dtemp = lapse_rate_correction(src_height)
    src_slt = src_temp.copy(data=src_temp.core_data() + src_dtemp.core_data())

    # Interpolate sea-level temperature to target grid
    target_slt = regrid(src_slt, target_grid=target_height, scheme='linear')

    # Convert sea-level temperature to new target elevation
    target_dtemp = lapse_rate_correction(target_height)
    target_temp = target_slt
    target_temp.data = target_slt.core_data() - target_dtemp.core_data()

    return target_temp


def tetens_derivative(temp):
    """Compute derivative of Teten's formula for saturated vapor pressure.

    Tetens formula (https://en.wikipedia.org/wiki/Tetens_equation) :=
    es(T) = e0 * exp(a * T / (T + b))

    Derivative (checked with Wolfram alpha)
    des / dT = a * b * e0 * exp(a * T / (b + T)) / (b + T)^2
    """
    # Ensure temperature is in degC
    temp.convert_units('degC')

    e0 = iris.coords.AuxCoord(
        6.112,
        long_name='Saturated vapour pressure at 273 Kelvin',
        units='hPa')
    emp_a = 17.67  # empirical constant a
    emp_b = iris.coords.AuxCoord(
        243.5,
        long_name='Empirical constant b in Tetens formula',
        units='degC')
    exponent = iris.analysis.maths.exp(emp_a * temp / (emp_b + temp))
    # return emp_a * emp_b * e0 * exponent / (emp_b + temp)**2
    # iris.exceptions.NotYetImplementedError: coord * coord (emp_b * e0)
    # workaround:
    tmp1 = emp_a * emp_b
    tmp2 = e0 * exponent / (emp_b + temp)**2
    return tmp1 * tmp2


def debruin_pet(tas, psl, rsds, rsdt):
    """Determine De Bruin (2016) reference evaporation.

    Implement equation 6 from De Bruin (10.1175/JHM-D-15-0006.1)
    """
    # Definition of constants
    rv = iris.coords.AuxCoord(
        461.51,
        long_name='Gas constant water vapour',
        # source='Wallace and Hobbs (2006), 2.6 equation 3.14',
        units='J K-1 kg-1')

    rd = iris.coords.AuxCoord(
        287.0,
        long_name='Gas constant dry air',
        # source='Wallace and Hobbs (2006), 2.6 equation 3.14',
        units='J K-1 kg-1')

    lambda_ = iris.coords.AuxCoord(
        2.5e6,
        long_name='Latent heat of vaporization',
        # source='Wallace and Hobbs 2006' divide by 86400 for seconds,
        # copy de Bruin method from Marmot model diag
        units='J kg-1')

    cp = iris.coords.AuxCoord(
        1004,
        long_name='Specific heat of dry air constant pressure',
        # source='Wallace and Hobbs 2006',
        units='J K-1 kg-1')

    beta = iris.coords.AuxCoord(
        20,
        long_name='Correction Constant',
        # source='De Bruin (2016), section 4a',
        units='W m-2')

    cs = iris.coords.AuxCoord(
        110,
        long_name='Empirical constant',
        # source = 'De Bruin (2016), section 4a',
        units='W m-2')

    # Unit checks:
    psl.convert_units('hPa')
    tas.convert_units('degC')

    # Variable derivation
    delta_svp = tetens_derivative(tas)
    # gamma = rv/rd * cp*msl/lambda_
    # iris.exceptions.NotYetImplementedError: coord / coord
    gamma = rv.points[0] / rd.points[0] * cp * psl / lambda_

    # Renaming for consistency with paper
    kdown = rsds
    kdown_ext = rsdt

    # Equation 6
    rad_term = (1 - 0.23) * kdown - cs * kdown / kdown_ext
    ref_evap = delta_svp / (delta_svp + gamma) * rad_term + beta

    pet = ref_evap / lambda_
    pet.var_name = 'pet'
    return pet


def main(cfg):
    """Process data for use as input to the wflow hydrological model."""
    input_metadata = cfg['input_data'].values()

    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        all_vars, provenance = get_input_cubes(metadata)
        # These keys are now available in all_vars:
        # > tas (air_temperature)
        # > pr (precipitation_flux)
        # > psl (air_pressure_at_mean_sea_level)
        # > rsds (surface_downwelling_shortwave_flux_in_air)
        # > rsdt (toa_incoming_shortwave_flux)
        # > orog (surface_altitude)

        # Interpolating precipitation to the target grid
        # Read the target cube, which contains target grid and target elevation
        dem_path = Path(cfg['auxiliary_data_dir']) / cfg['dem_file']
        dem = iris.load_cube(str(dem_path))
        dem = extract_region(dem, **cfg['region'])

        logger.info("Processing variable precipitation_flux")
        pr_dem = regrid(all_vars['pr'], target_grid=dem, scheme='linear')

        logger.info("Processing variable temperature")
        tas_dem = regrid_temperature(all_vars['tas'], all_vars['orog'], dem)

        logger.info("Processing variable Reference EvapoTranspiration (PET)")
        pet = debruin_pet(
            tas=all_vars['tas'],
            psl=all_vars['psl'],
            rsds=all_vars['rsds'],
            rsdt=all_vars['rsdt'],
        )
        pet_dem = regrid(pet, target_grid=dem, scheme='linear')

        logger.info("Converting units")
        pet_dem.units = pet_dem.units / 'kg m-3'
        pet_dem.data = pet_dem.core_data() / 1000.
        pet_dem.convert_units('mm day-1')

        pr_dem.units = pr_dem.units / 'kg m-3'
        pr_dem.data = pr_dem.core_data() / 1000.
        pr_dem.convert_units('mm day-1')

        tas_dem.convert_units('degC')

        cubes = iris.cube.CubeList([pr_dem, tas_dem, pet_dem])
        save(cubes, dataset, provenance, cfg)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
