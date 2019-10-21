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


def get_provenance_record(ancestor_file):
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
        'ancestors': [ancestor_file],
    }
    return record

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
        units='K/m')
    return height*gamma

def regrid_temperature(src_temp, src_height, target_height):
    """ Convert temperature to target grid with lapse rate correction """
    #TODO: Fix issue to get rid of workaround

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
    return target_temp

def main(cfg):
    """Process data for use as input to the wflow hydrological model """
    input_data = cfg['input_data'].values()
    logger.info(input_data)
    grouped_input_data = group_metadata(input_data,
                                        'standard_name',
                                        sort='dataset')

    # Make a dictionary containing all variables
    all_vars = {}
    for standard_name in grouped_input_data:
        logger.info("Processing variable %s", standard_name)

        # Combine multiple years into a singe iris cube
        all_years = []
        for attributes in grouped_input_data[standard_name]:
            input_file = attributes['filename']
            cube = iris.load_cube(input_file)
            all_years.append(cube)
        allyears = iris.cube.CubeList(all_years).concatenate_cube()

        all_vars[standard_name] = allyears

        # For now, save intermediate output
        output_file = get_diagnostic_filename(
            Path(input_file).stem + '_wflow_test', cfg)
        iris.save(allyears, output_file, fill_value=1.e20)


    # These keys are now available in all_vars:
    # print('###########3')
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
    # oro.coord('longitude').long_name="Longitude"
    # oro.coord('latitude').long_name="Latitude"

    ## Processing precipitation
    # Interpolate precipitation to the target grid and save new output
    pr = all_vars['precipitation_flux']
    pr_dem = preproc.regrid(pr, target_grid=dem, scheme='linear')
    iris.save(pr_dem, output_file, fill_value=1.e20) #TODO: change filename

    ## Processing temperature
    tas = all_vars['air_temperature']
    tas_dem = regrid_temperature(tas, oro, dem)

    ## Calculating
    # Save output
    iris.save(tas_dem, output_file, fill_value=1.e20)

    # TODO
    # - Add function(s) to calculate potential evapotranspiration
    # - Check whether the correct units are used
    # - See whether we can work with wflow pcraster .map files directly
    #   (currently, we use .nc dem files that Jerom converted externally)
    # - Compare output to prepared input during workshop
    # - Save the output files with correct variable- and file-names
    #   Output format: wflow_local_forcing_ERA5_Meuse_1990_2018.nc
    # - Add provencance stuff again in the diagnostic

def debruin_PET(t2m, msl, d2m, tp, ssrd, tisr):
    """ Determine De Bruin (2016) reference evaporation. """
    beta = 20 # empirical constant from De Bruin 2016 in W m-2
    Cs = 110 # empirical constant from De Bruin 2016 in W m-2
    pha = msl/100 # pressure in hPa ??? Replace by unit conversion/check
    cp = 1005 # specific heat of dry air in J kg-1 K-1 (??)
    epsilon = 0.622 # ratio of dry air over water vapour (mass/pressure) (WH06 eq 3.14)

    esat = saturated_vapor_pressure(t2m)
    delta_wvp = delta_wvp(t2m, esat)
    lambda_ = latent_heat(t2m)
    Gamma = (cp * pha)/(epsilon * lambda_)
    net_rad_down = (1-0.23)*ssrd - Cs * ssrd / (tisr + 0.00001)
    pot_evap = delta_wvp / (delta_wvp + Gamma) * net_rad_down + beta

    # What are we checking here?
    pot_evap = np.where(tisr == 0, 0, pot_evap) # no PET if there is no sunlight??
    pot_evap = np.where((pot_evap/lambda_)*1000*100 > 0,(pot_evap/lambda_)*1000*100,0)
    return pot_evap

def saturated_vapor_pressure(temperature):
    """ Calculate saturated vapor pressure according to Tetens' (??) formula """
    #TODO: Check formula
    #TODO: Add unit check (t2m should be in celcius??)
    e0 = iris.coords.AuxCoord(6.112,
        long_name='reference_vapor_pressure????',
        units='hPa')
    a = 17.67 # empirical constant 1  --> find source
    b = 243.5 # empirical constant 2  --> find source
    return e0*np.exp((a * temperature)/(temperature+b))

def delta_wvp(temperature, esat):
    """ Determine slope of water vapour pressure according to ..."""
    a = 17.269 # empirical constant 1 --> find source
    b = 243.5 # empirical constant 2  --> find source
    return esat*(a/(temperature+b))*(1-(temperature/temperature + b))

def latent_heat(temperature):
    """ Latent heat of vaporization ??? """
    #TODO: where does this formula come from?
    #TODO: Add units to these constants: J kg-1 K-1 ???
    return 2.502e6 - 2250 * t2m



if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
