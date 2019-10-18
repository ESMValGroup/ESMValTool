"""wflow diagnostic."""
import logging
from pathlib import Path
import os

# import dask.array as da
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
    # - air_temperature
    # - precipitation_flux

    # Interpolating precipitation to the target grid
    # Read the target cube, which contains target grid and target elevation
    dem_path = os.path.join(cfg['auxiliary_data_dir'],cfg['dem_file'])
    dem_cube = iris.load_cube(dem_path)

    # Read source orography (add era5 later) and try to make it cmor compatible
    era_orography_path = os.path.join(
        cfg['auxiliary_data_dir'], cfg['source_orography'])
    era_orography_cube = iris.load_cube(era_orography_path)
    oro.coord('longitude').long_name="Longitude"
    oro.coord('latitude').long_name="Latitude"

    ## Processing precipitation
    # Interpolate precipitation to the target grid and save new output
    pr = all_vars['precipitation_flux']
    pr_dem = preproc.regrid(pr, target_grid=dem_cube, scheme='linear')
    iris.save(pr_dem, output_file, fill_value=1.e20) #TODO: change filename

    ## Processing temperature
    tas = all_vars['air_temperature']

    # Convert geopotential to surface elevation
    g = iris.coords.AuxCoord(9.80665,
                             long_name='Acceleration due to gravity',
                             units='m s-2')
    era_height = era_orography_cube / g

    # Convert surface elevation to sea-level temperature correction
    gamma = iris.coords.AuxCoord(0.0065,
                                 long_name='Environmental lapse rate',
                                 units='K/m')

    # Convert surface (2m) temperature to sea-level temperature
    t_sea_level = tas + era_height*gamma
    #BUG this doesn't work:
    # ValueError: This operation cannot be performed as there are differing coordinates (time, longitude, latitude) remaining which cannot be ignored.
    # Related (?): https://github.com/SciTools/iris/issues/2765

    # Interpolate sea-level temperature to target grid
    tsl_target_grid = preproc.regrid(t_sea_level, target_grid=dem_cube, scheme='linear')

    # Convert sea-level temperature to new target elevation
    t_target_elevation = tsl_target_grid - dem_cube*gamma

    # Save output
    iris.save(t_target_elevation, output_file, fill_value=1.e20)

    # TODO
    # - Solve iris issue to make lapse rate temperature correction work
    # - Move variable manipulation to separate functions
    # - Add function(s) to calculate potential evapotranspiration
    # - Check whether the correct units are used
    # - See whether we can work with wflow pcraster .map files directly
    #   (currently, we use .nc dem files that Jerom converted externally)
    # - Compare output to prepared input during workshop
    # - Save the output files with correct variable- and file-names
    #   Output format: wflow_local_forcing_ERA5_Meuse_1990_2018.nc

if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
