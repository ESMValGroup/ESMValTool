"""wflow diagnostic."""
import logging
from pathlib import Path

import iris
import numpy as np
from osgeo import gdal

from esmvaltool.diag_scripts.hydrology.derive_evspsblpot import debruin_pet
from esmvaltool.diag_scripts.hydrology.lazy_regrid import lazy_regrid
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
        cube = iris.load_cube(filename)
        cube.attributes.clear()
        all_vars[short_name] = cube
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
        'wflow',
        dataset,
        cfg['basin'],
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
    gamma = iris.coords.AuxCoord(np.float32(0.0065),
                                 long_name='Environmental lapse rate',
                                 units='K m-1')
    return height * gamma


def regrid_temperature(src_temp, src_height, target_height, scheme):
    """Convert temperature to target grid with lapse rate correction."""
    # Convert 2m temperature to sea-level temperature (slt)
    src_dtemp = lapse_rate_correction(src_height)
    src_slt = src_temp.copy(data=src_temp.core_data() + src_dtemp.core_data())

    # Interpolate sea-level temperature to target grid
    target_slt = lazy_regrid(src_slt, target_height, scheme)

    # Convert sea-level temperature to new target elevation
    target_dtemp = lapse_rate_correction(target_height)
    target_temp = target_slt
    target_temp.data = target_slt.core_data() - target_dtemp.core_data()

    return target_temp


def load_dem(filename):
    """Load DEM into iris cube."""
    logger.info("Reading digital elevation model from %s", filename)
    if filename.suffix.lower() == '.nc':
        cube = iris.load_cube(str(filename))
    elif filename.suffix.lower() == '.map':
        cube = _load_pcraster_dem(filename)
    else:
        raise ValueError(f"Unknown file format {filename}. Supported formats "
                         "are '.nc' and '.map'.")
    for coord in 'longitude', 'latitude':
        if not cube.coord(coord).has_bounds():
            logger.warning("Guessing DEM %s bounds", coord)
            cube.coord(coord).guess_bounds()
    return cube


def _load_pcraster_dem(filename):
    """Load DEM from a PCRASTER .map file."""
    dataset = gdal.Open(str(filename))
    lon_offset, lon_step, _, lat_offset, _, lat_step = dataset.GetGeoTransform(
    )
    lon_size, lat_size = dataset.RasterXSize, dataset.RasterYSize
    data = dataset.ReadAsArray()
    data = np.ma.masked_less(data, -1e8)
    dataset = None

    lons = lon_offset + lon_step * (np.arange(lon_size) + 0.5)
    lats = lat_offset + lat_step * (np.arange(lat_size) + 0.5)

    lon_coord = iris.coords.DimCoord(
        lons,
        var_name='lon',
        standard_name='longitude',
        units='degrees',
    )
    lat_coord = iris.coords.DimCoord(
        lats,
        var_name='lat',
        standard_name='latitude',
        units='degrees',
    )

    cube = iris.cube.Cube(
        data,
        var_name='height',
        units='m',
        dim_coords_and_dims=[
            (lat_coord, 0),
            (lon_coord, 1),
        ],
    )
    return cube


def check_dem(dem, cube):
    """Check that the DEM and extract_region parameters match."""
    for coord in ('longitude', 'latitude'):
        start_dem_coord = dem.coord(coord).cell(0).point
        end_dem_coord = dem.coord(coord).cell(-1).point
        start_cube_coord = cube.coord(coord).cell(0).point
        end_cube_coord = cube.coord(coord).cell(-1).point
        if start_dem_coord < start_cube_coord:
            logger.warning(
                "Insufficient data available, input data starts at %s "
                "degrees %s, but should be at least one grid "
                "cell larger than the DEM start at %s degrees %s.",
                start_cube_coord, coord, start_dem_coord, coord)
        if end_dem_coord > end_cube_coord:
            logger.warning(
                "Insufficient data available, input data ends at %s "
                "degrees %s, but should be at least one grid "
                "cell larger than the DEM end at %s degrees %s.",
                end_cube_coord, coord, end_dem_coord, coord)


def shift_era5_time_coordinate(cube, shift=30):
    """Shift instantaneous variables (default = 30 minutes forward in time).

    After this shift, as an example:
    time format [1990, 1, 1, 11, 30, 0] will be [1990, 1, 1, 12, 0, 0].
    For aggregated variables, already time format is [1990, 1, 1, 12, 0, 0].
    """
    time = cube.coord(axis='T')
    time.points = time.points + shift / (24 * 60)
    time.bounds = None
    time.guess_bounds()
    return cube


def main(cfg):
    """Process data for use as input to the wflow hydrological model."""
    input_metadata = cfg['input_data'].values()

    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        all_vars, provenance = get_input_cubes(metadata)

        if dataset == 'ERA5':
            shift_era5_time_coordinate(all_vars['tas'])
            shift_era5_time_coordinate(all_vars['psl'])

        # Interpolating variables onto the dem grid
        # Read the target cube, which contains target grid and target elevation
        dem_path = Path(cfg['auxiliary_data_dir']) / cfg['dem_file']
        dem = load_dem(dem_path)
        check_dem(dem, all_vars['pr'])

        logger.info("Processing variable precipitation_flux")
        scheme = cfg['regrid']
        pr_dem = lazy_regrid(all_vars['pr'], dem, scheme)

        logger.info("Processing variable temperature")
        tas_dem = regrid_temperature(
            all_vars['tas'],
            all_vars['orog'],
            dem,
            scheme,
        )

        logger.info("Processing variable potential evapotranspiration")
        if 'evspsblpot' in all_vars:
            pet = all_vars['evspsblpot']
            pet_dem = lazy_regrid(pet, dem, scheme)
        else:
            logger.info("Potential evapotransporation not available, deriving")
            psl_dem = lazy_regrid(all_vars['psl'], dem, scheme)
            rsds_dem = lazy_regrid(all_vars['rsds'], dem, scheme)
            rsdt_dem = lazy_regrid(all_vars['rsdt'], dem, scheme)
            pet_dem = debruin_pet(
                tas=tas_dem,
                psl=psl_dem,
                rsds=rsds_dem,
                rsdt=rsdt_dem,
            )
        pet_dem.var_name = 'pet'

        logger.info("Converting units")
        pet_dem.units = pet_dem.units / 'kg m-3'
        pet_dem.data = pet_dem.core_data() / 1000.
        pet_dem.convert_units('mm day-1')

        pr_dem.units = pr_dem.units / 'kg m-3'
        pr_dem.data = pr_dem.core_data() / 1000.
        pr_dem.convert_units('mm day-1')

        tas_dem.convert_units('degC')

        # Adjust longitude coordinate to wflow convention
        for cube in [tas_dem, pet_dem, pr_dem]:
            cube.coord('longitude').points = (cube.coord('longitude').points +
                                              180) % 360 - 180

        cubes = iris.cube.CubeList([pr_dem, tas_dem, pet_dem])
        save(cubes, dataset, provenance, cfg)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
