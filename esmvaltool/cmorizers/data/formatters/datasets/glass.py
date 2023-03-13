"""ESMValTool CMORizer for GLASS data.

Tier
    Tier 2: other freely-available dataset.

Source
    http://www.glass.umd.edu/Download.html

Last access
    20230313

Download and processing instructions
    Use automatic download feature to get the data.

    By default, this dataset is regridded to a 0.5°x0.5° grid (original
    resolution is 0.05°). If you want to use the original resolution, remove
    the `regrid` section in the configuration file (`GLASS.yml`). Note that in
    this case, preprocessing the dataset with ESMValTool (i.e. every time you
    run the tool) can take a very long time (> 30 min).

"""

import logging
import re
from datetime import datetime, timedelta
from pathlib import Path
from xml.etree.ElementTree import parse as xml_parse

import dateutil.parser
import iris.coord_categorisation
import numpy as np
import pandas as pd
from cf_units import Unit
from esmvalcore.preprocessor import regrid
from iris.coords import DimCoord
from iris.cube import Cube, CubeList
from netCDF4 import Dataset

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


TIME_UNITS = Unit('days since 1950-01-01 00:00:00')
FILE_NAME_REGEX = (
    r'GLASS(?P<var_number>[0-9]{2})(?P<res_mark>[A-Z])'
    r'(?P<data_source>[0-9]{2})\.(?P<version>V[0-9]{2})\.'
    r'A(?P<year>[0-9]{4})(?P<doy>[0-9]{3})\.'
    r'(?P<prod_year>[0-9]{4})(?P<prod_doy>[0-9]{3})\.hdf'
)  # see http://www.glass.umd.edu/Overview.html


def _get_files(in_dir, filename, start_date, end_date):
    """Get input files."""
    all_files = sorted(Path(in_dir).glob(filename))
    files = {}
    for file in all_files:
        match = re.match(FILE_NAME_REGEX, file.name)
        year = int(match.group('year'))
        if start_date.year <= year <= end_date.year:
            files.setdefault(year, [])
            files[year].append(file)
    return files


def _get_daily_time_coord(metadata):
    """Get daily time coordinate from XML metadata."""
    one_day = timedelta(days=1)
    time_info = metadata.find('RangeDateTime')
    time_start = dateutil.parser.parse(
        f"{time_info.find('RangeBeginningDate').text} "
    )
    time_end = dateutil.parser.parse(
        f"{time_info.find('RangeEndingDate').text} "
    )
    time_bounds = TIME_UNITS.date2num(np.hstack([
        pd.date_range(
            time_start, time_end, freq='D'
        ).to_pydatetime().reshape(-1, 1),
        pd.date_range(
            time_start + one_day, time_end + one_day, freq='D'
        ).to_pydatetime().reshape(-1, 1),
    ]))
    time_points = np.mean(time_bounds, axis=1)
    time_coord = DimCoord(
        time_points,
        bounds=time_bounds,
        standard_name='time',
        var_name='time',
        long_name='time',
        units=TIME_UNITS,
    )
    return time_coord


def _get_horizontal_coords(metadata, n_lat, n_lon):
    """Get horizontal coords from XML metadata."""
    horizontal_info = (
        metadata
        .find('SpatialDomainContainer')
        .find('HorizontalSpatialDomainContainer')
        .find('GPolygon')
        .find('Boundary')
    )

    # Latitude (sorted 90.0 ... -90.0 in data)
    lat_min = float(horizontal_info.find('DL').find('PointLatitude').text)
    lat_max = float(horizontal_info.find('UR').find('PointLatitude').text)
    d_lat = (lat_max - lat_min) / n_lat
    lat_points = np.linspace(lat_max - d_lat / 2., lat_min + d_lat / 2., n_lat)
    lat_bounds = np.hstack([
        np.linspace(lat_max, lat_min + d_lat, n_lat).reshape(-1, 1),
        np.linspace(lat_max - d_lat, lat_min, n_lat).reshape(-1, 1),
    ])
    lat_coord = DimCoord(
        lat_points,
        bounds=lat_bounds,
        standard_name='latitude',
        var_name='lat',
        long_name='Latitude',
        units='degrees_north',
    )

    # Longitude (sorted -180 ... 180.0 in data)
    lon_min = float(horizontal_info.find('DL').find('PointLongitude').text)
    lon_max = float(horizontal_info.find('UR').find('PointLongitude').text)
    d_lon = (lon_max - lon_min) / n_lon
    lon_points = np.linspace(lon_min + d_lon / 2., lon_max - d_lon / 2., n_lon)
    lon_bounds = np.hstack([
        np.linspace(lon_min, lon_max - d_lon, n_lon).reshape(-1, 1),
        np.linspace(lon_min + d_lon, lon_max, n_lon).reshape(-1, 1),
    ])
    lon_coord = DimCoord(
        lon_points,
        bounds=lon_bounds,
        standard_name='longitude',
        var_name='lon',
        long_name='Longitude',
        units='degrees_east',
    )

    return (lat_coord, lon_coord)


def _get_dim_coords_and_dims(file, n_lat, n_lon):
    """Get coordinates for file from associated XML file."""
    xml_file = file.with_suffix('.hdf.xml')
    if not xml_file.is_file():
        raise ValueError(f"Metadata XML file {xml_file} not found")
    xml_tree = xml_parse(xml_file)
    metadata = xml_tree.getroot().find('GranuleURMetaData')

    daily_time_coord = _get_daily_time_coord(metadata)
    (lat_coord, lon_coord) = _get_horizontal_coords(metadata, n_lat, n_lon)

    return [(daily_time_coord, 0), (lat_coord, 1), (lon_coord, 2)]


def _fix_var_metadata(var_info, cmor_info, cube):
    """Fix variable metadata."""
    cube.convert_units(cmor_info.units)
    utils.fix_var_metadata(cube, cmor_info)


def _extract_variable(var_info, cmor_info, attrs, all_files, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    raw_var = var_info.get('raw_name', var)

    # Configure regriddind
    if var_info.get('regrid'):
        var_info['regrid'].setdefault('target_grid', '1x1')
        var_info['regrid'].setdefault('scheme', 'nearest')
        logger.info(
            "Final dataset will be regridded to %s grid using scheme '%s'",
            var_info['regrid']['target_grid'],
            var_info['regrid']['scheme'],
        )

    # Iterate over all years and save one cube per year
    for (year, files) in all_files.items():
        logger.info("Reading files for year %i", year)
        sub_cubes = CubeList([])
        for file in files:
            logger.debug("Reading file %s", file)
            dataset = Dataset(file, 'r')

            # Coordinates (daily time, lat, lon)
            n_lat = dataset.dimensions[var_info['raw_lat_name']].size
            n_lon = dataset.dimensions[var_info['raw_lon_name']].size
            dim_coords_and_dims = _get_dim_coords_and_dims(file, n_lat, n_lon)
            n_time = dim_coords_and_dims[0][0].shape[0]

            # Build daily cube (simply repeat data for each day in the given
            # interval, usually 8-days-long) to enable calculation of monthly
            # means later (necessary since original 8-day intervals may spread
            # over two months, so a simple mean aggregated_by using
            # month_number would be wrong.
            nc_var = dataset.variables[raw_var]
            cube_data = np.ma.array(
                np.broadcast_to(nc_var[:], (n_time, n_lat, n_lon)),
                mask=np.broadcast_to(nc_var[:].mask, (n_time, n_lat, n_lon)),
            )
            cube = Cube(
                cube_data,
                var_name=var,
                units=var_info.get('raw_units', nc_var.units),
                dim_coords_and_dims=dim_coords_and_dims,
            )
            sub_cubes.append(cube)
        cube = sub_cubes.concatenate_cube()

        # Regrid if desired
        if var_info.get('regrid'):
            cube = regrid(
                cube,
                var_info['regrid']['target_grid'],
                var_info['regrid']['scheme'],
            )

        # Calculate monthly means from daily data
        iris.coord_categorisation.add_month_number(cube, 'time')
        cube = cube.aggregated_by('month_number', iris.analysis.MEAN)
        cube.remove_coord('month_number')

        # Fix coords
        utils.flip_dim_coord(cube, 'latitude')
        cube = cube.intersection(longitude=(0.0, 360.0))

        # Fix metadata
        _fix_var_metadata(var_info, cmor_info, cube)
        utils.set_global_atts(cube, attrs)

        # Save variable
        utils.save_variable(
            cube,
            var,
            out_dir,
            attrs,
            unlimited_dimensions=['time'],
        )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    # Check time ranges
    if start_date is None:
        start_date = datetime(1981, 1, 1)
    if end_date is None:
        end_date = datetime(2023, 12, 31)
    logger.info("Considering years %i--%i", start_date.year, end_date.year)

    # Run the cmorization
    for (var, var_info) in cfg['variables'].items():
        logger.info(
            "CMORizing variable '%s' from files %s", var, var_info['filename']
        )
        files = _get_files(in_dir, var_info['filename'], start_date, end_date)
        glob_attrs['mip'] = var_info['mip']
        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        _extract_variable(var_info, cmor_info, glob_attrs, files, out_dir)
