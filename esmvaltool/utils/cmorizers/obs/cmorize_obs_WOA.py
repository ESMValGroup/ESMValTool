"""
# #############################################################################
# ESMValTool CMORizer for WOA data
# #############################################################################
#
# Tier
#    Tier 2: other freely-available dataset.
#
# Source
#    https://data.nodc.noaa.gov/woa/WOA13/DATAv2/
#
# Last access
#    20190131
#
# Download and processing instructions
#    Download the following files:
#      temperature/netcdf/decav81B0/1.00/woa13_decav81B0_t00_01.nc
#      salinity/netcdf/decav81B0/1.00/woa13_decav81B0_s00_01.nc
#      oxygen/netcdf/all/1.00/woa13_all_o00_01.nc
#      nitrate/netcdf/all/1.00/woa13_all_n00_01.nc
#      phosphate/netcdf/all/1.00/woa13_all_p00_01.nc
#      silicate/netcdf/all/1.00/woa13_all_i00_01.nc
#
# Modification history
#    20190131-A_pred_va: adapted to v2.
#    20190131-A_demo_le: written.
#
# #############################################################################
"""

import datetime
import logging
import os

import iris
from cf_units import Unit

from esmvaltool.utils.cmorizers.obs.utilities import (_add_metadata,
                                                      _convert_timeunits,
                                                      _save_variable)

logger = logging.getLogger(__name__)

timestamp = datetime.datetime.utcnow()
timestamp_format = "%Y-%m-%d %H:%M:%S"

# used vars
ALL_VARS = ['thetao', 'so', 'no3', 'po4', 'si', 'o2']

# project at hand
proj = {
    'dataset': 'WOA',
    'version': 'L3',
    'realm': 'clim',
    'field': 'TO3Y',
    'frequency': {
        'thetao': 'Omon',
        'so': 'Omon',
        'no3': 'Oyr',
        'po4': 'Oyr',
        'si': 'Oyr',
        'o2': 'Oyr'
    },
    'metadata_attributes': {
        'tier': '2',
        'source': 'https://data.nodc.noaa.gov/woa/WOA13/DATAv2/',
        'comment': 'cmorized for ESMValTool v2',
        'CMOR conventions': 'CF/CMOR3',
        'CMOR created': timestamp.strftime(timestamp_format)
    }
}

# all years to be analyzed
ALL_YEARS = [
    2000,
]

# specific CMOR nomenclature items
VAR_TO_FILENAME = {
    'thetao': 'woa13_decav81B0_t',
    'so': 'woa13_decav81B0_s',
    'o2': 'woa13_all_o',
    'no3': 'woa13_all_n',
    'po4': 'woa13_all_p',
    'si': 'woa13_all_i'
}

# Reference year
REF_YEAR = 2000

# specific fields names from raw obs files
FIELDS = {
    'si':
    'Objectively analyzed mean fields for moles_concentration_of_silicate_in_sea_water at standard depth levels.',
    'thetao':
    'Objectively analyzed mean fields for sea_water_temperature at standard depth levels.',
    'so':
    'Objectively analyzed mean fields for salinity at standard depth levels.',
    'po4':
    'Objectively analyzed mean fields for moles_concentration_of_phosphate_in_sea_water at standard depth levels.',
    'no3':
    'Objectively analyzed mean fields for moles_concentration_of_nitrate_in_sea_water at standard depth levels.',
    'o2':
    'Average of all unflagged interpolated values at each standard depth level for volume_fraction_of_oxygen_in_sea_water in each grid-square which contain at least one measurement.'
}

# cmor standard names
STANDARD_NAMES = {
    'si': 'mole_concentration_of_silicate_in_sea_water',
    'thetao': 'sea_water_potential_temperature',
    'so': 'sea_water_salinity',
    'po4': 'mole_concentration_of_phosphate_in_sea_water',
    'no3': 'mole_concentration_of_nitrate_in_sea_water',
    'o2': 'mole_concentration_of_dissolved_molecular_oxygen_in_sea_water'
}

# cmor long names
LONG_NAMES = {
    'si': 'Dissolved Silicate Concentration',
    'thetao': 'Sea Water Potential Temperature',
    'so': 'Sea Water Salinity',
    'po4': 'Dissolved Phosphate Concentration',
    'no3': 'Dissolved Nitrate Concentration',
    'o2': 'Dissolved Oxygen Concentration'
}


def _fix_coords(cube):
    """Fix the time units and values to something sensible."""
    # fix individual coords
    for cube_coord in cube.coords():
        # fix time
        if cube_coord.var_name == 'time':
            logger.info("Fixing time...")
            cube.coord('time').convert_units(
                Unit('days since 1950-1-1 00:00:00', calendar='gregorian'))
            if len(cube.coord('time').points) > 1:
                cube.coord('time').guess_bounds()
        # fix longitude
        if cube_coord.var_name in ['lon', 'longitude']:
            logger.info("Fixing longitude...")
            cube.coord('longitude').var_name = 'lon'
            cube.coord('longitude').long_name = 'longitude'
            if cube.coord('longitude').points[0] < 0.:
                cube.coord('longitude').points = \
                    cube.coord('longitude').points + 179.5
                if not cube.coord('longitude').has_bounds():
                    cube.coord('longitude').guess_bounds()
                else:
                    cube.coord('longitude').bounds = None
                    cube.coord('longitude').guess_bounds()
                cube.attributes['geospatial_lon_min'] = 0.
                cube.attributes['geospatial_lon_max'] = 360.
        # fix latitude
        if cube_coord.var_name in ['lat', 'latitude']:
            logger.info("Fixing latitude...")
            cube.coord('latitude').var_name = 'lat'
            cube.coord('latitude').long_name = 'latitude'
            if not cube.coord('latitude').has_bounds():
                cube.coord('latitude').guess_bounds()
        if cube_coord.var_name == 'depth':
            logger.info("Fixing depth...")
            cube.coord('depth').standard_name = 'depth'
            cube.coord('depth').long_name = 'ocean depth coordinate'
            cube.coord('depth').var_name = 'lev'

    # remove CS
    cube.coord('latitude').coord_system = None
    cube.coord('longitude').coord_system = None

    return cube


def _fix_metadata(cube, var):
    """Fix all aspects of metadata for different vars."""
    mol_m3 = ['si', 'po4', 'no3', 'o2']
    if var in mol_m3:
        cube.units = Unit('mol m-3')
    if var == 'thetao':
        cube.convert_units(Unit('kelvin'))
    if var == 'so':
        cube.units = Unit('Unknown')
    return cube


def _fix_data(cube, var):
    """Specific data fixes for different variables."""
    mll_to_mol = ['po4', 'si', 'no3']
    if var in mll_to_mol:
        cube.data = cube.data / 1000.  # Convert from ml/l to mol/m^3
    if var == 'o2':
        cube.data = cube.data * 44.661 / 1000.  # Convert from ml/l to mol/m^3
    return cube


def extract_variable(var, raw_file, out_dir, yr):
    """Extract to all vars."""
    cubes = iris.load(raw_file)
    field = FIELDS[var]
    for cube in cubes:
        if cube.long_name == field:
            cube.standard_name = STANDARD_NAMES[var]
            cube.long_name = LONG_NAMES[var]
            cube.var_name = var
            _convert_timeunits(cube, yr)
            _fix_coords(cube)
            _fix_data(cube, var)
            _fix_metadata(cube, var)
            _add_metadata(cube, proj)
            _save_variable(cube, var, out_dir, yr, proj)


def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    logger.info("Starting cmorization for WOA OBS files: Tier2")
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    # run the cmorization
    for var in ALL_VARS:
        if not os.path.exists(out_dir):
            os.path.makedirs(out_dir)
        for yr in ALL_YEARS:
            file_suffix = str(yr)[-2:] + '_' + str(yr + 1)[-2:] + '.nc'
            raw_file = os.path.join(in_dir, VAR_TO_FILENAME[var] + file_suffix)
            logger.info("CMORizing var %s in file %s", var, raw_file)
            extract_variable(var, raw_file, out_dir, yr)
