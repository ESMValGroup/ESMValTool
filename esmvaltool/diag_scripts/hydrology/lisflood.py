"""LISFLOOD diagnostic."""
import logging
from pathlib import Path

import iris
import numpy as np
import xarray as xr
from iris.analysis.maths import exp as iris_exp

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic)

logger = logging.getLogger(Path(__file__).name)


def get_provenance_record(ancestor_files):
    """Create a provenance record."""
    record = {
        'caption': "Forcings for the LISFLOOD hydrological model.",
        'domains': ['global'],
        'authors': [
            'verhoeven_stefan',
            'kalverla_peter',
            'camphuijsen_jaro',
        ],
        'projects': [
            'ewatercycle',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestor_files,
    }
    return record


def get_input_cubes(metadata):
    """Create a dict with all (preprocessed) input files."""
    inputs = {}
    ancestors = {}
    for attributes in metadata:
        short_name = attributes['short_name']
        if short_name in inputs:
            raise ValueError(f"Multiple input files found for variable "
                             f"'{short_name}'.")
        filename = attributes['filename']
        logger.info("Loading variable %s", short_name)
        cube = iris.load_cube(filename)
        cube.attributes.clear()
        inputs[short_name] = cube
        ancestors[short_name] = [filename]

    return inputs, ancestors


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


def compute_vapour_pressure(tdps):
    """Compute vapour pressure using tetens formula."""
    # taken from Eq. 3.21 of Goudriaan (1977;
    # https://library.wur.nl/WebQuery/wurpubs/70980)
    if tdps.units != 'degC':
        raise Exception('tdps should be in degC')
    esat = 6.10588 * iris_exp(17.32491 * tdps / (tdps + 238.102))
    esat.var_name = 'e'
    esat.long_name = 'Daily Actual Water Vapour Pressure'
    esat.standard_name = 'water_vapor_pressure'
    esat.units = 'hPa'
    esat.attributes['comment'] = ''.join(
        ('Actual water vapour pressure of air near the surface calculated',
         ' from tdps using Tetens formula'))
    return esat


def compute_windspeed(uas, vas):
    """Compute absolute wind speed from horizontal components."""
    sfc_wind = (uas**2 + vas**2)**.5
    sfc_wind.var_name = 'sfcWind'
    sfc_wind.long_name = 'Daily-Mean Near-Surface Wind Speed'
    sfc_wind.standard_name = 'wind_speed'
    comment = 'near-surface (usually, 10 meters) wind speed.'
    sfc_wind.attributes['comment'] = comment
    return sfc_wind


def save(xrds, var_name, dataset, cfg):
    """Save processed cube to a lisflood-compatible file."""
    start_year = int(xrds.time[0].dt.year)
    end_year = int(xrds.time[-1].dt.year)
    basename = '_'.join([
        'lisflood',
        dataset,
        cfg['catchment'],
        var_name,
        str(start_year),
        str(end_year),
    ])
    output_file = get_diagnostic_filename(basename, cfg)
    xrds.to_netcdf(output_file, encoding={var_name: {'_FillValue': 1.e20}})
    return output_file


def main(cfg):
    """Process data for use as input to the LISFLOOD hydrological model."""
    input_metadata = cfg['input_data'].values()
    logger.info(input_metadata)

    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        cubes, ancestors = get_input_cubes(metadata)

        if dataset == 'ERA5':
            shift_era5_time_coordinate(cubes['tas'])
            shift_era5_time_coordinate(cubes['tdps'])
            shift_era5_time_coordinate(cubes['uas'])
            shift_era5_time_coordinate(cubes['vas'])

        # Compute additional variables as input for lisvap
        tdps = cubes.pop('tdps')
        uas = cubes.pop('uas')
        vas = cubes.pop('vas')
        cubes['e'] = compute_vapour_pressure(tdps)
        ancestors['e'] = ancestors['tdps']
        cubes['sfcWind'] = compute_windspeed(uas, vas)
        ancestors['sfcWind'] = ancestors['uas'] + ancestors['vas']

        cubes['pr'].units = 'mm d-1'

        for var_name, cube in cubes.items():
            # Western emisphere longitudes should be negative
            points = cube.coord('longitude').points
            cube.coord('longitude').points = (points + 180) % 360 - 180
            # latitudes decreasing
            cube = cube[:, ::-1, ...]

            # convert to xarray dataset (xrds)
            # remove coordinate bounds drop extra coordinates and reorder
            xrds = xr.DataArray.from_iris(cube).to_dataset()
            ordered_coords = ['lon', 'lat', 'time']
            extra_coords = np.setdiff1d(xrds.coords, ordered_coords)
            xrds = xrds.drop(extra_coords)[ordered_coords + [var_name]]

            output_file = save(xrds, var_name, dataset, cfg)

            # Store provenance
            provenance_record = get_provenance_record(ancestors[var_name])
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(output_file, provenance_record)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
