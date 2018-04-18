#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
MJO assessment
"""
import argparse
import csv
import os
import os.path
import pickle
import tempfile

from datetime import datetime

# use Agg backend for running without X-server
import matplotlib as mpl
mpl.use('Agg')

import iris

import mjo_utils as mu
import mjo_plots as mjp
import reg_lat_long_grid

import diags_level1
import diags_level2
import diags_level3


AUTOASSESS_PARSER = argparse.ArgumentParser(description='Description')
AUTOASSESS_PARSER.add_argument('--ref-suite-id', required=True, help='Reference suite ID.')
AUTOASSESS_PARSER.add_argument('--out-dir', required=True, help='Output directory')
AUTOASSESS_PARSER.add_argument('--data-dir', required=True, help='Directory containing model data.')
AUTOASSESS_PARSER.add_argument('--start-date', required=True, help='Start date for assessment. Format YYYY/MM/DD')
AUTOASSESS_PARSER.add_argument('--end-date', required=True, help='End date for assessment. Format YYYY/MM/DD')
AUTOASSESS_PARSER.add_argument('--obs-dir', required=False)
AUTOASSESS_PARSER.add_argument('--tmp-dir', required=True)
AUTOASSESS_PARSER.add_argument('--ancil-dir', required=False)
AUTOASSESS_PARSER.add_argument('--ncpu', default=1, help='Number of available processors.')
AUTOASSESS_PARSER.add_argument('--no-netcdf', action='store_true', help='Do not write NetCDF files.')


def parse_args(parser):
    """Parse command line arguments, and check all paths are absolute."""
    args = parser.parse_args()

    for arg, val in vars(args).items():
        if '-dir' in arg and not val is None and not os.path.isabs(val):
            msg = ' '.join(['Cli argument ', str(arg), ' not absolute path:', str(val)])
            raise NotAbsolutePath(msg)
    return args


def unpickle_cubes(path):
    """Load cube list from path."""
    with open(path, 'r') as fh:
        cubes = pickle.load(fh)
    return cubes


class NotAbsolutePath(Exception):
    pass


def main():
    AREA_NAME = 'MJO'.lower()

    args = parse_args(AUTOASSESS_PARSER)

    if args.no_netcdf:
        def do_nothing(*args, **kwargs):
            pass
        iris.save = do_nothing

    start_date = datetime(*map(int, args.start_date.split('/')))
    end_date = datetime(*map(int, args.end_date.split('/')))

    model_run_dirs = [
        os.path.join(args.data_dir, dir_) for dir_ in os.walk(args.data_dir).next()[1]
    ]

    cubes_paths = [
        os.path.join(dir_, AREA_NAME, 'cubes.pickle') for dir_ in model_run_dirs
    ]

    cubes = iris.cube.CubeList()
    for cubes_path in cubes_paths:
        cube_list = unpickle_cubes(cubes_path)
        cubes.extend(cube_list)

    # create output directories for single run assessments
    suite_ids = list(set(
        [cube.attributes['MODEL_RUN_ID'] for cube in cubes]
    ))

    suite_ids.remove(args.ref_suite_id)
    suite_ids.append(args.ref_suite_id)
    AREA_OUT_DIR = os.path.join(args.out_dir, '_vs_'.join(suite_ids), AREA_NAME)

    # also run assessment on observational data
    suite_ids.append('OBS')

    for suite_id in suite_ids:
        out_dir = os.path.join(AREA_OUT_DIR, suite_id)
        os.makedirs(out_dir)

    # create tmp dir for large files
    tmp_dir = tempfile.mkdtemp(prefix='MJO_', dir=args.tmp_dir)


    ### Obs data
    # TODO obs data path
    # rename to same CF name as model data
    obs_dir = '/project/MJO_GCSS/eld523/data/local/hadpx/obs_data/MJO_OBS_DATA'

    precip_obs_cube = iris.load_cube(os.path.join(obs_dir, 'precip_daily_GPCP_1996_2009.nc'))
    precip_obs_cube.rename('precipitation_flux')
    precip_obs_cube.attributes['MODEL_RUN_ID'] = 'OBS'

    olr_obs_cube = iris.load_cube(os.path.join(obs_dir, 'olr_daily_NOAA_1989_2012.nc'))
    olr_obs_cube.rename('toa_outgoing_longwave_flux')
    olr_obs_cube.attributes['MODEL_RUN_ID'] = 'OBS'

    eastward_wind_850hPa_obs_cube = iris.load_cube(os.path.join(obs_dir, 'u850_daily_ERAINT_1989_2012.nc'))
    eastward_wind_850hPa_obs_cube.rename('x_wind')
    eastward_wind_850hPa_obs_cube.attributes['MODEL_RUN_ID'] = 'OBS'
    p_coord = iris.coords.DimCoord([850], long_name='pressure', units='hPa')
    eastward_wind_850hPa_obs_cube.add_aux_coord(p_coord)

    eastward_wind_200hPa_obs_cube = iris.load_cube(os.path.join(obs_dir, 'u200_daily_ERAINT_1989_2012.nc'))
    eastward_wind_200hPa_obs_cube.rename('x_wind')
    eastward_wind_200hPa_obs_cube.attributes['MODEL_RUN_ID'] = 'OBS'
    p_coord = iris.coords.DimCoord([200], long_name='pressure', units='hPa')
    eastward_wind_200hPa_obs_cube.add_aux_coord(p_coord)

    # TODO this is part of the pre-processing; do outside assessment

    # standardise units
    if precip_obs_cube.units == 'mm/day':
        precip_obs_cube.units = 'kg m-2 day-1'
    precip_obs_cube.convert_units('kg m-2 day-1')
    olr_obs_cube.convert_units('W m-2')
    eastward_wind_850hPa_obs_cube.convert_units('m s-1')
    eastward_wind_200hPa_obs_cube.convert_units('m s-1')

    # assign coord system (coord systems must be equal for regridding; or both
    # None)
    obs_cubes = iris.cube.CubeList(
        [precip_obs_cube, olr_obs_cube, eastward_wind_850hPa_obs_cube,
         eastward_wind_200hPa_obs_cube]
    )

    coord_sys = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
    _obs_cubes = []
    for cube in obs_cubes:
        if cube.coord_system() is None:
            cube.coord('latitude').coord_system = coord_sys
            cube.coord('longitude').coord_system = coord_sys
        _obs_cubes.append(cube)
    obs_cubes = _obs_cubes
    del _obs_cubes

    cubes.extend(obs_cubes)


    ### Assessment
    cubes = extract_assessment_period(cubes, start_date, end_date)
    cubes = extract_tropics(cubes)
    cubes = regrid_to_low_resolution(cubes)

    # assess single runs
    for suite_id in suite_ids:
        metrics = {}
        out_dir = os.path.join(AREA_OUT_DIR, suite_id)

        precip_cube, outgoing_longwave_cube, u_wind_200_cube, u_wind_850_cube = \
            None, None, None, None

        precip_cube = cubes.extract_strict(
            iris.Constraint('precipitation_flux') &
            iris.AttributeConstraint(MODEL_RUN_ID=suite_id)
        )
        precip_cube.convert_units('kg m-2 day-1')

        outgoing_longwave_cube = cubes.extract_strict(
            iris.Constraint('toa_outgoing_longwave_flux') &
            iris.AttributeConstraint(MODEL_RUN_ID=suite_id)
        )
        u_wind_200_cube = cubes.extract_strict(
            iris.Constraint('x_wind', pressure=200) &
            iris.AttributeConstraint(MODEL_RUN_ID=suite_id)
        )
        u_wind_850_cube = cubes.extract_strict(
            iris.Constraint('x_wind', pressure=850) &
            iris.AttributeConstraint(MODEL_RUN_ID=suite_id)
        )

        msg = ' '.join(
            ['No data for some cubes between', start_date.strftime('%Y-%m-%d'),
             end_date.strftime('%Y-%m-%d'), 'suite ID:', suite_id]
        )
        assert all([precip_cube, outgoing_longwave_cube, u_wind_200_cube,
                    u_wind_850_cube]), msg

        # Level 1 diagnostics
        # Mean, variance, filtered variance, filt variance/total variance
        for cube in [precip_cube, outgoing_longwave_cube, u_wind_200_cube,
                     u_wind_850_cube]:

            diags_level1.diagnos_level1(cube, out_dir, suite_id, tmp_dir)

        # Level 2 diagnostics
        # WK raw sym/antisym spectra, background spectra, (sym, antisym)/background
        for cube in [precip_cube, outgoing_longwave_cube, u_wind_200_cube,
                     u_wind_850_cube]:

            level2_metrics = diags_level2.diagnos_level2(cube, out_dir, suite_id)
            metrics.update(level2_metrics)

        # Level 3 diagnostics
        # Real-time multivariate MJO Index (RMM) calculations, and Summer/Winter
        # composites
        level3_metrics = diags_level3.diagnos_level3(
            outgoing_longwave_cube, u_wind_850_cube, u_wind_200_cube,
            precip_cube, suite_id, out_dir
        )
        metrics.update(level3_metrics)

        # write metrics to csv file for each suite_id
        with open(os.path.join(out_dir, 'metrics.csv'), 'w') as fh:
            writer = csv.writer(fh)
            for metric in metrics.items():
                writer.writerow(metric)


def extract_tropics(cubes):
    "Extract tropics region from cubes."
    # region (Tropics)
    lon = (0, 360)
    lat = (-30, 30)

    # Extract region
    region_constraint = iris.Constraint(
        longitude=lambda cell: lon[0] <= cell <= lon[1],
        latitude=lambda cell: lat[0] <= cell <= lat[1]
    )
    return cubes.extract(region_constraint)


def extract_assessment_period(cubes, start_date, end_date):
    "Extract assessment period from cubes."
    start = iris.Constraint(time=lambda cell: cell.point >= start_date)
    end = iris.Constraint(time=lambda cell: cell.point < end_date)
    with iris.FUTURE.context(cell_datetime_objects=True):
        cubes = cubes.extract(start & end)
    return cubes


def regrid_to_low_resolution(cubes):
    "Regrid cubes to regular lat-long 2.5 degree grid"
    ref_cube = reg_lat_long_grid.create_cube(latitudes=(-30, 30),
                                             longitudes=(0, 360), spacing=2.5)

    regridded_cubes = iris.cube.CubeList()
    for cube in cubes:
        if cube.coord('longitude').bounds is None:
            cube.coord('longitude').guess_bounds()
        if cube.coord('latitude').bounds is None:
            cube.coord('latitude').guess_bounds()
        # TODO is linear regridding sufficient?
        regridded_cubes.append(cube.regrid(ref_cube, iris.analysis.Linear()))
    return regridded_cubes


if __name__ == '__main__':
    main()
