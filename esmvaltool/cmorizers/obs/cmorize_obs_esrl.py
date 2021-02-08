"""ESMValTool CMORizer for ESRL data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://www.esrl.noaa.gov/gmd/dv/data/
    ftp://aftp.cmdl.noaa.gov/data/trace_gases/co2/


Last access
    20201126

Download and processing instructions
    Download the txt files from the ftp server/data interface

"""

import logging
import os
from datetime import datetime
from pprint import pformat
import glob
from ftplib import FTP

import iris
import numpy as np
import pandas as pd
import requests
from cf_units import Unit

from . import utilities as utils

logger = logging.getLogger(__name__)


def _download_files(in_dir, cfg, stations):
    """Download input files using FTP."""
    logger.info("Downloading data from FTP server %s", cfg['ftp_host'])
    files = {}
    for station in stations:
        # First look for baseline observatories
        if station.upper() in ['MLO', 'BRW', 'SMO', 'SPO']:
            files[station] = {'name': "co2_" + station.lower()
                                      + '_surface-insitu_1_ccgg_'
                                        'MonthlyData.txt',
                              'folder': cfg['data_dir'] + 'in-situ/surface/'
                                        + station.lower()}
        elif station.lower() == 'global':
            files[station] = {'name': 'co2_mm_gl.txt',
                              'folder': 'products/trends/co2/'}
        else:
            files[station] = {'name': 'co2_' + station.lower()
                                      + '_surface-flask_1_ccgg_month.txt',
                              'folder': cfg['data_dir'] + 'flask/surface/'}
    input_files = {}
    rm_stat = []
    with FTP(cfg['ftp_host']) as ftp_client:
        logger.info(ftp_client.getwelcome())
        ftp_client.login()
        for station in files:
            filename_full = os.path.join(files[station]["folder"],
                                         files[station]["name"])
            if filename_full in ftp_client.nlst(files[station]["folder"]):
                logger.info("Downloading %s", files[station]["name"])
                new_path = os.path.join(in_dir, files[station]["name"])
                with open(new_path, mode='wb') as outfile:
                    ftp_client.retrbinary(f'RETR {filename_full}',
                                          outfile.write)
                input_files[station] = [new_path]
            else:
                rm_stat.append(station)
    return input_files, rm_stat


def _get_cube(row, column_ind, fill_value, station_dict):
    """Create :class:`iris.cube.Cube` from :class:`pandas.Series`."""
    time_coord = _get_time_coord(int(row['year']),
                                 int(row['month']))
    lat_coord, lon_coord = _make_station_lat_lon_coord(station_dict)
    data = np.ma.masked_equal(float(row[column_ind[2]]), fill_value)
    cube = iris.cube.Cube(
        data.reshape((1, 1, 1)),
        dim_coords_and_dims=[(time_coord, 0), (lat_coord, 1), (lon_coord, 2)],
        units='ppm',
    )
    return cube


def _get_rows_and_fill_value(filepath):
    """Check which dataset type is present and return columns to use."""
    if 'insitu' in filepath:
        # Insitu tower monthly
        data_rows = [1, 2, 8]
        fill_v = -999.990
    elif 'month.' in filepath:
        # Monthly surface flask data, 1; year, 2: month, 3: data
        data_rows = [1, 2, 3]
        fill_v = -999.99  # not sure
    elif 'mm_gl' in filepath:
        data_rows = [0, 1, 3]
        fill_v = -999.99
    else:
        raise NotImplementedError("Unexpected number of columns, "
                                  "only monthly data from in situ or flask "
                                  "measurements currently supported")
    return data_rows, fill_v


def _get_station_dictionary():
    """Get station information from online table."""
    url = "https://www.esrl.noaa.gov/gmd/dv/site/?program=ccgg"
    stat_list = pd.read_html(requests.get(url).content)
    stats = stat_list[-1]
    # Remove asterisk from station names (flags inactive stations)
    stats['Code'] = stats['Code'].str.replace('*', '')
    stats.set_index("Code", drop=False, inplace=True)
    station_dict = stats.to_dict(orient="index")

    # Add entry for Global
    station_dict['GLOBAL'] = {'Latitude': 0.0, 'Longitude': 180.0,
                              'Elevation (meters)': 0, 'Code': 'GLOBAL'}
    return station_dict


def _get_time_coord(year, month):
    """Get time coordinate."""
    point = datetime(year=year, month=month, day=15)
    bound_low = datetime(year=year, month=month, day=1)
    if month == 12:
        month_bound_up = 1
        year_bound_up = year + 1
    else:
        month_bound_up = month + 1
        year_bound_up = year
    bound_up = datetime(year=year_bound_up, month=month_bound_up, day=1)
    time_units = Unit('days since 1950-01-01 00:00:00', calendar='standard')
    time_coord = iris.coords.DimCoord(
        time_units.date2num(point),
        bounds=time_units.date2num([bound_low, bound_up]),
        var_name='time',
        standard_name='time',
        long_name='time',
        units=time_units,
    )
    return time_coord


def _extract_variable(short_name, var, cfg, out_dir, station_dic):
    """Extract variable."""
    data = pd.read_csv(station_dic['filepath'], sep=' {1,}', comment='#',
                       engine='python', header=None)
    # Insitu tower monthly had uncommented header, remove
    if data.shape[1] == 17:
        data = data.drop(0)

    data_rows, fill_v = _get_rows_and_fill_value(station_dic['filepath'])

    # Resample data to monthly, pad with missing values as needed
    data[data_rows[2]] = pd.to_numeric(data[data_rows[2]])
    data = data.replace(fill_v, np.nan)
    data = data.rename(columns={data_rows[0]: "year", data_rows[1]: "month"})
    data['day'] = 15
    data['datetime'] = pd.to_datetime(data[['year', 'month', 'day']])
    data = data.resample('M', on="datetime").mean()
    data = data.fillna(fill_v)
    data['year'] = data.index.year
    data['month'] = data.index.month

    # Extract cube
    cubes = iris.cube.CubeList()
    for (_, row) in data.iterrows():
        cube = _get_cube(row, data_rows, fill_v, station_dic)
        cubes.append(cube)
    cube = cubes.concatenate_cube()
    cube.var_name = short_name

    # Fix metadata
    utils.convert_timeunits(cube, 1950)
    utils.fix_coords(cube)
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    cube.convert_units(cmor_info.units)
    attrs = cfg['attributes']
    attrs['version'] = station_dic['Code'].upper()
    attrs['mip'] = var['mip']
    attrs['altitude'] = station_dic['Elevation (meters)']
    attrs['altitude_units'] = 'm'
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def _make_station_lat_lon_coord(station_dic):
    """Make iris coordinates for given latitude and longitude."""
    lat = station_dic['Latitude']
    lon = station_dic['Longitude']
    if lon < 0:
        lon = lon + 360

    # Treat Global data differently
    if lat == 0.0 and lon == 180.0:
        lat_coord = iris.coords.DimCoord([0.0], bounds=[[-90.0, 90.0]],
                                         var_name='lat',
                                         standard_name='latitude',
                                         long_name='latitude',
                                         units=Unit('degrees_north'))
        lon_coord = iris.coords.DimCoord([180.0],
                                         bounds=[[0.0, 360.0]],
                                         var_name='lon',
                                         standard_name='longitude',
                                         long_name='longitude',
                                         units=Unit('degrees_east'))
    else:
        lat_coord = iris.coords.DimCoord([lat], var_name='lat',
                                         standard_name='latitude',
                                         long_name='latitude', units='degrees')
        lon_coord = iris.coords.DimCoord([lon], var_name='lon',
                                         standard_name='longitude',
                                         long_name='longitude',
                                         units='degrees')
    return lat_coord, lon_coord


def _get_filenames(stations, cfg, in_dir, all_stat):
    """Get filename given pattern and station name."""
    input_files = {}
    download_files = []
    for station in stations:
        if station.lower() == "global":
            st_filepattern = "co2_mm_gl.txt"
        else:
            # Replace first * with station name
            filename_pattern = cfg['input_filename_pattern']
            st_filepattern = filename_pattern.replace("*", station.lower(), 1)
        pattern = os.path.join(in_dir, st_filepattern)
        input_file = glob.glob(pattern)
        if not input_file:
            download_files.append(station)
        else:
            input_files[station] = input_file
    if len(download_files) > 0:
        if cfg['download']:
            input_files_dl, rm_stat = _download_files(in_dir, cfg,
                                                      download_files)
            input_files.update(input_files_dl)
            if len(rm_stat) > 0:
                if all_stat:
                    # When selecting "all", some stations may not have
                    # available data at the moment,
                    # so remove these from to process files
                    stations = [x for x in stations if x not in rm_stat]
                else:
                    raise ValueError("No data found for %s on the ftp server. "
                                     % rm_stat)
        else:
            if not all_stat:
                raise ValueError("No local data found for stations %s, "
                                 "consider turning on the download option."
                                 % download_files)
    logger.debug("Found input files:\n%s", pformat(input_files))
    return input_files, stations


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    # read station information
    station_dict = _get_station_dictionary()

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        # Read station names
        if 'all' in var['stations']:
            stations = station_dict.keys()
            all_stat = True
        else:
            stations = var['stations']
            all_stat = False
        # Check for wrong station names
        stat_upper = [element.upper() for element in stations]
        false_keys = np.setdiff1d(stat_upper, list(station_dict.keys()))
        if len(false_keys) == 0:
            filepath, stations = _get_filenames(stations, cfg, in_dir,
                                                all_stat)
            for station in stations:
                logger.info("Reading file '%s'", filepath[station][0])
                logger.info("CMORizing variable '%s' for station '%s'",
                            short_name, station)
                # Add filepath to station_dict
                station_dict[station.upper()]['filepath'] = \
                    filepath[station][0]
                _extract_variable(short_name, var, cfg, out_dir,
                                  station_dict[station.upper()])
        else:
            raise ValueError("Could not find the following station(s): %s. "
                             "Please double-check your spelling in the "
                             "cmor config file. The following is a list of "
                             "valid stations: %s."
                             % (false_keys, list(station_dict.keys())))
