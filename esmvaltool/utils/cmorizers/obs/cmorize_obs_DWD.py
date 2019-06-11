# pylint: disable=invalid-name
"""ESMValTool CMORizer for DWD CDC data.

Tier
   Tier 2: other freely-available dataset.

Source
   https://cdc.dwd.de/portal/

Last access
   20190605

Download and processing instructions
   Download any files:
      - tested on temperature/2m/hourly

Modification history
   20190605-A_muel_bn: written.

"""

import logging
import os

import iris
import dask.array as da
import dask.dataframe as dd
import pandas as pd
import numpy as np
import cf_units
from datetime import datetime

from .utilities import (constant_metadata, convert_timeunits, fix_coords,
                        fix_var_metadata, read_cmor_config, save_variable,
                        set_global_atts)

logger = logging.getLogger(__name__)

# read in CMOR configuration
CFG = read_cmor_config('DWD.yml')


def _fix_data(cube, var):
    """Specific data fixes for different variables."""
    logger.info("Fixing data ...")
    with constant_metadata(cube):
        mll_to_mol = ['po4', 'si', 'no3']
        if var in mll_to_mol:
            cube /= 1000.  # Convert from ml/l to mol/m^3
        elif var == 'thetao':
            cube += 273.15  # Convert to Kelvin
        elif var == 'o2':
            cube *= 44.661 / 1000.  # Convert from ml/l to mol/m^3
    return cube

def _fix_units(units):
    """Specific fixes for different units."""
    logger.info("Fixing unit ...")
    if units == "°C":
        return "deg_C"

def modify_cube(cube, var_info, raw_info, out_dir, attrs):
    """Extract to all vars."""
    logger.info(var_info)
    var = var_info.short_name
    rawvar = raw_info['name']

    if cube.var_name == rawvar:
        fix_var_metadata(cube, var_info)
#            convert_timeunits(cube, year)
        fix_coords(cube)
#            _fix_data(cube, var)
        set_global_atts(cube, attrs)
        save_variable(
            cube, var, out_dir, attrs, unlimited_dimensions=['time'])

def make_cube(vals, time, longitude, latitude,
              units = None, name = None, station_ID = None):
    """Produce cube from dataframe selections"""
    logger.info(time.flatten())
    logger.info(time.flatten().compute()[0])
#    t_unit = cf_units.Unit('{} since 1850-01-01 00:00:00'.format(timestep), calendar=cf_units.CALENDAR_STANDARD)
#    time = iris.coords.DimCoord(t_unit.date2num(time.flatten()),
#                                standard_name = "time")
    station_ID = iris.coords.DimCoord(np.array([station_ID]),
                                      standard_name = "station_ID")
    
    cube = iris.cube.Cube(
            vals,
            units=units,
            var_name = name['varname'],
            long_name = name['longname'],
            attributes=None,
            cell_methods=None,
            dim_coords_and_dims=[(time, 0),(station_ID, 1)],
            aux_coords_and_dims=None,
            aux_factories=None
            )
    logging.info(cube)
    
    logging.info("cube prepared")
    return cube

def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    
    cmor_table = CFG['cmor_table']
    glob_attrs = CFG['attributes']
    glob_vars = CFG['variables']
    
    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    data = dd.read_csv(in_dir + os.sep + 'data_*.csv',
                       dtype={'Zeitstempel': object})
    sdo = dd.read_csv(in_dir + os.sep + 'sdo_*.csv')
    prd = dd.read_csv(in_dir + os.sep + 'prd.csv')
    logger.info("Output will be written to: %s", out_dir)
    
    for var,vals in glob_vars.items():
        
        logger.info(prd.head())
        sub_prd = prd[prd['Produkt_Code'] == vals['raw']]
        units = (sub_prd['Einheit']).values.compute()
        description = sub_prd[
                ['Beschreibung_DWD','Beschreibung_Inspire']
                ].values.compute()[0]
        units = _fix_units(units)
        
        for _, [SDO_ID, name, longitude, latitude, height, link] in sdo.iterrows():
            logger.info("Processing variable %s for station %s", var, name)
            
            # extract data
            extracted_data = data[
                    (data['SDO_ID'] == SDO_ID) &
                    (data['Produkt_Code'] == vals['raw'])
                    ] # on might consider further selection based on Qualitaet_Niveau
            
            # convert time column              
            logger.info(extracted_data['Zeitstempel'].head())
            
            extracted_data['Zeitstempel'] = \
                extracted_data['Zeitstempel'].apply(
                        lambda x: datetime.strptime(x,'%Y%m%d%H%M'))
            
            # convert to dask arrays for cube
            values = extracted_data[['Wert']].to_dask_array(lengths=True)
            time = extracted_data[['Zeitstempel']].to_dask_array(lengths=True)

            # additional info
            var_info = cmor_table.get_variable(vals['mip'], var)
            raw_info = {'name': vals['raw'], 'file': link}
            glob_attrs['mip'] = vals['mip']
            glob_attrs['station'] = name 
            glob_attrs['elevation'] = height
            glob_attrs['description'] = description
            
            # produce cube
            cube = make_cube(values, time, longitude, latitude,
                             units = units, name = vals['raw'],
                             station_ID = SDO_ID)
            modify_cube(cube, var_info, raw_info, out_dir, glob_attrs)

#    logger.info("Input data from: %s", in_dir)
#    logger.info("Output will be written to: %s", out_dir)
#
#    # run the cmorization
#    for var, vals in CFG['variables'].items():
#        yr = None
#        for yr in CFG['custom']['years']:
#            file_suffix = str(yr)[-2:] + '_' + str(yr + 1)[-2:] + '.nc'
#            inpfile = os.path.join(in_dir, vals['file'] + file_suffix)
#            logger.info("CMORizing var %s from file %s", var, inpfile)
#            var_info = cmor_table.get_variable(vals['mip'], var)
#            raw_info = {'name': vals['raw'], 'file': inpfile}
#            glob_attrs['mip'] = vals['mip']
#            extract_variable(var_info, raw_info, out_dir, glob_attrs, yr)
