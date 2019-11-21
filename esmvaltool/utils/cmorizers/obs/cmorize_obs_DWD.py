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
import dask.dataframe as dd
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
    return cube

def _fix_units(units):
    """Specific fixes for different units."""
    logger.info("Fixing unit ...")
    if units == "°C":
        return "deg_C"

def modify_cube(cube, var_info, raw_info, out_dir, attrs):
    """process to cube."""
    var = var_info.short_name
    rawvar = raw_info['name']

    if cube.var_name == rawvar:
        fix_var_metadata(cube, var_info)
        fix_coords(cube)
        _fix_data(cube, var)
        set_global_atts(cube, attrs)
        logger.info(out_dir)
        save_variable(
            cube, var, out_dir, attrs, unlimited_dimensions=['time'])

def make_cube(vals, time, longitude, latitude,
              units = None, name = None, station_ID = None):
    """Produce cube from dataframe selections"""
#    apply(lambda x: datetime.strptime(x,'%Y%m%d%H%M'))
    time = np.array([datetime.strptime(str(ti),'%Y%m%d%H%M') 
        for ti in time.flatten().compute()])
    t_unit = cf_units.Unit('days since 1750-01-01 00:00:00',
                           calendar=cf_units.CALENDAR_STANDARD)
    time = iris.coords.DimCoord(t_unit.date2num(time),
                                standard_name = "time",
                                units = t_unit)
    
    station_ID = iris.coords.DimCoord(np.array([station_ID]),
                                      long_name = "SDO_ID")
    
    longitude = iris.coords.DimCoord(np.array([longitude]),
                                      standard_name = "longitude")
    
    latitude = iris.coords.DimCoord(np.array([latitude]),
                                      standard_name = "latitude")
    
    cube = iris.cube.Cube(
            vals,
            units=units,
            var_name = name['varname'],
            long_name = name['longname'],
            dim_coords_and_dims=[(time, 0),(station_ID, 1)],
            aux_coords_and_dims=[(latitude, 1),(longitude, 1)],
            )
    
    logger.info("cube prepared")
    
    return cube

def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    
    cmor_table = CFG['cmor_table']
    glob_attrs = CFG['attributes']
    glob_vars = CFG['variables']
    
    save_version = glob_attrs['version']
    
    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    data = dd.read_csv(in_dir + os.sep + 'data_*.csv',
                       dtype={'Zeitstempel': object})
    sdo = dd.read_csv(in_dir + os.sep + 'sdo_*.csv')
    prd = dd.read_csv(in_dir + os.sep + 'prd.csv')
    logger.info("Output will be written to: %s", out_dir)
    
    for var,vals in glob_vars.items():
        
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
            
            # convert to dask arrays for cube
            values = extracted_data[['Wert']].to_dask_array(lengths=True)
            time = extracted_data[['Zeitstempel']].to_dask_array(lengths=True)

            # additional info
            var_info = cmor_table.get_variable(vals['mip'], var)
            raw_info = {'name': vals['raw'], 'file': link}
            
            # produce cube
            cube_name = {'varname': vals['raw'], 'longname': name}
            cube = make_cube(values, time, longitude, latitude,
                             units = units, name = cube_name,
                             station_ID = SDO_ID)
            
            # global attributes
            glob_attrs['mip'] = vals['mip']
            glob_attrs['station'] = name 
            glob_attrs['elevation'] = height
            glob_attrs['description'] = " | ".join(description)
            glob_attrs['version'] = "-".join([save_version,
                      name.replace(os.sep,"-")])
            
            # modifying and saving cube
            modify_cube(cube, var_info, raw_info, out_dir, glob_attrs)
