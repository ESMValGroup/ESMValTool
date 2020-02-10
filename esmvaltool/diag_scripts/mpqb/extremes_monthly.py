#!/usr/bin/env python
import copy
import datetime
import logging
import os
import pandas as pd
import re
import warnings
from pprint import pformat

import esmvalcore.preprocessor as pp
import yaml
import skill_metrics as sm
import matplotlib.pyplot as plt
import iris
import numpy as np
from scipy.stats import pearsonr
from collections import OrderedDict
from iris.time import PartialDateTime

from esmvaltool.diag_scripts.shared.trend_mpqb_common.diag1d import rmsd1d
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import (ProvenanceLogger,
                                                  get_diagnostic_filename,
                                                  get_plot_filename)
from mpqb_plots import get_ecv_plot_config, mpqb_mapplot
from esmvaltool.diag_scripts.shared.trend_mpqb_common.sharedutils import parallel_apply_along_axis

logger = logging.getLogger(os.path.basename(__file__))

dataset_plotnames = {
  'ERA-Interim-Land' : 'ERA-Interim/Land',
  'CDS-SATELLITE-SOIL-MOISTURE' : 'ESA-CCI',
  'cds-era5-land-monthly' : 'ERA5-Land',
  'cds-era5-monthly' : 'ERA5',
  'MERRA2' : 'MERRA-2',
  'cds-satellite-lai-fapar' : 'SPOT-VGT',
}

def convert_human_readable_coords_to_iso(coord_in):
    '''
    # This function converts human readable single coords of latitude or longitude 
    # to iso-6709 standard
    
    parameters
    ----------
       
       coord_in : str
        the input coordinate (e.g. 53°N)
        
    returns
    -------
       coord_iso : float
        iso-6709 formatted coordinate
    '''
    val,compass = re.split('[deg]+', coord_in)
    if compass in ['N','E']:
        result = float(val)
    elif compass in ['S','W']:
        result = float(val)*-1.
    else:
        result = float(val)
    return result
 
def write_table_to_yaml(ex_table, filename):
    # First for regions
    name_mapping = {
        'Lat_from': 'start_latitude',
        'Lat_to': 'end_latitude',
        'Lon_from': 'start_longitude',
        'Lon_to': 'end_longitude',
    }
    ex_table = ex_table.rename(columns=name_mapping)
    extract_region = ex_table[name_mapping.values()]
    extract_region = extract_region.to_dict('index')
    for key in extract_region:
        extract_region[key] = {'extract_region' : extract_region[key]}
    # Now same for times
    ex_table['Time_start'].apply(lambda row: {'start_month' : row.month, 'start_year' : row.year})
    


    with open(filename, 'w') as handle:
       print(f"Writing to disk: {filename}")
       yaml.safe_dump(preprocessor_instructions, handle)

def read_extreme_event_catalogue():
    ''' This function reads the extreme event catalogue
    Returns
    -------
    (pandas.DataFrame,pandas.DataFrame)
         the extreme event catalogue (parsed), the extreme event catalogue (raw)
    '''
    ex_table_dir = './predef/extremes_catalogue/'
    #TODO move ex_table_file to the recipe (default value that can
    # be over written if specified)
    ex_table_file = 'latest_accessed20200206.csv'
    ex_table_loc = os.path.join(os.path.dirname(__file__),ex_table_dir,ex_table_file)

    ############################################
    ####  Reading in extreme event catalogue ###
    ############################################

    # First read complete table and parse into right data types
    ex_table = pd.read_csv(ex_table_loc,skiprows=1)
    raw_table = copy.deepcopy(ex_table)

    # Parameters related to formatting of the table
    fmt_dates = '%d/%m/%Y'
    
    # Convert the two datetime columns to datetime.datetime objects
    for colname in ['Time_start','Time_stop']:
        ex_table[colname] = ex_table[colname].apply(lambda x:\
                            datetime.datetime.strptime(x,fmt_dates))
    
    # Now split latitude start-end into two seperate columns
    lats = ex_table['Latitude extent (start, stop)'].str.split("-",n=1,expand=True)
    lons = ex_table['Longitude extent (start, stop)'].str.split("-",n=1,expand=True)
    
    ex_table['Lat_from'] = lats[0]
    ex_table['Lat_to'] = lats[1]
    
    ex_table['Lon_from'] = lons[0]
    ex_table['Lon_to'] = lons[1]
    
    # Now convert them to iso-6709
    for coord_key in ['Lat_from','Lat_to','Lon_from','Lon_to']:
        ex_table[coord_key] = ex_table[coord_key].apply(lambda x: convert_human_readable_coords_to_iso(x))

    # Now create event_id
    ex_table['extreme_event_id'] = ex_table[['Event Category','Region/Country', 'Year']].apply(lambda x: '_'.join(x),axis=1)
    # Now strip the whitespace
    ex_table['extreme_event_id'] = ex_table['extreme_event_id'].apply(lambda x: ''.join(x.split()))
    
    # Here we define the information that needs to end up in the dictionary 
    # for each event for further processing the event
    core_keys = ['extreme_event_id','Literature support',\
                 'Lat_from','Lat_to','Lon_from','Lon_to',\
                 'Time_start','Time_stop',\
                 'Estimated duration (1-3) [used for selection]',\
                 'Estimated scale  (1-3) [used for selection]',\
                 'Estimated impact (1-3) [used for selection]']
    ex_table = ex_table[core_keys]
    
    # Set extreme_event_id as table index
    ex_table.index = ex_table.pop('extreme_event_id')
    return ex_table,raw_table



def _extract_event(cube,ex_table,event_name):
    event = ex_table.loc[event_name]
    cube = pp.extract_region(cube,
                             event['Lon_from'],
                             event['Lon_to'],
                             event['Lat_from'],
                             event['Lat_to'])
    cube = pp.extract_time(cube,
                           event['Time_start'].year,
                           event['Time_start'].month,
                           event['Time_start'].day,
                           event['Time_stop'].year,
                           event['Time_stop'].month,
                           event['Time_stop'].day)
    return cube


def main(cfg):

    # Read all datasets that are provided.
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(
        input_data, 'dataset', sort='dataset')
    logger.info(
        "Starting MPQB extremes_monthly script."
    )

    ex_table = read_extreme_event_catalogue()[0]
#    event_names =  ['Drought-Heat_Europe_2013'] 
    event_names = [name for name in ex_table.index if ('Drought' in name)]

    # Initialize a pandas dataframe for saving the table of metrics
    df_metrics = pd.DataFrame(columns=[dataset_plotnames[key] for key in grouped_input_data.keys()],
                              index=event_names, dtype=float)

    aggregator = cfg.pop('aggregator', 'mean')
    aggregator = 'min' # TODO put back

    for dataset in grouped_input_data.keys():
        dataset_cfg = grouped_input_data[dataset]
        logger.info("Opening dataset: %s", dataset)

        # Opening the data
        cube = iris.load_cube(dataset_cfg[0]['filename'])
        for event_name in event_names:
#            try:
            if True:
                if dataset=='MERRA2':
                    cube.coord('time').bounds = None
                    logger.info(f"Guessing bounds for {dataset}")
                    cube.coord('time').guess_bounds()
                event_cube = _extract_event(cube, ex_table, event_name)
                event_cube = pp.area_statistics(event_cube, 'mean')
                logger.info(f"{event_cube.shape[0]} months for {event_name}")
                event_cube = pp.climate_statistics(event_cube, aggregator)
                meanvalue = float(event_cube.data)
                df_metrics.loc[event_name, dataset_plotnames[dataset]] = f"{meanvalue:0.2f}"
#            except ValueError:
#                logger.info(f"{event_name} not included in dataset (is the event too short?)")

    # Save to html with specified precision
    savename_html = os.path.join(cfg['plot_dir'],'event_table.html')
    html_metrics = df_metrics.style.set_precision(2).render()
    logger.info(
        "Saving metric html-table as: {0}".format(savename_html))
    with open(savename_html, mode='w+') as handle:
        handle.write(html_metrics)

    if cfg['write_plots']:
         print("No plotting implemented in this diagnostic")

    logger.info("Finished!")





if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
