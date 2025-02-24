"""Calculate Global Warming Level exceedance years."""
import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import iris 

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
    select_metadata,
    sorted_metadata
)

logger = logging.getLogger(Path(__file__).stem)


def log_provenance(filename, ancestors, caption, cfg):
    """Create a provenance record for the output file."""
    provenance = {
        'caption': caption,
        'domains': ['global'],
        'authors': ['swaminathan_ranjini'],
        'projects': ['ukesm'],
        'ancestors': ancestors,
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance)

def calculate_moving_average_centered(a, window_size):
    '''Calculate rolling means centred at a point'''
    
    if(window_size < 2):
        raise ValueError("Window size should be greater than 2")
    
    return pd.Series(a).rolling(window=window_size, center=True).mean().to_numpy()

def calculate_gwl_exceedance_years(input_data_sorted, gwls, window_size):
    '''Calculate GWL exceedance years for models and ensemble members'''
    
    col_names = ['Project', 'Exp', 'Model','Ens', 'GWL', 'Exceedance_Year']
    gwl_df = pd.DataFrame(columns = col_names)

    for data in input_data_sorted:
        project = data['project']
        exp = data['exp']
        dataset = data['dataset']
        ensemble = data['ensemble']
        filename = data['filename']
        logger.info('Processing %s %s %s %s ', project, exp, dataset)
        anomaly_ts_cube = iris.load_cube(data['filename'])
        years = anomaly_ts_cube.coord('year').points
        start_year = years[0]
        end_year = years[len(years) - 1]
        logger.info('Start and end years and averaging over  %s %s %s ', str(start_year), str(end_year), str(window_size))
        smoothed_ts = calculate_moving_average_centered(anomaly_ts_cube.data, window_size)

        for gwl in gwls:
            logger.info("GWL temperature %s", float(gwl))
            logger.info('Size of array %s ', str(smoothed_ts.shape))
            #only write exceedance_year if gwl is exceeded in the time series
            if(np.nanmax(smoothed_ts) > float(gwl)):                
                exceedance_year =  start_year + np.argmax(smoothed_ts > float(gwl))
                new_record_df = pd.DataFrame([{'Project': data['project'], 'Exp': data['exp'], 'Model': data['dataset'], 'Ens': data['ensemble'], 'GWL': gwl, 'Exceedance_Year': exceedance_year}])
                logger.info('Exceedance year %s ', exceedance_year)
                gwl_df = pd.concat([gwl_df, new_record_df], ignore_index = True)
    
    return gwl_df


def main(cfg):
    input_data = cfg['input_data'].values()

    #initialize pandas data frame to store GWL exceedance information
    window_size  = cfg['window_size']
    gwls = cfg['gwls']

  
    # group preproocessed input by project
    input_data_sorted = sorted_metadata(input_data, ['project','exp', 'dataset','ensemble'])

    # calculate GWL exceedance years and return in a dataframe
    gwl_df = calculate_gwl_exceedance_years(input_data_sorted, gwls, window_size)
    gwl_file = os.path.join(cfg['work_dir'], 'GWL_exceedance_years.csv')
    gwl_df.to_csv(gwl_file, sep = ',', index = False)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)