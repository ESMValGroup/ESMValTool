"""Python example diagnostic."""
import logging
import os
from pprint import pformat
import xarray as xr
import iris
import numpy as np 

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))

import surge_estimator
from surge_estimator import surge_estimator_main


def main(cfg):
    """Compute the surge time series for each input dataset."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(
        input_data, 'dataset')

    # Example of how to loop over variables/datasets in alphabetical order
    for dataset in grouped_input_data:
        logger.info("Processing dataset %s", dataset)
        # read data
        ua_info = select_metadata(grouped_input_data[dataset], short_name='ua')[0]
	ua_file = ua_info['filename']
        nc_ua   = xr.open_dataset(ua_file)
        ua_data = nc_ua.ua
        va_info = select_metadata(grouped_input_data[dataset], short_name='va')[0]
	va_file = va_info['filename']
        nc_va   = xr.open_dataset(va_file)
        va_data = nc_va.va
        psl_info = select_metadata(grouped_input_data[dataset], short_name='psl')[0]
	psl_file = psl_info['filename']
        nc_psl   = xr.open_dataset(psl_file)
        psl_data = nc_psl.psl
        # call surge estimator
        logger.info('Calling surge estimator')
        surge_estimator_main(psl_data, ua_data, va_data, cfg, dataset)
       

if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
