"""Python example diagnostic."""
import logging
import os
from pprint import pformat

import iris
import numpy as np
import xarray as xr

import surge_estimator
from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
from esmvaltool.diag_scripts.shared.plot import quickplot
from surge_estimator import surge_estimator_main

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Compute the surge time series for each input dataset."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(input_data, 'dataset')

    # Example of how to loop over variables/datasets in alphabetical order
    for dataset in grouped_input_data:
        logger.info("Processing dataset %s", dataset)
        # read data
        uas_info = select_metadata(
            grouped_input_data[dataset], short_name='ua')[0]
        uas_file = uas_info['filename']
        nc_ua = xr.open_dataset(uas_file)
        uas_data = nc_ua.uas
        vas_info = select_metadata(
            grouped_input_data[dataset], short_name='va')[0]
        vas_file = vas_info['filename']
        nc_va = xr.open_dataset(vas_file)
        vas_data = nc_va.vas
        psl_info = select_metadata(
            grouped_input_data[dataset], short_name='psl')[0]
        psl_file = psl_info['filename']
        nc_psl = xr.open_dataset(psl_file)
        psl_data = nc_psl.psl
        # call surge estimator
        logger.info('Calling surge estimator')
        srg_estim,dates = surge_estimator_main(psl_data, uas_data, vas_data, cfg, dataset)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
