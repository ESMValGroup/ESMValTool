"""
Look at this module for guidance how to write your own.

Read the README_PERSONAL_DIAGNOSTIC file associated with this example;

Module for personal diagnostics (example).
Internal imports from exmvaltool work e.g.:

from esmvalcore.preprocessor import regrid
from esmvaltool.diag_scripts.shared.supermeans import get_supermean

Pipe output through logger;

Please consult the documentation for help with esmvaltool's functionalities
and best coding practices.
"""
# place your module imports here:
# import gsw
import xarray as xr
import numpy as np

# import cmocean

import logging


# operating system manipulations (e.g. path constructions)
import os
import sys
from pathlib import Path

# to manipulate iris cubes
import iris
import matplotlib.pyplot as plt
from esmvalcore.preprocessor import area_statistics

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic, ProvenanceLogger

# reuse tools already developed for ocean diagnostics
# from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools

# This part sends debug statements to stdout
# logger = logging.getLogger(os.path.basename(__file__))
# logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

# This outputs to the run/transect/script1/log.txt file
logger = logging.getLogger(Path(__file__).stem)

# uses a class style following the example of eady_growth_rate.py.


class TransectDiagnostic():
    '''
        Class used to compute a transect of a velocity with
            a sigma2 contour overlaid.
    '''
    def __init__(self, config):
        self.config = config

    def _get_data(self):

        input_data = self.config["input_data"]

        logger.info("Input data: %s", input_data)     

        filenames = input_data.keys()
        logger.info("Filenames: %s", filenames)
        
        # groups = group_metadata(input_data, "variable_group", sort="dataset")
        # for group_name in groups:
        #     logger.info("Processing variable %s", group_name)
        #     for attributes in groups[group_name]:
        #         logger.info("Processing dataset %s", attributes["dataset"])
        #         input_file = attributes["filename"]
        #         logger.info("Loading data from %s", input_file)


        return input_data

    def _compute_sigma(self, ds):
        ''' 
            Compute the sigma2 contour from the input dataset.
        '''

        pass

    def _plot_transect(self):

        pass

    def call(self):
        input_data = self._get_data()
        return input_data

def main():
    """Run Eady Growth Rate diagnostic."""
    with run_diagnostic() as config:
        diagnostic = TransectDiagnostic(config).call()



if __name__ == "__main__":
    main()
