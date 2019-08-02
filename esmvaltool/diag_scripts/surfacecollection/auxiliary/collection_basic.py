#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:55:02 2019

@author: bmueller
"""

import os
import logging
import iris
import warnings
from pprint import pprint
from .libs import diagnostics as diags
from .libs import plots as plots
from .libs.utilities import cfg_checker, input_handler, usable_methods

warnings.simplefilter("ignore")

logger = logging.getLogger(os.path.basename(__file__))


class ecv_handler(object):
    """
    Basic class to implement any kind of diagnostic for the surface ecvs
    """    
    def __init__(self, **kwargs):
        """
        initializes the ecv_handler
        ---------------------------
        setting up the main attributes for the object
        """
        
        super(ecv_handler, self).__init__(**kwargs)
        
        self.cfg = {}
        self.input = input_handler()
        self.files_read = False
        
        logger.info("init completed")
        return
    
    def set_info(self, cfg):
        """
        handles the cfg information
        ---------------------------
        distributing the cfg entries to the object
        """
        
        self.cfg = cfg_checker(cfg)
        logger.setLevel(cfg['log_level'].upper())
        
        # handle input
        self.input.set_files(cfg['input_data'], read = False)
        
        logger.info("set_info completed")
        return
    
    def read(self):
        """
        reads the data sets from the filenames
        --------------------------------------
        overwrites the lists of names by a list of data
        """
        
        # read filenames into content
        self.input.read()
        
        self.files_read = True
        
        logger.info("read completed")
        return
    
    def run(self):
        """
        runs the diagnostics for data sets
        ----------------------------------
        the diagnostics defined by the cfg['requests'] are run
        """
        
        # report request options if all None
        if all([req is None for req in self.cfg["requests"]]):
            logger.info(usable_methods(diags))
        else:
            for diag in self.cfg["requests"]:
                results = getattr(diags, diag)(self.input,
                                 pthreshold=self.cfg["pthreshold"],
                                 percentiles=self.cfg["percentiles"],
                                 temporal_basis=self.cfg["temporal_basis"],
                                 )
                for r in results:
                    getattr(plots, diag)(r,
                           cmap=self.cfg["colormap"],
                           vminmax=self.cfg["vminmax"],
                           plotdir=self.cfg["plot_dir"],
                           fformat=self.cfg["output_file_type"])
            
#        logger.info('\n'.join("%s: %s" % item for item in vars(self).items()))
        
        logger.info("run completed")
        return