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

warnings.simplefilter("ignore")

logger = logging.getLogger(os.path.basename(__file__))

logger.info(dir())

import satan.diagnostics as diags
from satan.utilities import cfg_checker, input_handler


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
        self.input_dat = []
        self.input_ref = []
        self.files_read = False
        
        logger.info("init completed")
        return
    
    def set_info(self, cfg):
        """
        handles the cfg information
        ---------------------------
        distributing the cfg entries to the object
        """
        
        self.cfg = __cfg_checker__(cfg)
        
        # split inputs into ref and non ref
        self.input_ref = [ds for ds,ci in cfg['input_data'].items() if ci["dataset"]==ci["reference_dataset"]]
        self.input_dat = [ds for ds,ci in cfg['input_data'].items() if ci["dataset"]!=ci["reference_dataset"]]
        
        self.input.set_files(cfg['input_data'])
        
        logger.info("set_info completed")
        return
    
    def read(self):
        """
        reads the data sets from the filenames
        --------------------------------------
        overwrites the lists of names by a list of data
        """
        
        # read filenames into content
        self.input_ref = iris.cube.CubeList([iris.load_cube(iref) for iref in self.input_ref])
        self.input_dat = iris.cube.CubeList([iris.load_cube(idat) for idat in self.input_dat])
        
        self.input.read()
        
        logger.info(self.input)
        
        self.files_read = True
        
        logger.info("read completed")
        return
    
    def run(self):
        """
        runs the diagnostics for data sets
        ----------------------------------
        the diagnostics defined by the cfg['requests'] are run
        """
        
        logger.info('\n'.join("%s: %s" % item for item in vars(self).items()))
        
        logger.info("run completed")
        return