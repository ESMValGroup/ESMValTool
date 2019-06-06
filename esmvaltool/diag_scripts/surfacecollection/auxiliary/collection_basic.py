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

warnings.simplefilter("ignore")

logger = logging.getLogger(os.path.basename(__file__))

class ecv_handler(object):
    """
    Basic class to implement any kind of diagnostic for the surface ecvs
    """
    def __init__(self, **kwargs):
        
        self.cfg = None
        self.input_ref = []
        self.input_dat = []
        self.files_read = False
        
        logger.info("init completed")
        return
    
    def set_info(self, cfg):
        
        self.cfg = cfg
        
        self.input_ref = [ds for ds,ci in cfg['input_data'].items() if ci["dataset"]==ci["reference_dataset"]]
        
        self.input_dat = [ds for ds,ci in cfg['input_data'].items() if ci["dataset"]!=ci["reference_dataset"]]
        
        logger.info("set_info completed")
        return
    
    def read(self):
        
        self.input_ref = iris.cube.CubeList([iris.load_cube(iref) for iref in self.input_ref])
        
        self.input_dat = iris.cube.CubeList([iris.load_cube(idat) for idat in self.input_dat])
        
        self.files_read = True
        
        logger.info("read completed")
        return
    
    def run(self):
        
        logger.info('\n'.join("%s: %s" % item for item in vars(self).items()))
        
        logger.info("run completed")
        return