#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 13:16:54 2019

@author: bmueller
"""
import iris
import logging
import os

logger = logging.getLogger(os.path.basename(__file__))

def cfg_checker(cfg):
    """
    checks the cfg for relevant entries
    -----------------------------------
    adding the required cfg entries to the dictionary if missing
    """
    
    if "requests" not in cfg.keys():
#        cfg['requests'] = [None]
        cfg['requests'] = ["glob_temp_mean"]

    return cfg

class input_handler(object):
    """
    Basic class to handle the input lists
    """
    def __init__(self, **kwargs):
        """
        initializes the input_handler
        ---------------------------
        setting up the main attributes for the object
        """
        
        super(input_handler, self).__init__(**kwargs)
        self.files_read = False
        self.files_ref = []
        self.files_others = []
    
        return
    
    def set_files(self, names_list, read = True):
        """
        reads in the list of filenames
        ---------------------------
        preparing relevant input
        """
        
        self.files_ref = [ds for ds,ci in names_list.items() 
            if ci["dataset"]==ci["reference_dataset"]]
        self.files_others = [ds for ds,ci in names_list.items() 
            if ci["dataset"]!=ci["reference_dataset"]]
        
        if read:
            self.read()
        
        return
    
    def read(self):
        """
        reads in the files according to the list of filenames
        ---------------------------
        preparing relevant input
        """
        
        self.files_ref = iris.cube.CubeList([iris.load_cube(iref) 
            for iref in self.files_ref])
        self.input_dat = iris.cube.CubeList([iris.load_cube(idat) 
            for idat in self.input_dat])
        
        self.files_read = True
        
        return
    
    def get_ref(self):
        """
        return the reference files
        ---------------------------
        preparing relevant output
        """
        
        self.__unread_warning__()
    
        return self.files_ref

    def get_all(self):
        """
        return the all files
        ---------------------------
        preparing relevant output
        """
        
        self.__unread_warning__()
        
        return self.files_ref + self.files_others
    
    def get_others(self):
        """
        return the non-reference files
        ---------------------------
        preparing relevant output
        """
        
        self.__unread_warning__()
        
        return self.files_others
    
    def __unread_warning__(self):
        """
        returns a warning if files are not yet read in
        ---------------------------
        preparing relevant output
        """
        
        if not self.files_read:
            logger.warning("Files are not read in yet! " + 
                           "Only the list of filenames is returned!")
        
    
    