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

def __none_caster__(cfg_dict):
    """
    checks the cfg for textual Nones and casts to None
    --------------------------------------------------
    turns recipe Nones into actual Nones on the first level of the dict
    (lists and strings)
    """
    
    for name, content in cfg_dict.items():
        if isinstance(content, str):
            if content == "None":
                cfg_dict.update({name:None})
        if isinstance(content, list):
            if any([c == "None" for c in content]):
                content = [None for c in content if c == "None"]
                cfg_dict.update({name:content})
    

def cfg_checker(cfg):
    """
    checks the cfg for relevant entries
    -----------------------------------
    adding the required cfg entries to the dictionary if missing
    checking for textual Nones
    """
    
    if "requests" not in cfg.keys():
        cfg['requests'] = [None]
    if "colormap" not in cfg.keys():
        cfg['colormap'] = ["binary"]
    if "vminmax" not in cfg.keys():
        cfg['vminmax'] = [None, None]

    __none_caster__(cfg)

    logger.info(cfg)

    return cfg

class input_handler(object):
    """
    Basic class to handle the input lists
    """
    def __init__(self, **kwargs):
        """
        initializes the input_handler
        -----------------------------
        setting up the main attributes for the object
        """
        
        super(input_handler, self).__init__(**kwargs)
        self.files_read = False
        self.files_ref = []
        self.files_others = []
    
        return
        
    def __str__(self):
        """
        prepares the data to look like a dictionary when printed
        --------------------------------------------------------
        """
        return dict({"reference":self.files_ref,
                     "others":self.files_others}).__str__()
    
    def set_files(self, names_list, read = True):
        """
        reads in the list of filenames
        ------------------------------
        preparing relevant input
        """
        
        self.files_ref = [ds for ds,ci in names_list.items() 
            if ci["dataset"]==ci["reference_dataset"]]
        self.files_others = [ds for ds,ci in names_list.items() 
            if ci["dataset"]!=ci["reference_dataset"]]
        
        logger.info(self.files_ref)
        logger.info(self.files_others)
        
        if read:
            self.read()
        
        return
    
    def read(self):
        """
        reads in the files according to the list of filenames
        -----------------------------------------------------
        preparing relevant input
        """
        
        if self.files_read:
            logger.warning("Data was already read in! Nothing done.")
        else:
            self.files_ref = iris.cube.CubeList([iris.load_cube(iref) 
                for iref in self.files_ref])
            self.files_others = iris.cube.CubeList([iris.load_cube(idat) 
                for idat in self.files_others])
            self.files_read = True
            
        return
    
    def get_ref(self):
        """
        return the reference files
        --------------------------
        preparing relevant output
        """
        
        self.__unread_warning__()
    
        return self.files_ref

    def get_all(self):
        """
        return the all files
        --------------------
        preparing relevant output
        """
        
        self.__unread_warning__()
        
        return self.files_ref + self.files_others
    
    def get_others(self):
        """
        return the non-reference files
        ------------------------------
        preparing relevant output
        """
        
        self.__unread_warning__()
        
        return self.files_others
    
    def __unread_warning__(self):
        """
        returns a warning if files are not yet read in
        ----------------------------------------------
        preparing relevant output
        """
        
        if not self.files_read:
            logger.warning("Files are not read in yet! " + 
                           "Only the list of filenames is returned!")
        
    
    