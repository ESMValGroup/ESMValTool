#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 13:16:54 2019

@author: bmueller
"""
import iris
import logging
import os
import numpy as np
from matplotlib.colors import LinearSegmentedColormap as LS_cmap
from matplotlib.cm import get_cmap 
from copy import copy as copy
from scipy.stats.stats import kendalltau

logger = logging.getLogger(os.path.basename(__file__))

def set_metadata(c1, c2, fun):
    """
    combines the metadata of 2 cubes
    --------------------------------
    returns metadata with lists for each metadata attributes entry if distinct
    """
    
    c1m = c1.metadata
    c2m = c2.metadata
    
    sta = set([c1m.standard_name, c2m.standard_name])
    if len(sta)==1:
        sta = list(sta)[0]
        
    lon = set([c1m.long_name, c2m.long_name])
    if len(lon)==1:
        lon = list(lon)[0]
    
    var = set([c1m.var_name, c2m.var_name])
    if len(var)==1:
        var = list(var)[0]
        
    uni = set([c1m.units, c2m.units])
    if len(uni)==1:
        uni = list(uni)[0]
        
    att = __unify_attributes__(c1m.attributes, c2m.attributes)
    cel = fun
    
    meta = iris.cube.CubeMetadata(standard_name = sta,
                                  long_name = lon,
                                  var_name = var,
                                  units = uni,
                                  attributes = att,
                                  cell_methods = cel)
    
    return meta


def __unify_attributes__(att1, att2):
    """
    combines two attributes dictionaries from metadata
    --------------------------------------------------
    returns dictionary with lists for each attribute entry
    """
    
    attributes = dict()
    
    for key, entry in __common_entries__(att1, att2):
        if len(entry)==1:
            attributes.update({key:list(entry)[0]})
        elif isinstance(entry, str):
            attributes.update({key:entry})
        else:
            attributes.update({key:list(entry)})
    
    return attributes


def __flatten_to_strings__(strlist):
    """
    flattens a list of string lists or strings for 
    multiple levels of nesting
    --------------------------
    returns a flattened list of strings
    """
    result = []

    for sl in strlist:
        if isinstance(sl, str): # append if i is a string
            result.append(sl)
        else:                   # call this function recursively
            result.extend(__flatten_to_strings__(sl))
    return result


def __common_entries__(*dcts):
    """
    common entries with results for multiple dictionaries
    -----------------------------------------------------
    returns iterator
    """
    for i in set(dcts[0]).intersection(*dcts[1:]):
        try:
            yield (i,set([d[i] for d in dcts]))
        except TypeError:
            yield (i,set(__flatten_to_strings__([d[i] for d in dcts])))


def adjust_minmax(minmax, symmetric = False):
    """
    adjust the readability of calculated values
    -------------------------------------------
    returns adjusted values with an option for symmetric values
    """
    
    vmin = np.nanmin(minmax)
    vmax = np.nanmax(minmax)
    
    if symmetric:
        vmin = np.nanmin([vmin, -vmax])
        vmax = -vmin
        
    rounder = int(np.ceil(-np.log10(vmax - vmin)))
    vmin = np.floor(vmin * 10**rounder) / 10**rounder
    vmax = np.ceil(vmax * 10**rounder) / 10**rounder
    
    return vmin, vmax


def __return_methods__(module):
    """
    checks the module for available methods
    ---------------------------------------
    returns a list of available methods
    """
    
    return [method_name for method_name in dir(module) 
            if callable(getattr(module, method_name))]
    

def usable_methods(module):
    """
    nice printout for available methods in module
    ---------------------------------------------
    returns a well formatted string
    """
    
    splitter = "\n    - "
    methods = __return_methods__(module)
    meth_str = splitter.join(methods)
    
    return "Available diagnostics to request are:" + splitter + meth_str


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
        cfg['colormap'] = "binary"
    if "vminmax" not in cfg.keys():
        cfg['vminmax'] = [None, None]
    if "percentiles" not in cfg.keys():
        cfg['percentiles'] = [None]

    __none_caster__(cfg)

    return cfg


def cmap_switch(cmap, num):
    """
    checks the structure of the given cmap from cfg
    -----------------------------------------------
    switches between textual colormap or a list of named colors
    """
    
    if not isinstance(cmap, str):
        if isinstance(cmap, list):
            if all(isinstance(cm, str) for cm in cmap):
                cmap = LS_cmap.from_list("localcmap", cmap, N=num-1)
                
    else:
        cmap = get_cmap(cmap, num-1)
        
    return cmap


def checked_ref(refs, num=1):
    """
    checks the number of reference files and adjusts the attributes
    ---------------------------------------------------------------
    returns an error when too many ref files are available and adjusts
    the structure accordingly
    """
    
    if len(refs) != num:
        raise ValueError("reference data must have length {} " + 
                         "but is of length {}".format(num, len(refs)))
    for r in refs:
        r.metadata.attributes["source_file"] = "ref: {}".format(
            r.metadata.attributes["source_file"])
        
    if num == 1:
        refs = refs[0]
        
    return refs


def get_filenames(sourcefiles):
    """
    reduces the sourcefile attributes to readable filenames
    -------------------------------------------------------
    returns a list of strings
    """
    
    filenames = [sf.split(os.sep)[-1] for sf in sourcefiles]
    
    refs = [sf.split(" ") for sf in sourcefiles]
    refs = [r[0] if len(r)==2 else "" for r in refs]
    
    filenames = list(zip(refs,filenames))
    
    filenames = [" ".join(fn) for fn in filenames]
    
    logger.info(filenames)
    
    return filenames


def clean_filename(string):
    """
    changes all whitespace to underscores
    -------------------------------------
    returns a clean string
    """
    
    string = "_".join(string.split())
    string = string.replace(":","_")
    
    return string


def correlation(c1, c2):
    """
    calculates Kendall-Tau correlation for input cubes
    --------------------------------------------------
    returns correlation and p-value
    """
    
    kt = kendalltau(c1.core_data(), c2.core_data())
    
    return {"r": kt.correlation, "p-value": kt.pvalue}

def corr_extract(clist):
    """
    extracts a correlation pandas dataframe from a list of dictionaries
    -------------------------------------------------------------------
    returns pandas dataframe
    """
    
    corr_data = np.array([np.array([[c["r"] for c in v["corr"]] 
                    for v in d.values()][0]) for d in clist])
    
    corr_data = np.hstack((np.zeros((corr_data.shape[0], 1)), corr_data))
    corr_data[corr_data == 0] = np.nan
    
    return corr_data

def mean_std_txt(cube):
    """
    extracts mean and standard deviation from cube
    and returns print friendly text
    -------------------------------
    returns string
    """
    
    weights = iris.analysis.cartography.area_weights(cube)
    mean = cube.collapsed(["longitude", "latitude"],
                          iris.analysis.MEAN, weights = weights)
    std = cube.collapsed(["longitude", "latitude"],
                          iris.analysis.STD_DEV)
    
    txt = ("spatial statistics:\n" +
           "\u03BC: {:.3g}, " + 
           "\u03C3: {:.3g}").format(mean.data, std.data).expandtabs()
    
    return txt

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
        self.files_nonref = []
    
        return
        
    def __str__(self):
        """
        prepares the data to look like a dictionary when printed
        --------------------------------------------------------
        """
        
        return dict({"reference":self.files_ref,
                     "others":self.files_nonref}).__str__()
    
    def set_files(self, names_list, read = True):
        """
        reads in the list of filenames
        ------------------------------
        preparing relevant input
        """
        
        self.files_ref = [ds for ds,ci in names_list.items() 
            if ci["dataset"]==ci["reference_dataset"]]
        self.files_nonref = [ds for ds,ci in names_list.items() 
            if ci["dataset"]!=ci["reference_dataset"]]
        
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
            self.files_nonref = iris.cube.CubeList([iris.load_cube(idat) 
                for idat in self.files_nonref])
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
        
        return self.files_ref + self.files_nonref
    
    def get_nonref(self):
        """
        return the non-reference files
        ------------------------------
        preparing relevant output
        """
        
        self.__unread_warning__()
        
        return self.files_nonref
    
    def __unread_warning__(self):
        """
        returns a warning if files are not yet read in
        ----------------------------------------------
        preparing relevant output
        """
        
        if not self.files_read:
            logger.warning("Files are not read in yet! " + 
                           "Only the list of filenames is returned!")
            
    def ref_only(self):
        """
        returns a copy of the object with empty files_nonref
        ----------------------------------------------------
        preparing relevant output
        """
        
        ref = self.copy()
        
        ref.files_nonref = []
        
        return ref
        
    def nonref_only(self):
        """
        returns a copy of the object with empty files_ref
        -------------------------------------------------
        preparing relevant output
        """
        
        nonref = self.copy()
        
        nonref.files_ref = []
        
        return nonref
        
    def copy(self):
        """
        returns a copy of the object
        ----------------------------
        preparing relevant output
        """
        
        return copy(self)
    