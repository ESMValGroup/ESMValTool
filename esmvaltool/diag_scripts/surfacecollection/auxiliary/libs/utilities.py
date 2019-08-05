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
from scipy.stats import linregress
import cf_units

import functools
import time

def timer(func):
    """
    Print the runtime of the decorated function
    -------------------------------------------
    from https://realpython.com/primer-on-python-decorators/#simple-decorators
    """
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()    # 1
        value = func(*args, **kwargs)
        end_time = time.perf_counter()      # 2
        run_time = end_time - start_time    # 3
        logger.info(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value
    return wrapper_timer

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
        cfg['vminmax'] = [np.nan, np.nan]
    if "percentiles" not in cfg.keys():
        cfg['percentiles'] = [None]
    if "pthreshold" not in cfg.keys():
        cfg['pthreshold'] = None
    if "temporal_basis" not in cfg.keys():
        cfg['temporal_basis'] = [None]

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


def data_correlation(c1, c2):
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

def calculate_trend(cube, decadal = True):
    """
    calculates temporal trend and p-value from cube
    -----------------------------------------------
    returns cube dictionary for trend and p-value
    """
    
    slope_list = []
    pvalue_list = []
    
    # realization of data for speed
    cube.data 
    
    # deleting of aux_coords of data for speed
    delete_aux_coords(cube)
    
    for time_series in cube.slices("time"):
        slope = __make_dummy_cube__(time_series, ["latitude", "longitude"])
        pvalue = slope.copy()
        if not np.all(time_series.data.mask):
            ts_linmod = linregress(time_series.coord("time").points,
                                   time_series.data,
                                   )
            slope.data[0,0] = getattr(ts_linmod, "slope")
            pvalue.data[0,0] = getattr(ts_linmod, "pvalue")

        slope.units = cf_units.Unit(str(slope.units) +  " " +
                                    str(slope.coord(
                                            "time").units).split(" ")[0]
                                    + "-1" )
        pvalue.units = cf_units.Unit("-")
        
        slope_list.append(slope)
        pvalue_list.append(pvalue)
        
    slope_cube = iris.cube.CubeList(slope_list).concatenate_cube()
    pvalue_cube = iris.cube.CubeList(pvalue_list).concatenate_cube()
    
    if decadal:
        slope_cube.convert_units(cf_units.Unit(str(cube.units) + 
                                               " (10 years)-1"))
    else:
        slope_cube.convert_units(cf_units.Unit(str(cube.units) + 
                                               " year-1"))
        
    slope_cube.long_name = "Decadal Trend of " + cube.long_name
    pvalue_cube.long_name = "pvalue"
    
    return dict({"trend": slope_cube,
                 "p-value": pvalue_cube})
    
def calculate_correlation(cube, ref):
    """
    calculates temporal trend and p-value from cube
    -----------------------------------------------
    returns cube dictionary for trend and p-value
    """
    
    corr_list = []
    pvalue_list = []
    
    # realization of data for speed
    cube.data 
    
    # deleting of aux_coords of data for speed
    delete_aux_coords(cube)
    
    for time_series in cube.slices("time"):
        corr = __make_dummy_cube__(time_series, ["latitude", "longitude"])
        pvalue = corr.copy()
        ts_ref = ref.extract(iris.Constraint(
                longitude = time_series.coord("longitude").points))
        ts_ref = ts_ref.extract(iris.Constraint(
                latitude = time_series.coord("latitude").points))
        if not np.all(time_series.data.mask):
            ts_corr = data_correlation(ts_ref, time_series)
            corr.data[0,0] = ts_corr["r"]
            pvalue.data[0,0] = ts_corr["p-value"]

        corr_list.append(corr)
        pvalue_list.append(pvalue)
        
    corr_cube = iris.cube.CubeList(corr_list).concatenate_cube()
    pvalue_cube = iris.cube.CubeList(pvalue_list).concatenate_cube()
    
    corr_cube.metadata = set_metadata(corr_cube, ref, "correlation")
    pvalue_cube.metadata = set_metadata(pvalue_cube, ref, "pvalue")
    
    corr_cube.long_name = "Correlation of " + cube.long_name
    pvalue_cube.long_name = "pvalue"
    
    corr_cube.units = cf_units.Unit("-")
    pvalue_cube.units = cf_units.Unit("-")
    
    
    return dict({"correlation": corr_cube,
                 "p-value": pvalue_cube})

def __make_dummy_cube__(cube, dims="time"):
    """
    produces a dummy cube from a cube by calculating the means
    ----------------------------------------------------------
    returns cube with same dimcoords besides the aggregated ones
    """
    all_dims = [co.name() for co in cube.coords()
        if isinstance(co, iris.coords.DimCoord)]
    other_dims = all_dims
    [other_dims.remove(d) for d in dims]
    
    dummy = cube.collapsed(other_dims, iris.analysis.MEAN)
    
    dims.reverse()
    for d in dims:
        dummy = iris.util.new_axis(dummy, d)
    
    return dummy

def get_clim_categorisation(temporal_basis):
    """
    chooses climatology categorisation based on string
    --------------------------------------------------
    returns respective iris.coord_categorisation
    """
    
    if temporal_basis == "month":
        clim = iris.coord_categorisation.add_month_number
    elif temporal_basis == "day":
        clim = iris.coord_categorisation.add_day_of_year
    elif temporal_basis == "season":
        clim = iris.coord_categorisation.add_season
    else: 
        clim = iris.coord_categorisation.add_month_number
        logger.warning("No temporal_basis given (None)," +
                       " or temporal_basis unknown," +
                       " monthly climatology produced instead: " +
                       temporal_basis)
        logger.warning("Please use one of: day, month, season.")
        
    return clim

def __remove_all_aux_coords__(cube):
    """
    remove all auxiliary coordinates from cube
    ------------------------------------------
    returns cleaned cube
    """
    for dim in cube.coords():
        if isinstance(dim, iris.coords.AuxCoord):
             cube.remove_coord(dim)
            
    return
    
def get_differences_4_clim(cube, other):
    """
    calculate the differences between cube along the dimension clim
    ---------------------------------------------------------------
    returns difference cube
    """
    
    new_cube = []
    
    for map_slice in cube.slices(["latitude", "longitude"]):
        act_overlap_pos = map_slice.coord("clim").points[0]

        diff = map_slice - other.extract(
                iris.Constraint(clim = act_overlap_pos)
                )
        
        diff.add_aux_coord(map_slice.coord("time"))
        diff.remove_coord("clim")
        
        new_cube.append(diff)
    
    new_cube = iris.cube.CubeList(new_cube).merge_cube()
    
    diff.standard_name = cube.standard_name
    diff.long_name = cube.long_name
    diff.units = cube.units
    
    return new_cube

def change_long_name(cube, how, text):
    """
    wrapper to change the longname of the cube
    ------------------------------------------
    returns nothing (cube is adjusted)
    """
    
    new_name = cube.long_name
    
    if "overwrite" in how or "overwrite" == how:
        new_name = text
        
    elif "attach" in how:
        if "front" in how:
            new_name = text + new_name
        elif "end" in how:
            new_name = new_name + text
        else:
            ValueError("Command needed. " +
                       "Text can either be attached to 'front' or 'end'.")
    else:
        ValueError("Command needed. " +
                   "Text can be used to either 'overwrite' or 'attach'.")
    
    cube.long_name = new_name
    
    return 

def delete_aux_coords(cube):
    """
    deletes all aux coords from cube
    --------------------------------
    returns nothing (in place functions)
    """

    for ad in cube.coords():
        if isinstance(ad, iris.coords.AuxCoord):
            cube.remove_coord(ad)
            
    return

def copy_metadata(cube, other):
    """
    copies the metadata from other into cube
    ----------------------------------------
    returns nothing (in place functions)
    """
    cube.attributes = other.attributes.copy()
    
    return

def calculate_anomalies(data, temporal_basis):
    """
    calculates anomalies based on aggregation information
    -----------------------------------------------------
    returns anomalies
    """
            
    clim_fun = get_clim_categorisation(temporal_basis)
    
    added_clim = data.apply_iris_fun(clim_fun, ctype = "adjustment",
                                       coord="time", name="clim")

    clim_agg_data = added_clim.apply_iris_fun("aggregated_by",
                                              ctype = "method",
                                              coords="clim",
                                              aggregator = iris.analysis.MEAN)

    anomalies = added_clim.apply_iris_fun(get_differences_4_clim,
                                                 ctype = "elementwise",
                                                 other = clim_agg_data)
    
    anomalies.apply_iris_fun(change_long_name, ctype = "adjustment",
                             how = ["overwrite"],
                             text = "Anomalies of " + 
                                 list(set(
                                         [d.long_name for d in data.get_all()]
                                         ))[0])
    
    anomalies.apply_iris_fun(copy_metadata, ctype = "elementwise",
                             other = data)
    
    return anomalies

def calculate_climatology(data, temporal_basis):
    """
    calculates climatology based on aggregation information
    -----------------------------------------------------
    returns climatology
    """
            
    clim_fun = get_clim_categorisation(temporal_basis)
    
    added_clim = data.apply_iris_fun(clim_fun, ctype = "adjustment",
                                       coord="time", name="clim")

    clim_agg_data = added_clim.apply_iris_fun("aggregated_by",
                                              ctype = "method",
                                              coords="clim",
                                              aggregator = iris.analysis.MEAN)
    
    clim_agg_data.apply_iris_fun(change_long_name, ctype = "adjustment",
                                 how = ["overwrite"],
                                 text = "Climatology of " + 
                                     list(set(
                                         [d.long_name for d in data.get_all()]
                                         ))[0])
    
    clim_agg_data.apply_iris_fun(copy_metadata, ctype = "elementwise",
                                 other = data)
    
    return clim_agg_data

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
    
    def apply_iris_fun(self, fun, ctype, **kwargs):
        """
        applies functions to ref and nonref cubes
        -----------------------------------------
        returning a copy with changed data
        """
        
        copy_input = self.copy()
        
        if not self.files_read:
            copy_input.read()
        
        if ctype == "method":
            copy_input.files_ref = [getattr(fil,fun)(**kwargs) 
                for fil in copy_input.files_ref]
            copy_input.files_nonref = [getattr(fil,fun)(**kwargs)  
                for fil in copy_input.files_nonref]
        elif ctype == "adjustment":
            [fun(cube = fil, **kwargs) 
                for fil in copy_input.files_ref]
            [fun(cube = fil, **kwargs) 
                for fil in copy_input.files_nonref]
        elif ctype == "overwrite":
            copy_input.files_ref = [fun(cube = fil, **kwargs) 
                for fil in copy_input.files_ref]
            copy_input.files_nonref = [fun(cube = fil, **kwargs) 
                for fil in copy_input.files_nonref]
        elif ctype == "elementwise":
            if "other" not in kwargs.keys():
                NameError("Could not find 'other' variable for method input.")
            elif not isinstance(kwargs["other"], input_handler):
                TypeError("'Other' variable not of type input_handler.")
            else: 
                other = kwargs.pop("other")
                
            ref = []
            nonref = []
            
            for dat in zip(self.files_ref, other.files_ref):
                ref.append(fun(dat[0], dat[1], **kwargs))
            copy_input.files_ref = ref
            
            for dat in zip(self.files_nonref, other.files_nonref):
                nonref.append(fun(dat[0], dat[1], **kwargs))
            copy_input.files_nonref = nonref
            
        else:
            TypeError("Wrong ctype assigned. Options are: " + 
                      "method, adjustment, overwrite.")
        
        return copy_input
    
    def __sub__(self, data):
        """
        calculating difference
        ----------------------
        returning a copy with changed data
        """
        
        sub = self.copy()
        
        if not self.files_read:
            sub.read()
        
        sub_ref = []
        sub_nonref = []
        
        for (a, b) in zip(sub.files_ref, data.files_ref):
            sub_entry = a - b
            sub_entry.metadata = set_metadata(a, b, ["__sub__"])
            sub_ref.append(sub_entry)
        for (a, b) in zip(sub.files_nonref, data.files_nonref):
            sub_entry = a - b
            sub_entry.metadata = set_metadata(a, b, ["__sub__"])
            sub_nonref.append(sub_entry)
        
        sub.files_ref = sub_ref
        sub.files_nonref = sub_nonref
        
        return sub
    
    def __add__(self, data):
        """
        calculating sum
        ---------------
        returning a copy with changed data
        """
        
        add = self.copy()
        
        if not self.files_read:
            add.read()
        
        add_ref = []
        add_nonref = []
        
        for (a, b) in zip(add.files_ref, data.files_ref):
            add_entry = a + b
            add_entry.metadata = set_metadata(a, b, ["__add__"])
            add_ref.append(add_entry)
        for (a, b) in zip(add.files_nonref, data.files_nonref):
            add_entry = a + b
            add_entry.metadata = set_metadata(a, b, ["__add__"])
            add_nonref.append(add_entry)
        
        add.files_ref = add_ref
        add.files_nonref = add_nonref
        
        return add
    
    def __mul__(self, data):
        """
        calculating product
        -------------------
        returning a copy with changed data
        """
        
        mul = self.copy()
        
        if not self.files_read:
            mul.read()
        
        mul_ref = []
        mul_nonref = []
        
        for (a, b) in zip(mul.files_ref, data.files_ref):
            mul_entry = a * b
            mul_entry.metadata = set_metadata(a, b, ["__mul__"])
            mul_ref.append(mul_entry)
        for (a, b) in zip(mul.files_nonref, data.files_nonref):
            mul_entry = a * b
            mul_entry.metadata = set_metadata(a, b, ["__mul__"])
            mul_nonref.append(mul_entry)
        
        mul.files_ref = mul_ref
        mul.files_nonref = mul_nonref
        
        return mul
    
    def has_lazy_data(self):
        """
        checks all cubes in object on lazyness
        -------------------------------------
        returns a copy with boolean representation
        """
        
        hld = self.apply_iris_fun("has_lazy_data", "method")
        
        return hld
        
#    def realize_data(self):
#        """
#        realizes all iris cubes
#        -----------------------
#        returns nothing
#        """
#        
#        cd = self.apply_iris_fun("core_data", "method")
#        
#        hld = self.has_lazy_data()
#        
#        def __check_before_compute__(core_data, has_lazy_data):
#            """
#            checks if objects each are lazy 
#            ------------------------------
#            returns nothing
#            """
#            
#            if has_lazy_data:
#                core_data.compute()
#            
#            return
#        
#        cd.apply_iris_fun(__check_before_compute__,
#                          "elementwise",
#                          other = hld)
#        
#        return
