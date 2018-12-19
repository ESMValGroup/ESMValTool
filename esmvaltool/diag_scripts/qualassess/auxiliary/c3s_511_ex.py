#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 13:59:59 2018

@author: bmueller
"""

import iris
import os
import sys
import matplotlib.pyplot as plt
import datetime
import numpy as np

from .c3s_511_basic import Basic_Diagnostic_SP
from .libs.MD_old.ESMValMD import ESMValMD
from .libs.predef.ecv_lookup_table import ecv_lookup
from .plots.basicplot import \
    Plot2D, PlotHist, Plot2D_blank, Plot1D, PlotScales, plot_setup

class ex_Diagnostic_SP(Basic_Diagnostic_SP):
    """
    class to implement additional diagnostics
    """
    
    def set_info(self, **kwargs):
        
        super(ex_Diagnostic_SP, self).set_info(**kwargs)
        
        # add a region to the regions object
#        self.__regions__.update({
#            'MAR_region_Aug/Sept': {
#                'latitude': (-60, 50),
#                'longitude': (-60, 0),
#                'time': (datetime.datetime(2000, 8, 1),
#                         datetime.datetime(2000, 9, 30)
#                         )
#                }})
    
    def run_diagnostic(self):
#        self.sp_data = self.__spatiotemp_subsets__(self.sp_data)['Europe_2000']
        self.__do_extremes__()

        self.__do_full_report__()

    def __do_extremes__(self):
        
        this_function =  "extremes example"
        
        which_percentile = 90
        window_size = 5 # one directional 5 => 11
        
        # this the extremes example
        
        # handling of the cube restrictions copied from the basic diagnostic
        # needs to be set locally
        #################
        if not self.var3D:
            cube = self.sp_data
        else:
            cube = self.sp_data.extract(iris.Constraint(
                    coord_values={str(self.level_dim): lambda cell: cell == max(self.levels)}))

        # adjustment to ids and filenames
        if self.var3D is not None:
            basic_filename = self.__basic_filename__ + "_lev" + str(max(self.levels))
            dataset_id = [self.__dataset_id__[0], "at", "level", str(
                max(self.levels)), str(self.sp_data.coord(self.level_dim).units)]
        else:
            basic_filename = self.__basic_filename__
            dataset_id = [self.__dataset_id__[0]]
        ##################
        
        # we assume the regions to be spatio-temporally buffered.
        list_of_plots = []
        
        for r in self.__regions__:
            spat_r = self.__regions__[r]
            spat_r.pop("time")
            loc_cube = self.__spatiotemp_subsets__(cube,{r:spat_r})
            
            list_of_cubes = []
            
            for yx_slice in loc_cube[r].slices(['time']):
                agg_cube = yx_slice.aggregated_by("day_of_year",iris.analysis.MEAN) 
                list_of_sub_cubes=[]
                for doy in np.sort(agg_cube.coords("day_of_year")[0].points):
                    loc_slice = agg_cube.extract(iris.Constraint(day_of_year=doy))
                    tmin = (doy - 5) % 366
                    tmax = (doy + 5) % 366
                    if not tmin:
                        tmin = 366
                    if not tmax:
                        tmax = 366
                    if tmin > tmax:
                        window = lambda cell: 1 <= cell <= tmin or tmax <= cell <= 366
                    else:
                        window = lambda cell: tmin <= cell <= tmax
                    doy_sel = yx_slice.extract(iris.Constraint(day_of_year=doy)) # TODO: there is some issue here that "window" cannot be used
                    try:
                        perc = doy_sel.collapsed("time",iris.analysis.PERCENTILE,percent=[90]).data[0]
                    except:
                        perc = doy_sel.data
                     
                    loc_slice.data = perc
                    list_of_sub_cubes.append(loc_slice)
                self.__logger__.warning(iris.cube.CubeList(list_of_sub_cubes).merge()) # TODO: merge does not work properly
                list_of_cubes.append(iris.cube.CubeList(list_of_sub_cubes).merge())
                
            clim_cube = iris.cube.CubeList(list_of_cubes).merge()
            
            self.__logger__.warning(clim_cube)
        
        
        # produce report
        expected_input, found = \
            self.__file_anouncement__(subdir="c3s_511/single_extremes_input",
                                      expfile="_extremes.txt",
                                      protofile="empty.txt",
                                      function=this_function)

        if found:

            self.reporting_structure.update(
                    {"Extremes": 
                        {"plots": list_of_plots,
                         "freetext": expected_input}})
        else:

            self.reporting_structure.update(
                    {"Extremes": 
                        {"plots": list_of_plots}})