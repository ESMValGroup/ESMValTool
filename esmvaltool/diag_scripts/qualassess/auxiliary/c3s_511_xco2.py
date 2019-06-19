#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  16 10:15:00 2019
@author: bhassler
"""

import iris
import os
import sys
import matplotlib.pyplot as plt
import datetime

from .c3s_511_basic import Basic_Diagnostic_SP
from .libs.MD_old.ESMValMD import ESMValMD
from .libs.predef.ecv_lookup_table import ecv_lookup
from .plots.basicplot import \
    Plot2D, PlotHist, Plot2D_blank, Plot1D, PlotScales, plot_setup

class xco2_Diagnostic_SP(Basic_Diagnostic_SP):
    """
    class to implement additional diagnostics
    """
    
    def set_info(self, **kwargs):
        
        super(xco2_Diagnostic_SP, self).set_info(**kwargs)
        
        # add a region to the regions object
        self.__regions__ = {
            'MaunaLoa_region': {
                'latitude': (17, 20),
                'longitude': (203, 205)
                }}

    def __add_mean_var_procedures_2D__(self, cube=None, level=None):
        
        # handling of the cube restrictions copied from the basic diagnostic
        # needs to be set locally
        #################
        if cube is None:
            cube = self.sp_data

        # adjustment to ids and filenames
        if level is not None:
            basic_filename = self.__basic_filename__ + "_lev" + str(level)
            dataset_id = [self.__dataset_id__[0], "at", "level", str(
                level), str(self.sp_data.coord(self.level_dim).units)]
        else:
            basic_filename = self.__basic_filename__
            dataset_id = [self.__dataset_id__[0]]
        ##################
        
        list_of_plots = []
        
        
        cubes = self.__spatiotemp_subsets__(cube) 
        # if not further defined, self.__regions__ are taken.
        
        # go through all regional subsets 
        for c in cubes:
            
            loc_cube = cubes[c]
        
            # interquartiles range over time for the MAR area and for Europe
            # !!!this is not working, if there are all-Nan areas in the level!!!!
            try:
                cube_iqr = loc_cube.collapsed(
                                    iris.analysis.MEAN,
                                    ) 

            except:
                self.__logger__.warning("Weighted mean over time " + 
                                        "could not be calculated.")
                return list_of_plots
                        
            # define filename and set to plot list
            filename = self.__plot_dir__ + os.sep + basic_filename + \
                "_" + "iqr" + "_" + c +\
                "_lat_lon" + "." + self.__output_type__
            list_of_plots.append(filename)
    
            # produce plot
            try:
                print("in try loop")
                x = Plot1D(cube_iqr)
                print("after plot1d")
                caption = str("/".join(["lat","lon"]).title() +
                              ' ' +
                              "mean" +
                              ' time series of ' + c + " " +
                              ecv_lookup(self.__varname__) +
                              ' for the data set ' +
                              " ".join(dataset_id) + ' (' +
                              self.__time_period__ + ').')
    
                fig = plt.figure()
                print("after plt.figure")
                #(fig, ax, caption) = plot_setup(d="time",
                #                                fig=fig,
                #                                caption=caption)
                ax = [plt.subplot(1, 1, 1)]
                fig.set_figheight(1.2 * fig.get_figheight())
                
                x.plot(ax=ax,
                       colors=["blue"], 
                       title=["MeaN XCO2 at grid point closest to Mauna Loa, Hawaii"]
                       )
                fig.savefig(filename)
                plt.close(fig.number)
    
                ESMValMD("meta",
                         filename,
                         self.__basetags__ + 
                         ['DM_global', 'C3S_mean_var'],
                         caption + " NA-values are shown in grey.",
                         '#C3S' + "iqr" + "latlon" + \
                         self.__varname__,
                         self.__infile__,
                         self.diagname,
                         self.authors + ["me"])
    
            except Exception as e:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(
                    exc_tb.tb_frame.f_code.co_filename)[1]
                self.__logger__.error(
                    exc_type, fname, exc_tb.tb_lineno)
                self.__logger__.error("iqr")
                self.__logger__.error('Warning: blank figure!')
    
                print("in exception block")                

        return list_of_plots