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

from .c3s_511_basic import Basic_Diagnostic_SP
from .libs.MD_old.ESMValMD import ESMValMD
from .libs.predef.ecv_lookup_table import ecv_lookup
from .plots.basicplot import \
    Plot2D, PlotHist, Plot2D_blank, Plot1D, PlotScales, plot_setup

class ta_Diagnostic_SP(Basic_Diagnostic_SP):
    """
    class to implement additional diagnostics
    """
    
    def set_info(self, **kwargs):
        
        super(ta_Diagnostic_SP, self).set_info(**kwargs)
        
        # add a region to the regions object
        self.__regions__.update({
            'MAR_region': {
                'latitude': (-60, 50),
                'longitude': (-60, 0),
                # 'time': (datetime.datetime(2000, 1, 1),
                #          datetime.datetime(2000, 12, 31)
                #          )
                }})

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
                                    "time",
                                    iris.analysis.WPERCENTILE,
                                    percent=75,
                                    weights=iris.analysis.cartography.area_weights(
                                            loc_cube)) - \
                            loc_cube.collapsed(
                                    "time",
                                    iris.analysis.WPERCENTILE,
                                    percent=25,
                                    weights=iris.analysis.cartography.area_weights(
                                            loc_cube))
            except:
                self.__logger__.warning("Weighted percentiles over time " + 
                                        "could not be calculated.")
                return list_of_plots
                        
            # define filename and set to plot list
            filename = self.__plot_dir__ + os.sep + basic_filename + \
                "_" + "iqr" + "_" + c +\
                "_lat_lon" + "." + self.__output_type__
            list_of_plots.append(filename)
    
            # produce plot
            try:
                x = Plot2D(cube_iqr)
    
                caption = str("/".join(["lat","lon"]).title() +
                              ' ' +
                              "interquartile range (0.75 - 0.25)" +
                              ' maps of ' + c + " " +
                              ecv_lookup(self.__varname__) +
                              ' for the data set ' +
                              " ".join(dataset_id) + ' (' +
                              self.__time_period__ + ').')
    
                fig = plt.figure()
    
                (fig, ax, caption) = plot_setup(d="time",
                                                fig=fig,
                                                caption=caption)
    
                x.plot(ax=ax,
                       color=self.colormaps,
                       color_type="Sequential",
                       title=" ".join([self.__dataset_id__[indx] for
                                       indx in [0, 2, 1, 3]]) + \
                        " (" + self.__time_period__ + ")",
                       ext_cmap="both",
                       y_log=True,
                       contour_levels=[0,1,2],
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
    
                x = Plot2D_blank(loc_cube[0,:,:])
    
                fig = plt.figure()
    
                (fig, ax, caption) = plot_setup(d="time",
                                                fig=fig,
                                                caption=caption)
    
                x.plot(ax=ax,
                       color=self.colormaps,
                       color_type="Sequential",
                       title=" ".join([self.__dataset_id__[indx] for
                                       indx in [0, 2, 1, 3]]) + \
                        " (" + self.__time_period__ + ")",
                       ext_cmap="both")
                fig.savefig(filename)
                plt.close(fig.number)
    
                ESMValMD("meta",
                         filename,
                         self.__basetags__ + 
                         ['DM_global', 'C3S_mean_var'],
                         caption + ' Data can ' +
                         'not be displayed due to cartopy error!',
                         '#C3S' + "iqr" + "latlon" + \
                         self.__varname__,
                         self.__infile__,
                         self.diagname,
                         self.authors + ["me"])
                
            # relative interquartiles range over time
            try:
                cube_iqr_r = cube_iqr / \
                            loc_cube.collapsed(
                                    "time",
                                    iris.analysis.MEAN,
                                    weights=iris.analysis.cartography.area_weights(
                                            loc_cube))
            except:
                self.__logger__.warning("Relative IQR over time " + 
                                        "could not be calculated.")
                return list_of_plots
                        
            # define filename and set to plot list
            filename = self.__plot_dir__ + os.sep + basic_filename + \
                "_" + "iqr_r" + "_" + c + \
                "_lat_lon" + "." + self.__output_type__
            list_of_plots.append(filename)
    
            # produce plot
            try:
                x = Plot2D(cube_iqr_r)
    
                caption = str("/".join(["lat","lon"]).title() +
                              ' ' +
                              "relative interquartile range (0.75 - 0.25)" +
                              ' maps of ' + c + " " +
                              ecv_lookup(self.__varname__) +
                              ' for the data set ' +
                              " ".join(dataset_id) + ' (' +
                              self.__time_period__ + ').')
    
                fig = plt.figure()
    
                (fig, ax, caption) = plot_setup(d="time",
                                                fig=fig,
                                                caption=caption)
    
                x.plot(ax=ax,
                       color=self.colormaps,
                       color_type="Sequential",
                       title=" ".join([self.__dataset_id__[indx] for
                                       indx in [0, 2, 1, 3]]) + \
                        " (" + self.__time_period__ + ")",
                       ext_cmap="both")
                fig.savefig(filename)
                plt.close(fig.number)
    
                ESMValMD("meta",
                         filename,
                         self.__basetags__ + 
                         ['DM_global', 'C3S_mean_var'],
                         caption + " NA-values are shown in grey.",
                         '#C3S' + "iqrr" + "latlon" + \
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
    
                x = Plot2D_blank(cube[0,:,:])
    
                fig = plt.figure()
    
                (fig, ax, caption) = plot_setup(d="time",
                                                fig=fig,
                                                caption=caption)
    
                x.plot(ax=ax,
                       color=self.colormaps,
                       color_type="Sequential",
                       title=" ".join([self.__dataset_id__[indx] for
                                       indx in [0, 2, 1, 3]]) + \
                        " (" + self.__time_period__ + ")",
                       ext_cmap="both")
                fig.savefig(filename)
                plt.close(fig.number)
    
                ESMValMD("meta",
                         filename,
                         self.__basetags__ + 
                         ['DM_global', 'C3S_mean_var'],
                         caption + ' Data can ' +
                         'not be displayed due to cartopy error!',
                         '#C3S' + "iqrr" + "latlon" + \
                         self.__varname__,
                         self.__infile__,
                         self.diagname,
                         self.authors + ["me"])
        

        return list_of_plots
        