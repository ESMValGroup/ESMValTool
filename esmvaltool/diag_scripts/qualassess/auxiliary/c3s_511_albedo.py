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
import iris.coord_categorisation as cat

from .c3s_511_basic import Basic_Diagnostic_SP
from .libs.MD_old.ESMValMD import ESMValMD
from .libs.predef.ecv_lookup_table import ecv_lookup
from .plots.basicplot import \
    Plot2D, PlotHist, Plot2D_blank, Plot1D, PlotScales, plot_setup

class albedo_Diagnostic_SP(Basic_Diagnostic_SP):
    """
    class to implement additional diagnostics
    """
    
    def set_info(self, **kwargs):
        
        super(albedo_Diagnostic_SP, self).set_info(**kwargs)
        
        # add a region to the regions object
        self.__regions__ = {
            'SHem': {
                'latitude': (-90, 00),
                },
            'NHem': {
                'latitude': (90, 00),
                }}

    def __add_mean_var_procedures_2D__(self, cube=None, level=None):
        
        """
        - TODO: Histograms per land cover mask (distinguish albedo between classes) or
        - TODO: hemispheric land cover masked aggregated time series with standard deviation (find out reasonable natural variability and time series)
        - Seasonal (JJA, SON, DJF, MAM) maps
        """

        list_of_plots = []
        
#        hemispheres = self.__spatiotemp_subsets__(self.sp_data) 
#        
#        NH = hemispheres["NHemisphere"]
        
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + \
                "_" + "seasonal" + "_" +\
                "_lat_lon" + "." + self.__output_type__
        list_of_plots.append(filename)
        
        cat.add_season(cube, 'time', name='season')
        seasons = cube.aggregated_by('season', iris.analysis.MEAN)
        
        season_name = ['djf', 'mam', 'jja', 'son']
        
        season_dict = dict()
        for sn in season_name:
            season_dict.update({sn:seasons.extract(iris.Constraint(season=sn))})
            
        # produce plot
        try:
            x = Plot2D(list(season_dict.values()))
            self.__logger__.info(x.MAX_COLUMNS)
            x.MAX_COLUMNS = 2
            self.__logger__.info(dir(x))
            self.__logger__.info(x.MAX_COLUMNS)


            caption = str("/".join(["lat","lon"]).title() +
                          ' ' +
                          "mean seasonal" +
                          ' maps of ' + " " +
                          ecv_lookup(self.__varname__) +
                          ' for the data set ' +
                          " ".join(self.__dataset_id__) + ' (' +
                          self.__time_period__ + '). ' + 
                          'Subplots a) - f) are ' +
                          ", ".join(season_name).upper() +
                          ".")

            fig = plt.figure()

            (fig, ax, caption) = plot_setup(d="time",
                                            fig=fig,
                                            caption=caption)

            x.plot(ax=ax,
                   color=self.colormaps,
                   color_type="Data",
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
                     '#C3S' + "season" + "latlon" + \
                     self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors + ["bmueller"])

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(
                exc_tb.tb_frame.f_code.co_filename)[1]
            self.__logger__.error(
                exc_type, fname, exc_tb.tb_lineno)
            self.__logger__.error("season")
            self.__logger__.error('Warning: blank figure!')

            x = Plot2D_blank(season_dict.values()[0])

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
                     '#C3S' + "season" + "latlon" + \
                     self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors + ["bmueller"])
        

        return list_of_plots