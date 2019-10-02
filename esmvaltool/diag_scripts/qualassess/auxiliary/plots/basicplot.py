#/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
#from matplotlib.ticker import LogFormatterMathtext
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import cf_units as unit
import matplotlib.cm as mpl_cm
import matplotlib.colors as mpl_colors
import cartopy.crs as ccrs
from ..libs import c3s_511_util as utils
#import random
import sys
import string
from matplotlib.ticker import FuncFormatter
#import matplotlib.colors as colors
import logging
#from memory_profiler import profile
from dask import array as da

logger = logging.getLogger(os.path.basename(__file__))

def label_in_perc_multiple(x, pos=0):
    return '%1.1f%%' % (x)

def label_in_perc_single(x, pos=0):
    return '%1.1f%%' % (x * 100)

def label_in_exp(x,pos=0):
    return r"$10^{%.1f}$" % (x)

def label_in_exp_perc(x,pos=0):
    if x * 100 < 0.1:
        return "{}%".format(r"$10^{%d}$" % (np.log10(x*100)))
    else:
        return label_in_perc_single(x)


MPLSTYLE = os.path.dirname(
    os.path.realpath(__file__)
) + os.sep + 'default.mplstyle'


#class PlotHist2(object):
#    """
#    Description
#        Basic class for plotting histograms
#
#    Contents
#        method plot
#    """
#
################################################################################
#
#    def __init__(self, data):
#        """
#        Arguments
#            data : input data (iris cube or array)
#
#        Description
#            Initializes the class.
#
#        Modification history
#            20180209-A_schl_ma: written
#        """
#        self.logger = logging.getLogger(os.path.basename(__file__))
#        # Check arguments
#        if isinstance(data, iris.cube.Cube):
#            try:
#                self.name = data.long_name
#            except BaseException:
#                pass
#            self.units = ' [' + str(data.units) + ']'
#            self.data = np.ravel(data.data)
#        elif isinstance(data, np.ndarray):
#            self.name = 'data'
#            self.units = ''
#            self.data = data
#        elif isinstance(data, list):
#            self.name = []
#            self.units = []
#            self.data = []
#            for d in data:
#                if isinstance(d, iris.cube.Cube):
#                    try:
#                        self.name.append(d.long_name)
#                    except BaseException:
#                        pass
#                    self.units.append(' [' + str(d.units) + ']')
#                    self.data.append(np.ravel(d.data))
#                elif isinstance(d, np.ndarray):
#                    self.name.append('data')
#                    self.units.append('')
#                    self.data.append(d)
#
#        else:
#            raise TypeError(
#                "Invalid input: expected iris cube(s) or numpy array(s)")
#
#        # Matplotlib
#        plt.style.use(MPLSTYLE)
#        self.fig, self.ax = plt.subplots()
#
################################################################################
#
#    def plot(self, nbins=20, x_label=None, y_label=None, title=None,
#             color='#15b01a', alpha=1.):
#        """
#        Arguments
#            bins    : number of bins
#            x_label : label of x-axis
#            y_label : label of y-axis
#            title   : title of the plot
#            color   : color of the histogramm
#            alpha   : transparency of the historgramm
#
#        Returns
#            Matplotlib figure instance
#
#        Description
#            Actual plotting routine
#
#        Modification history
#            20180209-A_schl_ma: written
#        """
#
#        # Parse arguments
#        if (x_label is None):
#            x_label = self.name + self.units
#        if (y_label is None):
#            y_label = 'frequency'
#
#        # Create historgramm
##        self.ax.hist(self.data, bins, density=False, facecolor=color, alpha=alpha)
##        print self.ax.__dict__.keys()
#        if isinstance(self.data, list):
#            x_label = (" [").join(list(set(self.name)) +
#                                  list(set([su.split("[")[1]
#                                            for su in self.units]))) + "]"
#            cols = mpl_cm.gist_rainbow(np.linspace(0, 1, len(self.data)))
#            labels = title
#            vmin = np.inf
#            vmax = -np.inf
#            for d in self.data:
#                vmin = np.min([np.min(d), vmin])
#                vmax = np.max([np.max(d), vmax])
#            rounder = int(np.ceil(-np.log10((vmax + vmin) / 2) + 2))
##            vmin, vmax = np.round([vmin, vmax], rounder)
#            vmin = np.floor(vmin * 10**rounder) / 10**rounder
#            vmax = np.ceil(vmax * 10**rounder) / 10**rounder
#            levels = np.round(np.linspace(
#                vmin, vmax, num=nbins * len(self.data)), rounder)
#
#            for ind, _ in enumerate(self.data):
#                try:
#                    n, bins, patches = self.ax.hist(self.data[ind].
#                                                    compressed(),
#                                                    bins=levels,
#                                                    density=True,
#                                                    facecolor=cols[ind],
#                                                    alpha=alpha /
#                                                    len(self.data),
#                                                    label=labels[ind])
#                except BaseException:
#                    n, bins, patches = self.ax.hist(self.data[ind].data,
#                                                    bins=levels,
#                                                    density=True,
#                                                    facecolor=cols[ind],
#                                                    alpha=alpha /
#                                                    len(self.data),
#                                                    label=labels[ind])
#        else:
#            vmin = np.min(self.data)
#            vmax = np.max(self.data)
#            rounder = int(np.ceil(-np.log10((vmax + vmin) / 2) + 2))
#            vmin = np.floor(vmin * 10**rounder) / 10**rounder
#            vmax = np.ceil(vmax * 10**rounder) / 10**rounder
#            levels = np.round(np.linspace(vmin, vmax, num=nbins), rounder)
#
#            try:
#                n, bins, patches = self.ax.hist(
#                    self.data.compressed(),
#                    bins=levels, density=True, facecolor=color, alpha=alpha)
#
#            except BaseException:
#                n, bins, patches = self.ax.hist(
#                    self.data.data,
#                    bins=levels, density=True, facecolor=color, alpha=alpha)
#
#        self.ax.set_xlabel(x_label)
#        self.ax.set_ylabel(y_label)
#        if isinstance(self.data, list):
#            self.ax.yaxis.set_major_formatter(
#                FuncFormatter(label_in_perc_multiple))
#            box = self.ax.get_position()
#            self.ax.set_position([box.x0, box.y0,
#                                  box.width, box.y1 - box.height * 0.1])
#
#            # Put a legend above current axis
#            self.ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#                           ncol=3, mode="expand", borderaxespad=0.)
#        else:
#            self.ax.yaxis.set_major_formatter(
#                FuncFormatter(label_in_perc_single))
#            if (title is not None):
#                self.ax.set_title(title)
#
#        self.fig.tight_layout()
#        return self.fig


class PlotHist(object):
    """
    Description
        Basic class for plotting histograms

    Contents
        method plot
    """

###############################################################################

    def __init__(self, data):
        """
        Arguments
            data : input data (iris cube or array)

        Description
            Initializes the class.

        Modification history
            20180209-A_schl_ma: written
        """
        self.logger = logging.getLogger(os.path.basename(__file__))
        # Check arguments
        if isinstance(data, iris.cube.Cube):
            try:
                self.name = data.long_name
            except BaseException:
                pass
            self.units = ' [' + str(data.units) + ']'
            self.data = data
        elif isinstance(data, np.ndarray):
            self.name = 'data'
            self.units = ''
            self.data = data
        elif isinstance(data, list):
            self.name = []
            self.units = []
            self.data = []
            for d in data:
                if isinstance(d, iris.cube.Cube):
                    try:
                        self.name.append(d.long_name)
                    except BaseException:
                        pass
                    self.units.append(' [' + str(d.units) + ']')
                    self.data.append(np.ravel(d.data))
                elif isinstance(d, np.ndarray):
                    self.name.append('data')
                    self.units.append('')
                    self.data.append(d)

        else:
            raise TypeError(
                "Invalid input: expected iris cube(s) or numpy array(s)")

        # Matplotlib
        plt.style.use(MPLSTYLE)
        self.fig, self.ax = plt.subplots()

###############################################################################

    def plot(self, nbins=20, x_label=None, y_label=None, title=None,
             dat_log=False, color='#15b01a', alpha=1.):
        """
        Arguments
            bins    : number of bins
            x_label : label of x-axis
            y_label : label of y-axis
            title   : title of the plot
            color   : color of the histogramm
            alpha   : transparency of the historgramm

        Returns
            Matplotlib figure instance

        Description
            Actual plotting routine

        Modification history
            20180209-A_schl_ma: written
        """

        # Parse arguments
        if (x_label is None):
            x_label = self.name + self.units
        if (y_label is None):
            y_label = 'frequency'
            
        vmin = np.min(self.data.core_data().compute())
        vmax = np.max(self.data.core_data().compute())
        rounder = int(np.ceil(-np.log10((vmax - vmin) / 2) + 5))
        vmin = np.floor(vmin * 10**rounder) / 10**rounder
        vmax = np.ceil(vmax * 10**rounder) / 10**rounder
        levels = np.round(np.linspace(vmin, vmax, num=nbins), rounder)
        
        hist, bins = da.histogram(self.data.core_data(), bins=levels, range=[vmin, vmax])
        hist = hist.compute()
        hist = hist/np.sum(hist)
        
        x = 0.5 * (bins[1:] + bins[:-1])
        width = np.diff(bins)
        
        self.ax.bar(x,
                    hist,
                    width,
                    color = color)
        
        
        if dat_log:
            self.ax.set_yscale('log', nonposy='clip')
        
        self.ax.yaxis.set_major_formatter(FuncFormatter(label_in_exp_perc))#label_in_perc_single))
        if (title is not None):
            self.ax.set_title(title)
        self.ax.set_xlabel(x_label)
        self.ax.set_ylabel(y_label)
        self.ax.set_ylim(0,np.max(hist) * 1.1)

        self.fig.tight_layout()
        return self.fig


class PlotScatter(object):
    """
    Description
        Basic class for creating scatter plots

    Contents
        method plot
    """

###############################################################################

    def __init__(self, data1, data2):
        """
        Arguments
            data1 : single cube (y-axis)
            data2 : list of cubes which should be compared (x-axis)

        Description
            Initializes the class.

        Modification history:
            20180209-A_schl_ma: written
        """
        self.logging.getLogger(os.path.basename(__file__))
        # Check arguments
        if (not isinstance(data1, iris.cube.Cube)):
            raise TypeError("Invalid input: expected iris cube")
        if isinstance(data2, iris.cube.Cube):
            self.data2 = [data2]
        else:
            self.data2 = data2
        for d in self.data2:
            if (not isinstance(d, iris.cube.Cube)):
                raise TypeError("Invalid input: expected list of iris cubes " +
                                "or a single cube")
            if (data1.coords != d.coords):
                raise ValueError("Invalid input: all cubes need to have the " +
                                 "same dimensions")
        self.y = np.ravel(data1.data)

        # Matplotlib
        plt.style.use(MPLSTYLE)
        self.fig, self.ax = plt.subplots()

###############################################################################
#        
#        
#        counts_all = np.array(levels[1:]) * 0
#        
#        for subset in self.data.slices(["latitude","longitude"]):
#            counts = np.array(utils.count_cube_vals_for_levels(subset,levels))
#            counts_all = counts_all + counts
#        
#        plt.text(1.58,1.90,counts)
    def plot(self, x_label=None, y_label=None, title=None, add_info=None):
        """
        Arguments
            x_label  : label of x-axis
            y_label  : label of y-axis
            title    : title of the plot
            add_info : dictionay with additional text

        Returns
            Matplotlib figure instance

        Description
            Actual plotting routine

        Modification history
            20180209-A_schl_ma: written
        """

        # Parse arguments
        if (x_label is None):
            x_label = 'data 2'
        if (y_label is None):
            y_label = 'data 1'

        # Create scatter plot
        for d in self.data2:
            x = np.ravel(d.data)
            self.ax.plot(x, self.y, linestyle='none',
                         marker='o', markerfacecolor='none')

        # Plot appearance
        self.ax.set_xlabel(x_label)
        self.ax.set_ylabel(y_label)
        if (title is not None):
            self.ax.set_title(title)
        if (add_info is not None):
            text = ""
            for key in add_info:
                text += "{0} = {1}\n".format(key, add_info[key])
            self.ax.text(0.02, 0.8, text, transform=self.ax.transAxes)

        self.fig.tight_layout()
        return self.fig


class Plot2D(object):
    """
    Description
        Basic class for 2-dimensional plotting

    Contents
        method plot
    """

    # Class attributes
    LATS = ['latitude']     # accepted lat names
    LONS = ['longitude']    # accepted lon names
    TIME = ['time']         # accepted time names
    LEVS = ['air_pressure']  # accepted level names
    MAX_COLUMNS = 3
    TIME_LABEL = ''
    TIME_FORMAT = '%Y-%m-%d'
    VALSCALE = None

    def __init__(self, cubes):
        """
        Arguments
            cubes : iris cube or a list of iris cubes

        Description
            Initializes the class.

        TODO
            optional summary placement
            swap x/y axis

        Modification history
            20180207-A_muel_bn: copied Plot2D and adjusted
        """
        self.logger = logging.getLogger(os.path.basename(__file__))
        # Check arguments
        try:
            self.n_cubes = len(cubes)
        except TypeError:
            self.n_cubes = 1
            cubes = [cubes]
        for (idx, cube) in enumerate(cubes):
            if not isinstance(cube, iris.cube.Cube):
                raise TypeError("Invalid input: expected single iris cube or "
                                "list of iris cubes")
            cubes[idx] = iris.util.squeeze(cube)
            if cubes[idx].ndim != 2:
                raise TypeError("Invalid input: expected 2-dimensional iris "
                                "cubes at list index {}".format(idx))
        self.cubes = cubes
        if self.n_cubes < 1:
            raise TypeError("Invalid input: expected at least one cube, not "
                            "an empty list")

        elif self.n_cubes > 21:
            # TODO break images instead of just plotting the last 21
            self.cubes = cubes[-21:]
            self.n_cubes = 21

        # Get properties of the cubes
        self.names = [None] * self.n_cubes
        self.units = [None] * self.n_cubes
        self.time_vars = [None] * self.n_cubes
        self.lev_vars = [None] * self.n_cubes
        self.lat_vars = [None] * self.n_cubes
        self.lon_vars = [None] * self.n_cubes
        for (idx, cube) in enumerate(self.cubes):

            # Name
            try:
                self.names[idx] = cube.long_name
            except BaseException:
                pass

            # Units
            self.units[idx] = (' [' + str(cube.units) + ']')

            # Dimension names
            dim_names = [dim.standard_name for dim in cube.dim_coords]
            for dim in dim_names:
                if (dim in self.__class__.TIME):
                    self.time_vars[idx] = dim
                elif (dim in self.__class__.LEVS):
                    self.lev_vars[idx] = dim
                elif (dim in self.__class__.LATS):
                    self.lat_vars[idx] = dim
                elif (dim in self.__class__.LONS):
                    self.lon_vars[idx] = dim

        # Lat/lon plot
        if (all(lat is not None for lat in self.lat_vars) and
                all(lon is not None for lon in self.lon_vars)):
            self.plot_type = 'latlon'
            self.x_axis = 'longitude'
            self.y_axis = 'latitude'
            self.data_transposed = False

        # Lat/time plot
        elif (all(lat is not None for lat in self.lat_vars) and
              all(time is not None for time in self.time_vars)):
            self.plot_type = 'lattime'
            self.x_axis = 'time'
            self.y_axis = 'latitude'
            self.data_transposed = True

        # Lon/time plot
        elif (all(lon is not None for lon in self.lon_vars) and
              all(time is not None for time in self.time_vars)):
            self.plot_type = 'lontime'
            self.x_axis = 'longitude'
            self.y_axis = 'time'
            self.data_transposed = False

        # Lat/lev plot
        elif (all(lat is not None for lat in self.lat_vars) and
              all(lev is not None for lev in self.lev_vars)):
            self.plot_type = 'latlev'
            self.x_axis = 'latitude'
            self.y_axis = 'air_pressure'
            self.data_transposed = False

        # Lon/lev plot
        elif (all(lon is not None for lon in self.lon_vars) and
              all(lev is not None for lev in self.lev_vars)):
            self.plot_type = 'lonlev'
            self.x_axis = 'longitude'
            self.y_axis = 'air_pressure'
            self.data_transposed = False

        # time/lev plot
        elif (all(time is not None for time in self.time_vars) and
              all(lev is not None for lev in self.lev_vars)):
            self.plot_type = 'timelev'
            self.x_axis = 'time'
            self.y_axis = 'air_pressure'
            self.data_transposed = True

        # Default case
        else:
            raise TypeError("Invalid input: Not all cubes have the same "
                            "dimensions or at least on of them has an "
                            "unknown dimension name")

        # Setup matplotlib
        plt.style.use(MPLSTYLE)

###############################################################################

    def plot(
        self,
        summary_plot=False,
        title=None,
        ax=None,
        vminmax=None,
        color=None,
        color_type=None,
        color_reverse=False,
        ext_cmap='neither',
        y_log=False,
        dat_log=False,
        contour_levels=None,
    ):
        """
        Arguments
            summary_plot   : Add summary line plot
            colorbar_ticks : Ticks of the colorbar
            x_label        : label of x-axis
            y_label        : label of y-axis
            title          : title of the plot

        Returns
            Matplotlib figure

        Description
            Actual plotting routine

        Modification history
            20180515-A_muel_bn: copied Plot2D and adjusted
        """
        # Preprocessing cube information
        if self.n_cubes == 1:
            self.cubes[0].rename(title)
        for (idx, cube) in enumerate(self.cubes):
            try:
                self.cubes[idx] = cube.intersection(longitude=(-180, 180))
            except BaseException:
                pass

        # Color
        if color is None:
            brewer_cmap = mpl_cm.get_cmap('brewer_Spectral_11')
        else:
            if color_type is None or color_type not in color.keys():
                if color_reverse:
                    col_save = color["default"]
                    if color["default"][-2:] == "_r":
                        color["default"] = color["default"][:-2]
                    else:
                        color["default"] = color["default"] + "_r"
                brewer_cmap = mpl_cm.get_cmap(color["default"], lut=11)
            else:
                if color_reverse:
                    col_save = color[color_type]
                    if color[color_type][-2:] == "_r":
                        color[color_type] = color[color_type][:-2]
                    else:
                        color[color_type] = color[color_type] + "_r"
                brewer_cmap = mpl_cm.get_cmap(color[color_type], lut=11)

        brewer_cmap.set_bad("grey", 0.1)

        # Vminmax
        if vminmax is None:
            vmin, vmax = np.inf, -np.inf
            for cube in self.cubes:
                try:
                    if cube.data is np.ma.masked:
                        temp_min, temp_max = np.nanpercentile(
                            cube.data.compressed(), [5.0, 95.0])
                    else:
                        temp_min, temp_max = np.nanpercentile(
                            cube.data, [5.0, 95.0])
                    if (temp_min < vmin):
                        vmin = temp_min
                    if (temp_max > vmax):
                        vmax = temp_max
                except BaseException:
                    pass
            if vmin > vmax:
                vmin = None
                vmax = None
                levels = None
        else:
            if len(vminmax) == 2:
                vmin, vmax = vminmax
            else:
                vmin = vmax = vminmax[0]

        if vmin is not None:
            if vmax == vmin:
                vmin = vmin - 1.
                vmax = vmax + 1.
            elif vmin > vmax:
                vmin, vmax = [vmax, vmin]

            rounder = int(np.ceil(-np.log10(vmax - vmin) + 1))
            vmin = np.floor(vmin * 10**rounder) / 10**rounder
            vmax = np.ceil(vmax * 10**rounder) / 10**rounder
            levels = np.round(np.linspace(vmin, vmax, num=11), rounder)

        # Logarithmic y axis
        y_logarithmic = False
        if y_log:
            if 'lev' in self.plot_type:
                y_logarithmic = True
            else:
                self.logger.warning("Logarithmic y axis is only supported for "
                                    "level plots")

        # Check axes
        if ax is None:
            ax = [plt.gca()]
        try:
            if len(ax) == 3 and not summary_plot:
                summary_plot = True
                self.logger.warning("Too many plot axes (3)! Summary plot "
                                    "will automatically be added.")
            elif len(ax) > 3 or len(ax) < 2:
                raise ValueError("Invalid input: axes should not be more "
                                 "than 3 or less than 2!")
        except TypeError:
            ax = [ax]

        # Plot summary if necessary
        if summary_plot:
            if len(ax) < 3:
                raise ValueError("Invalid input: Need 3 axes for summary plot")
            plt.sca(ax[1])

            # Latlon plot and multiple cubes not supported yet
            if self.n_cubes > 1:
                raise ValueError("Invalid input: summary plot for multiple "
                                 "cubes is not supported")
            cube = self.cubes[0]
            time_var = self.time_vars[0]
            lev_var = self.lev_vars[0]
            lat_var = self.lat_vars[0]
            lon_var = self.lon_vars[0]

            # Get variable to collapse
            if (self.plot_type == 'latlon'):
                raise ValueError("Invalid input: latlon should not have a "
                                 "summary plot")
            elif (self.plot_type == 'lattime'):
                collapse = time_var
            elif (self.plot_type == 'lontime'):
                collapse = time_var
            else:
                collapse = [c.name() for c in cube.coords() if c.name() !=
                            lev_var and len(c.points) > 1][0]
            if (collapse == lat_var):
                grid_areas = iris.analysis.cartography.area_weights(cube)
            else:
                grid_areas = None
            cube_line = cube.collapsed(collapse, iris.analysis.MEAN,
                                       weights=grid_areas)

            # Lattime
            if self.plot_type == 'lattime':
                x = cube_line
                y = cube_line.coord(lat_var)
                latrange = cube_line.coords(lat_var).pop()
                latrange = (np.min(latrange.points), np.max(latrange.points))
                plt.gca().set_xlim(vmin, vmax)
                plt.gca().set_ylim(latrange)

            # Lontime
            elif self.plot_type == 'lontime':
                x = cube_line.coord(lon_var)
                y = cube_line
                lonrange = cube_line.coords(lon_var).pop()
                lonrange = (np.min(lonrange.points), np.max(lonrange.points))
                plt.gca().set_ylim(vmin, vmax)
                plt.gca().set_xlim(lonrange)
            # Lev plots
            else:
                if y_logarithmic:
                    cube_line.coord(lev_var).points = np.log10(
                        cube_line.coord(lev_var).points)
                x = cube_line
                y = cube_line.coord(lev_var)
                levrange = cube_line.coords(lev_var).pop()
                levrange = (np.min(levrange.points), np.max(levrange.points))
                plt.gca().set_xlim(vmin, vmax)
                plt.gca().set_ylim(levrange)

            # Plot
            iplt.plot(x, y)
            if y_logarithmic:
#                locs = plt.gca().get_yticks()
#                (bottom, top) = plt.gca().get_ylim()
#                labels = np.round(10**locs, decimals=1)
#                plt.gca().set_yticklabels(labels)
                plt.gca().set_ylim(levrange[::-1])
                plt.gca().yaxis.set_major_formatter(
                        FuncFormatter(label_in_exp))


        # Setup axes for plotting multiple cubes
        n_columns = min([self.n_cubes, self.MAX_COLUMNS])
        rel_plots_cbar = 4
        n_rows = rel_plots_cbar * \
            ((self.n_cubes - 1) // self.MAX_COLUMNS + 1) + 1
        gs = gridspec.GridSpec(n_rows, n_columns)

        if dat_log:
        
            vmin = np.inf
            
            for cube in self.cubes:
                
                loc_cube = cube.copy()
                loc_data = loc_cube.core_data()
                if isinstance(np.ma.getmask(loc_data), bool):
                    if np.ma.getmask(loc_data) == np.ma.nomask:
                        loc_data = da.ma.masked_where(loc_data<=0.0, loc_data)
                else:
                    loc_data = da.ma.masked_where(da.logical_or(np.ma.getmask(loc_data),
                                                  loc_data<=0.0), loc_data)
                loc_cube.data = loc_data.compute()
                temp_min= da.min(loc_data)
                temp_min=10**np.floor(np.log10(temp_min))
                if temp_min<vmin:
                    vmin = temp_min

        # Iterate over cubes
        for (idx, cube) in enumerate(self.cubes):
            column = idx % self.MAX_COLUMNS
            row = idx // self.MAX_COLUMNS
            if self.n_cubes > 1:
                curplot = gs[
                    (rel_plots_cbar * row):(rel_plots_cbar * (row + 1)),
                    column
                ]
                plt.sca(plt.subplot(curplot, transform='robinson'))
            else:
                plt.sca(ax[0])

            # Plot map
            # (this needs to be done due to an error in cartopy)
#            try:
            for i in [1]:
                if y_logarithmic:
                    cube.coord(lev_var).points = np.log10(
                        cube.coord(lev_var).points)
                if contour_levels:
                    options = {
                        'colors': 'black',
                        'linewidths': 2,
                    }
                    if 'time' not in self.plot_type:
                        cs = iplt.contour(cube,
                                          contour_levels,
                                          **options)
                    else:
                        if self.data_transposed:
                            data = cube.data.transpose()
                        else:
                            data = cube.data
                        cs = plt.contour(cube.coord(self.x_axis).points,
                                         cube.coord(self.y_axis).points,
                                         data,
                                         contour_levels,
                                         **options)
                    plt.clabel(cs, cs.levels)

                if dat_log:
#                    loc_cube = cube.copy()
#                    loc_data = loc_cube.core_data()
#                    if np.ma.getmask(loc_data) == np.ma.nomask:
#                        loc_data = np.ma.MaskedArray(loc_data,mask=loc_data<=0.0)
#                    else:
#                        loc_data.mask = np.logical_or(loc_data.mask,loc_data<=0.0)
#                    loc_cube.data = loc_data
                    
                    logticks=np.linspace(np.log10(vmin),
                                            np.log10(vmax),
                                            len(levels)+1)
                    ticks = 10**logticks
                    
                    pcm = iplt.pcolormesh(loc_cube,
                                    cmap=brewer_cmap,
                                    vmin=vmin,
                                    vmax=vmax,
                                    norm=mpl_colors.LogNorm())
                else:
                    pcm = iplt.pcolormesh(cube,
                                          cmap=brewer_cmap,
                                          vmin=vmin,
                                          vmax=vmax)
                if y_logarithmic:
                    plt.ylim(levrange[::-1])
                    plt.gca().yaxis.set_major_formatter(
                            FuncFormatter(label_in_exp))
                if 'time' in self.plot_type:
                    time_coord = cube.coord(self.time_vars[idx])
                    locs = [tp for tp in time_coord.points if
                            np.isclose(tp % 365.2425, 0,
                                       atol=np.mean(
                                           np.diff(time_coord.points)))]
                    while len(locs) > 7:
                        locs = [locs[lind] for (lind, _) in
                                enumerate(locs) if not lind % 2]
                    locs = locs + [time_coord.points[-1]]
                    labels = unit.num2date(locs,
                                           time_coord.units.name,
                                           time_coord.units.calendar)
                    for (idx, _) in enumerate(labels):
                        labels[idx] = labels[idx].strftime(
                            self.__class__.TIME_FORMAT)
                    if self.plot_type == 'lontime':
                        plt.yticks(locs, labels, rotation=25)
                        plt.ylabel(self.__class__.TIME_LABEL)
                    else:
                        plt.xticks(locs, labels, rotation=25)
                        plt.xlabel(self.__class__.TIME_LABEL)
#            except Exception as e:
#                self.logger.exception(e)
#                exc_type, exc_obj, exc_tb = sys.exc_info()
#                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
#                self.logger.debug(exc_type, fname, exc_tb.tb_lineno)
#                pcm = qplt.pcolormesh(cube, cmap=brewer_cmap,
#                                      vmin=vmin, vmax=vmax, norm=None)
#                plt.text(0.5, 0.5, 'Data cannot be displayed as intended due '
#                         'to cartopy issues with the data cube!',
#                         horizontalalignment='center',
#                         verticalalignment='center',
#                         transform=plt.gca().transAxes)
                        
            if self.plot_type == 'latlon':
                mean = cube.collapsed(
                    [coord.name() for coord in cube.coords()],
                    iris.analysis.MEAN,
                    weights=iris.analysis.cartography.area_weights(cube)
                        ).data
                std = utils.weighted_STD_DEV(
                    cube,
                    [coord.name() for coord in cube.coords()],
                    weights=iris.analysis.cartography.area_weights(cube)
                        ).data
                plt.gca().coastlines()
                plt.gca().gridlines(crs=ccrs.Geodetic(), color="k",
                                    linestyle=':')
                plt.gca().text(0,
                       - 0.1 * ((self.n_cubes // n_columns) * 0.7 if self.n_cubes>n_columns else 1),
                       r'mean: {:.2g} $\pm$ {:.2g} '.format(
                           mean,#str(np.round(mean, 2)),
                           std),#str(np.round(std, 2))),
                       transform=plt.gca().transAxes)

            # Label plots (supports at most 26 figures at the moment)
            if self.n_cubes > 1:
                letter = '(' + __my_string_ascii_lc__(idx) + ')'
                plt.text(0.01, 0.99, letter,
                         horizontalalignment='left',
                         verticalalignment='top',
                         transform=plt.gca().transAxes)

        if self.n_cubes > 1:
            cax = plt.subplot(gs[n_rows - 1, :])
        else:
            cax = ax[len(ax) - 1]
        
        if dat_log:
            cbar = plt.colorbar(pcm,
                cax=cax,
                orientation='horizontal',
                spacing='proportional',
                fraction=1.,
                extend=ext_cmap,
                ticks = ticks,
#                boundaries=(levels)
                )
            cbar.ax.set_xticklabels([('{:.'+str(2)+'g}').format(x) for x in ticks])
#            pass
        else:
            plt.colorbar(pcm,
                cax=cax,
                orientation='horizontal',
                fraction=1.,
                extend=ext_cmap,
                boundaries=levels)
        cax.set_xlabel((list(set(self.names))[0] if list(
            set(self.names))[0] else "") + list(set(self.units))[0])

        # Colors
        if color_type is None or color_type not in color.keys():
            if color_reverse:
                color["default"] = col_save
        else:
            if color_reverse:
                color[color_type] = col_save

        # Set position of summary plot
        if summary_plot:
            if self.plot_type not in ['lontime', 'latlon']:
                bb = ax[1].get_position()
                bb.y0 = ax[0].get_position().y0
                ax[1].set_position(bb)

        try:
            plt.tight_layout()
        except BaseException:
            pass
        return


class Plot2D_blank(Plot2D):
    """
    Description
        Blank class for 2-dimensional plotting

    Contents
        method plot
    """

    def __init__(self, cubes):
        """
        Arguments
            cube : iris cube

        Description
            Initializes the class.

        TODO
            optional summary placement
            swap x/y axis

        Modification history
            20180515-A_muel_bn: written
        """
        super(Plot2D_blank, self).__init__(cubes)
        # erase all data
        for cube in self.cubes:
            cube.data = cube.data * 0.


class Plot1D(object):
    """
    Description
        Basic class for 1-dimensional plotting

    Contents
        method plot
    """

    # Class attributes
    LATS = ['latitude']     # accepted lat names
    LONS = ['longitude']    # accepted lon names
    TIME = ['time']         # accepted time names
    TIME_LABEL = ''
    TIME_FORMAT = '%Y-%m-%d' #'%Y-%m-%d %H:%M:%S'

    def __init__(self, cube):
        """
        Arguments
            cube : iris cube

        Description
            Initializes the class.


        Modification history
            20180527-A_muel_bn: copied Plot2D and adjusted
        """
        self.logger = logging.getLogger(os.path.basename(__file__))
        # Check arguments
        if (not (isinstance(cube, iris.cube.Cube) or
                 isinstance(cube, iris.cube.CubeList))):
            raise TypeError("Invalid input: expected iris cube(s)")

        if isinstance(cube, iris.cube.Cube):
            cube = iris.cube.CubeList([cube])

        self.cube = iris.cube.CubeList([iris.util.squeeze(c) for c in cube])
        if not all([c.ndim == 1] for c in self.cube):
            raise TypeError("Invalid input: expected 1-dimensional iris cube")
        try:
            self.name = list(filter(None,list(set(c.long_name for c in cube))))[0]
        except BaseException:
            pass
        self.units = ' [' + str(list(set(c.units for c in cube))[0]) + ']'

        # Get dimension names
        dim_names = [list(
                set(dim.standard_name for dim in c.dim_coords)
                        )[0] for c in cube]
        for dim in dim_names:
            if (dim in self.__class__.LATS):
                self.lat_var = dim
                break
            else:
                self.lat_var = None
        for dim in dim_names:
            if (dim in self.__class__.LONS):
                self.lon_var = dim
                break
            else:
                self.lon_var = None
        for dim in dim_names:
            if (dim in self.__class__.TIME):
                self.time_var = dim
                break
            else:
                self.time_var = None

        # Lat/lon plot
        if (self.lat_var is not None):
            self.plot_type = 'lat'

        # Lat/time plot
        elif (self.time_var is not None):
            self.plot_type = 'time'

        # Lon/time plot
        elif (self.lon_var is not None):
            self.plot_type = 'lon'

        # Default case
        else:
            raise TypeError("Invalid input: cube does not contain supported " +
                            "dimensions")

        # Setup matplotlib
        plt.style.use(MPLSTYLE)

###############################################################################

    def plot(self, title=None, ax=None, colors=None, dat_log=True):
        """
        Arguments
            title          : title of the plot
            ax             : ax to put the plot in

        Returns
            Matplotlib figure instance

        Description
            Actual plotting routine

        Modification history
            20180527-A_muel_bn: copied Plot2D_2 and adjusted
        """
        #brewer_cmap = mpl_cm.get_cmap('brewer_Spectral_11')

        if not (len(colors) == len(self.cube) 
                and len(title) == len(self.cube)):
            raise ValueError("Invalid input: color or title length not" +
                                 " compatibplotcubes_std_p,le:" + str(len(self.cube))) 


        # check axes
        if ax is None:
            ax = [plt.gca()]
        try:
            if len(ax) >= 2:
                raise ValueError("Invalid input: axes should not be more " + 
                                 "than 1!")
        except TypeError:
            ax = [ax]

        if colors is None:
            cols = mpl_cm.gist_rainbow(np.linspace(0, 1, len(self.cube)))
        else:
            cols = colors
            
        plt.sca(ax[0])
        
        ymin = ymax = np.nan
        
        plotlist = []

        # plot line
        try:
            for ind, c in enumerate(self.cube):
                c.data = 10**(c.data / 70)
                linplot = plt.plot(
                                c.coords("time")[0].points,
                                c.data,
                                color=cols[ind],
                                label=title[ind],
                                )
                if dat_log:
                    plt.gca().set_yscale("log")
                plotlist.append(linplot.copy())
                ymin=np.nanmin([ymin,np.nanmin(c.data)])
                ymax=np.nanmax([ymax,np.nanmax(c.data)])
                buffer = 0.1 * (ymax - ymin)
            plt.gca().set_ylabel(self.name + " " + str(self.units),
                                 rotation=90)

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            self.logger.debug(exc_type, fname, exc_tb.tb_lineno)
            self.logger.exception(e)
            self.logger.debug('We did not expect this to fail!')
            plt.plot()
            
#        if dat_log:
#            seq_y = 10**np.linspace((np.floor(np.log10(ymin - buffer)*100)/100), (np.ceil(np.log10(ymax + buffer)*100)/100),10)
#            plt.yticks(seq_y, seq_y)
#            self.logger.info(seq_y)
            
        if 'time' == self.plot_type:
            pointsum = list(
                set([np.sum(c.coords("time")[0].points) for c in self.cube]))
            pointlen = list(
                set([len(c.coords("time")[0].points) for c in self.cube]))
            cal = list(set([np.sum(c.coords("time")[0].units)
                            for c in self.cube]))
            if not (len(pointsum) == 1 and len(
                    pointlen) == 1 and len(cal) == 1):
                raise ValueError("Time coordinates are not consistent!")
            time_coord = self.cube[0].coords("time")[0]
            (locs, _) = plt.xticks()

            locs = [tp for tp in time_coord.points if
                    np.isclose(tp % 365.2425, 0,
                               atol=np.mean(
                                   np.diff(time_coord.points)))]
            while len(locs) > 7:
                locs = [locs[lind] for (lind, _) in
                        enumerate(locs) if not lind % 2]
            locs = locs + [time_coord.points[-1]]
            labels = unit.num2date(locs,
                                   time_coord.units.name,
                                   time_coord.units.calendar)
            for (idx, _) in enumerate(labels):
                labels[idx] = labels[idx].strftime(
                    self.__class__.TIME_FORMAT)

            plt.xticks(locs, labels, rotation=25)
            plt.xlabel(self.__class__.TIME_LABEL)
            
            plt.ylim(ymin - buffer, ymax + buffer)
            plt.grid()  
            

    
            if len(self.cube) > 1:
                box = plt.gca().get_position()
                plt.gca().set_position([box.x0, box.y0 + 0.05 * box.height,
                                        box.width, box.y1 - box.height * 0.15])
                handles, labels = plt.gca().get_legend_handles_labels()
                order = [labels.index(t) for t in title]
                handles = [handles[o] for o in order]
                
                # Put a legend above current axis
                plt.gca().legend(handles, title, 
                       bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                       ncol=7, mode="expand", borderaxespad=0.)
            
#        # removes ylabel instead of whitespace?!
#        try:
#            plt.tight_layout()
#        except BaseException:
#            pass

        return


def plot_setup(d="time", m="module", numfigs=1, fig=plt.figure(), caption=''):

    if d in ["longitude", "levels"]:
        gs = gridspec.GridSpec(10, 5)
        ax = np.array([plt.subplot(gs[:-2, :-1]),
                       plt.subplot(gs[:-2, -1]), plt.subplot(gs[-2:, :])])
        fig.set_figwidth(1.7 * fig.get_figwidth())
        fig.set_figheight(1.2 * fig.get_figheight())
    elif "time" == d:
        gs = gridspec.GridSpec(9, 1)
        ax = np.array([plt.subplot(gs[:-2, 0]), plt.subplot(gs[-1, 0])])
        fig.set_figheight(1.2 * fig.get_figheight())
    elif "latitude" == d:
        gs = gridspec.GridSpec(16, 1)
        ax = np.array([plt.subplot(gs[:-5, 0]),
                       plt.subplot(gs[-5:-2, 0]), plt.subplot(gs[-2:, 0])])
        fig.set_figheight(1.7 * fig.get_figheight())
    if "climatology" == m:
        fig.set_figheight(1.7 * fig.get_figheight())
        caption = caption + ' Subplots a) - l) show months January - December.'
    if "percentiles" == m:
        fig.set_figheight(1.3 * fig.get_figheight())
        caption = caption + ' Subplots a) - g) show percentiles ' + \
            '1%, 5%, 10%, 50%, 90%, 95%, and 99%.'
    if "anomalies" in m:
        fig.set_figheight(np.ceil(numfigs / 9.) * fig.get_figheight())
        caption = caption + ' Subplots ' + __my_string_ascii_lc__(0) + \
            ') - ' + __my_string_ascii_lc__(numfigs - 1) +  \
            ') show single years (max. last 21 years only).'

#    caption = caption + " : " + m + " / " + d

    return fig, ax, caption


def __my_string_ascii_lc__(n):
    if n > 701:
        raise ValueError(
            "You are trying to get more than 702 plots into one multiple plot."
            "This is not possible due to limited plot numbering.")
    numlet = int(n / 26)
    numrest = n % 26

    if n < 26:
        return string.ascii_lowercase[n]
    else:
        return string.ascii_lowercase[numlet - 1] + \
    string.ascii_lowercase[numrest]


class PlotScales(object):
    """
    Description
        Basic class for plotting redolution comparisons

    Contents
        method plot
    """

###############################################################################

    def __init__(self, data):
        """
        Arguments
            data : input data (dictionary of resolution infos)

        Description
            Initializes the class.

        Modification history
            20181016-A_muel_bn: written
        """

        if not isinstance(data, dict):
            raise ValueError(
                "Wrong data format for input, should be dict: " +
                str(type(data)))

        if "names" not in data.keys():
            raise KeyError("Dict is missing 'names' key: " + str(data.keys()))

        if not len(list(set(len(data[k]) for k in data.keys()))) == 1:
            raise ValueError(
                "Length of dictionary elements is not equal: " +
                str(data))

        self.data = data

###############################################################################

    def plot(self, ax=plt.gca()):

        def shinescale_0_1(l):
            l0 = np.array(l)

            vmin = np.min(l0)
            vmax = np.max(l0)
            diff = vmax - vmin
            if diff == 0:
                diff = 1.
            rounder = int(np.ceil(-np.log10(diff)))
            vmin = (np.floor(vmin * 10**rounder) - 1) / 10**rounder
            vmax = (np.ceil(vmax * 10**rounder) + 1) / 10**rounder
            levels = np.round(np.linspace(vmin, vmax, num=4), rounder)

            l0 = (l0 - vmin) / (vmax - vmin)

            labels = list(levels)
            levels = (levels - np.min(levels)) / \
                (np.max(levels) - np.min(levels))

            return((list(l0), levels, labels))

        ax.plot()
        plt.quiver([0] * len(self.data["names"]),
                   range(len(self.data["names"])),
                   [1.15] * len(self.data["names"]),
                   [0] * len(self.data["names"]),
                   angles='xy', scale_units='xy', scale=0.95)
        plt.quiver([1.] * len(self.data["names"]),
                   range(len(self.data["names"])),
                   [-1.15] * len(self.data["names"]),
                   [0] * len(self.data["names"]),
                   angles='xy', scale_units='xy', scale=0.95)
        plt.xlim(-0.6, 1.25)
        plt.ylim(-1, len(self.data["names"]))

        coords = self.data.copy()
        coords.pop("names")
        coords.pop("units")
        arrow = 0.
        for c in coords.keys():
            plt.text(-0.6, arrow, c + " [" + self.data["units"][c] + "]")
            plot_info = shinescale_0_1(self.data[c])
            COLS = mpl_cm.gist_rainbow(np.linspace(0, 1, len(plot_info[0])))
            plt.scatter(plot_info[1],
                        [arrow] * len(plot_info[1]),
                        color="black",
                        marker='|',
                        s=1000,
                        alpha=1.)
            for l, v in enumerate(plot_info[1]):
                plt.text(v, arrow + 0.3, str(plot_info[2][l]), ha="center")
            plt.scatter(plot_info[0],
                        [arrow - 0.15] * len(plot_info[0]),
                        color=COLS,
                        marker='^',
                        s=200,
                        alpha=0.5)
            arrow += 1

        for ind, lab in enumerate(self.data["names"]):
            plt.plot(-100, -100, color=COLS[ind], marker="^",
                     markersize=10, linewidth=0, alpha=0.5, label=lab)

        ax.axis("off")

        box = ax.get_position()
        ax.set_position([box.x0, box.y0,
                         box.width, box.y1 - box.height * 0.1])

        # Put a legend below current axis
        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=3, mode="expand", borderaxespad=0.)

        plt.tight_layout()
        return
