#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 11:32:36 2017

@author: bmueller

Will be included in next geoval version

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from geoval.core.data import GeoData


class easyhov(object):
    """
    Basic class to implement simple hovmoeller plots from GeoData or
        3d-Matrix Objects
    """

    def __init__(self, **kwargs):
        super(easyhov, self).__init__(**kwargs)
        """
        Default values to experiment with the plot
        """
        self.__data__ = None
        self.__lat__ = True
        self.__data_hov__ = None

    def setdata(self, l_GeoData):
        """
        set the data
        """
        # TODO: check if data is GeoData and initialize securing issues
        self.__data__ = l_GeoData
        self.__calculate__()

    def getdata(self):
        """
        get the data
        """
        return self.__data__

    def __setHOVdata__(self, HOVdata, l_type):
        """
        set the hovmoeller data
        """
        self.__data_hov__ = HOVdata
        self.settype(l_type)

    def getHOVdata(self):
        """
        get the hovmoeller data
        """
        return self.__data_hov__

    def settype(self, l_type):
        """
        set the type
        """
        self.__lat__ = l_type == "lat"
        if self.__data__ is not None:
            self.__calculate__()

    def gettype(self):
        """
        get the type
        """
        return "lat" if self.__lat__ else "lon"

    def __init_trial__(self):
        """
        Initiates trial Geodata (For debugging purpose only)
        """
        self.setdata(GeoData("path/to/data.nc",
                             "variable", read=True))

    def __calculate__(self):
        """
        calculates the hovmoeller matrix (for climatologies only!)
        """
        # TODO: check if self.__data__ is GeoData
        if self.__data__ is not None:
            self.__data_hov__ = self.__data__.get_climatology(
                    return_object=True)
            self.__data_hov__ = self.__data_hov__.get_zonal_mean(
                return_object=False, lat=self.__lat__)
        else:
            assert False, "no data available!"

    def __set0mid__(self, Data):
        """
        adjusts GeoData to reach from -180° to 180° instead of 0° to 360°
        """
        def coord_reshape(D):
            coord = D.shape[1]
            D = np.hstack((D[:, np.arange(0, coord/2)+coord/2-1],
                           D[:, np.arange(0, coord/2)]))
            return(D)

        d_hov = coord_reshape(Data.data.copy())
        m_hov = coord_reshape(Data.mask.copy())

        Data = np.ma.array(d_hov, mask=m_hov)

        return(Data)

    def plot(self, minhov=None, maxhov=None, unit="",
             cmap='YlGn_r', title=None, ydesc=True, xdesc=True):
        """
        plots the hovmoeller matrix
        """
        if minhov is None:
            minhov = self.__data_hov__.min()
        if maxhov is None:
            maxhov = self.__data_hov__.max()

        if self.__lat__:
            Data = np.flipud(np.rot90(self.__data_hov__))
        else:
            Data = self.__set0mid__(self.__data_hov__)

        the_shape = Data.shape[0]

        imgplot = plt.imshow(Data,
                             cmap=cm.get_cmap(cmap, 10))
        ax = plt.gca()

        if self.__lat__:
            the_shape = Data.shape[0]
            if ydesc:
                plt.ylabel('latitude [deg]')
            if xdesc:
                plt.xlabel('climatology [months]')
            plt.xticks(range(0, 12), range(1, 13))
            plt.yticks([-0.5, the_shape*1/6, the_shape/3, the_shape/2,
                        the_shape*2/3, the_shape*5/6, the_shape-0.5],
                       ["90 N", "60 N", "30 N", "0", "30 S", "60 S", "90 S"])
            ax.set_xticks(np.arange(-.5, 12, 1), minor=True)

        else:
            the_shape = Data.shape[1]
            if ydesc:
                plt.ylabel('climatology [months]')
            plt.yticks(range(0, 12), range(1, 13))
            if xdesc:
                plt.xlabel('longitude [deg]')
            plt.xticks([-0.5, the_shape/4, the_shape/2,
                        3*the_shape/4, the_shape-0.5],
                       ["180 W", "90 W", "0", "90 E", "180 E"])
            ax.set_yticks(np.arange(-.5, 12, 1), minor=True)

        cbar = plt.colorbar(ticks=np.linspace(minhov,
                                              maxhov,
                                              num=11))
        cbar.ax.set_title(unit, ha="right", x=1.1)

        imgplot.set_clim((minhov, maxhov))

        plt.axis('tight')
        ax.grid(which='minor', color='black', linestyle='-', linewidth=1)
        plt.title(title, ha="left", x=-0.1, y=1.05)

        return(imgplot)


class easyhov_diff(object):
    """
    Class to implement simple hovmoeller difference plots from GeoData or
        3d-Matrix Objects
    """
    def __init__(self, **kwargs):
        super(easyhov_diff, self).__init__(**kwargs)
        """
        Default values to experiment with the plot
        """
        self.__data__ = [None, None]
        self.__lat__ = True
        self.__data_hov__ = [None, None, None]

    def setdata(self, l_GeoDatas):
        """
        set the data
        """
        # TODO: check size and if data is GeoDatas
        # and initialize securing issues
        A = easyhov()
        B = easyhov()
        A.setdata(l_GeoDatas[0])
        B.setdata(l_GeoDatas[1])
        A.settype("lat" if self.__lat__ else "lon")
        B.settype("lat" if self.__lat__ else "lon")
        self.__data__ = [A, B]
        self.__calculate__()

    def __update_type__(self):
        """
        update the type
        """
        for x in self.__data__:
            if x is not None:
                x.settype("lat" if self.__lat__ else "lon")

    def settype(self, l_type):
        """
        set the type
        """
        self.__lat__ = l_type == "lat"
        self.__update_type__()
        if all([x is not None for x in self.__data__]):
            self.__calculate__()

    def gettype(self):
        """
        get the type
        """
        return "lat" if self.__lat__ else "lon"

    def getdata(self):
        """
        get the data
        """
        return self.__data__

    def getHOVdata(self):
        """
        get the hovmoeller data
        """
        return self.__data_hov__

    def __init_trial__(self):
        """
        Initiates trial Geodatas (For debugging purpose only)
        """

        X = GeoData("/path/to/data.nc",
                    "variable", read=True)

        Y = X.copy()
        Y.data = Y.data * 2

        self.setdata([X, Y])

    def __calculate__(self):
        """
        calculates the hovmoeller matrices and differences
        (for climatologies only!)
        """

        self.__data_hov__ = []

        for x in self.__data__:
            if x is not None:
                x.__calculate__()
                self.__data_hov__.append(x.getHOVdata())
            else:
                assert False, "no data available!"

        self.__data_hov__.append(self.__data_hov__[0]-self.__data_hov__[1])

        self.__data_hov__ = dict(zip(["A", "B", "A-B"],
                                     self.__data_hov__))

    def plot(self, minhov=None, maxhov=None, unit="",
             mindiff=None, maxdiff=None,
             cmap='YlGn_r', cmapdiff='RdBu', titles=["A", "A-B", "B"]):
        """
        plots the hovmoeller matrices and differences
        """
        # TODO: check all the sizes

        if minhov is None:
            minhov = min([X.getHOVdata().min() for X in self.__data__])
        if maxhov is None:
            maxhov = max([X.getHOVdata().max() for X in self.__data__])

        diff_span = maxhov - minhov
        diff_span = [diff_span * sg for sg in [-1, 1]]

        if mindiff is None:
            mindiff = min(diff_span)
        if maxdiff is None:
            maxdiff = max(diff_span)

        HOV = plt.figure()

        HOV.add_subplot(1, 3, 1)
        here = "A"
        A = easyhov()
        A.__setHOVdata__(self.__data_hov__[here],
                         "lat" if self.__lat__ else "lon")
        A.plot(title=titles[0], ydesc=True, xdesc=False,
               minhov=minhov, maxhov=maxhov, cmap=cmap, unit=unit)

        HOV.add_subplot(1, 3, 2)
        here = "A-B"
        B = easyhov()
        B.__setHOVdata__(self.__data_hov__[here],
                         "lat" if self.__lat__ else "lon")
        B.plot(title=titles[1], ydesc=False, xdesc=True,
               minhov=mindiff, maxhov=maxdiff, cmap=cmapdiff)

        HOV.add_subplot(1, 3, 3)
        here = "B"
        C = easyhov()
        C.__setHOVdata__(self.__data_hov__[here],
                         "lat" if self.__lat__ else "lon")
        C.plot(title=titles[2], ydesc=False, xdesc=False,
               minhov=minhov, maxhov=maxhov, cmap=cmap, unit=unit)

        return(HOV)
