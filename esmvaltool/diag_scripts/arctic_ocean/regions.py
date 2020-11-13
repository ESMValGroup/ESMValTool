# -*- coding: utf-8 -*-
"""Part of the ESMValTool Arctic Ocean diagnostics.

This module contains functions with definitions of regions.
"""
import logging
import os
import numpy as np
from scipy.interpolate import interp1d

logger = logging.getLogger(os.path.basename(__file__))


def hofm_regions(region, lon2d, lat2d):
    """Define regions for data selection.

    Parameters
    ----------
    region: str
        the name of the region
    lon2d: 2d numpy array
    lat2d: 2d numpy array

    Returns
    -------
    indexesi: 1d numpy array
        i indexes of the selected points
    indexesj: 1d numpy array
        j indexes of the selected points
    """
    if region == 'EB':
        # Eurasian Basin of the Arctic Ocean
        indi, indj = np.where((lon2d > 300) & (lat2d > 80))
        indi2, indj2 = np.where((lon2d < 100) & (lat2d > 80))
        indi3, indj3 = np.where((lon2d > 100) & (lon2d < 140) & (lat2d > 66))

        indexesi = np.hstack((indi, indi2, indi3))
        indexesj = np.hstack((indj, indj2, indj3))
    elif region == 'AB':
        # Amerasian basin of the Arctic Ocean
        indi, indj = np.where((lon2d >= 260) & (lon2d <= 300) & (lat2d >= 80))
        indi2, indj2 = np.where((lon2d >= 140) & (lon2d < 260) & (lat2d > 66))

        indexesi = np.hstack((indi, indi2))
        indexesj = np.hstack((indj, indj2))
    elif region == 'Barents_sea':

        indi, indj = np.where((lon2d >= 20) & (lon2d <= 55) & (lat2d >= 70)
                              & (lat2d <= 80))

        indexesi = indi
        indexesj = indj
    elif region == 'North_sea':
        # Amerasian basin of the Arctic Ocean
        indi, indj = np.where((lon2d >= 355) & (lon2d <= 360) & (lat2d >= 50)
                              & (lat2d <= 60))
        indi2, indj2 = np.where((lon2d >= 0) & (lon2d <= 10) & (lat2d >= 50)
                                & (lat2d <= 60))

        indexesi = np.hstack((indi, indi2))
        indexesj = np.hstack((indj, indj2))
    else:
        print('Region {} is not recognized'.format(region))
    return indexesi, indexesj


def transect_points(transect, mult=2):
    """Return a collection of points for transect.

    Parameters
    ----------
    transect: str
        Name of the predefined transect
    mult: int
        multiplicator that allow to increase
        the number of points by `mult` times

    Returns
    -------
    lon_s4new: 1d numpy array
        longitude points of the transect
    lat_s4new: 1d numpy array
        latitude points of the transect
    """
    if transect == 'AWpath':
        lon_s4 = np.array([
            17.6, 16.5, 16.05, 15.6, 15.1, 14.1, 13.0, 12.0, 10.0, 8.0, 4.0,
            4.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0,
            110.0, 120.0, 130.0, 140.0
        ])
        lat_s4 = np.array([
            69.0, 70.6, 71.3, 72.02, 72.8, 73.8, 75.0, 76.0, 77.0, 78.0, 79.0,
            80.0, 81.0, 81.8, 81.8, 82.6, 83.0, 83.2, 83.1, 82.8, 82.5, 81.8,
            79.7, 78.2, 78.7, 79.7
        ])
    elif transect == 'Fram':
        lon_s4 = np.array([
            -22.3732, -21.4796, -19.6479, -18.1074, -16.7828, -15.504,
            -14.2042, -12.9771, -11.6642, -10.1892, -8.7414, -7.719, -6.3646,
            -4.4803, -3.4232, -2.4435, -1.615, -0.6752, 0.343, 1.6947, 2.7157,
            3.7374, 4.6099, 5.5097, 6.3754, 7.2394, 7.9238, 8.7029, 9.7338,
            10.4462, 11.0559, 12.0102, 13.3313
        ])
        lat_s4 = np.array([
            78.9373, 78.9276, 78.9183, 78.9356, 78.9346, 78.9334, 78.9425,
            78.9434, 78.9274, 78.9392, 78.9287, 78.9262, 78.9392, 78.95,
            78.9405, 78.9347, 78.9334, 78.922, 78.9287, 78.9131, 78.919,
            78.9215, 78.9242, 78.909, 78.8995, 78.8874, 78.8865, 78.9026,
            78.8992, 78.8841, 78.8793, 78.8715, 78.9012
        ])

    else:
        print('Transect {} is not recognized'.format(transect))

    point_number = np.linspace(1,
                               lon_s4.shape[0],
                               num=lon_s4.shape[0],
                               endpoint=True)
    f_lons = interp1d(point_number, lon_s4)
    g_lats = interp1d(point_number, lat_s4)
    xnew = np.linspace(1,
                       lon_s4.shape[0],
                       num=mult * lon_s4.shape[0],
                       endpoint=True)

    lon_s4new = f_lons(xnew)
    lat_s4new = g_lats(xnew)
    return lon_s4new, lat_s4new
