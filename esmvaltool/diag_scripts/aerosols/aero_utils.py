"""Part of the ESMValTool Aerosol diagnostics.

This module contains utility functions commonly used by aerosol
assessment routines.
"""

import iris
import numpy as np


class AeroAnsError(Exception):
    """Exception class for errors raised when model data is checked in the
      extract_pt module.
    """
    pass


def add_bounds(cube):
    """Add bounds to a cubes latitude and longitude coordinates.

    Parameters
    ----------
    cube : Iris cube
        Iris cube with latitude and longitude coordinates.

    Returns
    -------
    cube : cube
        Iris cube with bounds added to the latitude and longitude coordinates.
    """

    if not cube.coord('latitude').has_bounds():
        cube.coord('latitude').guess_bounds()
    if not cube.coord('longitude').has_bounds():
        cube.coord('longitude').guess_bounds()

    return cube


def extract_pt(icube, pt_lat, pt_lon, **kwargs):
    """Extracts given location(s) (3-D) from a cube.

    Method
    ------
    Uses Iris module Analysis.Interpolate to extract values
    first based on horizontal coordinates and later based on
    height, if specified.

    If height ('altitude') is requested, checks if cube heights
    include orography, i.e. HybridHeights have been derived.

    Parameters
    ----------
    icube : Iris cube
    pt_lat, pt_lon : latitude and longitude coordinates of desired points
    kwargs: (optional)
        height  : Altitude (above geoid) of point.
        level   : Model/ pseudo/ tile number, default = all.
        nearest : Boolean. Specify whether to use 'nearest neighbour', instead
                  of 'linear' method while extracting data. Default is False.

    Returns
    -------
    data_out : List
        List of single point values, corresponding to each point specified

    Raises
    ------
    AeroAnsError : If the numbers of latitude and longitude points are
        mismatched. OR if both and level and height are passed as kwargs.
        OR if the cude contains a time coordinate. OR if a pseudo level
        coordinate is requested, but not present in the cube. OR if the numbers
        of latitude/longitude and height points are mismatched. OR if height
        is requested but the cube does not contain and altitude coordinate.
    """

    # Check that input data is a (single) cube
    if not isinstance(icube, iris.cube.Cube):
        raise AeroAnsError('Extract_pt:First argument must be a single cube')

    # Check that equal number of lat/lon pairs are passed in point coordinates
    # convert arguments to lists as easier to process later
    pt_lat1 = []
    pt_lon1 = []

    if not isinstance(pt_lat, list):
        pt_lat1.append(pt_lat)
        pt_lon1.append(pt_lon)

    else:
        for n in np.arange(len(pt_lat)):
            pt_lat1.append(pt_lat[n])
            pt_lon1.append(pt_lon[n])

    if len(pt_lat1) != len(pt_lon1):
        raise AeroAnsError('Extract_pt:Mismatch in number of lat/long values')

    if 'level' in kwargs and 'height' in kwargs:
        raise AeroAnsError('Extract_pt: Both Level and Height requested')

    if icube.coords()[0].name() == 'time':  # only lat,lon,lev accepted
        raise AeroAnsError(
            'Extract_pt:Cannot handle time dimension at present')

    if 'level' in kwargs and not icube.coord(
            'model_level_number') and not icube.coord('pseudo_level'):
        raise AeroAnsError('Extract_pt:Level requested, but not found in cube')

    if 'height' in kwargs:
        pt_hgt = []

        if not isinstance(kwargs['height'], list):
            pt_hgt.append(kwargs['height'])

        else:
            pt_hgt.extend(kwargs['height'])

        if len(pt_lat1) != len(pt_hgt):
            raise AeroAnsError(
                'Extract_pt:Mismatch in number of points for lat/long/height')

        # Check that heights have been merged with orography.
        if not icube.coords('altitude'):
            raise AeroAnsError(
                'Extract_pt:Height requested but input data does not contain \
                "Altitude" coordinate')

        # Store the min and max altitudes from cube data so that user
        # cannot request points located below/ above that.
        # Will extract =min/max if beyond limits.
        hgt_min = icube.coord('altitude').points.min()
        hgt_max = icube.coord('altitude').points.max()

    # ---------- Finished checks -- begin processing -------------------------

    # If level specified, extract slice first
    if 'level' in kwargs:

        if icube.coord('model_level_number'):
            icube = icube.extract(
                iris.Constraint(model_level_number=kwargs['level']))

        else:

            if icube.coord('pseudo_level'):

                icube = icube.extract(
                    iris.Constraint(pseudo_level=kwargs['level']))

    # Extract values for specified points lat/lon

    # ******** Does not seem to handle multiple points if 3-D *************

    data_out = []

    for n in np.arange(len(pt_lat1)):
        pt_coords = [('latitude', pt_lat1[n]), ('longitude', pt_lon1[n])]

        if 'nearest' in kwargs and kwargs['nearest']:
            tcube = icube.interpolate(pt_coords, iris.analysis.Nearest())

        else:
            tcube = icube.interpolate(pt_coords, iris.analysis.Linear())

        # Extract at requested height

        if 'height' in kwargs:
            pt = max(pt_hgt[n], hgt_min)
            pt = min(pt_hgt[n], hgt_max)
            pt_coords = [('altitude', pt)]

            if 'nearest' in kwargs and kwargs['nearest']:
                tcube = tcube.interpolate(pt_coords, iris.analysis.Nearest())

            else:
                tcube = tcube.interpolate(pt_coords, iris.analysis.Linear())

        data_out.append(tcube.data)

    return data_out
