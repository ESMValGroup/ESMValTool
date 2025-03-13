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


def add_bounds(cube):
    """Add bounds to a cube's latitude and longitude coordinates.

    Parameters
    ----------
    cube : Iris cube
        Iris cube with latitude and longitude coordinates.

    Returns
    -------
    cube : Iris cube.
        Iris cube with bounds added to the latitude and longitude coordinates.
    """

    if not cube.coord("latitude").has_bounds():
        cube.coord("latitude").guess_bounds()
    if not cube.coord("longitude").has_bounds():
        cube.coord("longitude").guess_bounds()

    return cube


def extract_pt(icube, pt_lat, pt_lon, height=None, level=None, nearest=False):
    """Extracts given location(s) (3-D) from a cube.

    Method
    ------
    Uses Iris module Analysis.Interpolate to extract values,
    initially based on horizontal coordinates, and then based on
    height, if specified.

    If height ('altitude') is requested, checks if cube heights
    include orography, i.e. HybridHeights have been derived.

    Parameters
    ----------
    icube : Iris cube
    pt_lat, pt_lon : Float or list/array of floats. Latitude and longitude
                     coordinates of desired points.
    args:
        height  : Float or list/array of floats. Altitude (above geoid) of
                  point. Initialized to None.
        level   : Integer . Model level or pseudo level or tile number.
                  Initialized to None, meaning that all available levels in
                  the cube are used.
        nearest : Boolean. Specify whether to use 'nearest neighbour', instead
                  of 'linear' method while extracting data. Default is False.

    Returns
    -------
    data_out : List
        List of single point values, corresponding to each point specified.

    Raises
    ------
    AeroAnsError : If the number of latitude and longitude points are
        mismatched. OR if both and level and height are passed as args.
        OR if the cube contains a time coordinate. OR if a pseudo level
        coordinate is requested, but not present in the cube. OR if the numbers
        of latitude/longitude and height points are mismatched. OR if height
        is requested but the cube does not contain an altitude coordinate.
    """

    # Check that input data is a (single) cube
    if not isinstance(icube, iris.cube.Cube):
        raise AeroAnsError("Extract_pt:First argument must be a single cube")

    # Check if the cube contains a time dimension, which is
    # currently unsupported.
    if icube.coords()[0].name() == "time":
        raise AeroAnsError(
            "Extract_pt:Cannot handle time dimension at present"
        )

    # Check that equal number of lat/lon pairs are passed in point coordinates.
    # Convert arguments to lists for easier processing if necessary.
    pt_lat1 = []
    pt_lon1 = []

    if not isinstance(pt_lat, list):
        pt_lat1.append(pt_lat)
        pt_lon1.append(pt_lon)

    else:
        for n_lat in np.arange(len(pt_lat)):
            pt_lat1.append(pt_lat[n_lat])
            pt_lon1.append(pt_lon[n_lat])

    if len(pt_lat1) != len(pt_lon1):
        raise AeroAnsError("Extract_pt:Mismatch in number of lat/long values")

    # Check that both level and height haven't been requested.
    if level is not None and height is not None:
        raise AeroAnsError("Extract_pt: Both Level and Height requested")

    # Check that the cube has a level coordinate if level has been requested.
    if (
        level is not None
        and not icube.coord("model_level_number")
        and not icube.coord("pseudo_level")
    ):
        raise AeroAnsError("Extract_pt:Level requested, but not found in cube")

    # Check that the number of height points is equal to the number of
    # lat/lon pairs. Convert the argument to a list for easier
    # processing if necessary.
    if height is not None:
        pt_hgt = []

        #        if isinstance(height, list):
        #            pt_hgt.extend(height)
        #        else:
        #            pt_hgt.append(height)
        pt_hgt.extend(height) if isinstance(height, list) else pt_hgt.append(
            height
        )

        if len(pt_lat1) != len(pt_hgt):
            raise AeroAnsError(
                "Extract_pt:Mismatch in number of points for lat/long/height"
            )

        # Check that heights have been merged with orography.
        if not icube.coords("altitude"):
            raise AeroAnsError(
                'Extract_pt:Height requested but input data does not contain \
                "Altitude" coordinate'
            )

        # Store the min and max altitudes from cube data so that user
        # cannot request points located below/ above that.
        # Will extract =min/max if beyond limits.
        hgt_min = icube.coord("altitude").points.min()
        hgt_max = icube.coord("altitude").points.max()

    # ---------- Finished checks -- begin processing -------------------------

    # If level specified, extract slice first
    if level is not None:
        try:
            icube = icube.extract(iris.Constraint(model_level_number=level))

        except Exception:
            print("Model level number not available. Use pseudo level")

        else:
            icube = icube.extract(iris.Constraint(pseudo_level=level))

    # Extract values for specified points lat/lon
    # NOTE: Does not seem to handle multiple points if 3-D
    data_out = []

    # Set lat/lon coordinates for model grid cell interpolation
    for n_lat1 in np.arange(len(pt_lat1)):
        latlon_coords = [
            ("latitude", pt_lat1[n_lat1]),
            ("longitude", pt_lon1[n_lat1]),
        ]

        if nearest:
            tcube = icube.interpolate(latlon_coords, iris.analysis.Nearest())
        else:
            tcube = icube.interpolate(latlon_coords, iris.analysis.Linear())

        # If height specified, interpolate to requested height
        if height is not None:
            # Set vertical coordinates for model grid cell interpolation
            point = max(pt_hgt[n_lat1], hgt_min)
            point = min(pt_hgt[n_lat1], hgt_max)
            hgt_coords = [("altitude", point)]

            if nearest:
                tcube = tcube.interpolate(hgt_coords, iris.analysis.Nearest())
            else:
                tcube = tcube.interpolate(hgt_coords, iris.analysis.Linear())

        # Append processed data point
        data_out.append(tcube.data)

    return data_out
