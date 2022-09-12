"""Part of the ESMValTool Aerosol diagnostics.

This module contains utility functions commonly used by aerosol
assessment routines.
"""

import os

import iris
import numpy as np


class AeroAnsError(Exception):
    pass


# Global variables

# STASH - name mapping file
PATH_AERO = os.path.abspath(os.path.dirname(__file__))

seasons = ["djf", "mam", "jja", "son"]
nseas = len(seasons)

mm_dust = 0.100  # Mol mass dust(kg mol-1) - see UKCA_MODE_SETUP: 'mm' arrays
earth_g = 9.80665  # m/s^2


def string_stc(stash):
    """Convert stashcode to mXXsYYiZZZ format.

    Parameters
    ----------
    stash : str
        UM stashcode with format YYZZZ.

    Returns
    -------
    : str
        UM stashcode with format mXXsYYiZZZ.
    """

    stash = int(stash)
    if stash < 1 or stash > 99999:
        raise AeroAnsError("STR_STC: Invalid stashcode " + str(stash))
    return "m01s{0:02d}i{1:03d}".format(int(stash / 1000), int(stash % 1000))


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

    if not cube.coord("latitude").has_bounds():
        cube.coord("latitude").guess_bounds()
    if not cube.coord("longitude").has_bounds():
        cube.coord("longitude").guess_bounds()
    return cube


def get_wavel_index(wavel):
    """Return the pseudo-level index for given wavelength in AOD diags.

    Parameters
    ----------
    wavel : str
        A given wavelength.

    Returns
    -------
    pseudo-level index : str
        The pseudo-level index corresponding to the given wavelength.
    """
    # Wavelengths (nm)
    cwavel = ["380", "440", "500", "675", "870", "1020"]
    try:
        return cwavel.index(str(wavel)) + 1  # for 1: indexing in files
    except Exception:
        raise AeroAnsError("Get_Wavelen_Index: requested Wavelength not found")


def extract_pt(icube, pt_lat, pt_lon, **keys):
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
    keys: (optional)
        height  : altitude (above geoid) of point
        level   : model/ pseudo/ tile number, default = all
        nearest : (optional) Whether to use 'nearest neighbour', instead
                  of 'linear' method while extracting data.

    Returns
    -------
    data_out : List
        List of single point values, corresponding to each point specified.
    """

    # Check that input data is a (single) cube
    if not isinstance(icube, iris.cube.Cube):
        raise AeroAnsError("Extract_pt:First argument must be a single cube")

    # Convert arguments to lists as easier to process later
    pt_lat1 = []
    pt_lon1 = []

    # Check that equal number of lat/lon pairs are passed in point coordinates
    if not isinstance(pt_lat, list):
        pt_lat1.append(pt_lat)
        pt_lon1.append(pt_lon)

    else:
        for n in np.arange(len(pt_lat)):
            pt_lat1.append(pt_lat[n])
            pt_lon1.append(pt_lon[n])

    if len(pt_lat1) != len(pt_lon1):
        raise AeroAnsError("Extract_pt:Mismatch in number of lat/long values")

    if "level" in keys and "height" in keys:
        raise AeroAnsError("Extract_pt: Both Level and Height requested")

    if icube.coords()[0].name() == "time":  # only lat,lon,lev accepted
        raise AeroAnsError(
            "Extract_pt:Cannot handle time dimension at present")

    if ("level" in keys and not icube.coord("model_level_number")
            and not icube.coord("pseudo_level")):
        raise AeroAnsError("Extract_pt:Level requested, but not found in cube")

    if "height" in keys:
        pt_hgt = []
        if not isinstance(keys["height"], list):
            pt_hgt.append(keys["height"])
        else:
            for n in np.arange(len(keys["height"])):
                pt_hgt.append(keys["height"][n])

        if len(pt_lat1) != len(pt_hgt):
            raise AeroAnsError("Extract_pt:Mismatch in number of points for\
           lat/long/height")

        # Check that heights have been merged with orography.
        if not icube.coords("altitude"):
            raise AeroAnsError(
                'Extract_pt:Height requested but input data does\
           not contain "Altitude" coordinate')
        else:
            # Store the min and max altitudes from cube data so that user
            # cannot request points located below/ above that.
            # Will extract =min/max if beyond limits.
            hgt_min = icube.coord("altitude").points.min()
            hgt_max = icube.coord("altitude").points.max()

    # ---------- Finished checks -- begin processing --------------------------
    # If level specified, extract slice first
    if "level" in keys:
        if icube.coord("model_level_number"):
            icube = icube.extract(
                iris.Constraint(model_level_number=keys["level"]))
        else:
            if icube.coord("pseudo_level"):
                icube = icube.extract(
                    iris.Constraint(pseudo_level=keys["level"]))

    # Extract values for specified points lat/lon
    #   ******** Does not seem to handle multiple points if 3-D *************
    data_out = []
    for n in np.arange(len(pt_lat1)):
        pt_coords = [("latitude", pt_lat1[n]), ("longitude", pt_lon1[n])]
        if "nearest" in keys and keys["nearest"]:
            tcube = icube.interpolate(pt_coords, iris.analysis.Nearest())
        else:
            tcube = icube.interpolate(pt_coords, iris.analysis.Linear())

        # Extract at requested height
        if "height" in keys:
            pt = max(pt_hgt[n], hgt_min)
            pt = min(pt_hgt[n], hgt_max)
            pt_coords = [("altitude", pt)]
            if "nearest" in keys and keys["nearest"]:
                tcube = tcube.interpolate(pt_coords, iris.analysis.Nearest())
            else:
                tcube = tcube.interpolate(pt_coords, iris.analysis.Linear())

        data_out.append(tcube.data)

    return data_out


def calc_mass(pstar, prho, ptheta_top, mmr, nlev):
    """Calculates mass in model layer from mass mixing ratios.

    Method
    ------
    Simply calculates mass in layer using mmr*rho*dz = mmr*g*(-dp)
       - at lowest level uses dp = p2 - pstar
       - at highest level uses dp = p_theta(top) - p_rho(top_rho_lev)
    Orig/IDL: Stephanie Woodward, 2003

    Parameters
    ----------
    pstar :
    prho :
    ptheta :
    mmr : Iris cube
       Iris cube containing mass mixing ratio data.
    nlev : int (optional)

    Returns
    -------
    tmp : Iris cube
        Iris cube containing mass data.
    """

    tmp = mmr.copy()
    tmp.data[:, :, :] = 0.0
    tmp[0, :, :] = mmr[0, :, :] * (pstar - prho[1, :, :]) / earth_g

    # Loop over model levels. Range will iterate up to nlev-2
    for i in range(1, nlev - 1):
        tmp[i, :, :] = mmr[i, :, :] * (prho[i, :, :] -
                                       prho[i + 1, :, :]) / earth_g

    tmp[nlev - 1, :, :] = (mmr[nlev - 1, :, :] *
                           (prho[nlev - 1, :, :] - ptheta_top) / earth_g)

    return tmp
