"""
Module for getting land fraction on correct resolution for
global hydrocycle metrics
"""

import csv
import os.path
import sys

import numpy as np

import iris
import iris.coord_categorisation


def read_info_file(infofile):
    """
    Reads the info file and returns the data held in it as a single dictionary
    """

    # Initialise variables
    info_dict = {}

    # Read the valorder file as a csv file
    with open(infofile, 'r') as inf:
        f = csv.reader(inf, delimiter=':')
        for line in f:
            # Ignore any blank lines or those beginning with a hash
            if len(line) > 0:
                if line[0][0] != '#':
                    info_dict[line[0].strip()] = line[1].strip()

    return info_dict


def fix_lat_coord_bounds(lat_coord):
    """
    function to fix the bounds on a latitude coordinate
     - in case bounds go across poles
    """

    LATITUDE_VALID_RANGE = np.array((-90.0, 90.0))
    b = lat_coord.bounds.copy()
    inds = np.where((b > LATITUDE_VALID_RANGE.max()) |
                    (b < LATITUDE_VALID_RANGE.min()))

    if len(inds) == 0:
        return
    if len(inds) != 2:
        raise Exception("Unexpected number of bounds on a cell")
    for i, j in zip(inds[0], inds[1]):
        b[i, j] = max(LATITUDE_VALID_RANGE.min(), b[i, j])
        b[i, j] = min(LATITUDE_VALID_RANGE.max(), b[i, j])
    lat_coord.bounds = b


def guess_dim(cube):
    for co in cube.coords():
        if "lat" in co.name():
            lat_coord = co.name()
        if "lon" in co.name():
            lon_coord = co.name()
    return lat_coord, lon_coord


def guessbounds(cube):
    tnam = guess_dim(cube)
    for c in [tnam[0], tnam[1]]:
        cube.coord(c).guess_bounds()
        fix_lat_coord_bounds(cube.coord('latitude'))


def get_landfrac(incube):
    """To get the right land fraction file for the run resolution"""
    # TODO local file paths in to land_frac in 'extras_file.dat'.
    # TODO This does also not work for other models than the UM
    extras_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'extras_file.dat')
    extras_dict = read_info_file(extras_file)

    land_fraction_files_dict = {}
    for key in extras_dict.iterkeys():
        land_fraction_files_dict[key] = iris.load_cube(extras_dict[key])

    cube = incube.copy()
    guessbounds(cube)

    # What is the resolution?
    resolution = cube.coords(axis='x')[0].shape[0]/2.0

    # Load the correct resolution land fraction
    if resolution < 72:
        land_frac = land_fraction_files_dict['n48_land_frac']
    elif resolution >= 72 and resolution < 156:
        land_frac = land_fraction_files_dict['n96_land_frac']
    elif resolution >= 156 and resolution < 384:
        land_frac = land_fraction_files_dict['n216_land_frac']
    else:
        land_frac = land_fraction_files_dict['n512_land_frac']

    guessbounds(land_frac)

    incube_str = str(land_frac.shape) + \
        str(land_frac.coords(axis='x')[0].points[0]) + \
        str(land_frac.coords(axis='y')[0].points[0])
    newgrid_str = str(cube.shape) + \
        str(cube.coords(axis='x')[0].points[0]) + \
        str(cube.coords(axis='y')[0].points[0])

    # Regrid it if needed
    if incube_str != newgrid_str:
        scheme = iris.analysis.AreaWeighted(mdtol=0.5)
        regridder = scheme.regridder(land_frac, cube)
        land_frac_r = regridder(land_frac)
    else:
        land_frac_r = land_frac

    return land_frac_r

