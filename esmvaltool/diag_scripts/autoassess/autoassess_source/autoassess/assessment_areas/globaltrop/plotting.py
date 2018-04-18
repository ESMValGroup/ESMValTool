'''
Module for functions useful for image manipulation and plotting
'''

import os
import subprocess as subp

import matplotlib.colors as mcol
import matplotlib.pyplot as mplt
import numpy as np
import numpy.ma as ma
import scipy.interpolate as spint

import iris.analysis
import iris.analysis.cartography


def segment2list(cmap, numcol, reverse=False):
    """
    Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        numcol: Number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = segment2list(cm.jet, 5)
        imshow(x, cmap=djet)

    Routine takes a LinearSegmentColormap (continuous colormapping) and
    converts it to a ListColormap (discrete color list)
    This is used so that the use of this with BoundaryNorm mimics default
    behaviour for linear colour scales, and actually get required behaviour
    for non-linear colour scales.
    """

    cdict = cmap._segmentdata.copy()
    # numcol colors
    colors_i = np.linspace(0, 1., numcol+1)
    # Set array to keep colors in
    colors = np.zeros((numcol+1, 3), float)
    for (i, key) in enumerate(('red', 'green', 'blue')):
        # Find the numcol colors
        if callable(cdict[key]):
            print ('ERR>> Colour map ({}) defined by continuous functions,' +
                   ' cannot yet convert to list').format(cmap.name)
            return None
            # This doesn't work as end up with negative values in colors[]
            # for j,icol in enumerate(colors_i):
            #     colors[j, i]=cdict[key](icol)
        else:
            xy_int = np.array(cdict[key])
            interpolator = spint.interp1d(xy_int[:, 0], xy_int[:, 1])
            colors[:, i] = interpolator(colors_i)

    if reverse:
        colors = colors[::-1, :]

    # Return colormap object.
    colormap = mcol.ListedColormap(colors[1:-1, :], name='colormap')

    # Use zeroth color as negative extreme
    colormap.set_under(color=colors[0, :])
    # Use last color as positive extreme
    colormap.set_over(color=colors[-1, :])

    return colormap


def has_ll(coord_group, coords):
    """
    Routine used to check if any coordinate names in coords are in
    coord_group
    """
    # Does 'any' need to be 'all' in cases where coords do not match up?
    return any(co.name() in coords for co in coord_group if co is not None)


def regridcube(icube, gcube, debug=False):
    """
    Function to regrid a cube to that given by example cube
    - Suitable for both maps and zonal x-sections
    - Will only regrid latitude and longitude, what if other coords in
      compare['not_equal']
    -- Can I extend to any horizontal coordinate?
    """

    # Compare field1 and field2 cubes
    compare = iris.analysis.coord_comparison(icube, gcube)
    if debug:
        print compare['not_equal']

    # Keep copy of original fields in case no interpolation required
    rcube = icube.copy()

    # If longitude or latitude do not match then interpolate field
    ll_coords = ['latitude', 'longitude']
    if any(has_ll(grp, ll_coords) for grp in compare['not_equal']):

        # Set empty list to hold new lat/lon points
        samplepoints = []

        # Add coordinate values for coords that need regridding
        if any(has_ll(grp, ['longitude']) for grp in compare['not_equal']):
            # If scalar then just replace icube coordinate with cube coordinate
            if gcube.coord('longitude').points.size == 1:
                icube.replace_coord(gcube.coord('longitude'))
            else:
                samplepoints.append(('longitude',
                                     gcube.coord('longitude').points))
        if any(has_ll(grp, ['latitude']) for grp in compare['not_equal']):
            # If scalar then just replace icube coordinate with cube coordinate
            if gcube.coord('latitude').points.size == 1:
                icube.replace_coord(gcube.coord('latitude'))
            else:
                samplepoints.append(('latitude',
                                     gcube.coord('latitude').points))

        if samplepoints:
            # Perform interpolation
            # - Note that linear interpolation does not respect MDI (Need to
            #   continously test this)
            # - Note that iris.analysis.interpolate.regrid only does horizontal
            #   maps
            rcube = iris.analysis.interpolate.linear(icube, samplepoints)

            # Set any values close to MDI (1e20) to be MDI
            if ma.isMaskedArray(icube.data):
                rcube.data = ma.masked_values(rcube.data,
                                              icube.data.fill_value,
                                              rtol=0.9)
        else:
            rcube = icube

    return rcube


def flddiff(fld1, fld2, rev=False, res_check=True, error=False, debug=False):
    """
    Routine to calculate field differences
    - Will cope with fields on different horizontal grids, may need to deal
      with vertical grids?
    - Do I want to try to determine which field to interpolate? e.g. Always
      ND onto EG
    -- Current test ignores bounded coordinates
    -- Tries to avoid extrapolation
    -- May need to protect this in case user wants finer control of what is
       going on
    """

    # Keep copy of original fields in case no interpolation required
    newf1 = fld1.copy()
    newf2 = fld2.copy()

    if res_check:
        # Need only test on latitude as longitude is a circular coordinate
        if fld1.coords('latitude') and fld2.coords('latitude'):
            if not fld1.coord('latitude').has_bounds() and \
               not fld2.coord('latitude').has_bounds():

                # First check if extreme latitudinal points are the same, if
                # so then check on number of grid points. Lower resolution
                # should be target grid(?)
                if ((fld1.coord('latitude').points[0] ==
                     fld2.coord('latitude').points[0]) and
                    (fld1.coord('latitude').points[-1] ==
                     fld2.coord('latitude').points[-1])):
                    # Default assumes fld2 lower resolution, change if fld1 is
                    # actually lower resolution
                    if fld1.coord('latitude').shape[0] < \
                       fld2.coord('latitude').shape[0]:
                        rev = True
                else:
                    # This test asks if extreme latitude points for fld1 are
                    # inside fld2, if so then need to reverse regrid in order
                    # to avoid extrapolation
                    # - Assumes increasing latitude coordinate
                    if ((fld1.coord('latitude').points[0] >=
                         fld2.coord('latitude').points[0]) and
                        (fld1.coord('latitude').points[-1] <=
                         fld2.coord('latitude').points[-1])):
                        rev = True

    # Regrid fields
    if rev:
        newf2 = regridcube(fld2, fld1, debug=debug)
    else:
        newf1 = regridcube(fld1, fld2, debug=debug)

    # Do difference
    # - some reason will not ignore forecast_reference_time differences
    # - sometimes this will be used and forecast_reference_time will be
    #   different (e.g. Forecast Errors)
    # - will need to test this is still the case from time to time
    if newf1.coords('forecast_reference_time') and \
       newf2.coords('forecast_reference_time'):
        newf2.coord('forecast_reference_time').points = \
            newf1.coord('forecast_reference_time').points
    diff = newf1 - newf2

    # Standard name is not maintained, and attributes are lost in subtraction
    # above. Need to be added
    # - We may want to do this in a better manner, e.g. is there a simple
    #   copy metadata from one cube to another?
    diff.rename(fld1.name())
    diff.attributes = fld1.attributes

    # Need to re-add forecast_period coordinate since it is 'ignored' above
    # - Use add_aux_coord as dimension is scalar, even though a DimCoord
    #   Don't have to use data_dim keyword in this case
    if not diff.coords('forecast_period'):
        if fld1.coords('forecast_period'):
            diff.add_aux_coord(fld1.coord('forecast_period'))

    # Need to re-add height coordinate since it is 'ignored' above
    # - Use add_aux_coord as dimension is scalar, even though a DimCoord
    #   Don't have to use data_dim keyword in this case
    if fld1.coords('height') and not diff.coords('height'):
        diff.add_aux_coord(fld1.coord('height'))

    return diff
