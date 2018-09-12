"""
Plot tool for stratosphere.

Module for functions useful for image manipulation and plotting.
"""

import matplotlib.colors as mcol
import numpy as np
import scipy.interpolate as spint


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
    colors_i = np.linspace(0, 1., numcol + 1)
    # Set array to keep colors in
    colors = np.zeros((numcol + 1, 3), float)
    for (i, key) in enumerate(('red', 'green', 'blue')):
        # Find the numcol colors
        if callable(cdict[key]):
            print('ERR>> Colour map ({}) defined by continuous functions,' +
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
