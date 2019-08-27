import numpy as np

from iris.analysis import SUM


def low_pass_weights(window, freq):
    """
    Calculate weights for a low pass Lanczos filter.

    Arguments:

    * window: int
        The length of the filter window in days.

    * freq: str
        The frequency of the data.

    Notes: (MSM)
        From "http://www.scitools.org.uk/iris/docs/latest/examples/"
        "graphics/SOI_filtering.html"
        Similar to LN's LanczosWeights function, but arrays are two elements
        shorter.
        LN's LanczosWeights returns very small numbers for the extra elements
        so I'm happy to stick with this implementation.
        Fractional differences in the resulting values are of the order
        of 1e-16
    """
    if freq.endswith('hr'):
        hours = int(freq[:-2])
        window = window * 24 // hours

    cutoff = 1. / window

    order = ((window - 1) // 2) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
    w[n-1:0:-1] = firstfactor * sigma
    w[n+1:-1] = firstfactor * sigma

    weights = w[1:-1]
    weights /= weights.sum()

    return weights


def lanczos_filter(cube, weights):
    """
    Perform Lanczos filtering

    Also fix bounds on time coordinates to remove overlapping bounds

    Arguments:
     * cube:
        iris cube to filter

     * weights:
        numpy 1D array of weights to do the filtering with.
    """
    cube.coord('time').bounds = None
    cube_filtered = cube.rolling_window(
        'time', SUM, len(weights), weights=weights
    )
    for coord in ['time', 'forecast_period']:
        try:
            cube_filtered.coord(coord).bounds = None
        except Exception:
            pass
    return cube_filtered
