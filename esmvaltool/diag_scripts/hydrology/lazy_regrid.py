"""Lazy regridding, because this is not supported by iris (yet).

Iris issue requesting the feature:
https://github.com/SciTools/iris/issues/3808
"""
import numpy as np
import iris

HORIZONTAL_SCHEMES = {
    'linear': iris.analysis.Linear(extrapolation_mode='mask'),
    'linear_extrapolate':
    iris.analysis.Linear(extrapolation_mode='extrapolate'),
    'nearest': iris.analysis.Nearest(extrapolation_mode='mask'),
    'area_weighted': iris.analysis.AreaWeighted(),
}
"""Supported horizontal regridding schemes."""


def _compute_chunks(src, tgt):
    """Compute the chunk sizes needed to regrid src to tgt."""
    block_bytes = 50 * (1 << 20)  # 50 MB block size

    if src.dtype == np.float32:
        dtype_bytes = 4  # size of float32 in bytes
    else:
        dtype_bytes = 8  # size of float64 in bytes

    ntime = src.coord('time').shape[0]
    tgt_nlat = tgt.coord('latitude').shape[0]
    tgt_nlon = tgt.coord('longitude').shape[0]

    # Define blocks along the time dimension
    min_nblocks = int(ntime * tgt_nlat * tgt_nlon * dtype_bytes / block_bytes)
    min_nblocks = max(min_nblocks, 1)
    timefull = ntime // min_nblocks
    timepart = ntime % timefull

    nfullblocks = ntime // timefull
    npartblocks = int(timepart > 0)

    time_chunks = (timefull, ) * nfullblocks + (timepart, ) * npartblocks
    src_chunks = (
        time_chunks,
        (src.coord('latitude').shape[0], ),
        (src.coord('longitude').shape[0], ),
    )

    return src_chunks


def lazy_regrid(src, tgt, scheme):
    """Regrid cube src onto the grid of cube tgt."""
    src_chunks = _compute_chunks(src, tgt)

    if scheme not in HORIZONTAL_SCHEMES:
        raise ValueError(f"Regridding scheme {scheme} not supported, "
                         f"choose from {HORIZONTAL_SCHEMES.keys()}.")
    regridder = HORIZONTAL_SCHEMES[scheme].regridder(src, tgt)
    src.data = src.lazy_data().rechunk(src_chunks)
    return regridder(src)
