"""Lazy regridding, because this is not supported by iris (yet).

Iris issue requesting the feature:
https://github.com/SciTools/iris/issues/3700
"""
import copy

import iris
import numpy as np


def _compute_chunks(src, tgt):
    """Compute the chunk sizes needed to regrid src to tgt."""
    block_bytes = 10 * (1 << 20)  # 10 MB block size

    if src.dtype == np.float32:
        dtype_bytes = 4  # size of float32 in bytes
    else:
        dtype_bytes = 8  # size of float64 in bytes

    ntime = src.coord('time').shape[0]
    tgt_nlat = tgt.coord('latitude').shape[0]
    tgt_nlon = tgt.coord('longitude').shape[0]

    # Define blocks along the time dimension
    min_nblocks = int(ntime * tgt_nlat * tgt_nlon * dtype_bytes / block_bytes)
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
    tgt_chunks = (
        time_chunks,
        (tgt_nlat, ),
        (tgt_nlon, ),
    )

    return src_chunks, tgt_chunks


def _regrid_data(src, tgt):
    """Regrid data from cube src onto grid of cube tgt."""
    src_chunks, tgt_chunks = _compute_chunks(src, tgt)

    # Define the block regrid function
    regridder = iris.analysis.Linear().regridder(src, tgt)

    def regrid(block):
        tlen = block.shape[0]
        cube = src[:tlen].copy(block)
        return regridder(cube).core_data()

    # Regrid
    data = src.core_data().rechunk(src_chunks).map_blocks(regrid,
                                                          dtype=src.dtype,
                                                          chunks=tgt_chunks)

    return data


def lazy_linear_regrid(src, tgt):
    """Regrid cube src onto the grid of cube tgt."""
    data = _regrid_data(src, tgt)

    result = iris.cube.Cube(data)
    result.metadata = copy.deepcopy(src.metadata)

    def copy_coords(src_coords, add_method):
        for coord in src_coords:
            dims = src.coord_dims(coord)
            if coord == src.coord('longitude'):
                coord = tgt.coord('longitude')
            elif coord == src.coord('latitude'):
                coord = tgt.coord('latitude')
            result_coord = coord.copy()
            add_method(result_coord, dims)

    copy_coords(src.dim_coords, result.add_dim_coord)
    copy_coords(src.aux_coords, result.add_aux_coord)

    return result
