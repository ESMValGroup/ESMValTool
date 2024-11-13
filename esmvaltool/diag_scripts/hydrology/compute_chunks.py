"""Re-chunk the time dimension, to be used by the regrid processor.

For large cubes, regridding to a high resolution grid increases the size
of the data. To reduce memory use, we re-chunk the time dimension.

Related iris issue:
https://github.com/SciTools/iris/issues/3808
"""
import numpy as np


def compute_chunks(src, tgt):
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
