"""ESMValTool CMORizer for Earth System Data Cube data.

<We will add some useful info here later>
"""
import logging
from esmvaltool.cmorizers.data import utilities as utils
from pathlib import Path
import xarray as xr

logger = logging.getLogger(__name__)

def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize the dataset."""

    # This is where you'll add the cmorization code
    # 1. find the input data
    logger.info("in_dir: '%s'", in_dir)
    logger.info("cfg: '%s'", cfg)

    cfg['attributes']['cube'] = "esdc-8d-0.25deg-1x720x1440"

    version = cfg['attributes']['version']
    zarr_cube = f"{cfg['attributes']['cube']}-{version}.zarr"
    version = f"v{version}"

    ds = xr.open_zarr(Path(in_dir, version, zarr_cube), consolidated=True)

    logger.info(ds)

    # 2. apply the necessary fixes
    
    # 3. store the data with the correct filename

    