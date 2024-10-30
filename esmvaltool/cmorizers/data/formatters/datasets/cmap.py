"""ESMValTool CMORizer for CMAP (CPC Merged Analysis of Precipitation) data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://psl.noaa.gov/data/gridded/data.cmap.html

Last access
    20240909

Download and processing instructions
    To facilitate the download, the links to the ftp server are provided.

    https://downloads.psl.noaa.gov/Datasets/cmap/enh/
        precip.mon.mean.nc

Caveats

"""

import logging
import re
from copy import deepcopy
from pathlib import Path

import iris

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, cfg, raw_filepath, out_dir):
    cmor_info = cfg["cmor_table"].get_variable(var["mip"], short_name)
    attributes = deepcopy(cfg["attributes"])
    attributes["mip"] = var["mip"]

    cubes = iris.load(raw_filepath)
    for cube in cubes:
        assert cube.units == "mm/day", f"unknown units:{cube.units}"
        # convert data from mm/day to kg m-2 s-1
        # mm/day ~ density_water * mm/day
        # = 1000 kg m-3 * 1/(1000*86400) m s-1 = 1/86400 kg m-2 s-1
        cube = cube / 86400
        cube.units = "kg m-2 s-1"

        utils.fix_var_metadata(cube, cmor_info)
        cube = utils.fix_coords(cube)
        utils.set_global_atts(cube, attributes)

        logger.info("Saving file")
        utils.save_variable(cube, short_name, out_dir, attributes,
                            unlimited_dimensions=["time"])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    for short_name, var in cfg["variables"].items():
        logger.info("CMORizing variable '%s'", short_name)
        short_name = var["short_name"]
        raw_filenames = Path(in_dir).rglob("*.nc")
        filenames = []
        for raw_filename in raw_filenames:
            if re.search(var["file"], str(raw_filename)) is not None:
                filenames.append(raw_filename)

        for filename in sorted(filenames):
            _extract_variable(short_name, var, cfg, filename, out_dir)
