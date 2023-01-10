"""ESMValTool CMORizer for GLOBMAP LAI data.
Tier
   Tier 2: other freely-available dataset.
Source
   https://zenodo.org/record/4700264
Download and processing instructions
   Original data files are hdf4 files
"""

import datetime
import logging
#from calendar import monthrange
# hdf readers 

import iris

from ...utilities import fix_coords, save_variable

logger = logging.getLogger(__name__)


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
