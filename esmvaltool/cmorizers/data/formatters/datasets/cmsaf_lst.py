"""ESMValTool CMORizer for ESACCI-LST data.
Tier
   Tier 2: other freely-available dataset.
Source
   CMSAF website www.cmsaf.eu
Download and processing instructions
   Order .tar files from CMSAF website
   Untar to get netcdf files, one for each hour of every day
   To fix non-CF complient naming of semiminor and semimajor axes in the
   netcdf files run this shell script in the directory:
   for file in ls *.nc; 
   do
   ncrename -O -h -a grid_mapping@semiminoraxis,semi_minor_axis \
                  -a grid_mapping@semimajoraxis,semi_major_axis ${file};
   done

   This means that Iris will load the files sucessfully, otherwise it
   will raise ValueError: No ellipsoid specified
"""

import datetime
import logging
#from calendar import monthrange

import iris

from ...utilities import fix_coords, save_variable

logger = logging.getLogger(__name__)


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""


