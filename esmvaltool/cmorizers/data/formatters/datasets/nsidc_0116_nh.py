"""ESMValTool CMORizer for NSIDC-0116 data.

Tier
   Tier 3: restricted dataset.

Source
   https://nsidc.org/data/NSIDC-0116

Last access
   20190513

Download and processing instructions
    Download daily data from:
    https://nsidc.org/data/NSIDC-0116

    Login required for download, but requires citation only to use
"""
from esmvaltool.cmorizers.data.formatters.nsidc_common import cmorize


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    cmorize(cfg, 'nh', in_dir, out_dir)
