# pylint: disable=invalid-name
"""ESMValTool CMORizer for OSI-SAF data.

Tier
   Tier 2: other freely-available dataset.

Source
   http://osisaf.met.no/p/ice/

Last access
   20190502

Download and processing instructions
    Download the desired years from the following ftp:
        ftp://osisaf.met.no/reprocessed/ice/conc/v2p0
    Please, keep folder structure.

"""
from .osi_common import cmorize_osi


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    cmorize_osi(in_dir, out_dir, cfg, 'nh')
