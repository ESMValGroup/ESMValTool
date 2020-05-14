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

    If you want to use only the nh data download only the nh files,
    using, e.g., wget -r -A '*_nh*.nc'.

    If you also want to cmorize the sh, download everything and create a link
    for OSI-450-sh pointing to the data folder. Both cmorizers will ignore
    files belonging to the other hemisphere

"""
from esmvaltool.cmorizers.obs.formatters.osi_common import OSICmorizer


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    cmorizer = OSICmorizer(in_dir, out_dir, cfg, 'nh')
    cmorizer.cmorize()
