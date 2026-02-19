"""ESMValTool CMORizer for NOAA GML surface flask CH4 data.

Tier
    Tier 2: freely available dataset.

Source
    https://gml.noaa.gov/

Last access
    20240730

Download and processing instructions
    Download the following file:
    https://gml.noaa.gov/aftp/data/trace_gases/ch4/flask/surface/ch4_surface-flask_ccgg_text.tar.gz
"""

from esmvaltool.cmorizers.data.formatters.datasets import (
    noaa_gml_surface_flask,
)


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    noaa_gml_surface_flask.cmorization(
        in_dir,
        out_dir,
        cfg,
        cfg_user,
        start_date,
        end_date,
    )
