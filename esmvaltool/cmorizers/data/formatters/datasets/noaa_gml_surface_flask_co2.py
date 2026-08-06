"""ESMValTool CMORizer for NOAA GML surface flask CO2 data.

Tier
    Tier 2: freely available dataset.

Source
    https://gml.noaa.gov/

Last access
    20240730

Download and processing instructions
    Download the following file:
    https://gml.noaa.gov/aftp/data/trace_gases/co2/flask/surface/co2_surface-flask_ccgg_text.tar.gz

NOTE this formatter is currently not working for Python 3.14 since pys2index is
not supported under Python 3.14! You will have to use Python<3.14, and install pys2index
from conda-forge.
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
