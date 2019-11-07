"""ESMValTool CMORizer for ERA-Interim-Land data.

Tier
    Tier 3: restricted datasets (i.e., dataset which requires a registration
 to be retrieved or provided upon request to the respective contact or PI).

Source
    https://apps.ecmwf.int/datasets/data/interim-land/type=fc/

Last access
    20191104

Download and processing instructions
     See script cmorize_obs_era_interim.py

"""

from .cmorize_obs_era_interim import cmorization

__all__ = ['cmorization']
