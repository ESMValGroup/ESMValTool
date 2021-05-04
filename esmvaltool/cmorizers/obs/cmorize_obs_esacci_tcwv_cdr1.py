"""ESMValTool CMORizer for ESACCI-TCWV-CDR1 data.

Tier
   Tier 3: currently still restricted because preliminary.

Source
   Marc Schröder, ftp.brockmann-consult.de

Last access
   20212903

Download and processing instructions
   TBD

Modification history
   20210408-weigel_katja: written.

"""

import logging

from esmvaltool.cmorizers.obs.cmorize_obs_esacci_tcwv import (cmorization,
                                                              extract_variable)

logger = logging.getLogger(__name__)
logger.info("Load def from cmorize_obs_esacci_tcwv.py")
logger.info(extract_variable)
logger.info(cmorization)
