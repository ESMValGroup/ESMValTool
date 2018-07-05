"""
Execute

List the functions which are to be called by auto-assess
There are two types of functions:
 1. those which compute metrics (and optionally make plots) for a single run
 2. those which make plots comparing multiple runs


Each metrics function (type 1) has the following interface:
    Arguments:
        run - dictionary containing metadata for a single model run
              (see auto_assess.model_run for description of
              the contents of this dictionary)

    Returns:
        metrics - dictionary of metrics names and values
        optionally also writes image files to the current working dir


Each "multi" function (type 2) has the following interface:
    Arguments:
        runs - list of run dictionaries.  Each dictionary contains
               metadata for a single model run.  The first dictionary
               in this list is the control experiment.
               (see auto_assess.model_run for description of
               the contents of this dictionary)

    Returns:
        doesn't return any objects - it only writes image files to the
        current working dir

"""

# local modules
from . import atmospheric_mass as am
from . import energy_budget as eb
from . import GC_water_conservation as gcwc

metrics_functions = [am.global_atmos_mass_conservation,
                     eb.atmos_energy_budget,
                     gcwc.global_freshwater_fluxes_over_various_GC_cubmodels]
