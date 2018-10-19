"""
Fixes for EC-Earth3-LR PRIMAVERA project data.

Until now, all fixes affect  HR and LR versions of EC-Earth, so all classes in
this module are derived from the HR versions
"""
from . import EC_Earth3_HR


class allvars(EC_Earth3_HR.allvars):
    """Fixes common to all variables."""
    pass


class siconc(EC_Earth3_HR.siconc):
    """Fixes for siconc."""
    pass


class zg(EC_Earth3_HR.zg):
    """Fixes for zg."""
    pass


class tas(EC_Earth3_HR.tas):
    """Fixes for tas."""
    pass
