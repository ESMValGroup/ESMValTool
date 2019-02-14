"""
Zonal-mean annular mode cleaning routine.

Author: Federico Serva (ISAC-CNR & ISMAR-CNR, Italy)
Copernicus C3S 34a lot 2 (MAGIC)
"""

import subprocess


def zmnam_clean():
    """Dispose of temporary files."""
    subprocess.check_call('rm tmp_gh_da_an_zm_hem.nc', shell=True)
    subprocess.check_call('rm tmp_gh_mo_an_hem.nc', shell=True)
