"""

Zonal-mean annular mode preproc routine.

Author: Federico Serva (ISAC-CNR & ISMAR-CNR, Italy)
Copernicus C3S 34a lot 2 (MAGIC)

"""

import subprocess


def zmnam_preproc(ifile):
    """Preprocessing of the input dataset files."""
    # Delete leap day, if any.
    subprocess.check_call('cdo -delete,month=2,day=29 ' + ifile +
                          ' tmp_full_da_nl.nc', shell=True)

    # Compute anomalies from the daily/monthly means.
    subprocess.check_call('cdo ydaymean tmp_full_da_nl.nc tmp_gh_da_dm.nc',
                          shell=True)
    subprocess.check_call('cdo sub tmp_full_da_nl.nc tmp_gh_da_dm.nc '
                          'tmp_gh_da_an.nc', shell=True)
    subprocess.check_call('cdo zonmean tmp_gh_da_an.nc tmp_gh_da_an_zm_hem.nc',
                          shell=True)

    subprocess.check_call('cdo monmean tmp_full_da_nl.nc tmp_gh_mo.nc',
                          shell=True)
    subprocess.check_call('cdo ymonmean tmp_gh_mo.nc tmp_gh_mo_mm.nc',
                          shell=True)
    subprocess.check_call('cdo sub tmp_gh_mo.nc tmp_gh_mo_mm.nc '
                          'tmp_gh_mo_an_hem.nc', shell=True)

    # Cleanup unnecessary files. Retain hemispheric only.
    subprocess.check_call('rm tmp_full_da_nl.nc tmp_gh_da_dm.nc '
                          'tmp_gh_da_an.nc', shell=True)
    subprocess.check_call('rm tmp_gh_mo.nc tmp_gh_mo_mm.nc', shell=True)
