"""

Zonal-mean annular mode preproc routine.

Author: Federico Serva (ISAC-CNR & ISMAR-CNR, Italy)
Copernicus C3S 34a lot 2 (MAGIC)

"""

import os


def zmnam_preproc(ifile):
    """

    cdo preprocessing of the input dataset files

    """

    os.system('cdo -delete,month=2,day=29 ' + ifile + ' tmp_full_da_nl.nc')

    # Compute anomalies from the daily/monthly means. Regrid.
    os.system('cdo ydaymean tmp_full_da_nl.nc tmp_gh_da_nl_dm.nc')
    os.system('cdo sub tmp_full_da_nl.nc tmp_gh_da_nl_dm.nc tmp_gh_da_an.nc')
    os.system('cdo zonmean tmp_gh_da_an.nc tmp_gh_da_an_zm_hem.nc')

    os.system('cdo monmean tmp_full_da_nl.nc tmp_gh_mo.nc')
    os.system('cdo ymonmean tmp_gh_mo.nc tmp_gh_mo_mm.nc')
    os.system('cdo sub tmp_gh_mo.nc tmp_gh_mo_mm.nc tmp_gh_mo_an_hem.nc')

    # Cleanup unnecessary files. Retain hemispheric only.
    os.system('rm tmp_full_da_nl.nc tmp_gh_da_nl_dm.nc ' +
              'tmp_gh_da_an.nc tmp_gh_da_an_zm.nc')

    os.system('rm tmp_gh_mo.nc tmp_gh_mo_mm.nc tmp_gh_mo_an.nc ')
