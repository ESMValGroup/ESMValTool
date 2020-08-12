"""
Zonal-mean annular mode preproc routine.

Author: Federico Serva (ISAC-CNR & ISMAR-CNR, Italy)
Copernicus C3S 34a lot 2 (MAGIC)
"""

import cdo as cd


def zmnam_preproc(ifile):
    """Preprocessing of the input dataset files."""
    cdo = cd.Cdo()
    # Delete leap day, if any.
    full_da_nl = cdo.delete('month=2,day=29', input=ifile)

    # Fill missing values with bilinear method.
    full_da_nl_nn = cdo.fillmiss(input=full_da_nl)

    # Detrend with memory-efficient method.
    full_da_nl_a, full_da_nl_b = cdo.trend(input=full_da_nl_nn)
    full_da_nl_detr = cdo.subtrend(input=full_da_nl_nn + ' ' + full_da_nl_a +
                                   ' ' + full_da_nl_b)

    # Compute anomalies from the daily/monthly means.
    gh_da_dm = cdo.ydaymean(input=full_da_nl_detr)
    gh_da_an = cdo.sub(input=full_da_nl_detr + ' ' + gh_da_dm)
    gh_da_an_zm = cdo.zonmean(input=gh_da_an)

    # Compute monthly mean anomalies.
    gh_mo = cdo.monmean(input=full_da_nl_detr)
    gh_mo_mm = cdo.ymonmean(input=gh_mo)
    gh_mo_an = cdo.sub(input=gh_mo + ' ' + gh_mo_mm)

    return (gh_da_an_zm, gh_mo_an)
