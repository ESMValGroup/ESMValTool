'''
South Asian Monsoon Metrics

Split into several modules:
  SAM_area_metrics
  SAM_monsoon_indices
  SAM_jet_metrics
'''

import iris
import numpy as np

from ..utility.area_utils import area_average
from . import monsoon_general as mgl


def sam_area_metrics(temp_jjas, temp_djf, ppn_jjas,
                     u850_jjas, u200_jjas, u10m_son,
                     vmfl_jjas=None, land_jjas=None):
    'South Asian Monsoon metrics area averaged'

    metrics = dict()

    if land_jjas:
        land_jjas.data = np.asarray(land_jjas.data, dtype=np.float64)
        land_jjas_avg = land_jjas.collapsed('time', iris.analysis.MEAN)
        land_jjas_avg.data = np.asarray(land_jjas_avg.data, dtype=np.float32)
    ppn_jjas_avg = ppn_jjas.collapsed('time', iris.analysis.MEAN)
    ppn_jjas_std = ppn_jjas.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    temp_jjas_avg = temp_jjas.collapsed('time', iris.analysis.MEAN)
    temp_jjas_std = temp_jjas.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    temp_djf_avg = temp_djf.collapsed('time', iris.analysis.MEAN)
    temp_djf_std = temp_djf.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    u850_jjas_avg = u850_jjas.collapsed('time', iris.analysis.MEAN)
    u850_jjas_std = u850_jjas.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    u200_jjas_avg = u200_jjas.collapsed('time', iris.analysis.MEAN)
    u200_jjas_std = u200_jjas.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    u10m_son_avg = u10m_son.collapsed('time', iris.analysis.MEAN)
    u10m_son_std = u10m_son.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    if vmfl_jjas:
        vmfl_jjas_avg = vmfl_jjas.collapsed('time', iris.analysis.MEAN)
        vmfl_jjas_std = vmfl_jjas.collapsed('time',
                                            iris.analysis.STD_DEV,
                                            ddof=0)

    # Precipitation

    if land_jjas:
        reg = area_average(ppn_jjas_avg, mask=land_jjas_avg,
                           **mgl.regll(70, 85, 5, 30))
        metrics['Mean Precip: Indian Land (JJAS)'] = float(reg.data)
    else:
        metrics['Mean Precip: Indian Land (JJAS)'] = -10000.
    reg = area_average(ppn_jjas_avg, **mgl.regll(70, 85, 5, 30))
    metrics['Mean Precip: India (JJAS)'] = float(reg.data)
    reg = area_average(ppn_jjas_avg, **mgl.regll(70, 75, 7.5, 20))
    metrics['Mean Precip: SE ArSea+Ghats (JJAS)'] = float(reg.data)
    reg = area_average(ppn_jjas_avg, **mgl.regll(60, 80, -5, 7.5))
    metrics['Mean Precip: W Indian Oc. (JJAS)'] = float(reg.data)
    reg = area_average(ppn_jjas_avg, **mgl.regll(72, 80, 7.5, 15))
    metrics['Mean Precip: S India (JJAS)'] = float(reg.data)
    reg = area_average(ppn_jjas_avg, **mgl.regll(70, 85, 15, 22.5))
    metrics['Mean Precip: C India (JJAS)'] = float(reg.data)
    reg = area_average(ppn_jjas_avg, **mgl.regll(70, 90, 22.5, 30))
    metrics['Mean Precip: N India (JJAS)'] = float(reg.data)
    reg = area_average(ppn_jjas_avg, **mgl.regll(80, 95, 7.5, 20))
    metrics['Mean Precip: B.o.Bengal (JJAS)'] = float(reg.data)
    if land_jjas:
        reg = area_average(ppn_jjas_std, mask=land_jjas_avg,
                           **mgl.regll(70, 85, 5, 30))
        metrics['St dev Precip: Indian Land (JJAS)'] = float(reg.data)
    else:
        metrics['St dev Precip: Indian Land (JJAS)'] = -10000.
    reg = area_average(ppn_jjas_std, **mgl.regll(70, 85, 5, 30))
    metrics['St dev Precip: India (JJAS)'] = float(reg.data)
    reg = area_average(ppn_jjas_std, **mgl.regll(70, 75, 7.5, 20))
    metrics['St dev Precip: SE ArSea+Ghats (JJAS)'] = float(reg.data)
    reg = area_average(ppn_jjas_std, **mgl.regll(60, 80, -5, 7.5))
    metrics['St dev Precip: W Indian Oc. (JJAS)'] = float(reg.data)
    reg = area_average(ppn_jjas_std, **mgl.regll(72, 80, 7.5, 15))
    metrics['St dev Precip: S India (JJAS)'] = float(reg.data)
    reg = area_average(ppn_jjas_std, **mgl.regll(70, 85, 15, 22.5))
    metrics['St dev Precip: C India (JJAS)'] = float(reg.data)
    reg = area_average(ppn_jjas_std, **mgl.regll(70, 90, 22.5, 30))
    metrics['St dev Precip: N India (JJAS)'] = float(reg.data)
    reg = area_average(ppn_jjas_std, **mgl.regll(80, 95, 7.5, 20))
    metrics['St dev Precip: B.o.Bengal (JJAS)'] = float(reg.data)

    # Temperature

    if land_jjas:
        reg = area_average(temp_jjas_avg, mask=land_jjas_avg,
                           **mgl.regll(70, 85, 5, 30))
        metrics['Mean Indian land 1.5mT (JJAS)'] = float(reg.data)
        reg = area_average(temp_jjas_std, mask=land_jjas_avg,
                           **mgl.regll(70, 85, 5, 30))
        metrics['St dev Indian land 1.5mT (JJAS)'] = float(reg.data)
    else:
        metrics['Mean Indian land 1.5mT (JJAS)'] = -10000.
        metrics['St dev Indian land 1.5mT (JJAS)'] = -10000.
    reg = area_average(temp_djf_avg, **mgl.regll(40, 60, 20, 40))
    metrics['Mean Arab Penins-Iran Plat 1.5mT (DJF)'] = float(reg.data)
    reg = area_average(temp_djf_std, **mgl.regll(40, 60, 20, 40))
    metrics['St dev Arab Penins-Iran Plat 1.5mT (DJF)'] = float(reg.data)
    reg = area_average(temp_djf_avg, **mgl.regll(60, 100, 30, 40))
    metrics['Mean Tibet Plat 1.5mT (DJF)'] = float(reg.data)

    # 850hPa U Wind

    reg = area_average(u850_jjas_avg, **mgl.regll(70, 85, 5, 30))
    metrics['Mean 850hPa U: India (JJAS)'] = float(reg.data)
    reg = area_average(u850_jjas_std, **mgl.regll(70, 85, 5, 30))
    metrics['St dev 850hPa U: India (JJAS)'] = float(reg.data)

    # 200hPa U Wind

    reg = area_average(u200_jjas_avg, **mgl.regll(70, 85, 5, 30))
    metrics['Mean 200hPa U: India (JJAS)'] = float(reg.data)
    reg = area_average(u200_jjas_std, **mgl.regll(70, 85, 5, 30))
    metrics['St dev 200hPa U: India (JJAS)'] = float(reg.data)

    # 10m U Wind

    reg = area_average(u10m_son_avg, **mgl.regll(70, 90, -5, 5))
    metrics['Mean 10mZonal winds EQ Indian Ocean(SON)'] = float(reg.data)
    reg = area_average(u10m_son_std, **mgl.regll(70, 90, -5, 5))
    metrics['St dev 10mZonal winds EQ Indian Oc(SON)'] = float(reg.data)

    # Vertically integrated moisture flux
    # - at the moment there is no observational equivalent of this data 30428

    if vmfl_jjas:
        reg = area_average(vmfl_jjas_avg, **mgl.regll(50, 75, 7.5, 30))
        metrics['Mean Arabian Sea zonal VINT MFL (JJAS)'] = float(reg.data)
        reg = area_average(vmfl_jjas_std, **mgl.regll(50, 75, 7.5, 30))
        metrics['St dev Arabian Sea zonal VINT MFL (JJAS)'] = float(reg.data)
    else:
        metrics['Mean Arabian Sea zonal VINT MFL (JJAS)'] = -10000.
        metrics['St dev Arabian Sea zonal VINT MFL (JJAS)'] = -10000.

    return metrics


def sam_monsoon_indices(u850_jjas, u200_jjas, v850_jjas, v200_jjas):
    '''Output the SAM Seasonally and area averaged metrics'''

    metrics = dict()

    # Mean Dynamical Monsoon Index (Webster-Yang Index)
    #  = U850(40-110E,EQ-20N)- U200(40-110E,EQ-20N)
    lower = area_average(u850_jjas, **mgl.regll(40, 110, 0, 20))
    upper = area_average(u200_jjas, **mgl.regll(40, 110, 0, 20))
    index = lower - upper

    index_avg = index.collapsed('time', iris.analysis.MEAN)
    metrics['Mean Dynamical Monsoon Index'] = float(index_avg.data)
    lower_avg = lower.collapsed('time', iris.analysis.MEAN)
    metrics['Mean Dynamical Monsoon Index L'] = float(lower_avg.data)
    upper_avg = upper.collapsed('time', iris.analysis.MEAN)
    metrics['Mean Dynamical Monsoon Index U'] = float(upper_avg.data)
    index_std = index.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    metrics['St dev Dynamical Monsoon Index'] = float(index_std.data)
    lower_std = lower.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    metrics['St dev Dynamical Monsoon Index L'] = float(lower_std.data)
    upper_std = upper.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    metrics['St dev Dynamical Monsoon Index U'] = float(upper_std.data)

    # Mean Goswami merid.sh Index = meridional wind-shear anomaly (V850-V200)
    #  for JJAS averaged over [70-110E, 10-30N]
    lower = area_average(v850_jjas, **mgl.regll(70, 110, 10, 30))
    upper = area_average(v200_jjas, **mgl.regll(70, 110, 10, 30))
    index = lower - upper

    index_avg = index.collapsed('time', iris.analysis.MEAN)
    metrics['Mean Goswami merid.sh Index'] = float(index_avg.data)
    lower_avg = lower.collapsed('time', iris.analysis.MEAN)
    metrics['Mean Goswami merid.sh IndexL'] = float(lower_avg.data)
    upper_avg = upper.collapsed('time', iris.analysis.MEAN)
    metrics['Mean Goswami merid.sh IndexU'] = float(upper_avg.data)
    index_std = index.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    metrics['St dev Goswami merid.sh Index'] = float(index_std.data)
    lower_std = lower.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    metrics['St dev Goswami merid.sh IndexL'] = float(lower_std.data)
    upper_std = upper.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    metrics['St dev Goswami merid.sh IndexU'] = float(upper_std.data)

    # Mean WFan monsoon Index = U850(40E-80E,5N-15N)-U850(70E-90E,20N-30N)
    #  using a north box of [60,20,90,30] to be comparable with old code
    north = area_average(u850_jjas, **mgl.regll(60, 90, 20, 30))
    south = area_average(u850_jjas, **mgl.regll(40, 80, 5, 15))
    index = south - north

    index_avg = index.collapsed('time', iris.analysis.MEAN)
    metrics['Mean WFan monsoon Index'] = float(index_avg.data)
    north_avg = north.collapsed('time', iris.analysis.MEAN)
    metrics['Mean WFan monsoon IndexN'] = float(north_avg.data)
    south_avg = south.collapsed('time', iris.analysis.MEAN)
    metrics['Mean WFan monsoon IndexS'] = float(south_avg.data)
    index_std = index.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    metrics['St dev WFan monsoon Index'] = float(index_std.data)
    north_std = north.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    metrics['St dev WFan monsoon IndexN'] = float(north_std.data)
    south_std = south.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    metrics['St dev WFan monsoon IndexS'] = float(south_std.data)

    # Return the dictionary at the end
    return metrics


def sam_jet_metrics(u850_as):
    '''
    Output the SAM metrics which are specifically related to the size,
    position and shape of the jet.
    '''

    metrics = dict()

    data = u850_as.collapsed('time', iris.analysis.MEAN)
    reg = area_average(data, **mgl.regll(70, 85, 5, 30))
    metrics['Mean 850hPa U: jet mean strength (AS)'] = float(reg.data)

    reg = area_average(u850_as, **mgl.regll(70, 85, 5, 30))
    dmax = reg.collapsed('time', iris.analysis.MAX)
    metrics['850hPa monsoon jet strength - max'] = float(dmax.data)

    # Return the dictionary at the end
    return metrics


def sam_other_metrics(run):
    '''
    Output the SAM metrics which are neither area averages or monsoon metrics

    Note that these are not currently calculated
    '''

    metrics = dict()
    region = mgl.regll(40, 160, -20, 50)

    # RMSE metrics

    # Precip JJAS RMSE Region = [40,-20,160,50]
    metrics['Precip JJAS RMSE 40-160E 20S-50N'] = \
        mgl.calc_rms_error(run, 'precip', 'jjas', region)

    # UVWind JJAS RMSE Region = [40,-20,160,50] - Not sure what level????

    # Pattern Correlation metrics

    metrics['Precip JJAS Pattern Corr 40-160E 20S-50N'] = \
        mgl.calc_pattern_corr(run, 'precip', 'jjas', region)

    return metrics
