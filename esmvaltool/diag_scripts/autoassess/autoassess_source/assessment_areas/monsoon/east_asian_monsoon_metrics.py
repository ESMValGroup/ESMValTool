'''
East Asian Monsoon Metrics

Split into several modules:
  EAM_area_metrics
  EAM_monsoon_indices
  EAM_other_metrics
'''

import numpy as np

import iris

from ..utility.area_utils import area_average

from . import monsoon_general as mgl


def eam_area_metrics(temp_jja, temp_djf, ppn_jja, ppn_djf,
                     mslp_jja, mslp_djf, z500_jja, z500_djf,
                     v850_jja, v850_djf, u200_jja, u200_djf,
                     land_jja=None, land_djf=None):
    'East Asian Monsoon metrics area averaged'

    metrics = dict()

    if land_jja:
        land_jja.data = np.asarray(land_jja.data, dtype=np.float64)
        land_jja_avg = land_jja.collapsed('time', iris.analysis.MEAN)
        land_jja_avg.data = np.asarray(land_jja_avg.data, dtype=np.float32)
    if land_djf:
        land_djf.data = np.asarray(land_djf.data, dtype=np.float64)
        land_djf_avg = land_djf.collapsed('time', iris.analysis.MEAN)
        land_djf_avg.data = np.asarray(land_djf_avg.data, dtype=np.float32)
    temp_jja_avg = temp_jja.collapsed('time', iris.analysis.MEAN)
    temp_djf_avg = temp_djf.collapsed('time', iris.analysis.MEAN)
    ppn_jja_avg = ppn_jja.collapsed('time', iris.analysis.MEAN)
    ppn_jja_std = ppn_jja.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    ppn_djf_avg = ppn_djf.collapsed('time', iris.analysis.MEAN)
    ppn_djf_std = ppn_djf.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    mslp_jja_avg = mslp_jja.collapsed('time', iris.analysis.MEAN)
    mslp_jja_std = mslp_jja.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    mslp_djf_std = mslp_djf.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    z500_jja_avg = z500_jja.collapsed('time', iris.analysis.MEAN)
    z500_jja_std = z500_jja.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    z500_djf_avg = z500_djf.collapsed('time', iris.analysis.MEAN)
    z500_djf_std = z500_djf.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    v850_jja_avg = v850_jja.collapsed('time', iris.analysis.MEAN)
    v850_jja_std = v850_jja.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    v850_djf_avg = v850_djf.collapsed('time', iris.analysis.MEAN)
    v850_djf_std = v850_djf.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    u200_jja_avg = u200_jja.collapsed('time', iris.analysis.MEAN)
    u200_jja_std = u200_jja.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    u200_djf_avg = u200_djf.collapsed('time', iris.analysis.MEAN)
    u200_djf_std = u200_djf.collapsed('time', iris.analysis.STD_DEV, ddof=0)

    # Temperature

    reg = area_average(temp_jja_avg, **mgl.regll(100, 145, 20, 50))
    metrics['Temp JJA: 100-145E 20-50N'] = float(reg.data)
    reg = area_average(temp_jja_avg, **mgl.regll(100, 145, 20, 35))
    metrics['Temp JJA: 100-145E 20-35N'] = float(reg.data)
    reg = area_average(temp_jja_avg, **mgl.regll(100, 145, 35, 50))
    metrics['Temp JJA: 100-145E 35-50N'] = float(reg.data)

    reg = area_average(temp_djf_avg, **mgl.regll(100, 145, 20, 50))
    metrics['Temp DJF: 100-145E 20-50N'] = float(reg.data)
    reg = area_average(temp_djf_avg, **mgl.regll(100, 145, 20, 35))
    metrics['Temp DJF: 100-145E 20-35N'] = float(reg.data)
    reg = area_average(temp_djf_avg, **mgl.regll(100, 145, 35, 50))
    metrics['Temp DJF: 100-145E 35-50N'] = float(reg.data)

    # Precipitation

    if land_jja:
        reg = area_average(ppn_jja_avg, mask=land_jja_avg,
                               **mgl.regll(100, 145, 20, 50))
        metrics['Prcp(land) JJA: 100-145E 20-50N'] = float(reg.data)
    else:
        metrics['Prcp(land) JJA: 100-145E 20-50N'] = -10000.
    reg = area_average(ppn_jja_avg, **mgl.regll(100, 145, 20, 50))
    metrics['Precip JJA: 100-145E 20-50N'] = float(reg.data)
    reg = area_average(ppn_jja_avg, **mgl.regll(100, 120, 20, 40))
    metrics['Precip JJA: 100-120E 20-40N'] = float(reg.data)
    reg = area_average(ppn_jja_avg, **mgl.regll(120, 145, 20, 25))
    metrics['Precip JJA: 120-145E 20-25N'] = float(reg.data)
    reg = area_average(ppn_jja_avg, **mgl.regll(120, 145, 30, 40))
    metrics['Precip JJA: 120-145E 30-40N'] = float(reg.data)
    reg = area_average(ppn_jja_std, **mgl.regll(100, 145, 20, 50))
    metrics['Prcp-std JJA: 100-145E 20-50N'] = float(reg.data)

    if land_djf:
        reg = area_average(ppn_djf_avg, mask=land_djf_avg,
                               **mgl.regll(100, 145, 20, 50))
        metrics['Prcp(land) DJF: 100-145E 20-50N'] = float(reg.data)
    else:
        metrics['Prcp(land) DJF: 100-145E 20-50N'] = -10000.
    reg = area_average(ppn_djf_avg, **mgl.regll(100, 145, 20, 50))
    metrics['Precip DJF: 100-145E 20-50N'] = float(reg.data)
    reg = area_average(ppn_djf_avg, **mgl.regll(100, 120, 20, 40))
    metrics['Precip DJF: 100-120E 20-40N'] = float(reg.data)
    reg = area_average(ppn_djf_avg, **mgl.regll(100, 120, 30, 40))
    metrics['Precip DJF: 100-120E 30-40N'] = float(reg.data)
    reg = area_average(ppn_djf_avg, **mgl.regll(120, 145, 20, 25))
    metrics['Precip DJF: 120-145E 20-25N'] = float(reg.data)
    reg = area_average(ppn_djf_avg, **mgl.regll(120, 145, 30, 40))
    metrics['Precip DJF: 120-145E 30-40N'] = float(reg.data)
    reg = area_average(ppn_djf_std, **mgl.regll(100, 145, 20, 50))
    metrics['Prcp-std DJF: 100-145E 20-50N'] = float(reg.data)

    # Mean Sea Level Pressure

    reg = area_average(mslp_jja_avg, **mgl.regll(120, 150, 15, 30))
    metrics['MSLP JJA: 120-150E 15-30N'] = float(reg.data)
    reg = area_average(mslp_jja_std, **mgl.regll(120, 150, 15, 30))
    metrics['MSLP-std JJA: 120-150E 15-30N'] = float(reg.data)

    reg = area_average(mslp_djf_std, **mgl.regll(120, 150, 15, 30))
    metrics['MSLP-std DJF: 120-150E 15-30N'] = float(reg.data)

    # 500 hPa Geopotential height

    reg = area_average(z500_jja_avg, **mgl.regll(120, 140, 35, 55))
    metrics['Z500 JJA: 120-140E 35-55N'] = float(reg.data)
    reg = area_average(z500_jja_std, **mgl.regll(120, 140, 35, 55))
    metrics['Z500-std JJA: 120-140E 35-55N'] = float(reg.data)

    reg = area_average(z500_djf_avg, **mgl.regll(120, 140, 35, 55))
    metrics['Z500 DJF: 120-140E 35-55N'] = float(reg.data)
    reg = area_average(z500_djf_std, **mgl.regll(120, 140, 35, 55))
    metrics['Z500-std DJF: 120-140E 35-55N'] = float(reg.data)

    # 850hPa V Wind

    reg = area_average(v850_jja_avg, **mgl.regll(110, 130, 15, 30))
    metrics['V850 JJA: 110-130E 15-30N'] = float(reg.data)
    reg = area_average(v850_jja_avg, **mgl.regll(110, 120, 20, 35))
    metrics['V850 JJA: 110-120E 20-35N'] = float(reg.data)
    reg = area_average(v850_jja_avg, **mgl.regll(120, 130, 15, 30))
    metrics['V850 JJA: 120-130E 15-30N'] = float(reg.data)
    reg = area_average(v850_jja_std, **mgl.regll(120, 140, 30, 50))
    metrics['V850-std JJA: 120-140E 30-50N'] = float(reg.data)

    reg = area_average(v850_djf_avg, **mgl.regll(120, 140, 30, 50))
    metrics['V850 DJF: 120-140E 30-50N'] = float(reg.data)
    reg = area_average(v850_djf_std, **mgl.regll(120, 140, 15, 35))
    metrics['V850-std DJF: 120-140E 15-35N'] = float(reg.data)

    # 200hPa U Wind

    reg = area_average(u200_jja_avg, **mgl.regll(100, 145, 20, 50))
    metrics['U200 JJA: 100-145E 20-50N'] = float(reg.data)
    reg = area_average(u200_jja_std, **mgl.regll(90, 150, 25, 55))
    metrics['U200-std JJA 90-150E 25-55N'] = float(reg.data)

    reg = area_average(u200_djf_avg, **mgl.regll(120, 145, 20, 50))
    metrics['U200 DJF: 100-145E 20-50N'] = float(reg.data)
    reg = area_average(u200_djf_std, **mgl.regll(90, 150, 25, 55))
    metrics['U200-std DJF 90-150E 25-55N'] = float(reg.data)

    return metrics


def eam_monsoon_indices(v850_jja, v850_djf, mslp_jja, mslp_djf,
                        u850_jja, u300_djf):
    '''Outputs the EAM Monsoon Index values for several different indices'''

    metrics = dict()

    # Wang JJA 110-140E 20-40N 850 hPa V winds (ea_sub2.pro)
    region = area_average(v850_jja, **mgl.regll(110, 140, 20, 40))
    index_avg = region.collapsed('time', iris.analysis.MEAN)
    index_std = region.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    name = 'Wang JJA 110-140E 20-40N 850 hPa V winds'
    metrics[name] = float(index_avg.data)
    name = 'Wang Stddev JJA 110-140E 20-40N 850 hPa V winds'
    metrics[name] = float(index_std.data)

    # Wang DJF 110-140E 20-40N 850 hPa V winds (ea_sub2.pro) ???
    region = area_average(v850_djf, **mgl.regll(110, 140, 20, 40))
    index_avg = region.collapsed('time', iris.analysis.MEAN)
    index_std = region.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    name = 'Wang DJF 110-140E 20-40N 850 hPa V winds'
    metrics[name] = float(index_avg.data)
    name = 'Wang Stddev DJF 110-140E 20-40N 850 hPa V winds'
    metrics[name] = float(index_std.data)

    # Guo-std JJA MSLP south=[110, 10, 120, 50] & north=[150,10,160,50]
    # ea_sub2.pro
    west = area_average(mslp_jja, **mgl.regll(110, 120, 10, 50))
    east = area_average(mslp_jja, **mgl.regll(150, 160, 10, 50))
    index = west - east
    index_std = index.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    metrics['Guo-std JJA'] = float(index_std.data)

    # Guo-std DJF
    west = area_average(mslp_djf, **mgl.regll(110, 120, 10, 50))
    east = area_average(mslp_djf, **mgl.regll(150, 160, 10, 50))
    index = west - east
    index_std = index.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    metrics['Guo-std DJF'] = float(index_std.data)

    # WNPMI-std JJA added from Kung-On's updates
    # U850 (5-15N,100-130E - 23-33N,110-140E) (20-30N is also seen)
    # This index is the inverse of the EASM index below
    south = area_average(u850_jja, **mgl.regll(100, 130, 5, 15))
    north = area_average(u850_jja, **mgl.regll(110, 140, 23, 33))
    index = south - north
    index_avg = index.collapsed('time', iris.analysis.MEAN)
    name = 'WNPMI-std JJA U850 (5-15N 100-130E - 23-33N 110-140E)'
    metrics[name] = float(index_avg.data)

    # EASM index (taken from IAP via Malcom Robert's copy of slides)
    # Added by Amanda
    # U850 (22.5-32.5N,110-140E - 5-15N,90-130E)
    north = area_average(u850_jja, **mgl.regll(110, 140, 22.5, 32.5))
    south = area_average(u850_jja, **mgl.regll(90, 130, 5, 15))
    index = north - south
    index_avg = index.collapsed('time', iris.analysis.MEAN)
    name = 'EASM Index JJA U850 (22.5-32.5N 110-140E - 5-15N 90-130E)'
    metrics[name] = float(index_avg.data)

    # Jhun and Lee (2004) EAWMI - New index added by Amanda.
    south = area_average(u300_djf, **mgl.regll(110, 170, 27.5, 37.5))
    north = area_average(u300_djf, **mgl.regll(80, 140, 50, 60))
    index = south - north
    index_avg = index.collapsed('time', iris.analysis.MEAN)
    index_std = index.collapsed('time', iris.analysis.STD_DEV, ddof=0)
    metrics['Jhun and Lee EAWMI mean'] = float(index_avg.data)
    metrics['Jhun and Lee EAWMI std'] = float(index_std.data)

    return metrics


def eam_other_metrics(mslp_jja, mslp_djf, u850_jja, u850_djf,
                      u200_jja, u200_djf, v850_jja, v850_djf):
    '''
    Output the EAM metrics which are neither area averages or monsoon metrics
    '''

    metrics = dict()

    # Gradient/Sheer Calculations
    # MSLP EW gradJJA (new in ea_sub1a / ea_sub1 in Kyung-On's updates)
    # [120,20,130,35] - [105,20,115,35]
    west = area_average(mslp_jja, **mgl.regll(105, 115, 20, 35))
    east = area_average(mslp_jja, **mgl.regll(120, 130, 20, 35))
    index = west - east
    index_avg = index.collapsed('time', iris.analysis.MEAN)
    metrics['MSLP EW grad JJA'] = float(index_avg.data)

    # MSLP EW gradDJF 16222 (from ea_sub1a / ea_sub1)
    # [100,35,130,60] - [140,35,170,60]
    west = area_average(mslp_djf, **mgl.regll(100, 130, 35, 60))
    east = area_average(mslp_djf, **mgl.regll(140, 170, 35, 60))
    index = west - east
    index_avg = index.collapsed('time', iris.analysis.MEAN)
    metrics['MSLP EW grad DJF'] = float(index_avg.data)

    # dudz JJA: 100-145E 20-50N (u200-u850 (30201)) cf ea_sub1b.pro and ea_sub1
    lo = area_average(u850_jja, **mgl.regll(100, 145, 20, 50))
    hi = area_average(u200_jja, **mgl.regll(100, 145, 20, 50))
    index = hi - lo
    index_avg = index.collapsed('time', iris.analysis.MEAN)
    metrics['dudz JJA: 100-145E 20-50N'] = float(index_avg.data)

    # dudz DJF: 100-145E 20-50N (u200-u850 (30201)) cf ea_sub1b.pro and ea_sub1
    lo = area_average(u850_djf, **mgl.regll(100, 145, 20, 50))
    hi = area_average(u200_djf, **mgl.regll(100, 145, 20, 50))
    index = hi - lo
    index_avg = index.collapsed('time', iris.analysis.MEAN)
    metrics['dudz DJF: 100-145E 20-50N'] = float(index_avg.data)

    # Magnitude calculations

    # UV850 JJA: 110-130E 15-30N
    index = mgl.mag(u850_jja, v850_jja, **mgl.regll(110, 130, 15, 30))
    index_avg = index.collapsed('time', iris.analysis.MEAN)
    metrics['UV850 JJA: 110-130E 15-30N'] = float(index_avg.data)

    # UV850 DJF: 120-140E 20-50N
    index = mgl.mag(u850_djf, v850_djf, **mgl.regll(120, 140, 20, 50))
    index_avg = index.collapsed('time', iris.analysis.MEAN)
    metrics['UV850 DJF: 120-140E 20-50N'] = float(index_avg.data)

    return metrics


def eam_rmse_metrics(run):
    '''
    Output the EAM metrics which are RMS Errors

    Note that these are not currently calculated
    '''

    metrics = dict()
    region = mgl.regll(100, 145, 20, 50)

    # RMSE calculations
    # - variation compared to observational climatology

    # RMS(Temp:land)JJA 100-145E 20-50N
    metrics['RMS(Temp:land)JJA 100-145E 20-50N'] = \
        mgl.calc_rms_error(run, 'temp', 'jja', region, mask=True)

    # RMS(Temp:land)DJF 100-145E 20-50N
    metrics['RMS(Temp:land)DJF 100-145E 20-50N'] = \
        mgl.calc_rms_error(run, 'temp', 'djf', region, mask=True)

    # RMS(MSLP) JJA 100-145E 20-50N
    metrics['RMS(MSLP)JJA 100-145E 20-50N'] = \
        mgl.calc_rms_error(run, 'mslp', 'jja', region)

    # RMS(MSLP) DJF 100-145E 20-50N
    metrics['RMS(MSLP)DJF 100-145E 20-50N'] = \
        mgl.calc_rms_error(run, 'mslp', 'djf', region)

    # RMS(Prcp:land) JJA 100-145E 20-50N
    metrics['RMS(Prcp:land)JJA 100-145E 20-50N'] = \
        mgl.calc_rms_error(run, 'precip', 'jja', region, mask=True)

    # RMS(Prcp:land) DJF 100-145E 20-50N
    metrics['RMS(Prcp:land)DJF 100-145E 20-50N'] = \
        mgl.calc_rms_error(run, 'precip', 'djf', region, mask=True)

    # RMS(Prcp) JJA 100-145E 20-50N
    metrics['RMS(Prcp)JJA 100-145E 20-50N'] = \
        mgl.calc_rms_error(run, 'precip', 'jja', region)

    # RMS(Prcp) DJF 100-145E 20-50N
    metrics['RMS(Prcp)DJF 100-145E 20-50N'] = \
        mgl.calc_rms_error(run, 'precip', 'djf', region)

    # Pattern Correlation metrics
    # - These are currently outputting a two by two matrix.

    # Prcp JJA pattern Corr [100,20,145,50]
    metrics['Prcp JJA pattern Corr 100-145E 20-50N'] = \
        mgl.calc_pattern_corr(run, 'precip', 'jja', region)

    # Prcp DJF pattern Corr [100,20,145,50]
    metrics['Prcp DJF pattern Corr 100-145E 20-50N'] = \
        mgl.calc_pattern_corr(run, 'precip', 'djf', region)

    return metrics
