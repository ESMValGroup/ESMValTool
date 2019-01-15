"""Module for LEC computation in climate models.

This module contains all the instructions to compute the atmospheric
Lorenz Energy Cycle in spectral coordinates.

@author: Valerio Lembo, University of Hamburg
"""

from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
from fluxogram import Fluxogram
import sys
import math
import os
import warnings
import numpy as np
warnings.filterwarnings("ignore")

G = 9.81
R = 287.00
CP = 1003.5
AA = 6.371E6
PS = 101100.0
NW_1 = 3
NW_2 = 9
NW_3 = 21


class LorenzCycle():
    """PROGRAM FOR LEC COMPUTATION.

    The class consists of the following functions:
        - lorenz: it is the main program, controlling the file input,
                  separating the real from imaginary part of the Fourier
                  coefficients, reordering the latitudinal dimension
                  (from N to S), interpolating on a reference sigma coordinate,
                  then computing the reservoirs and conversion terms, storing
                  them separately in NetCDF files and providing a flux diagram
                  and a table outputs, the latter separately for the two
                  hemispheres;
        - bsslzr: it contains the coefficients for the conversion from regular
                  lonlat grid to Gaussian grid;
        - diagram: it is the interface between the main program and a
                   class "Fluxogram", producing the flux diagram;
        - gauaw: it uses the coefficients provided in bsslzr for the lonlat to
                 Gaussian grid conversion;
        - globall_cg: it computes the global and hemispheric means at each
                      timestep;
        - makek: computes the KE reservoirs;
        - makea: computes the APE reservoirs;
        - mka2k: computes the APE->KE conversion terms;
        - mkaeaz: computes the zonal APE - eddy APE conversion terms;
        - mkkekz: computes the zonal KE - eddy KE conversion terms;
        - mkatas: computes the stationay eddy - transient eddy APE conversions;
        - mkktks: computes the stationay eddy - transient eddy KE conversions;
        - ncdump: provides the attributes of a given variable in a Nc dataset;
        - pr_output: prints a single component of the LEC computations to a
                     single Nc file;
        - removeif: removes a file if it exists;
        - stabil: calculates the stability parameter;
        - table: prints the global and hemispheric mean values of
                 the reservoirs;
        - table_conv: prints the global and hemispheric mean values of the
                      conversion terms;
        - varatts: prints the attributes of a variable in a Nc file;

    Constants:
        G: gravitational acceleration;
        R: gas constant;
        CP: specific heat of dry air;
        AA: planet radius;
        PS: mean surface pressure;
        NW_1: wavenumber limit for long waves;
        NW_2: wavenumber limit for synoptic waves;
        NW_3: wavenumber limit for short waves;
        
    References:
        Ulbrich P. and P. Speth (1991) The global energy cycle of stationary
        and transient atmospheric waves: Results from ECMWF analyses, Met.
 
    Authors:
        Frank Lunkeit, Meteorology Department, University of Hamburg
        Valerio Lembo, Meteorology Department, University of Hamburg
        
        Contact author: valerio.lembo@uni-hamburg.de.
    """
        
    from lorenz_cycle import LorenzCycle
        
    def lorenz(self, outpath, model, year, filenc, plotfile, logfile):
        """Main script, managing input and output fields and calling functions.

        Receive fields t,u,v,w as input fields in Fourier
        coefficients (time,level,wave,lon) and compute the LEC.
        
        Arguments:
            - outpath: ath where otput fields are stored (as NetCDF fields);
            - model: name of the model that is analysed;
            - year: year that is considered;
            - filenc: name of the file containing the input fields;
            - plotfile: name of the file that will contain the flux diagram;
            - logfile: name of the file containing the table as a .txt file.
        """
        lorenz = LorenzCycle()
        log = open(logfile, 'w')
        log.write('########################################################\n')
        log.write('#                                                      #\n')
        log.write('#      LORENZ     ENERGY    CYCLE                      #\n')
        log.write('#                                                      #\n')
        log.write('########################################################\n')
        filep = filenc
        dataset0 = Dataset(filenc)
        t_a = dataset0.variables['ta'][:, :, :, :]
        u_a = dataset0.variables['ua'][:, :, :, :]
        v_a = dataset0.variables['va'][:, :, :, :]
        wap = dataset0.variables['wap'][:, :, :, :]
        nlat = np.shape(t_a)[2]
        nfc = np.shape(t_a)[3]
        lev = dataset0.variables['plev'][:]
        nlev = len(lev)
        time = dataset0.variables['time'][:]
        ntime = len(time)
        lat = dataset0.variables['lat'][:]
        nlat = len(lat)
        if max(lev) < 1000:
            lev = lev * 100
            wap = wap * 100
        t_a = np.transpose(t_a, (1, 0, 2, 3))
        ta_r = t_a[:, :, :, 0::2]
        ta_i = t_a[:, :, :, 1::2]
        u_a = np.transpose(u_a, (1, 0, 2, 3))
        ua_r = u_a[:, :, :, 0::2]
        ua_i = u_a[:, :, :, 1::2]
        v_a = np.transpose(v_a, (1, 0, 2, 3))
        va_r = v_a[:, :, :, 0::2]
        va_i = v_a[:, :, :, 1::2]
        wap = np.transpose(wap, (1, 0, 2, 3))
        wap_r = wap[:, :, :, 0::2]
        wap_i = wap[:, :, :, 1::2]
        ta_c = ta_r + 1j * ta_i
        ua_c = ua_r + 1j * ua_i
        va_c = va_r + 1j * va_i
        wap_c = wap_r + 1j * wap_i
        log.write(' \n')
        log.write(' \n')
        log.write('INPUT DATA:\n')
        log.write('-----------\n')
        log.write(' \n')
        log.write('SPECTRAL RESOLUTION : {}\n'.format(nfc))
        log.write('NUMBER OF LATITUDES : {}\n'.format(nlat))
        log.write('NUMBER OF LEVEL : {}'.format(nlev))
        log.write('LEVEL  : {} Pa\n'.format(lev))
        log.write(' \n')
        ntp = nfc / 2 + 1
        log.write('WAVES:\n')
        log.write(' \n')
        log.write('(1) : 1 - {}\n'.format(NW_1))
        log.write('(2) : {} - {}\n'.format(NW_1, NW_2))
        log.write('(3) : {} - {}\n'.format(NW_2, NW_3))
    # Compute sigma level and dsigma
        sig = np.zeros(len(lev))
        for j_l in range(nlev):
            sig[j_l] = lev[j_l] / PS
        d_s = np.zeros(len(lev))
        for j_l in range(1, nlev - 1, 1):
            d_s[j_l] = 0.5 * abs(sig[j_l + 1] - sig[j_l - 1])
        d_s[0] = sig[0] + 0.5 * abs(sig[1] - sig[0])
        d_s[nlev - 1] = 1 - sig[nlev - 1] + 0.5 * abs(sig[nlev - 1] -
                                                      sig[nlev - 2])
    # Compute Gaussian weights
        g_w = np.zeros(nlat)
        y_l = np.zeros(nlat)
        for j_y in range(nlat):
            y_l[j_y] = np.deg2rad(lat[j_y])
            g_w[j_y] = np.cos(y_l[j_y])
        log.write(' \n')
        log.write('GLOBAL DIAGNOSTIC: \n')
        log.write('  \n')
        log.write('                            I GLOBAL I NORTH I SOUTH I\n')
        log.write('------------------------------------------------------\n')
    # Compute time mean
        ta_tmn = np.nanmean(ta_c, axis=1)
        ua_tmn = np.nanmean(ua_c, axis=1)
        va_tmn = np.nanmean(va_c, axis=1)
        wap_tmn = np.nanmean(wap_c, axis=1)
    # Compute zonal mean of time means
        ta_ztmn = np.squeeze(np.real(ta_tmn[:, :, 0]))
        ua_ztmn = np.squeeze(np.real(ua_tmn[:, :, 0]))
        va_ztmn = np.squeeze(np.real(va_tmn[:, :, 0]))
        wap_ztmn = np.squeeze(np.real(wap_tmn[:, :, 0]))
    # Compute global mean of time means
        ta_gmn = np.zeros(len(lev))
        ua_gmn = np.zeros(len(lev))
        va_gmn = np.zeros(len(lev))
        wap_gmn = np.zeros(len(lev))
        for j in range(nlev):
            ta_gmn[j] = np.nansum(ta_ztmn[j, :] * g_w) / np.nansum(g_w)
            ua_gmn[j] = np.nansum(ua_ztmn[j, :] * g_w) / np.nansum(g_w)
            va_gmn[j] = np.nansum(va_ztmn[j, :] * g_w) / np.nansum(g_w)
            wap_gmn[j] = np.nansum(wap_ztmn[j, :] * g_w) / np.nansum(g_w)
    # Compute stability parameter
        gam_ztmn = np.zeros([nlev, nlat])
        for l_l in range(nlat):
            gam_ztmn[:, l_l] = lorenz.stabil(ta_ztmn[:, l_l], lev, nlev)
        gam_tmn = lorenz.stabil(ta_gmn, lev, nlev)
        ek = np.zeros([nlev, ntime, nlat, ntp-1])
        ape = np.zeros([nlev, ntime, nlat, ntp-1])
        a2k = np.zeros([nlev, ntime, nlat, ntp-1])
        ae2az = np.zeros([nlev, ntime, nlat, ntp-1])
        ke2kz = np.zeros([nlev, ntime, nlat, ntp-1])
        at2as = np.zeros([nlev, ntime, nlat, ntp-1])
        kt2ks = np.zeros([nlev, ntime, nlat, ntp-1])
        for t_t in range(ntime):
            ta_t = ta_c[:, t_t, :, :]
            ua_t = ua_c[:, t_t, :, :]
            va_t = va_c[:, t_t, :, :]
            wap_t = wap_c[:, t_t, :, :]
            ta_tan = ta_c[:, t_t, :, :] - ta_tmn
            ua_tan = ua_c[:, t_t, :, :] - ua_tmn
            va_tan = va_c[:, t_t, :, :] - va_tmn
            wap_tan = wap_c[:, t_t, :, :] - wap_tmn
        # Compute zonal means
            ta_tzmn = np.squeeze(np.real(ta_t[:, :, 0]))
            ua_tzmn = np.squeeze(np.real(ua_t[:, :, 0]))
            va_tzmn = np.squeeze(np.real(va_t[:, :, 0]))
            wap_tzmn = np.squeeze(np.real(wap_t[:, :, 0]))
            ta_tzan = np.squeeze(np.real(ta_tan[:, :, 0]))
            ua_tzan = np.squeeze(np.real(ua_tan[:, :, 0]))
            va_tzan = np.squeeze(np.real(va_tan[:, :, 0]))
            wap_tzan = np.squeeze(np.real(wap_tan[:, :, 0]))
        # Compute global means as a function of levels
            ta_tgmn = np.zeros(len(lev))
            ua_tgmn = np.zeros(len(lev))
            va_tgmn = np.zeros(len(lev))
            wap_tgmn = np.zeros(len(lev))
            ta_tgan = np.zeros(len(lev))
            ua_tgan = np.zeros(len(lev))
            va_tgan = np.zeros(len(lev))
            wap_tgan = np.zeros(len(lev))
            for j in range(nlev):
                ta_tgmn[j] = np.nansum(ta_tzmn[j, :] * g_w) / np.nansum(g_w)
                ua_tgmn[j] = np.nansum(ua_tzmn[j, :] * g_w) / np.nansum(g_w)
                va_tgmn[j] = np.nansum(va_tzmn[j, :] * g_w) / np.nansum(g_w)
                wap_tgmn[j] = np.nansum(wap_tzmn[j, :] * g_w) / np.nansum(g_w)
                ta_tgan[j] = np.nansum(ta_tzan[j, :] * g_w) / np.nansum(g_w)
                ua_tgan[j] = np.nansum(ua_tzan[j, :] * g_w) / np.nansum(g_w)
                va_tgan[j] = np.nansum(va_tzan[j, :] * g_w) / np.nansum(g_w)
                wap_tgan[j] = np.nansum(wap_tzan[j, :] * g_w) / np.nansum(g_w)
         # Compute kinetic energy   
            ek[:, t_t, :, :] = lorenz.makek(ua_tan, va_tan, nlat, ntp, nlev)
         # Compute available potential energy
            ape[:, t_t, :, :] = lorenz.makea(ta_tan, ta_tgan, gam_tmn, nlat,
                                             ntp, nlev)
         # Compute conversion between kin.en. and pot.en.
            a2k[:, t_t, :, :] = lorenz.mka2k(wap_tan, ta_tan, wap_tgan,
                                             ta_tgan, lev, nlat, ntp, nlev)
         # Compute conversion between zonal and eddy APE
            ae2az[:, t_t, :, :] = lorenz.mkaeaz(va_tan, wap_tan, ta_tan,
                                                ta_tmn, ta_gmn, lev, y_l,
                                                gam_tmn, nlat, ntp, nlev)
         # Compute conversion between zonal and eddy KE
            ke2kz[:, t_t, :, :] = lorenz.mkkekz(ua_tan, va_tan, wap_tan,
                                                ua_tmn, va_tmn, lev, y_l,
                                                nlat, ntp, nlev)
         # Compute conversion between stationary and transient eddy APE
            at2as[:, t_t, :, :] = lorenz.mkatas(ua_tan, va_tan, wap_tan,
                                                ta_tan, ta_ztmn, gam_ztmn,
                                                lev, y_l, nlat, ntp, nlev)
         # Compute conversion between stationary and transient eddy KE
            kt2ks[:, t_t, :, :] = lorenz.mkktks(ua_tan, va_tan, wap_tan,
                                                ua_tmn, va_tmn, wap_tmn,
                                                lev, y_l, nlat, ntp, nlev)
        ek_tmn = np.nanmean(ek, axis=1)
        ek_tgmn = lorenz.globall_cg(ek_tmn, g_w, d_s, nlat, ntp, nlev)
        lorenz.table(ek_tgmn, ntp, 'TOT. KIN. EN.    ', log)
        ape_tmn = np.nanmean(ape, axis=1)
        ape_tgmn = lorenz.globall_cg(ape_tmn, g_w, d_s, nlat, ntp, nlev)
        lorenz.table(ape_tgmn, ntp, 'TOT. POT. EN.   ', log)
        a2k_tmn = np.nanmean(a2k, axis=1)
        a2k_tgmn = lorenz.globall_cg(a2k_tmn, g_w, d_s, nlat, ntp, nlev)
        lorenz.table_conv(a2k_tgmn, ntp, 'KE -> APE (trans) ', log)
        ae2az_tmn = np.nanmean(ae2az,axis=1)
        ae2az_tgmn = lorenz.globall_cg(ae2az_tmn, g_w, d_s, nlat, ntp, nlev)
        lorenz.table_conv(ae2az_tgmn, ntp, 'AZ <-> AE (trans) ', log)
        ke2kz_tmn = np.nanmean(ke2kz,axis=1)
        ke2kz_tgmn = lorenz.globall_cg(ke2kz_tmn, g_w, d_s, nlat, ntp, nlev)
        lorenz.table_conv(ke2kz_tgmn, ntp, 'KZ <-> KE (trans) ', log)
        at2as_tmn = np.nanmean(at2as, axis=1)
        at2as_tgmn = lorenz.globall_cg(at2as_tmn, g_w, d_s, nlat, ntp, nlev)
        lorenz.table_conv(at2as_tgmn, ntp, 'ASE  <->  ATE   ', log)
        kt2ks_tmn = np.nanmean(kt2ks, axis=1)
        kt2ks_tgmn = lorenz.globall_cg(kt2ks_tmn, g_w, d_s, nlat, ntp, nlev)
        lorenz.table_conv(kt2ks_tgmn, ntp, 'KSE  <->  KTE   ', log)
        ek_st = lorenz.makek(ua_tmn, va_tmn, nlat, ntp, nlev)
        ek_stgmn = lorenz.globall_cg(ek_st, g_w, d_s, nlat, ntp, nlev)
        lorenz.table(ek_stgmn, ntp, 'STAT. KIN. EN.    ', log)
        ape_st = lorenz.makea(ta_tmn,ta_gmn,gam_tmn,nlat,ntp,nlev)
        ape_stgmn = lorenz.globall_cg(ape_st, g_w, d_s, nlat, ntp, nlev)
        lorenz.table(ape_stgmn, ntp, 'STAT. POT. EN.    ', log)
        a2k_st = lorenz.mka2k(wap_tmn, ta_tmn, wap_gmn, ta_gmn, lev, nlat,
                              ntp, nlev)
        a2k_stgmn = lorenz.globall_cg(a2k_st, g_w, d_s, nlat, ntp, nlev)
        lorenz.table_conv(a2k_stgmn, ntp, 'KE -> APE (stat)', log)
        ae2az_st = lorenz.mkaeaz(va_tmn, wap_tmn, ta_tmn, ta_tmn, ta_gmn, lev,
                                 y_l, gam_tmn, nlat, ntp, nlev)
        ae2az_stgmn = lorenz.globall_cg(ae2az_st, g_w, d_s, nlat, ntp, nlev)
        lorenz.table_conv(ae2az_stgmn, ntp, 'AZ <-> AE (stat)', log)
        ke2kz_st = lorenz.mkkekz(ua_tmn, va_tmn, wap_tmn, ua_tmn, va_tmn, lev,
                                 y_l, nlat, ntp, nlev)
        ke2kz_stgmn = lorenz.globall_cg(ke2kz_st, g_w, d_s, nlat, ntp, nlev)
        lorenz.table_conv(ke2kz_stgmn, ntp, 'KZ <-> KE (stat)', log)
        apz = '{:.2f}'.format(ape_tgmn[0, 0] + ape_stgmn[0, 0])
        az2kz = '{:.2f}'.format(-1e5 * (a2k_tgmn[0, 0]))
        az2at = '{:.2f}'.format(-1e5 * np.nansum(ae2az_tgmn[0, 1:ntp - 1]))
        aps = '{:.2f}'.format(np.nansum(ape_stgmn[0, 1:ntp - 1]))
        as2ks = '{:.2f}'.format(1e5 * np.nansum(a2k_stgmn[0, 1:ntp - 1]))
        apt = '{:.2f}'.format(np.nansum(ape_tgmn[0, 1:ntp - 1]))
        at2kt = '{:.2f}'.format(1e5 * np.nansum(a2k_tgmn[0, 1:ntp - 1]))
        az2as = '{:.2f}'.format(-1e5 * np.nansum(ae2az_stgmn[0, 1:ntp - 1]))
        as2at = '{:.2f}'.format(1e5 * np.nansum(at2as_tgmn[0, 1:ntp - 1]))
        azin = '{:.2f}'.format((float(az2at) + float(az2as) - float(az2kz)))
        asein = '{:.2f}'.format((float(as2ks) + float(as2at) - float(az2as)))
        atein = '{:.2f}'.format(float(at2kt) - float(az2at) - float(as2at))
        kz = '{:.2f}'.format(ek_tgmn[0, 0] + ek_stgmn[0, 0])
        kte = '{:.2f}'.format(np.nansum(ek_tgmn[0, 1:ntp - 1]))
        kse = '{:.2f}'.format(np.nansum(ek_stgmn[0, 1:ntp - 1]))
        kt2kz = '{:.2f}'.format(1e5 * np.nansum(ke2kz_tgmn[0, 1:ntp - 1]))
        kt2ks = '{:.2f}'.format(-1e5 * np.nansum(kt2ks_tgmn[0, 1:ntp - 1]))
        ks2kz = '{:.2f}'.format(1e5 * np.nansum(ke2kz_stgmn[0, 1:ntp - 1]))
        kteout = '{:.2f}'.format(float(at2kt) - float(kt2ks) - float(kt2kz))
        kseout = '{:.2f}'.format(float(kt2ks) + float(as2ks) - float(ks2kz))
        kzout = '{:.2f}'.format(float(kt2kz) + float(ks2kz) - float(az2kz))
        lorenz.diagram(plotfile, azin, apz, asein, aps, atein, apt, as2ks,
                       at2kt, kteout, kte, kseout, kse, kzout, kz, az2kz,
                       az2at, az2as, as2at, kt2kz, kt2ks, ks2kz)
        lec_strength = float(kteout) + float(kseout) + float(kzout)
        ek_aux = np.zeros([nlev, nlat, ntp-1])
        ape_aux = np.zeros([nlev, nlat, ntp-1])
        a2k_aux = np.zeros([nlev, nlat, ntp-1])
        ae2az_aux = np.zeros([nlev, nlat, ntp-1])
        ke2kz_aux = np.zeros([nlev, nlat, ntp-1])
        for l in range(nlev):
            ek_aux[l, :, :] = ek_tmn[l, :, :] * d_s[l]
            ape_aux[l, :, :] = ape_tmn[l, :, :] * d_s[l]
            a2k_aux[l, :, :] = a2k_tmn[l, :, :] * d_s[l]
            ae2az_aux[l, :, :] = ae2az_tmn[l, :, :] * d_s[l]
            ke2kz_aux[l, :, :] = ke2kz_tmn[l, :, :] * d_s[l]
        ek_vmn = np.nansum(ek_aux, axis=0) / np.nansum(d_s)
        ape_vmn = np.nansum(ape_aux, axis=0) / np.nansum(d_s)
        a2k_vmn = np.nansum(a2k_aux, axis=0) / np.nansum(d_s)
        ae2az_vmn = np.nansum(ae2az_aux, axis=0) / np.nansum(d_s)
        ke2kz_vmn = np.nansum(ke2kz_aux, axis=0) / np.nansum(d_s)
        nc_f = outpath + '/ek_tmap_{}_{}.nc'.format(model,year)
        lorenz.removeif(nc_f)
        lorenz.pr_output(ek_vmn, 'ek', filep, nc_f, 1, verb=True)
        nc_f = outpath + '/ape_tmap_{}_{}.nc'.format(model, year)
        lorenz.removeif(nc_f)
        lorenz.pr_output(ape_vmn, 'ape', filep, nc_f, 1, verb=True)
        nc_f = outpath + '/a2k_tmap_{}_{}.nc'.format(model, year)
        lorenz.removeif(nc_f)
        lorenz.pr_output(a2k_vmn, 'a2k', filep, nc_f, 1, verb=True)
        nc_f = outpath + '/ae2az_tmap_{}_{}.nc'.format(model,year)
        lorenz.removeif(nc_f)
        lorenz.pr_output(ae2az_vmn, 'ae2az', filep, nc_f, 1, verb=True)
        nc_f = outpath + '/ke2kz_tmap_{}_{}.nc'.format(model,year)
        lorenz.removeif(nc_f)
        lorenz.pr_output(ke2kz_vmn, 'ke2kz', filep, nc_f, 1, verb=True)
        log.close()
        return lec_strength

    def bsslzr(self, kdim):
            
        NDIM = 50
        
        PI = math.pi
      
        zbes = [2.4048255577, 5.5200781103, 8.6537279129,  11.7915344391,
                14.9309177086, 18.0710639679, 21.2116366299, 24.3524715308,
                27.4934791320, 30.6346064684, 33.7758202136, 36.9170983537,
                40.0584257646, 43.1997917132, 46.3411883717, 49.4826098974,
                52.6240518411, 55.7655107550, 58.9069839261, 62.0484691902,
                65.1899648002, 68.3314693299, 71.4729816036, 74.6145006437,
                77.7560256304, 80.8975558711, 84.0390907769, 87.1806298436,
                90.3221726372, 93.4637187819, 96.6052679510, 99.7468198587,
                102.8883742542, 106.0299309165,109.1714896498, 112.3130502805,
                115.4546126537, 118.5961766309, 121.7377420880, 124.8793089132,
                128.0208770059, 131.1624462752, 134.3040166383, 137.4455880203,
                140.5871603528, 143.7287335737, 146.8703076258, 150.0118824570,
                153.1534580192, 156.2950342685]
        pbes = np.zeros(kdim)
        idim = min([kdim, NDIM])
        pbes[0:idim] = zbes[0:idim]
        for j in range(idim, kdim - 1, 1):
            pbes[j] = pbes[j - 1] + PI
        return(pbes)

    def diagram(self, filen, azin, apz, asein, aps, atein, apt, as2ks, at2kt, 
                kteout, kte, kseout, kse, kzout, kz, az2kz, az2at, az2as, 
                as2at, kt2kz, kt2ks, ks2kz):       
        """Diagram interface script.

        Call the class fluxogram, serving as        
        interface between the main script and the class for flux 
        diagrams design.
        
        @author: Valerio Lembo
        """
        FL = Fluxogram(1000, 1000, grid_size=20)
        FL.add_storage("AZ", 600, 0, 0)
        FL.add_storage("ASE", 600, 0.75, 0.25)
        FL.add_storage("ATE", 600, 1.5, 0)
        FL.add_storage("KTE", 600, 1.5, 1.5)
        FL.add_storage("KSE", 600, 0.75,1.25)
        FL.add_storage("KZ", 600, 0, 1.5)
        FL.add_storage("AZ+", 0, 0, -1)
        FL.add_storage("ASE+", 0, 0.75,-1)
        FL.add_storage("ATE+", 0, 1.5, -1)
        FL.add_storage("KTE-", 0, 1.5, 2.5)
        FL.add_storage("KSE-", 0, 0.75, 2.5)
        FL.add_storage("KZ-", 0, 0, 2.5)
        FL.add_flux("A2KZ", FL.storages[5], FL.storages[0], 100)
        FL.add_flux("AE2AZ", FL.storages[0], FL.storages[2], 150)
        FL.add_flux("AE2AS", FL.storages[0], FL.storages[1], 60)
        FL.add_flux("AE2AT", FL.storages[1], FL.storages[2], 60)
        FL.add_flux("A2KS", FL.storages[1], FL.storages[4], 60)
        FL.add_flux("A2KT", FL.storages[2], FL.storages[3], 100)
        FL.add_flux("KE2KS", FL.storages[3], FL.storages[4], 60)
        FL.add_flux("KS2KZ", FL.storages[4], FL.storages[5], 60)
        FL.add_flux("KE2KZ", FL.storages[3], FL.storages[5], 150)
        FL.add_flux("AZ+", FL.storages[6], FL.storages[0], 60)
        FL.add_flux("ASE+", FL.storages[7], FL.storages[1], 60)
        FL.add_flux("ATE+", FL.storages[8], FL.storages[2], 60)
        FL.add_flux("KTE-", FL.storages[3], FL.storages[9], 60)
        FL.add_flux("KSE-", FL.storages[4], FL.storages[10], 60)
        FL.add_flux("KZ-", FL.storages[5], FL.storages[11], 60)
        FL.draw(filen, azin, apz, asein, aps, atein, apt, as2ks, at2kt, kteout, 
                kte, kseout, kse, kzout, kz, az2kz, az2at, az2as, as2at, 
                kt2kz, kt2ks, ks2kz) 
         
    def gauaw(self, ny):
        """Compute the Gaussian coefficients for the Gaussian grid conversion.
        
        @author: Valerio Lembo
        """       
        lorenz = LorenzCycle()        
        c = (1 - (2 / math.pi) ** 2) / 4
        eps = 0.00000000000001    
        kk = ny/2
        pa = np.zeros(ny)
        pa[0:kk] = lorenz.bsslzr(kk)
        pw = np.zeros(ny)
        for i in range(kk):
            xz = np.cos(pa[i] / math.sqrt((ny + 0.5) ** 2 + c))
            iterr = 0.
            zsp = 1.0
            while (abs(zsp) > eps and iterr <= 10):
                pkm1 = xz
                pkm2 = 1.0
                for n in range(2, ny, 1):
                    pk = ((n * 2 - 1.0) * xz * pkm1 - (n - 1.0) * pkm2) / n
                    pkm2 = pkm1
                    pkm1 = pk
                pkm1 = pkm2
                pkmrk = (ny * (pkm1 - xz * pk)) / (1.0 - xz ** 2)
                zsp = pk / pkmrk
                xz = xz - zsp
                iterr = iterr + 1
            if iterr > 15:
                sys.exit("*** no convergence in gauaw ***")
            pa[i] = xz
            pw[i] = (2.0 * (1.0 - xz ** 2)) / ((ny ** 2) * (pkm1 ** 2))
            pa[ny - 1 - i] = - pa[i]
            pw[ny - 1 - i] = pw[i]
        psi = pa
        pgw = pw
        return psi, pgw

    def globall_cg(self, d3v, g_w, d_s, nlat, ntp, nlev):
        """Compute the global and hemispheric averages.
        
        @author: Valerio Lembo
        """    
        gmn = np.zeros([3, ntp-1])
        aux1 = np.zeros([nlev, nlat / 2, ntp - 1])
        aux2 = np.zeros([nlev, nlat / 2, ntp - 1])
        aux1v = np.zeros([nlev, ntp - 1])
        aux2v = np.zeros([nlev, ntp - 1])
        nhem = nlat / 2
        fac = 1 / G * PS / 1e5
        for l in range(nlev):
            for i in range(nhem):
                aux1[l, i, :] = fac * np.real(d3v[l, i, :]) * g_w[i]
                aux2[l, i, :] = fac * np.real(d3v[l, i + nhem - 1, :]) *\
                g_w[i + nhem - 1]
            aux1v[l, :] = np.nansum(aux1[l, :, :],
                                    axis = 0) / np.nansum(g_w[0:nhem]) * d_s[l]
            aux2v[l, :] = np.nansum(aux2[l, :, :],
                                    axis = 0) / np.nansum(g_w[0:nhem]) * d_s[l]
        gmn[1, :] = (np.nansum(aux1v, axis=0) / np.nansum(d_s))
        gmn[2, :] = (np.nansum(aux2v, axis=0) / np.nansum(d_s))
        gmn[0, :] = 0.5 * (gmn[1, :] + gmn[2, :])
        return(gmn)

    def makek(self, u, v, nlat, ntp, nlev):
        """Compute the kinetic energy reservoirs from u and v.
        
        @author: Valerio Lembo
        """  
        ek = np.zeros([nlev, nlat, ntp-1])
        ck1 = u * np.conj(u)
        ck2 = v * np.conj(v)
        ek[:, :, 0] = 0.5 * np.real(u[:, :, 0] * u[:, :, 0] +
                                    v[:, :, 0] * v[:, :, 0])
        ek = np.real(ck1 + ck2)
        ek[:, :, 0] = 0.5 * np.real(u[:, :, 0] * u[:, :, 0] +
                                    v[:, :, 0] * v[:, :, 0])
        return(ek)
    
        
    def makea(self, t, tg, gam, nlat, ntp, nlev):
        """Compute the kinetic energy reservoirs from t.
        
        @author: Valerio Lembo
        """  
        a = gam[:, np.newaxis, np.newaxis] * np.real(t * np.conj(t))
        a[:, :, 0] = gam[:, np.newaxis] * 0.5 * \
        np.real((t[:, :, 0] - tg[:, np.newaxis]) *
                (t[:, :, 0]-tg[:, np.newaxis]))        
        return(a)

    def mka2k(self, wap,t,wg,tg,p,nlat,ntp,nlev):
        """Compute the KE to APE energy conversions from t and w.
        
        @author: Valerio Lembo
        """  
        a2k = -R / p[:, np.newaxis, np.newaxis] *\
        (t * np.conj(wap) + np.conj(t) * wap)
        a2k[:, :, 0]=-R/p[:, np.newaxis]*(t[:, :, 0] -\
        tg[:, np.newaxis])*(wap[:, :, 0] - wg[:, np.newaxis])
        return(a2k)   
    
    def mkaeaz(self, v,wap, t, tt, ttg, p, lat, gam, nlat, ntp, nlev):
        """Compute the zonal mean - eddy APE conversions from t and v.
        
        @author: Valerio Lembo
        """  
        ae2az = np.zeros([nlev, nlat, ntp-1])
        dtdp = np.zeros([nlev, nlat])
        dtdy = np.zeros([nlev, nlat])
        for l in np.arange(nlev):
            if l == 0:
                t1 = np.real(tt[l, :, 0]) - ttg[l]
                t2 = np.real(tt[l + 1, :, 0]) - ttg[l + 1]
                dtdp[l, :] = (t2 - t1)/(p[l+1] - p[l])
            elif l == nlev-1:
                t1 = np.real(tt[l - 1, :, 0]) - ttg[l - 1]
                t2 = np.real(tt[l, :, 0]) - ttg[l]
                dtdp[l, :] = (t2 - t1)/(p[l] - p[l - 1])
            else:
                t1 = np.real(tt[l, :, 0])-ttg[l]
                t2 = np.real(tt[l + 1, :, 0])-ttg[l + 1]
                dtdp1 = (t2 - t1) / (p[l + 1]-p[l])
                t2 = t1
                t1 = np.real(tt[l - 1, :, 0])-ttg[l - 1]
                dtdp2 = (t2 - t1)/(p[l]-p[l - 1])
                dtdp[l, :] = (dtdp1*(p[l] - p[l - 1]) + \
                    dtdp2*(p[l + 1]-p[l]))/(p[l + 1]-p[l - 1])
            dtdp[l, :] = dtdp[l, :]-R / (CP * p[l])*(tt[l, :, 0] - ttg[l])
        for i in np.arange(nlat):
            if i == 0:
                t1 = np.real(tt[:, i, 0])
                t2 = np.real(tt[:, i + 1, 0])
                dtdy[:, i] = (t2 - t1)/(lat[i + 1] - lat[i])
            elif i == nlat-1:
                t1 = np.real(tt[:, i - 1,0])
                t2 = np.real(tt[:, i, 0])
                dtdy[:, i] = (t2 - t1) / (lat[i] - lat[i - 1])
            else:
                t1 = np.real(tt[:, i - 1, 0])
                t2 = np.real(tt[:, i + 1, 0])
                dtdy[:, i] = (t2 - t1) / (lat[i + 1] - lat[i - 1])
        dtdy = dtdy / AA
        c1 = np.real(v * np.conj(t) + t * np.conj(v))
        c2 = np.real(wap * np.conj(t) + t * np.conj(wap))
        ae2az = gam[:, np.newaxis, np.newaxis]\
        * (dtdy[:, :, np.newaxis] * c1 + dtdp[:, :, np.newaxis] * c2)      
        ae2az[:, :, 0] = 0.
        return(ae2az)

    def mkkekz(self, u, v, wap, ut, vt, p, lat, nlat, ntp, nlev):
        """Compute the zonal mean - eddy KE conversions from u and v.
        
        @author: Valerio Lembo
        """  
        dudp = np.zeros([nlev, nlat])
        dvdp = np.zeros([nlev, nlat])
        dudy = np.zeros([nlev, nlat])
        dvdy = np.zeros([nlev, nlat])
        for l in np.arange(nlev):
            if l == 0:
                dudp[l, :] = (np.real(ut[l + 1, :, 0] - ut[l, :, 0]))\
                / (p[l + 1] - p[l])
                dvdp[l, :] = (np.real(vt[l + 1, :, 0] - vt[l, :, 0]))\
                / (p[l + 1] - p[l])
            elif l == nlev-1:
                dudp[l, :] = (np.real(ut[l, :, 0] - ut[l - 1, :, 0]))\
                / (p[l]-p[l - 1])
                dvdp[l, :] = (np.real(vt[l, :, 0] - vt[l - 1, :, 0]))\
                / (p[l]-p[l - 1])
            else:
                dudp1 = (np.real(ut[l + 1, :, 0] - ut[l, :, 0]))\
                / (p[l + 1] - p[l])
                dvdp1 = (np.real(vt[l + 1, :, 0] - vt[l, :, 0]))\
                / (p[l + 1] - p[l])
                dudp2 = (np.real(ut[l, :, 0] - ut[l - 1, :, 0]))\
                / (p[l] - p[l - 1])
                dvdp2 = (np.real(vt[l, :, 0] - vt[l - 1, :, 0]))\
                / (p[l] - p[l - 1])
                dudp[l, :] = (dudp1 * (p[l] - p[l - 1]) +\
                dudp2 * (p[l + 1] - p[l])) / (p[l + 1] - p[l - 1])
                dvdp[l, :] = (dvdp1 * (p[l]-p[l - 1]) +\
                dvdp2 * (p[l + 1] - p[l])) / (p[l + 1] - p[l - 1])         
        for i in np.arange(nlat):    
            if i == 0:
                dudy[:, i]  = (np.real(ut[:, i + 1, 0] - ut[:, i, 0]))\
                / (lat[i + 1] - lat[i])
                dvdy[:, i]  = (np.real(vt[:, i + 1, 0] - vt[:, i, 0]))\
                / (lat[i + 1]-lat[i])
            elif i == nlat - 1:
                dudy[:, i]  = (np.real(ut[:, i, 0] - ut[:, i - 1, 0]))\
                / (lat[i] - lat[i - 1])
                dvdy[:, i]  = (np.real(vt[:, i, 0] - vt[:, i - 1, 0]))\
                / (lat[i] - lat[i - 1])
            else:
                dudy[:, i]  = (np.real(ut[:, i + 1, 0] - ut[:, i-1, 0]))\
                / (lat[i + 1] - lat[i - 1])
                dvdy[:, i]  = (np.real(vt[:, i+1, 0] - vt[:, i - 1, 0]))\
                / (lat[i + 1] - lat[i - 1])
        dudy  = dudy / AA
        dvdy  = dvdy / AA       
        c1 = np.zeros([nlev, nlat, ntp - 1])
        c2 = np.zeros([nlev, nlat, ntp - 1])
        c3 = np.zeros([nlev, nlat, ntp - 1])
        c4 = np.zeros([nlev, nlat, ntp - 1])
        c5 = np.zeros([nlev, nlat, ntp - 1])
        c6 = np.zeros([nlev, nlat, ntp - 1])
        ke2kz = np.zeros([nlev, nlat, ntp - 1])
        uu = u * np.conj(u) + u * np.conj(u)
        uv = u * np.conj(v) + v * np.conj(u)
        vv = v * np.conj(v) + v * np.conj(v)
        uw = u * np.conj(wap) + wap * np.conj(u)
        vw = v * np.conj(wap) + wap * np.conj(v)
        for i in np.arange(nlat):
            c1[:, i, :] = dudy[:, i][:, np.newaxis] * uv[:, i, :]
            c2[:, i, :] = dvdy[:, i][:, np.newaxis] * vv[:, i, :]
            c5[:, i, :] = np.tan(lat[i]) / AA \
            * np.real(ut[:, i, 0])[:, np.newaxis] * (uv[:, i, :])
            c6[:, i, :] = - np.tan(lat[i]) / AA \
            * np.real(vt[:, i, 0])[:, np.newaxis] * (uu[:, i, :])    
        for l in np.arange(nlev):
            c3[l, :, :] = dudp[l, :][:, np.newaxis] * uw[l, :, :]
            c4[l, :, :] = dvdp[l, :][:, np.newaxis] * vw[l, :, :]
        ke2kz = (c1 + c2 + c3 + c4 + c5 + c6)
        ke2kz[:, :, 0] = 0.
        return(ke2kz)

    def mkatas(self, u, v, wap, t, tt, g_w, p, lat, nlat, ntp, nlev):
        """Compute the stat.-trans. eddy APE conversions from u, v, wap and t.
        
        @author: Valerio Lembo
        """  
        tr = np.fft.ifft(t, axis=2)
        ur = np.fft.ifft(u, axis=2)
        vr = np.fft.ifft(v, axis=2)
        wr = np.fft.ifft(wap, axis=2)
        tur = tr * ur
        tvr = tr * vr
        twr = tr * wr
        tu = np.fft.fft(tur, axis=2)
        tv = np.fft.fft(tvr, axis=2)
        tw = np.fft.fft(twr, axis=2)
        c1 = np.zeros([nlev, nlat, ntp-1])
        c6 = np.zeros([nlev, nlat, ntp-1])
        c1 = tu * np.conj(tt[:, :, np.newaxis]) \
        - tt[:, :, np.newaxis] * np.conj(tu)
        c6 = tw * np.conj(tt[:, :, np.newaxis]) \
        - tt[:, :, np.newaxis] * np.conj(tw)
        c2 = np.zeros([nlev, nlat, ntp - 1])
        c3 = np.zeros([nlev, nlat, ntp - 1])
        c5 = np.zeros([nlev, nlat, ntp - 1])
        for i in range(nlat):
                if i == 0:
                    c2[:, i, :] = tv[:, i, :] / (AA * (lat[i + 1]-lat[i]))\
                    * np.conj(tt[:, i + 1, np.newaxis] - tt[:, i, np.newaxis]) 
                    c3[:, i, :] = np.conj(tv[:, i, :]) /\
                    (AA * (lat[i + 1] - lat[i])) * (tt[:, i + 1, np.newaxis] - 
                                                    tt[:, i, np.newaxis])
                elif i == nlat - 1:
                    c2[:, i, :] = tv[:, i, :] / (AA * (lat[i] - lat[i - 1])) \
                    * np.conj(tt[:, i, np.newaxis] - tt[:, i - 1, np.newaxis])
                    c3[:, i, :] = np.conj(tv[:, i, :]) /\
                    (AA * (lat[i] - lat[i - 1])) * (tt[:, i, np.newaxis] -
                                                    tt[:, i - 1, np.newaxis])
                else:
                    c2[:, i, :] = tv[:, i, :] / (AA * (lat[i + 1] - lat[i - 1]))\
                    * np.conj(tt[:, i + 1, np.newaxis] - tt[:, i - 1, np.newaxis])
                    c3[:, i, :] = np.conj(tv[:, i, :]) /\
                    (AA * (lat[i + 1] - lat[i - 1])) * (tt[:, i + 1, np.newaxis] -\
                                                        tt[:, i - 1, np.newaxis])       
        for l in range(nlev):
            if l == 0:
                c5[l, :, :] = (tt[l + 1, :, np.newaxis] 
                               - tt[l, :, np.newaxis]) / (p[l + 1] - p[l])
            elif l == nlev-1:
                c5[l, :, :] = (tt[l, :, np.newaxis] 
                               - tt[l - 1, :, np.newaxis]) / (p[l] - p[l - 1])
            else:
                c51 = (tt[l + 1, :, np.newaxis] - tt[l, :, np.newaxis])\
                / (p[l + 1] - p[l])
                c52 = (tt[l, :, np.newaxis] - tt[l - 1, :, np.newaxis])\
                / (p[l] - p[l - 1])
                c5[l, :, :]  = (c51 * (p[l] - p[l - 1])
                                + c52 * (p[l + 1] - p[l]))\
                / (p[l + 1] - p[l - 1])         
        K = np.arange(0, ntp - 1)            
        at2as = ((K - 1)[np.newaxis, np.newaxis, :] * np.imag(c1) \
                 / (AA * np.cos(lat[np.newaxis, :, np.newaxis])) \
                 + np.real(tw * np.conj(c5) + np.conj(tw) * c5) \
                 + np.real(c2 + c3) + R / (CP * p[:, np.newaxis, np.newaxis]) \
                 * np.real(c6)) * g_w[:, :, np.newaxis]        
        at2as[:, :, 0] = 0.
        return(at2as)    

    def mkktks(self, u, v, wap, ut, vt, wt, p, lat, nlat, ntp, nlev):    
        """Compute the stat.-trans. eddy KE conversions from u, v, wap and t.
        
        @author: Valerio Lembo
        """
        kt2ks = np.zeros([nlev, nlat, ntp - 1])
        c1 = np.zeros([nlev, nlat, ntp - 1])
        c21 = np.zeros([nlev, nlat, ntp - 1])
        c22 = np.zeros([nlev, nlat, ntp - 1])
        c3 = np.zeros([nlev, nlat, ntp - 1])
        c41 = np.zeros([nlev, nlat, ntp - 1])
        c42 = np.zeros([nlev, nlat, ntp - 1])
        c5 = np.zeros([nlev, nlat, ntp - 1])
        c6 = np.zeros([nlev, nlat, ntp - 1])
        dut = np.zeros([nlev, nlat, ntp - 1])
        dvt = np.zeros([nlev, nlat, ntp - 1])
        dlat =np.zeros([nlat])
        ur = np.fft.irfft(u, axis=2)
        vr = np.fft.irfft(v, axis=2)
        uur = ur * ur
        uvr = ur * vr
        vvr = vr * vr
        uu = np.fft.rfft(uur, axis=2)
        vv = np.fft.rfft(vvr, axis=2)
        uv = np.fft.rfft(uvr, axis=2)
        c1 = uu * np.conj(ut) - ut * np.conj(uu)
        c3 = uv * np.conj(ut) + ut * np.conj(uv)
        c5 = uu * np.conj(vt) + vt * np.conj(uu)
        c6 = uv * np.conj(vt) - vt * np.conj(uv)
        for i in range(nlat):
            if i == 0:
                dut[:, i, :]=(ut[:, i + 1, :] - ut[:, i, :])
                dvt[:, i, :]=(vt[:, i + 1, :] - vt[:, i, :])
                dlat[i] = (lat[i + 1] - lat[i])
            elif i == nlat-1:
                dut[:, i, :] = (ut[:, i, :] - ut[:, i - 1, :])
                dvt[:, i, :] = (vt[:, i, :] - vt[:, i - 1, :])
                dlat[i] = (lat[i] - lat[i - 1])
            else:
                dut[:, i, :] = (ut[:, i + 1, :] - ut[:, i - 1, :])
                dvt[:, i, :] = (vt[:, i + 1, :] - vt[:, i - 1, :])
                dlat[i] = (lat[i + 1] - lat[i - 1])           
        c21 = np.conj(uu) * dut / dlat[np.newaxis, :, np.newaxis]
        c22 = uu * np.conj(dut) / dlat[np.newaxis, :, np.newaxis]
        c41 = np.conj(vv) * dvt / dlat[np.newaxis, :, np.newaxis]
        c42 = vv * np.conj(dvt) / dlat[np.newaxis, :, np.newaxis]
        K = np.arange(0, ntp - 1)      
        kt2ks =  np.real(c21 + c22 + c41 + c42) / AA \
        + np.tan(lat)[np.newaxis, :, np.newaxis] * np.real(c1 - c5) / AA \
        + np.imag(c1 + c6) * (K - 1)[np.newaxis, np.newaxis, :] \
        / (AA * np.cos(lat)[np.newaxis, :, np.newaxis])
        kt2ks[:, :, 0] = 0
        return(kt2ks)   

    def ncdump(self, nc_fid, key):
        """Print the NetCDF file attributes for a given key.

        Arguments:
        - nc_fid: the ID of a NetCDF file containing variable 'key';
        - key: the name of a variable to obtain the attributes from.
        """
        nc_attrs = nc_fid.ncattrs()
        nc_dims = [dim for dim in nc_fid.dimensions]
        nc_vars = [var for var in nc_fid.variables]
        return nc_attrs, nc_dims, nc_vars

    def pr_output(self, varo, varname, filep, nc_f, opt, verb=True):
        """Print outputs to NetCDF.

        Save fields to NetCDF, retrieving information from an existing
        NetCDF file. Metadata are transferred from the existing file to the
        new one.
        Arguments:
            - varo: the field to be stored;
            - varname: the name of the variables to be saved;
            - filep: the existing dataset, containing the metadata;
            - nc_f: the name of the output file;
            - opt: depends on the shape of the output variable to be saved;
            
        PROGRAMMER(S)
            Chris Slocum (2014), modified by Valerio Lembo (2018).
        """
        lorenz = LorenzCycle()
        nc_fid = Dataset(filep, 'r')
        nc_attrs, nc_dims, nc_vars = lorenz.ncdump(nc_fid, 'var130', verb)
        # Extract data from NetCDF file
        lats = nc_fid.variables['lat'][:]  # extract/copy the data
        time = nc_fid.variables['time'][:]
        wave = nc_fid.variables['wave'][:]
        ntp = len(wave) / 2
        # Writing NetCDF files
        w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
        w_nc_fid.description = "Outputs of LEC program"
        if opt == 0:
            w_nc_fid.createDimension('time', None)
            w_nc_dim = w_nc_fid.createVariable('time',
                                               nc_fid.variables['time'].dtype,
                                               ('time',))
            for ncattr in nc_fid.variables['time'].ncattrs():
                w_nc_dim.setncattr(ncattr,
                                   nc_fid.variables['time'].getncattr(ncattr))
            w_nc_fid.variables['time'][:] = time
            w_nc_fid.createDimension('lat', len(lats))
            w_nc_dim = w_nc_fid.createVariable('lat',
                                               nc_fid.variables['lat'].dtype,\
                                               ('lat',))
            for ncattr in nc_fid.variables['lat'].ncattrs():
                w_nc_dim.setncattr(ncattr,
                                   nc_fid.variables['lat'].getncattr(ncattr))
            w_nc_fid.variables['lat'][:] = lats
            w_nc_fid.createDimension('wave', ntp)
            w_nc_dim = w_nc_fid.createVariable('wave',
                                               nc_fid.variables['wave'].dtype,\
                                               ('wave',))
            w_nc_fid.variables['wave'][:] = wave[0:ntp]
            w_nc_var = w_nc_fid.createVariable(varname, 'f8',
                                               ('time', 'lat', 'wave'))
            lorenz.varatts(w_nc_var,varname,0,0)
            w_nc_fid.variables[varname][:] = varo
            w_nc_fid.close()
        elif opt == 1:
            w_nc_fid.createDimension('lat', len(lats))
            w_nc_dim = w_nc_fid.createVariable('lat',
                                               nc_fid.variables['lat'].dtype,\
                                               ('lat',))
            for ncattr in nc_fid.variables['lat'].ncattrs():
                w_nc_dim.setncattr(ncattr,
                                   nc_fid.variables['lat'].getncattr(ncattr))
            w_nc_fid.variables['lat'][:] = lats
            w_nc_fid.createDimension('wave', ntp)
            w_nc_dim = w_nc_fid.createVariable('wave',
                                               nc_fid.variables['wave'].dtype,\
                                               ('wave',))
            for ncattr in nc_fid.variables['wave'].ncattrs():
                w_nc_dim.setncattr(ncattr,
                                   nc_fid.variables['wave'].getncattr(ncattr))
            w_nc_fid.variables['wave'][:] = wave[0:ntp]
            w_nc_var = w_nc_fid.createVariable(varname, 'f8', ('lat', 'wave'))
            lorenz.varatts(w_nc_var,varname, 1, 0)
            w_nc_fid.variables[varname][:] = varo
            w_nc_fid.close()
        elif opt == 2:
            w_nc_fid.createDimension('lat', len(lats))
            w_nc_dim = w_nc_fid.createVariable('lat',
                                               nc_fid.variables['lat'].dtype,\
                                               ('lat',))
            for ncattr in nc_fid.variables['lat'].ncattrs():
                w_nc_dim.setncattr(ncattr,
                                   nc_fid.variables['lat'].getncattr(ncattr))
            w_nc_fid.variables['lat'][:] = lats
            w_nc_var = w_nc_fid.createVariable(varname, 'f8', ('lat'))
            lorenz.varatts(w_nc_var,varname, 1, 0)
            w_nc_fid.variables[varname][:] = varo
            w_nc_fid.close()
        elif opt == 3:
            w_nc_fid.createDimension('hem', 3)
            w_nc_dim = w_nc_fid.createVariable('hem',
                                               nc_fid.variables['hem'].dtype,\
                                               ('hem',))
            w_nc_fid.variables['hem'][:] = [0, 1, 2]
            w_nc_fid.createDimension('time', None)
            w_nc_dim = w_nc_fid.createVariable('time',
                                               nc_fid.variables['time'].dtype,\
                                               ('time',))
            for ncattr in nc_fid.variables['time'].ncattrs():
                w_nc_dim.setncattr(ncattr,
                                   nc_fid.variables['time'].getncattr(ncattr))
            w_nc_fid.variables['time'][:] = time
            w_nc_var = w_nc_fid.createVariable(varname, 'f8', ('hem', 'time'))
            lorenz.varatts(w_nc_var,varname,0,0)
            w_nc_fid.variables[varname][:] = varo
            w_nc_fid.close()  # close the new file
        elif opt == 4:
            w_nc_fid.createDimension('hem', 3)
            w_nc_dim = w_nc_fid.createVariable('hem',
                                               nc_fid.variables['hem'].dtype,\
                                               ('hem',))
            w_nc_fid.variables['hem'][:] = [0, 1, 2]
            w_nc_fid.createDimension('time', None)
            w_nc_dim = w_nc_fid.createVariable('time',
                                               nc_fid.variables['time'].dtype,\
                                               ('time',))
            for ncattr in nc_fid.variables['time'].ncattrs():
                w_nc_dim.setncattr(ncattr,
                                   nc_fid.variables['time'].getncattr(ncattr))
            w_nc_fid.variables['time'][:] = time
            w_nc_fid.createDimension('wave', ntp)
            w_nc_dim = w_nc_fid.createVariable('wave',
                                               nc_fid.variables['wave'].dtype,\
                                               ('wave',))
            for ncattr in nc_fid.variables['wave'].ncattrs():
                w_nc_dim.setncattr(ncattr,
                                   nc_fid.variables['wave'].getncattr(ncattr))
            w_nc_fid.variables['wave'][:] = wave[0:ntp]
            w_nc_var = w_nc_fid.createVariable(varname, 'f8',
                                               ('hem', 'time', 'wave'))
            lorenz.varatts(w_nc_var,varname, 0, 0)
            w_nc_fid.variables[varname][:] = varo
            w_nc_fid.close()
        nc_fid.close()
        
    def removeif(self, filename):
        """Remove filename if it exists."""
        try:
            os.remove(filename)
        except OSError:
            pass
  
    def stabil(self, ta_gmn,p,nlev):
        """Compute the stability parameter from temp. and pressure levels.
        
        @author: Valerio Lembo
        """
        cpdr = CP / R
        t = ta_gmn
        gs = np.zeros(nlev)
        for i in range(nlev):
            if i == 0:
                dtdp=(t[i + 1] - t[i]) / (p[i + 1] - p[i])
            elif i == nlev - 1:
                dtdp = (t[i] - t[i - 1]) / (p[i] - p[i-1])
            else:
                dtdp1 = (t[i + 1] - t[i]) / (p[i + 1] - p[i])
                dtdp2 = (t[i] - t[i - 1]) / (p[i] - p[i - 1])
                dtdp = (dtdp1 * (p[i] - p[i - 1]) + dtdp2 * (p[i + 1] - p[i]))\
                / (p[i + 1] - p[i - 1])
            gs[i] = CP / (t[i] - p[i] * dtdp * cpdr)
        return gs

    def table(self, varin, ntp, name, log):
        """Write global and hem. storage terms to .txt table.
        
        @author: Valerio Lembo
        """
        varzon = varin[:,0]
        vared = np.nansum(varin[:, 1:ntp - 1],axis = 1)
        vared1 = np.nansum(varin[:, 1:NW_1-1],axis = 1)
        vared2 = np.nansum(varin[:, NW_1:NW_2-  1],axis = 1)
        vared3 = np.nansum(varin[:, NW_2:NW_3 - 1],axis = 1)
        vartot = varzon + vared
        log.write(' {} TOTAL    {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,
                  vartot[0], vartot[1], vartot[2]))
        log.write('--------------------------------------\n')
        log.write(' {} ZONAL    {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,
                  varzon[0], varzon[1], varzon[2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY     {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,
                  vared[0], vared[1], vared[2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY(LW) {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,
                  vared1[0], vared1[1], vared1[2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY(SW) {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,
                  vared2[0], vared2[1], vared2[2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY(KW) {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,
                  vared3[0], vared3[1], vared3[2]))
        log.write('--------------------------------------\n')
        
    
    def table_conv(self, varin, ntp, name, log):
        """Write global and hem. conversion terms to .txt table.
        
        @author: Valerio Lembo
        """
        fac = 1e5
        varin = fac * varin
        varzon = varin[:,0]
        vared  = np.nansum(varin[:, 1:ntp - 1], axis=1)
        vared1 = np.nansum(varin[:, 1:NW_1 - 1], axis=1)
        vared2 = np.nansum(varin[:, NW_1:NW_2 - 1], axis=1)
        vared3 = np.nansum(varin[:, NW_2:NW_3 - 1], axis=1)
        vartot = varzon + vared
        log.write(' {} TOTAL    {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,
                  vartot[0], vartot[1], vartot[2]))
        log.write('--------------------------------------\n')
        log.write(' {} ZONAL    {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,
                  varzon[0], varzon[1], varzon[2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY     {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,
                  vared[0], vared[1], vared[2]))
        log.write('-------------------------------------\n')
        log.write(' {} EDDY(LW) {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,
                  vared1[0], vared1[1], vared1[2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY(SW) {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,
                  vared2[0], vared2[1], vared2[2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY(KW) {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,
                  vared3[0], vared3[1], vared3[2]))
        log.write('--------------------------------------\n')   
    
    def varatts(self, w_nc_var, varname, tres, vres):
        """Add attibutes to the variables, depending on name and time res.

        Arguments:
        - w_nc_var: a variable object;
        - varname: the name of the variable, among ta, ua, va and wap;
        - tres: the time resolution (daily or annual);
        - vres: the vertical resolution (pressure levels or vert. integr.).
        """
        if tres == 0:
            tatt = u"Daily\nM"
        elif tres == 1:
            tatt = u"Annual mean\nM"
        if vres == 0:
            vatt = u"Pressure levels\n"
        elif vres == 1:
            vatt = u"Vertically integrated\n"
        if varname == 'a':
            w_nc_var.setncatts({'long_name': u"Available Potential Energy",
                                'units': u"W m-2", 'level_desc': vatt,
                                'var_desc': u"APE -> KE",
                                'statistic': tatt})
        elif varname == 'ek':
            w_nc_var.setncatts({'long_name': u"Kinetic Energy",
                                'units': u"W m-2", 'level_desc': vatt,
                                'var_desc': u"APE -> KE",
                                'statistic': tatt})
        elif varname == 'a2k':
            w_nc_var.setncatts({'long_name': u"Conversion between APE and KE",
                                'units': u"W m-2", 'level_desc': vatt,
                                'var_desc': u"APE <-> KE",
                                'statistic': tatt})
        elif varname == 'k':
            w_nc_var.setncatts({'long_name': u"Kinetic Energy",
                                'units': u"W m-2", 'level_desc': vatt,
                                'var_desc': u"APE -> KE",
                                'statistic': tatt})
