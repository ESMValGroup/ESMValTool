"""PROGRAM FOR LEC COMPUTATION.

The module contains the following functions:
    - lorenz: it is the main program, calling functions that compute the
              reservoirs and conversion terms, storing them separately in
              NetCDF files and providing a flux diagram and a table outputs,
              the latter separately for the two hemispheres;
    - averages: a script computing time, global and zonal averages;
    - averages_comp: a script computing global mean of the output fields;
    - bsslzr: it contains the coefficients for the conversion from regular
              lonlat grid to Gaussian grid;
    - diagram: it is the interface between the main program and a
               class "Fluxogram", producing the flux diagram;
    - gauaw: it uses the coefficients provided in bsslzr for the lonlat to
             Gaussian grid conversion;
    - globall_cg: it computes the global and hemispheric means at each
                  timestep;
    - init: initializes the table and ingests input fields;
    - makek: computes the KE reservoirs;
    - makea: computes the APE reservoirs;
    - mka2k: computes the APE->KE conversion terms;
    - mkaeaz: computes the zonal APE - eddy APE conversion terms;
    - mkkekz: computes the zonal KE - eddy KE conversion terms;
    - mkatas: computes the stationay eddy - transient eddy APE conversions;
    - mkktks: computes the stationay eddy - transient eddy KE conversions;
    - output: compute vertical integrals and print NC output;
    - preprocess_lec: a script handling the input files, separating the real
                      from imaginary part of the Fourier coefficients,
                      reordering the latitudinal dimension (from N to S),
                      interpolating on a reference sigma coordinate,
    - pr_output: prints a single component of the LEC computations to a
                 single Nc file;
    - removeif: removes a file if it exists;
    - stabil: calculates the stability parameter;
    - table: prints the global and hemispheric mean values of
             the reservoirs;
    - table_conv: prints the global and hemispheric mean values of the
                  conversion terms;
    - varatts: prints the attributes of a variable in a Nc file;
    - weights: computes the weights for vertical integrations and meridional
               averages;
    - write_to_tab: a script for writing global and hemispheric means to table;

References.
    Ulbrich P. and P. Speth (1991) The global energy cycle of stationary
    and transient atmospheric waves: Results from ECMWF analyses, Met.

@author: valerio.lembo@uni-hamburg.de, Valerio Lembo, Hamburg University, 2018.
"""

import math
import os
import sys

import numpy as np
from cdo import Cdo
from netCDF4 import Dataset

import esmvaltool.diag_scripts.shared as e
from esmvaltool.diag_scripts.thermodyn_diagtool import (fluxogram,
                                                        fourier_coefficients)

G = 9.81
R = 287.00
CP = 1003.5
AA = 6.371E6
PS = 101100.0
NW_1 = 3
NW_2 = 9
NW_3 = 21


def lorenz(outpath, model, year, filenc, plotfile, logfile):
    """Manage input and output fields and calling functions.

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
    ta_c, ua_c, va_c, wap_c, dims, lev, lat, log = init(logfile, filenc)
    nlev = int(dims[0])
    ntime = int(dims[1])
    nlat = int(dims[2])
    ntp = int(dims[3])
    d_s, y_l, g_w = weights(lev, nlev, lat)
    # Compute time mean
    ta_tmn = np.nanmean(ta_c, axis=1)
    ta_ztmn, ta_gmn = averages(ta_tmn, g_w)
    ua_tmn = np.nanmean(ua_c, axis=1)
    va_tmn = np.nanmean(va_c, axis=1)
    wap_tmn = np.nanmean(wap_c, axis=1)
    _, wap_gmn = averages(wap_tmn, g_w)
    # Compute stability parameter
    gam_ztmn = np.zeros([nlev, nlat])
    for l_l in range(nlat):
        gam_ztmn[:, l_l] = stabil(ta_ztmn[:, l_l], lev, nlev)
    gam_tmn = stabil(ta_gmn, lev, nlev)
    e_k = np.zeros([nlev, ntime, nlat, ntp - 1])
    ape = np.zeros([nlev, ntime, nlat, ntp - 1])
    a2k = np.zeros([nlev, ntime, nlat, ntp - 1])
    ae2az = np.zeros([nlev, ntime, nlat, ntp - 1])
    ke2kz = np.zeros([nlev, ntime, nlat, ntp - 1])
    at2as = np.zeros([nlev, ntime, nlat, ntp - 1])
    kt2ks = np.zeros([nlev, ntime, nlat, ntp - 1])
    for t_t in range(ntime):
        ta_tan = ta_c[:, t_t, :, :] - ta_tmn
        ua_tan = ua_c[:, t_t, :, :] - ua_tmn
        va_tan = va_c[:, t_t, :, :] - va_tmn
        wap_tan = wap_c[:, t_t, :, :] - wap_tmn
        # Compute zonal means
        _, ta_tgan = averages(ta_tan, g_w)
        _, wap_tgan = averages(wap_tan, g_w)
        # Compute kinetic energy
        e_k[:, t_t, :, :] = makek(ua_tan, va_tan)
        # Compute available potential energy
        ape[:, t_t, :, :] = makea(ta_tan, ta_tgan, gam_tmn)
        # Compute conversion between kin.en. and pot.en.
        a2k[:, t_t, :, :] = mka2k(wap_tan, ta_tan, wap_tgan, ta_tgan, lev)
        # Compute conversion between zonal and eddy APE
        ae2az[:, t_t, :, :] = mkaeaz(va_tan, wap_tan, ta_tan, ta_tmn, ta_gmn,
                                     lev, y_l, gam_tmn, nlat, nlev)
        # Compute conversion between zonal and eddy KE
        ke2kz[:, t_t, :, :] = mkkekz(ua_tan, va_tan, wap_tan, ua_tmn, va_tmn,
                                     lev, y_l, nlat, ntp, nlev)
        # Compute conversion between stationary and transient eddy APE
        at2as[:, t_t, :, :] = mkatas(ua_tan, va_tan, wap_tan, ta_tan, ta_ztmn,
                                     gam_ztmn, lev, y_l, nlat, ntp, nlev)
        # Compute conversion between stationary and transient eddy KE
        kt2ks[:, t_t, :, :] = mkktks(ua_tan, va_tan, ua_tmn, va_tmn, y_l, nlat,
                                     ntp, nlev)
    ek_tgmn = averages_comp(e_k, g_w, d_s, dims)
    table(ek_tgmn, ntp, 'TOT. KIN. EN.    ', logfile, flag=0)
    ape_tgmn = averages_comp(ape, g_w, d_s, dims)
    table(ape_tgmn, ntp, 'TOT. POT. EN.   ', logfile, flag=0)
    a2k_tgmn = averages_comp(a2k, g_w, d_s, dims)
    table(a2k_tgmn, ntp, 'KE -> APE (trans) ', logfile, flag=1)
    ae2az_tgmn = averages_comp(ae2az, g_w, d_s, dims)
    table(ae2az_tgmn, ntp, 'AZ <-> AE (trans) ', logfile, flag=1)
    ke2kz_tgmn = averages_comp(ke2kz, g_w, d_s, dims)
    table(ke2kz_tgmn, ntp, 'KZ <-> KE (trans) ', logfile, flag=1)
    at2as_tgmn = averages_comp(at2as, g_w, d_s, dims)
    table(at2as_tgmn, ntp, 'ASE  <->  ATE   ', logfile, flag=1)
    kt2ks_tgmn = averages_comp(kt2ks, g_w, d_s, dims)
    table(kt2ks_tgmn, ntp, 'KSE  <->  KTE   ', logfile, flag=1)
    ek_st = makek(ua_tmn, va_tmn)
    ek_stgmn = globall_cg(ek_st, g_w, d_s, dims)
    table(ek_stgmn, ntp, 'STAT. KIN. EN.    ', logfile, flag=0)
    ape_st = makea(ta_tmn, ta_gmn, gam_tmn)
    ape_stgmn = globall_cg(ape_st, g_w, d_s, dims)
    table(ape_stgmn, ntp, 'STAT. POT. EN.    ', logfile, flag=0)
    a2k_st = mka2k(wap_tmn, ta_tmn, wap_gmn, ta_gmn, lev)
    a2k_stgmn = globall_cg(a2k_st, g_w, d_s, dims)
    table(a2k_stgmn, ntp, 'KE -> APE (stat)', logfile, flag=1)
    ae2az_st = mkaeaz(va_tmn, wap_tmn, ta_tmn, ta_tmn, ta_gmn, lev, y_l,
                      gam_tmn, nlat, nlev)
    ae2az_stgmn = globall_cg(ae2az_st, g_w, d_s, dims)
    table(ae2az_stgmn, ntp, 'AZ <-> AE (stat)', logfile, flag=1)
    ke2kz_st = mkkekz(ua_tmn, va_tmn, wap_tmn, ua_tmn, va_tmn, lev, y_l, nlat,
                      ntp, nlev)
    ke2kz_stgmn = globall_cg(ke2kz_st, g_w, d_s, dims)
    table(ke2kz_stgmn, ntp, 'KZ <-> KE (stat)', logfile, flag=1)
    list_diag = [
        ape_tgmn, ape_stgmn, ek_tgmn, ek_stgmn, ae2az_tgmn, ae2az_stgmn,
        a2k_tgmn, a2k_stgmn, at2as_tgmn, kt2ks_tgmn, ke2kz_tgmn, ke2kz_stgmn
    ]
    lec_strength = diagram(plotfile, list_diag, dims)
    nc_f = outpath + '/ek_tmap_{}_{}.nc'.format(model, year)
    output(e_k, d_s, filenc, 'ek', nc_f)
    nc_f = outpath + '/ape_tmap_{}_{}.nc'.format(model, year)
    output(ape, d_s, filenc, 'ape', nc_f)
    nc_f = outpath + '/a2k_tmap_{}_{}.nc'.format(model, year)
    output(a2k, d_s, filenc, 'a2k', nc_f)
    nc_f = outpath + '/ae2az_tmap_{}_{}.nc'.format(model, year)
    output(ae2az, d_s, filenc, 'ae2az', nc_f)
    nc_f = outpath + '/ke2kz_tmap_{}_{}.nc'.format(model, year)
    output(ke2kz, d_s, filenc, 'ke2kz', nc_f)
    log.close()
    return lec_strength


def averages(x_c, g_w):
    """Compute time, zonal and global mean averages of initial fields.

    Arguments:
    - x_c: the input field as (lev, lat, wave);
    - g_w: the Gaussian weights for meridional averaging;
    """
    xc_ztmn = np.squeeze(np.real(x_c[:, :, 0]))
    xc_gmn = np.nansum(xc_ztmn * g_w[np.newaxis, :], axis=1) / np.nansum(g_w)
    return xc_ztmn, xc_gmn


def averages_comp(fld, g_w, d_s, dims):
    """Compute the global mean averages of reservoirs and conversion terms.

    Arguments:
    - fld: the component of the LEC (time, lev, lat, wave);
    - g_w: the Gaussian weights for meridional averaging;
    - d_s: the Delta sigma of the sigma levels;
    - dims: a list containing the dimensions length0;
    """
    fld_tmn = np.nanmean(fld, axis=1)
    fld_tgmn = globall_cg(fld_tmn, g_w, d_s, dims)
    return fld_tgmn


def bsslzr(kdim):
    """Obtain parameters for the Gaussian coefficients.

    @author: Valerio Lembo
    """
    ndim = 50
    p_i = math.pi
    zbes = [
        2.4048255577, 5.5200781103, 8.6537279129, 11.7915344391, 14.9309177086,
        18.0710639679, 21.2116366299, 24.3524715308, 27.4934791320,
        30.6346064684, 33.7758202136, 36.9170983537, 40.0584257646,
        43.1997917132, 46.3411883717, 49.4826098974, 52.6240518411,
        55.7655107550, 58.9069839261, 62.0484691902, 65.1899648002,
        68.3314693299, 71.4729816036, 74.6145006437, 77.7560256304,
        80.8975558711, 84.0390907769, 87.1806298436, 90.3221726372,
        93.4637187819, 96.6052679510, 99.7468198587, 102.8883742542,
        106.0299309165, 109.1714896498, 112.3130502805, 115.4546126537,
        118.5961766309, 121.7377420880, 124.8793089132, 128.0208770059,
        131.1624462752, 134.3040166383, 137.4455880203, 140.5871603528,
        143.7287335737, 146.8703076258, 150.0118824570, 153.1534580192,
        156.2950342685
    ]
    pbes = np.zeros(kdim)
    idim = min([kdim, ndim])
    pbes[0:idim] = zbes[0:idim]
    for j in range(idim, kdim - 1, 1):
        pbes[j] = pbes[j - 1] + p_i
    return pbes


def diagram(filen, listf, dims):
    """Diagram interface script.

    Call the class fluxogram, serving as
    interface between the main script and the class for flux
    diagrams design.

    Arguments:
    - filen: the filename of the diagram flux;
    - listf: a list containing the fluxes and storages;
    - dims: the dimensions of the variables;
    """
    ntp = int(dims[3])
    apet = listf[0]
    apes = listf[1]
    ekt = listf[2]
    eks = listf[3]
    ae2azt = listf[4]
    ae2azs = listf[5]
    a2kt = listf[6]
    a2ks = listf[7]
    at2as = listf[8]
    kt2ks = listf[9]
    ke2kzt = listf[10]
    ke2kzs = listf[11]
    apz = '{:.2f}'.format(apet[0, 0] + apes[0, 0])
    az2kz = '{:.2f}'.format(-1e5 * (a2kt[0, 0]))
    az2at = '{:.2f}'.format(-1e5 * np.nansum(ae2azt[0, 1:ntp - 1]))
    aps = '{:.2f}'.format(np.nansum(apes[0, 1:ntp - 1]))
    as2ks = '{:.2f}'.format(1e5 * np.nansum(a2ks[0, 1:ntp - 1]))
    apt = '{:.2f}'.format(np.nansum(apet[0, 1:ntp - 1]))
    at2kt = '{:.2f}'.format(1e5 * np.nansum(a2kt[0, 1:ntp - 1]))
    az2as = '{:.2f}'.format(-1e5 * np.nansum(ae2azs[0, 1:ntp - 1]))
    as2at = '{:.2f}'.format(1e5 * np.nansum(at2as[0, 1:ntp - 1]))
    azin = '{:.2f}'.format((float(az2at) + float(az2as) - float(az2kz)))
    asein = '{:.2f}'.format((float(as2ks) + float(as2at) - float(az2as)))
    atein = '{:.2f}'.format(float(at2kt) - float(az2at) - float(as2at))
    k_z = '{:.2f}'.format(ekt[0, 0] + eks[0, 0])
    kte = '{:.2f}'.format(np.nansum(ekt[0, 1:ntp - 1]))
    kse = '{:.2f}'.format(np.nansum(eks[0, 1:ntp - 1]))
    kt2kz = '{:.2f}'.format(1e5 * np.nansum(ke2kzt[0, 1:ntp - 1]))
    ks2kt = '{:.2f}'.format(-1e5 * np.nansum(kt2ks[0, 1:ntp - 1]))
    ks2kz = '{:.2f}'.format(1e5 * np.nansum(ke2kzs[0, 1:ntp - 1]))
    kteout = '{:.2f}'.format(float(at2kt) - float(ks2kt) - float(kt2kz))
    kseout = '{:.2f}'.format(float(ks2kt) + float(as2ks) - float(ks2kz))
    kzout = '{:.2f}'.format(float(kt2kz) + float(ks2kz) - float(az2kz))
    list_lorenz = [
        azin, apz, asein, aps, atein, apt, as2ks, at2kt, kteout, kte, kseout,
        kse, kzout, k_z, az2kz, az2at, az2as, as2at, kt2kz, ks2kt, ks2kz
    ]
    flux = fluxogram.Fluxogram(1000, 1000)
    flux.add_storage("AZ", 600, 0, 0)
    flux.add_storage("ASE", 600, 0.75, 0.25)
    flux.add_storage("ATE", 600, 1.5, 0)
    flux.add_storage("KTE", 600, 1.5, 1.5)
    flux.add_storage("KSE", 600, 0.75, 1.25)
    flux.add_storage("KZ", 600, 0, 1.5)
    flux.add_storage("AZ+", 0, 0, -1)
    flux.add_storage("ASE+", 0, 0.75, -1)
    flux.add_storage("ATE+", 0, 1.5, -1)
    flux.add_storage("KTE-", 0, 1.5, 2.5)
    flux.add_storage("KSE-", 0, 0.75, 2.5)
    flux.add_storage("KZ-", 0, 0, 2.5)
    flux.add_flux("A2KZ", flux.storages[5], flux.storages[0], 100)
    flux.add_flux("AE2AZ", flux.storages[0], flux.storages[2], 150)
    flux.add_flux("AE2AS", flux.storages[0], flux.storages[1], 60)
    flux.add_flux("AE2AT", flux.storages[1], flux.storages[2], 60)
    flux.add_flux("A2KS", flux.storages[1], flux.storages[4], 60)
    flux.add_flux("A2KT", flux.storages[2], flux.storages[3], 100)
    flux.add_flux("KE2KS", flux.storages[3], flux.storages[4], 60)
    flux.add_flux("KS2KZ", flux.storages[4], flux.storages[5], 60)
    flux.add_flux("KE2KZ", flux.storages[3], flux.storages[5], 150)
    flux.add_flux("AZ+", flux.storages[6], flux.storages[0], 60)
    flux.add_flux("ASE+", flux.storages[7], flux.storages[1], 60)
    flux.add_flux("ATE+", flux.storages[8], flux.storages[2], 60)
    flux.add_flux("KTE-", flux.storages[3], flux.storages[9], 60)
    flux.add_flux("KSE-", flux.storages[4], flux.storages[10], 60)
    flux.add_flux("KZ-", flux.storages[5], flux.storages[11], 60)
    flux.draw(filen, list_lorenz)
    lec = float(kteout) + float(kseout) + float(kzout)
    return lec


def gauaw(n_y):
    """Compute the Gaussian coefficients for the Gaussian grid conversion.

    Arguments:
    - n_y: the latitude dimension;
    """
    c_c = (1 - (2 / math.pi)**2) / 4
    eps = 0.00000000000001
    k_k = n_y / 2
    p_a = np.zeros(n_y)
    p_a[0:k_k] = bsslzr(k_k)
    p_w = np.zeros(n_y)
    for i_l in range(k_k):
        x_z = np.cos(p_a[i_l] / math.sqrt((n_y + 0.5)**2 + c_c))
        iterr = 0.
        zsp = 1.0
        while (abs(zsp) > eps and iterr <= 10):
            pkm1 = x_z
            pkm2 = 1.0
            for n_n in range(2, n_y, 1):
                p_k = ((n_n * 2 - 1.0) * x_z * pkm1 - (n_n - 1.0) * pkm2) / n_n
                pkm2 = pkm1
                pkm1 = p_k
            pkm1 = pkm2
            pkmrk = (n_y * (pkm1 - x_z * p_k)) / (1.0 - x_z**2)
            zsp = p_k / pkmrk
            x_z = x_z - zsp
            iterr = iterr + 1
        if iterr > 15:
            sys.exit("*** no convergence in gauaw ***")
        p_a[i_l] = x_z
        p_w[i_l] = (2.0 * (1.0 - x_z**2)) / ((n_y**2) * (pkm1**2))
        p_a[n_y - 1 - i_l] = -p_a[i_l]
        p_w[n_y - 1 - i_l] = p_w[i_l]
    psi = p_a
    pgw = p_w
    return psi, pgw


def globall_cg(d3v, g_w, d_s, dims):
    """Compute the global and hemispheric averages.

    Arguments:
    - d3v: the 3D dataset to be averaged;
    - g_w: the gaussian weights;
    - d_s: the vertical levels;
    - dims: a list containing the sizes of the dimensions;
    """
    nlev = int(dims[0])
    nlat = int(dims[2])
    ntp = int(dims[3])
    gmn = np.zeros([3, ntp - 1])
    aux1 = np.zeros([nlev, int(nlat / 2), ntp - 1])
    aux2 = np.zeros([nlev, int(nlat / 2), ntp - 1])
    aux1v = np.zeros([nlev, ntp - 1])
    aux2v = np.zeros([nlev, ntp - 1])
    nhem = int(nlat / 2)
    fac = 1 / G * PS / 1e5
    for l_l in range(nlev):
        for i_h in range(nhem):
            aux1[l_l, i_h, :] = fac * np.real(d3v[l_l, i_h, :]) * g_w[i_h]
            aux2[l_l, i_h, :] = (fac * np.real(d3v[l_l, i_h + nhem - 1, :]) *
                                 g_w[i_h + nhem - 1])
        aux1v[l_l, :] = (np.nansum(aux1[l_l, :, :], axis=0) /
                         np.nansum(g_w[0:nhem]) * d_s[l_l])
        aux2v[l_l, :] = (np.nansum(aux2[l_l, :, :], axis=0) /
                         np.nansum(g_w[0:nhem]) * d_s[l_l])
    gmn[1, :] = (np.nansum(aux1v, axis=0) / np.nansum(d_s))
    gmn[2, :] = (np.nansum(aux2v, axis=0) / np.nansum(d_s))
    gmn[0, :] = 0.5 * (gmn[1, :] + gmn[2, :])
    return gmn


def init(logfile, filep):
    """Ingest input fields as complex fields and initialise tables.

    Receive fields t,u,v,w as input fields in Fourier
    coefficients  (time,level,wave,lon), with real as even and imaginary parts
    as odd. Convert them to complex fields for Python.

    Arguments:
        - filenc: name of the file containing the input fields;
        - logfile: name of the file containing the table as a .txt file.
    """
    with open(logfile, 'w') as log:
        log.write('########################################################\n')
        log.write('#                                                      #\n')
        log.write('#      LORENZ     ENERGY    CYCLE                      #\n')
        log.write('#                                                      #\n')
        log.write('########################################################\n')
    with Dataset(filep) as dataset0:
        t_a = dataset0.variables['ta'][:, :, :, :]
        u_a = dataset0.variables['ua'][:, :, :, :]
        v_a = dataset0.variables['va'][:, :, :, :]
        wap = dataset0.variables['wap'][:, :, :, :]
        lev = dataset0.variables['plev'][:]
        time = dataset0.variables['time'][:]
        lat = dataset0.variables['lat'][:]
    nfc = np.shape(t_a)[3]
    nlev = len(lev)
    ntime = len(time)
    nlat = len(lat)
    ntp = nfc / 2 + 1
    dims = [nlev, ntime, nlat, ntp]
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
    with open(logfile, 'w') as log:
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
        log.write('WAVES:\n')
        log.write(' \n')
        log.write('(1) : 1 - {}\n'.format(NW_1))
        log.write('(2) : {} - {}\n'.format(NW_1, NW_2))
        log.write('(3) : {} - {}\n'.format(NW_2, NW_3))
        log.write(' \n')
        log.write('GLOBAL DIAGNOSTIC: \n')
        log.write('  \n')
        log.write('                            I GLOBAL I NORTH I SOUTH I\n')
        log.write('------------------------------------------------------\n')
    return ta_c, ua_c, va_c, wap_c, dims, lev, lat, log


def makek(u_t, v_t):
    """Compute the kinetic energy reservoirs from u and v.

    Arguments:
    - u_t: a 3D zonal velocity field;
    - v_t: a 3D meridional velocity field;
    """
    ck1 = u_t * np.conj(u_t)
    ck2 = v_t * np.conj(v_t)
    e_k = np.real(ck1 + ck2)
    e_k[:, :, 0] = 0.5 * np.real(u_t[:, :, 0] * u_t[:, :, 0] +
                                 v_t[:, :, 0] * v_t[:, :, 0])
    return e_k


def makea(t_t, t_g, gam):
    """Compute the kinetic energy reservoirs from t.

    Arguments:
    - t_t_ a 3D temperature field;
    - t_g: a temperature vertical profile;
    - gam: a vertical profile of the stability parameter;
    """
    ape = gam[:, np.newaxis, np.newaxis] * np.real(t_t * np.conj(t_t))
    ape[:, :, 0] = (gam[:, np.newaxis] * 0.5 * np.real(
        (t_t[:, :, 0] - t_g[:, np.newaxis]) *
        (t_t[:, :, 0] - t_g[:, np.newaxis])))
    return ape


def mka2k(wap, t_t, w_g, t_g, p_l):
    """Compute the KE to APE energy conversions from t and w.

    Arguments:
    - wap: a 3D vertical velocity field;
    - t_t: a 3D temperature field;
    - w_g: a vertical velocity vertical profile;
    - t_g: a temperature vertical profile;
    - p_l: the pressure levels;
    """
    a2k = -(R / p_l[:, np.newaxis, np.newaxis] *
            (t_t * np.conj(wap) + np.conj(t_t) * wap))
    a2k[:, :, 0] = -(R / p_l[:, np.newaxis] *
                     (t_t[:, :, 0] - t_g[:, np.newaxis]) *
                     (wap[:, :, 0] - w_g[:, np.newaxis]))
    return a2k


def mkaeaz(v_t, wap, t_t, ttt, ttg, p_l, lat, gam, nlat, nlev):
    """Compute the zonal mean - eddy APE conversions from t and v.

    Arguments:
    - v_t: a 3D meridional velocity field;
    - wap: a 3D vertical velocity field;
    - t_t: a 3D temperature field;
    - ttt: a climatological mean 3D temperature field;
    - p_l: the pressure levels;
    - lat: the latudinal dimension;
    - gam: a vertical profile of the stability parameter;
    - nlat: the number of latitudes;
    - nlev: the number of levels;
    """
    dtdp = np.zeros([nlev, nlat])
    dtdy = np.zeros([nlev, nlat])
    for l_l in np.arange(nlev):
        if l_l == 0:
            t_1 = np.real(ttt[l_l, :, 0]) - ttg[l_l]
            t_2 = np.real(ttt[l_l + 1, :, 0]) - ttg[l_l + 1]
            dtdp[l_l, :] = (t_2 - t_1) / (p_l[l_l + 1] - p_l[l_l])
        elif l_l == nlev - 1:
            t_1 = np.real(ttt[l_l - 1, :, 0]) - ttg[l_l - 1]
            t_2 = np.real(ttt[l_l, :, 0]) - ttg[l_l]
            dtdp[l_l, :] = (t_2 - t_1) / (p_l[l_l] - p_l[l_l - 1])
        else:
            t_1 = np.real(ttt[l_l, :, 0]) - ttg[l_l]
            t_2 = np.real(ttt[l_l + 1, :, 0]) - ttg[l_l + 1]
            dtdp1 = (t_2 - t_1) / (p_l[l_l + 1] - p_l[l_l])
            t_2 = t_1
            t_1 = np.real(ttt[l_l - 1, :, 0]) - ttg[l_l - 1]
            dtdp2 = (t_2 - t_1) / (p_l[l_l] - p_l[l_l - 1])
            dtdp[l_l, :] = ((dtdp1 * (p_l[l_l] - p_l[l_l - 1]) + dtdp2 *
                             (p_l[l_l + 1] - p_l[l_l])) /
                            (p_l[l_l + 1] - p_l[l_l - 1]))
        dtdp[l_l, :] = dtdp[l_l, :] - (R / (CP * p_l[l_l]) *
                                       (ttt[l_l, :, 0] - ttg[l_l]))
    for i_l in np.arange(nlat):
        if i_l == 0:
            t_1 = np.real(ttt[:, i_l, 0])
            t_2 = np.real(ttt[:, i_l + 1, 0])
            dtdy[:, i_l] = (t_2 - t_1) / (lat[i_l + 1] - lat[i_l])
        elif i_l == nlat - 1:
            t_1 = np.real(ttt[:, i_l - 1, 0])
            t_2 = np.real(ttt[:, i_l, 0])
            dtdy[:, i_l] = (t_2 - t_1) / (lat[i_l] - lat[i_l - 1])
        else:
            t_1 = np.real(ttt[:, i_l - 1, 0])
            t_2 = np.real(ttt[:, i_l + 1, 0])
            dtdy[:, i_l] = (t_2 - t_1) / (lat[i_l + 1] - lat[i_l - 1])
    dtdy = dtdy / AA
    c_1 = np.real(v_t * np.conj(t_t) + t_t * np.conj(v_t))
    c_2 = np.real(wap * np.conj(t_t) + t_t * np.conj(wap))
    ae2az = (gam[:, np.newaxis, np.newaxis] *
             (dtdy[:, :, np.newaxis] * c_1 + dtdp[:, :, np.newaxis] * c_2))
    ae2az[:, :, 0] = 0.
    return ae2az


def mkkekz(u_t, v_t, wap, utt, vtt, p_l, lat, nlat, ntp, nlev):
    """Compute the zonal mean - eddy KE conversions from u and v.

    Arguments:
    - u_t: a 3D zonal velocity field;
    - v_t: a 3D meridional velocity field;
    - wap: a 3D vertical velocity field;
    - utt: a climatological mean 3D zonal velocity field;
    - vtt: a climatological mean 3D meridional velocity field;
    - p_l: the pressure levels;
    - lat: the latitude dimension;
    - nlat: the number of latitudes;
    - ntp: the number of wavenumbers;
    - nlev: the number of vertical levels;
    """
    dudp = np.zeros([nlev, nlat])
    dvdp = np.zeros([nlev, nlat])
    dudy = np.zeros([nlev, nlat])
    dvdy = np.zeros([nlev, nlat])
    for l_l in np.arange(nlev):
        if l_l == 0:
            dudp[l_l, :] = ((np.real(utt[l_l + 1, :, 0] - utt[l_l, :, 0])) /
                            (p_l[l_l + 1] - p_l[l_l]))
            dvdp[l_l, :] = ((np.real(vtt[l_l + 1, :, 0] - vtt[l_l, :, 0])) /
                            (p_l[l_l + 1] - p_l[l_l]))
        elif l_l == nlev - 1:
            dudp[l_l, :] = ((np.real(utt[l_l, :, 0] - utt[l_l - 1, :, 0])) /
                            (p_l[l_l] - p_l[l_l - 1]))
            dvdp[l_l, :] = ((np.real(vtt[l_l, :, 0] - vtt[l_l - 1, :, 0])) /
                            (p_l[l_l] - p_l[l_l - 1]))
        else:
            dudp1 = ((np.real(utt[l_l + 1, :, 0] - utt[l_l, :, 0])) /
                     (p_l[l_l + 1] - p_l[l_l]))
            dvdp1 = ((np.real(vtt[l_l + 1, :, 0] - vtt[l_l, :, 0])) /
                     (p_l[l_l + 1] - p_l[l_l]))
            dudp2 = ((np.real(utt[l_l, :, 0] - utt[l_l - 1, :, 0])) /
                     (p_l[l_l] - p_l[l_l - 1]))
            dvdp2 = ((np.real(vtt[l_l, :, 0] - vtt[l_l - 1, :, 0])) /
                     (p_l[l_l] - p_l[l_l - 1]))
            dudp[l_l, :] = ((dudp1 * (p_l[l_l] - p_l[l_l - 1]) + dudp2 *
                             (p_l[l_l + 1] - p_l[l_l])) /
                            (p_l[l_l + 1] - p_l[l_l - 1]))
            dvdp[l_l, :] = ((dvdp1 * (p_l[l_l] - p_l[l_l - 1]) + dvdp2 *
                             (p_l[l_l + 1] - p_l[l_l])) /
                            (p_l[l_l + 1] - p_l[l_l - 1]))
    for i_l in np.arange(nlat):
        if i_l == 0:
            dudy[:, i_l] = ((np.real(utt[:, i_l + 1, 0] - utt[:, i_l, 0])) /
                            (lat[i_l + 1] - lat[i_l]))
            dvdy[:, i_l] = ((np.real(vtt[:, i_l + 1, 0] - vtt[:, i_l, 0])) /
                            (lat[i_l + 1] - lat[i_l]))
        elif i_l == nlat - 1:
            dudy[:, i_l] = ((np.real(utt[:, i_l, 0] - utt[:, i_l - 1, 0])) /
                            (lat[i_l] - lat[i_l - 1]))
            dvdy[:, i_l] = ((np.real(vtt[:, i_l, 0] - vtt[:, i_l - 1, 0])) /
                            (lat[i_l] - lat[i_l - 1]))
        else:
            dudy[:, i_l] = (
                (np.real(utt[:, i_l + 1, 0] - utt[:, i_l - 1, 0])) /
                (lat[i_l + 1] - lat[i_l - 1]))
            dvdy[:, i_l] = (
                (np.real(vtt[:, i_l + 1, 0] - vtt[:, i_l - 1, 0])) /
                (lat[i_l + 1] - lat[i_l - 1]))
    dudy = dudy / AA
    dvdy = dvdy / AA
    c_1 = np.zeros([nlev, nlat, ntp - 1])
    c_2 = np.zeros([nlev, nlat, ntp - 1])
    c_3 = np.zeros([nlev, nlat, ntp - 1])
    c_4 = np.zeros([nlev, nlat, ntp - 1])
    c_5 = np.zeros([nlev, nlat, ntp - 1])
    c_6 = np.zeros([nlev, nlat, ntp - 1])
    u_u = u_t * np.conj(u_t) + u_t * np.conj(u_t)
    u_v = u_t * np.conj(v_t) + v_t * np.conj(u_t)
    v_v = v_t * np.conj(v_t) + v_t * np.conj(v_t)
    u_w = u_t * np.conj(wap) + wap * np.conj(u_t)
    v_w = v_t * np.conj(wap) + wap * np.conj(v_t)
    for i_l in np.arange(nlat):
        c_1[:, i_l, :] = dudy[:, i_l][:, np.newaxis] * u_v[:, i_l, :]
        c_2[:, i_l, :] = dvdy[:, i_l][:, np.newaxis] * v_v[:, i_l, :]
        c_5[:, i_l, :] = (np.tan(lat[i_l]) / AA *
                          np.real(utt[:, i_l, 0])[:, np.newaxis] *
                          (u_v[:, i_l, :]))
        c_6[:, i_l, :] = -(np.tan(lat[i_l]) / AA *
                           np.real(vtt[:, i_l, 0])[:, np.newaxis] *
                           (u_u[:, i_l, :]))
    for l_l in np.arange(nlev):
        c_3[l_l, :, :] = dudp[l_l, :][:, np.newaxis] * u_w[l_l, :, :]
        c_4[l_l, :, :] = dvdp[l_l, :][:, np.newaxis] * v_w[l_l, :, :]
    ke2kz = (c_1 + c_2 + c_3 + c_4 + c_5 + c_6)
    ke2kz[:, :, 0] = 0.
    return ke2kz


def mkatas(u_t, v_t, wap, t_t, ttt, g_w, p_l, lat, nlat, ntp, nlev):
    """Compute the stat.-trans. eddy APE conversions from u, v, wap and t.

    Arguments:
    - u_t: a 3D zonal velocity field;
    - v_t: a 3D meridional velocity field;
    - wap: a 3D vertical velocity field;
    - t_t: a 3D temperature field;
    - ttt: a climatological mean 3D temperature field;
    - g_w: the gaussian weights;
    - p_l: the pressure levels;
    - lat: the latitude dimension;
    - nlat: the number of latitudes;
    - ntp: the number of wavenumbers;
    - nlev: the number of vertical levels;
    """
    t_r = np.fft.ifft(t_t, axis=2)
    u_r = np.fft.ifft(u_t, axis=2)
    v_r = np.fft.ifft(v_t, axis=2)
    w_r = np.fft.ifft(wap, axis=2)
    tur = t_r * u_r
    tvr = t_r * v_r
    twr = t_r * w_r
    t_u = np.fft.fft(tur, axis=2)
    t_v = np.fft.fft(tvr, axis=2)
    t_w = np.fft.fft(twr, axis=2)
    c_1 = (t_u * np.conj(ttt[:, :, np.newaxis]) -
           ttt[:, :, np.newaxis] * np.conj(t_u))
    c_6 = (t_w * np.conj(ttt[:, :, np.newaxis]) -
           ttt[:, :, np.newaxis] * np.conj(t_w))
    c_2 = np.zeros([nlev, nlat, ntp - 1])
    c_3 = np.zeros([nlev, nlat, ntp - 1])
    c_5 = np.zeros([nlev, nlat, ntp - 1])
    for i_l in range(nlat):
        if i_l == 0:
            c_2[:, i_l, :] = (
                t_v[:, i_l, :] / (AA * (lat[i_l + 1] - lat[i_l])) *
                np.conj(ttt[:, i_l + 1, np.newaxis] - ttt[:, i_l, np.newaxis]))
            c_3[:, i_l, :] = (
                np.conj(t_v[:, i_l, :]) / (AA * (lat[i_l + 1] - lat[i_l])) *
                (ttt[:, i_l + 1, np.newaxis] - ttt[:, i_l, np.newaxis]))
        elif i_l == nlat - 1:
            c_2[:, i_l, :] = (
                t_v[:, i_l, :] / (AA * (lat[i_l] - lat[i_l - 1])) *
                np.conj(ttt[:, i_l, np.newaxis] - ttt[:, i_l - 1, np.newaxis]))
            c_3[:, i_l, :] = (
                np.conj(t_v[:, i_l, :]) / (AA * (lat[i_l] - lat[i_l - 1])) *
                (ttt[:, i_l, np.newaxis] - ttt[:, i_l - 1, np.newaxis]))
        else:
            c_2[:, i_l, :] = (t_v[:, i_l, :] /
                              (AA * (lat[i_l + 1] - lat[i_l - 1])) *
                              np.conj(ttt[:, i_l + 1, np.newaxis] -
                                      ttt[:, i_l - 1, np.newaxis]))
            c_3[:, i_l, :] = (
                np.conj(t_v[:, i_l, :]) / (AA *
                                           (lat[i_l + 1] - lat[i_l - 1])) *
                (ttt[:, i_l + 1, np.newaxis] - ttt[:, i_l - 1, np.newaxis]))
    for l_l in range(nlev):
        if l_l == 0:
            c_5[l_l, :, :] = (
                (ttt[l_l + 1, :, np.newaxis] - ttt[l_l, :, np.newaxis]) /
                (p_l[l_l + 1] - p_l[l_l]))
        elif l_l == nlev - 1:
            c_5[l_l, :, :] = (
                (ttt[l_l, :, np.newaxis] - ttt[l_l - 1, :, np.newaxis]) /
                (p_l[l_l] - p_l[l_l - 1]))
        else:
            c51 = ((ttt[l_l + 1, :, np.newaxis] - ttt[l_l, :, np.newaxis]) /
                   (p_l[l_l + 1] - p_l[l_l]))
            c52 = ((ttt[l_l, :, np.newaxis] - ttt[l_l - 1, :, np.newaxis]) /
                   (p_l[l_l] - p_l[l_l - 1]))
            c_5[l_l, :, :] = ((c51 * (p_l[l_l] - p_l[l_l - 1]) + c52 *
                               (p_l[l_l + 1] - p_l[l_l])) /
                              (p_l[l_l + 1] - p_l[l_l - 1]))
    k_k = np.arange(0, ntp - 1)
    at2as = (((k_k - 1)[np.newaxis, np.newaxis, :] * np.imag(c_1) /
              (AA * np.cos(lat[np.newaxis, :, np.newaxis])) +
              np.real(t_w * np.conj(c_5) + np.conj(t_w) * c_5) +
              np.real(c_2 + c_3) + R /
              (CP * p_l[:, np.newaxis, np.newaxis]) * np.real(c_6)) *
             g_w[:, :, np.newaxis])
    at2as[:, :, 0] = 0.
    return at2as


def mkktks(u_t, v_t, utt, vtt, lat, nlat, ntp, nlev):
    """Compute the stat.-trans. eddy KE conversions from u, v and t.

    Arguments:
    - u_t: a 3D zonal velocity field;
    - v_t: a 3D meridional velocity field;
    - utt: a climatological mean 3D zonal velocity field;
    - vtt: a climatological mean 3D meridional velocity field;
    - lat: the latitude dimension;
    - nlat: the number of latitudes;
    - ntp: the number of wavenumbers;
    - nlev: the number of vertical levels;
    """
    dut = np.zeros([nlev, nlat, ntp - 1])
    dvt = np.zeros([nlev, nlat, ntp - 1])
    dlat = np.zeros([nlat])
    u_r = np.fft.irfft(u_t, axis=2)
    v_r = np.fft.irfft(v_t, axis=2)
    uur = u_r * u_r
    uvr = u_r * v_r
    vvr = v_r * v_r
    u_u = np.fft.rfft(uur, axis=2)
    v_v = np.fft.rfft(vvr, axis=2)
    u_v = np.fft.rfft(uvr, axis=2)
    c_1 = u_u * np.conj(u_t) - u_t * np.conj(u_u)
    # c_3 = u_v * np.conj(u_t) + u_t * np.conj(u_v)
    c_5 = u_u * np.conj(v_t) + v_t * np.conj(u_u)
    c_6 = u_v * np.conj(v_t) - v_t * np.conj(u_v)
    for i_l in range(nlat):
        if i_l == 0:
            dut[:, i_l, :] = (utt[:, i_l + 1, :] - utt[:, i_l, :])
            dvt[:, i_l, :] = (vtt[:, i_l + 1, :] - vtt[:, i_l, :])
            dlat[i_l] = (lat[i_l + 1] - lat[i_l])
        elif i_l == nlat - 1:
            dut[:, i_l, :] = (utt[:, i_l, :] - utt[:, i_l - 1, :])
            dvt[:, i_l, :] = (vtt[:, i_l, :] - vtt[:, i_l - 1, :])
            dlat[i_l] = (lat[i_l] - lat[i_l - 1])
        else:
            dut[:, i_l, :] = (utt[:, i_l + 1, :] - utt[:, i_l - 1, :])
            dvt[:, i_l, :] = (vtt[:, i_l + 1, :] - vtt[:, i_l - 1, :])
            dlat[i_l] = (lat[i_l + 1] - lat[i_l - 1])
    c21 = np.conj(u_u) * dut / dlat[np.newaxis, :, np.newaxis]
    c22 = u_u * np.conj(dut) / dlat[np.newaxis, :, np.newaxis]
    c41 = np.conj(v_v) * dvt / dlat[np.newaxis, :, np.newaxis]
    c42 = v_v * np.conj(dvt) / dlat[np.newaxis, :, np.newaxis]
    k_k = np.arange(0, ntp - 1)
    kt2ks = (np.real(c21 + c22 + c41 + c42) / AA +
             np.tan(lat)[np.newaxis, :, np.newaxis] * np.real(c_1 - c_5) / AA +
             np.imag(c_1 + c_6) * (k_k - 1)[np.newaxis, np.newaxis, :] /
             (AA * np.cos(lat)[np.newaxis, :, np.newaxis]))
    kt2ks[:, :, 0] = 0
    return kt2ks


def output(fld, d_s, filenc, name, nc_f):
    """Compute vertical integrals and print (time,lat,ntp) to NC output.

    Arguments:
    - fld: the annual mean fields (lev, lat, wave);
    - d_s: Delta sigma;
    - filenc: the input file containing the Fourier coefficients of t,u,v,w;
    - name: the variable name;
    - nc_f: the name of the output file (with path)
    """
    fld_tmn = np.nanmean(fld, axis=1)
    fld_aux = fld_tmn * d_s[:, np.newaxis, np.newaxis]
    fld_vmn = np.nansum(fld_aux, axis=0) / np.nansum(d_s)
    removeif(nc_f)
    pr_output(fld_vmn, name, filenc, nc_f)


def pr_output(varo, varname, filep, nc_f):
    """Print outputs to NetCDF.

    Save fields to NetCDF, retrieving information from an existing
    NetCDF file. Metadata are transferred from the existing file to the
    new one.
    Arguments:
    - varo: the field to be stored;
    - varname: the name of the variables to be saved;
    - filep: the existing dataset, containing the metadata;
    - nc_f: the name of the output file;

    PROGRAMMER(S)
        Chris Slocum (2014), modified by Valerio Lembo (2018).
    """
    fourc = fourier_coefficients
    with Dataset(nc_f, 'w', format='NETCDF4') as w_nc_fid:
        w_nc_fid.description = "Outputs of LEC program"
        with Dataset(filep, 'r') as nc_fid:
            # Extract data from NetCDF file
            wave = nc_fid.variables['wave'][:]
            ntp = int(len(wave) / 2)
            # Writing NetCDF files
            fourc.extr_lat(nc_fid, w_nc_fid, 'lat')
            w_nc_fid.createDimension('wave', ntp)
            w_nc_dim = w_nc_fid.createVariable('wave',
                                               nc_fid.variables['wave'].dtype,
                                               ('wave', ))
            for ncattr in nc_fid.variables['wave'].ncattrs():
                w_nc_dim.setncattr(ncattr,
                                   nc_fid.variables['wave'].getncattr(ncattr))
        w_nc_fid.variables['wave'][:] = wave[0:ntp]
        w_nc_var = w_nc_fid.createVariable(varname, 'f8', ('lat', 'wave'))
        varatts(w_nc_var, varname, 1, 0)
        w_nc_fid.variables[varname][:] = varo


def preproc_lec(model, wdir, pdir, input_data):
    """Preprocess fields for LEC computations and send it to lorenz program.

    This function computes the interpolation of ta, ua, va, wap daily fields to
    fill gaps using near-surface data, then computes the Fourier coefficients
    and performs the LEC computations. For every year, (lev,lat,wave) fields,
    global and hemispheric time series of each conversion and reservoir term
    of the LEC is provided.

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - pdir: a new directory is created as a sub-directory of the plot directory
      to store tables of conversion/reservoir terms and the flux diagram for
      year;
    - filelist: a list of file names containing the input fields;
    """
    cdo = Cdo()
    fourc = fourier_coefficients
    ta_file = e.select_metadata(input_data, short_name='ta',
                                dataset=model)[0]['filename']
    tas_file = e.select_metadata(input_data, short_name='tas',
                                 dataset=model)[0]['filename']
    ua_file = e.select_metadata(input_data, short_name='ua',
                                dataset=model)[0]['filename']
    uas_file = e.select_metadata(input_data, short_name='uas',
                                 dataset=model)[0]['filename']
    va_file = e.select_metadata(input_data, short_name='va',
                                dataset=model)[0]['filename']
    vas_file = e.select_metadata(input_data, short_name='vas',
                                 dataset=model)[0]['filename']
    wap_file = e.select_metadata(input_data, short_name='wap',
                                 dataset=model)[0]['filename']
    ldir = os.path.join(pdir, 'LEC_results')
    os.makedirs(ldir)
    maskorog = wdir + '/orog.nc'
    ua_file_mask = wdir + '/ua_fill.nc'
    va_file_mask = wdir + '/va_fill.nc'
    energy3_file = wdir + '/energy_short.nc'
    cdo.setmisstoc('0',
                   input='-setmisstoc,1 -sub {0} {0}'.format(ua_file),
                   options='-b F32',
                   output=maskorog)
    cdo.add(input=('-setmisstoc,0 -selvar,ua {} '
                   '-setmisstoc,0 -mul {} -selvar,ua {}').format(
                       ua_file, uas_file, maskorog),
            options='-b F32',
            output=ua_file_mask)
    cdo.add(input=('-setmisstoc,0 -selvar,va {} '
                   '-setmisstoc,0 -mul {} -selvar,ua {}').format(
                       va_file, vas_file, maskorog),
            options='-b F32',
            output=va_file_mask)
    cdo.setmisstoc('0',
                   input=('-invertlat -sellevel,10000/90000 '
                          '-merge {} {} {} {}').format(ta_file, ua_file_mask,
                                                       va_file_mask, wap_file),
                   options='-b F32',
                   output=energy3_file)
    yrs = cdo.showyear(input=energy3_file)
    yrs = str(yrs)
    yrs2 = yrs.split()
    y_i = 0
    lect = np.zeros(len(yrs2))
    for y_r in yrs2:
        y_rl = [y_n for y_n in y_r]
        y_ro = ''
        for e_l in y_rl:
            e_l = str(e_l)
            if e_l.isdigit() is True:
                y_ro += e_l
        # print(filter(str.isdigit, str(y_r)))
        enfile_yr = wdir + '/inputen.nc'
        tasfile_yr = wdir + '/tas_yr.nc'
        tadiag_file = wdir + '/ta_filled.nc'
        ncfile = wdir + '/fourier_coeff.nc'
        cdo.selyear(y_ro,
                    input=energy3_file,
                    options='-b F32',
                    output=enfile_yr)
        cdo.selyear(y_ro, input=tas_file, options='-b F32', output=tasfile_yr)
        fourc.fourier_coeff(tadiag_file, ncfile, enfile_yr, tasfile_yr)
        diagfile = (ldir + '/{}_{}_lec_diagram.png'.format(model, y_ro))
        logfile = (ldir + '/{}_{}_lec_table.txt'.format(model, y_ro))
        lect[y_i] = lorenz(wdir, model, y_ro, ncfile, diagfile, logfile)
        y_i = y_i + 1
        os.remove(enfile_yr)
        os.remove(tasfile_yr)
        os.remove(tadiag_file)
        os.remove(ncfile)
    os.remove(maskorog)
    os.remove(ua_file_mask)
    os.remove(va_file_mask)
    os.remove(energy3_file)
    return lect


def removeif(filename):
    """Remove filename if it exists."""
    try:
        os.remove(filename)
    except OSError:
        pass


def stabil(ta_gmn, p_l, nlev):
    """Compute the stability parameter from temp. and pressure levels.

    Arguments
    - ta_gmn: a temperature vertical profile;
    - p_l: the vertical levels;
    - nlev: the number of vertical levels;
    """
    cpdr = CP / R
    t_g = ta_gmn
    g_s = np.zeros(nlev)
    for i_l in range(nlev):
        if i_l == 0:
            dtdp = (t_g[i_l + 1] - t_g[i_l]) / (p_l[i_l + 1] - p_l[i_l])
        elif i_l == nlev - 1:
            dtdp = (t_g[i_l] - t_g[i_l - 1]) / (p_l[i_l] - p_l[i_l - 1])
        else:
            dtdp1 = (t_g[i_l + 1] - t_g[i_l]) / (p_l[i_l + 1] - p_l[i_l])
            dtdp2 = (t_g[i_l] - t_g[i_l - 1]) / (p_l[i_l] - p_l[i_l - 1])
            dtdp = ((dtdp1 * (p_l[i_l] - p_l[i_l - 1]) + dtdp2 *
                     (p_l[i_l + 1] - p_l[i_l])) /
                    (p_l[i_l + 1] - p_l[i_l - 1]))
        g_s[i_l] = CP / (t_g[i_l] - p_l[i_l] * dtdp * cpdr)
    return g_s


def table(varin, ntp, name, logfile, flag):
    """Write global and hem. storage terms to .txt table.

    Arguments:
    - varin: the variable to be printed out;
    - ntp: the number of wavenumbers;
    - name: the name of the variable to be printed out;
    - logfile: the filename of the .txt where the variable is printed out;
    - flag: a flag for NH, SH, global;
    """
    if flag is True:
        fac = 1e5
        varin = fac * varin
    varzon = varin[:, 0]
    vared = np.nansum(varin[:, 1:ntp - 1], axis=1)
    vared1 = np.nansum(varin[:, 1:NW_1 - 1], axis=1)
    vared2 = np.nansum(varin[:, NW_1:NW_2 - 1], axis=1)
    vared3 = np.nansum(varin[:, NW_2:NW_3 - 1], axis=1)
    vared_tog = [vared, vared1, vared2, vared3]
    write_to_tab(logfile, name, vared_tog, varzon)


def varatts(w_nc_var, varname, tres, vres):
    """Add attibutes to the variables, depending on name and time res.

    Arguments:
    - w_nc_var: a variable object;
    - varname: the name of the variable, among ta, ua, va and wap;
    - tres: the time resolution (daily or annual);
    - vres: the vertical resolution (pressure levels or vert. integr.).

    @author: Chris Slocum (2014), modified by Valerio Lembo (2018).
    """
    if tres == 0:
        tatt = "Daily\nM"
    elif tres == 1:
        tatt = "Annual mean\nM"
    if vres == 0:
        vatt = "Pressure levels\n"
    elif vres == 1:
        vatt = "Vertically integrated\n"
    if varname == 'a':
        w_nc_var.setncatts({
            'long_name': "Available Potential Energy",
            'units': "W m-2",
            'level_desc': vatt,
            'var_desc': "APE -> KE",
            'statistic': tatt
        })
    elif varname == 'ek':
        w_nc_var.setncatts({
            'long_name': "Kinetic Energy",
            'units': "W m-2",
            'level_desc': vatt,
            'var_desc': "APE -> KE",
            'statistic': tatt
        })
    elif varname == 'a2k':
        w_nc_var.setncatts({
            'long_name': "Conversion between APE and KE",
            'units': "W m-2",
            'level_desc': vatt,
            'var_desc': "APE <-> KE",
            'statistic': tatt
        })
    elif varname == 'k':
        w_nc_var.setncatts({
            'long_name': "Kinetic Energy",
            'units': "W m-2",
            'level_desc': vatt,
            'var_desc': "APE -> KE",
            'statistic': tatt
        })


def weights(lev, nlev, lat):
    """Compute weigths for vertical integration and meridional averages.

    Arguments:
    - lev: the pressure levels;
    - nlev: the number of pressure levels;
    - lat: the latitudes in degrees;
    - nlat: the number of latitudinal gridsteps;
    """
    # Compute sigma level and dsigma
    sig = lev / PS
    d_s = np.zeros(nlev)
    for j_l in range(1, nlev - 1, 1):
        d_s[j_l] = 0.5 * abs(sig[j_l + 1] - sig[j_l - 1])
    d_s[0] = sig[0] + 0.5 * abs(sig[1] - sig[0])
    d_s[nlev -
        1] = 1 - sig[nlev - 1] + 0.5 * abs(sig[nlev - 1] - sig[nlev - 2])
    # Compute Gaussian weights
    y_l = np.zeros(lat.shape)
    np.deg2rad(lat, out=y_l)
    g_w = np.cos(y_l)
    return d_s, y_l, g_w


def write_to_tab(logfile, name, vared, varzon):
    """Specify the formats for table entries.

    Arguments:
    - log: the logfile where the entries must be written;
    - name: the name of the variable;
    - vared: a list of arrays containing the overall eddy components, the LW,
      the SW and the KW components;
    - varzon: an array containing the zonal mean component;
    """
    vartot = varzon + vared[0]
    with open(logfile, 'w') as log:
        log.write(' {} TOTAL    {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(
            name, vartot[0], vartot[1], vartot[2]))
        log.write('--------------------------------------\n')
        log.write(' {} ZONAL    {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(
            name, varzon[0], varzon[1], varzon[2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY     {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(
            name, vared[0][0], vared[0][1], vared[0][2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY(LW) {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(
            name, vared[1][0], vared[1][1], vared[1][2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY(SW) {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(
            name, vared[2][0], vared[2][1], vared[2][2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY(KW) {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(
            name, vared[3][0], vared[3][1], vared[3][2]))
        log.write('--------------------------------------\n')
