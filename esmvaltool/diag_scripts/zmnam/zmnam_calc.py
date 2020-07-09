"""
Zonal-mean annular mode calculation routine.

Author: Federico Serva (ISAC-CNR & ISMAR-CNR, Italy)
Copernicus C3S 34a lot 2 (MAGIC)
"""

import netCDF4
import numpy as np
from scipy import signal


def butter_filter(data, freq, lowcut=None, order=2):
    """Function to perform time filtering."""
    if lowcut is not None:
        filttype = 'lowpass'

    # Sampling determines Nyquist frequency
    nyq = 0.5 * freq

    if filttype == 'lowpass':
        low = lowcut / nyq
        freqs = low

    bpoly, apoly = signal.butter(order, freqs, btype=filttype)
    #    _, h = signal.freqs(b, a)
    ysig = signal.filtfilt(bpoly, apoly, data, axis=0)

    return ysig


def zmnam_calc(da_fname, outdir, src_props):
    """Function to do EOF/PC decomposition of zg field."""
    deg_to_r = np.pi / 180.
    lat_weighting = True
    outfiles = []

    # Note: daily/monthly means have been
    # already subtracted from daily/monthly files

    # Open daily data

    with netCDF4.Dataset(da_fname, "r") as in_file:
        time_dim = in_file.variables['time'][:]
        time_lnam = getattr(in_file.variables['time'], 'long_name', '')
        time_snam = getattr(in_file.variables['time'], 'standard_name', '')
        time_uni = in_file.variables['time'].units
        time_cal = in_file.variables['time'].calendar
        time = np.array(time_dim[:], dtype='d')
        date = netCDF4.num2date(time, in_file.variables['time'].units,
                                in_file.variables['time'].calendar)

        lev = np.array(in_file.variables['plev'][:], dtype='d')
        lev_lnam = getattr(in_file.variables['plev'], 'long_name', '')
        lev_snam = getattr(in_file.variables['plev'], 'standard_name', '')
        lev_uni = in_file.variables['plev'].units
        lev_pos = in_file.variables['plev'].positive
        lev_axi = in_file.variables['plev'].axis

        lat = np.array(in_file.variables['lat'][:], dtype='d')
        # lat_nam = in_file.variables['lat'].long_name
        lat_uni = in_file.variables['lat'].units
        lat_axi = in_file.variables['lat'].axis

        lon = np.array(in_file.variables['lon'][:], dtype='d')
        # lon_nam = in_file.variables['lon'].long_name
        lon_uni = in_file.variables['lon'].units
        lon_axi = in_file.variables['lon'].axis

        zg_da = np.squeeze(np.array(in_file.variables['zg'][:], dtype='d'))

    n_tim = len(time_dim)
    print('end infile close')

    # Start zmNAM index calculation

    # Lowpass filter
    zg_da_lp = butter_filter(zg_da, 1, lowcut=1. / 90, order=2)

    # Outputs: stored by level
    # EOFs, eigenvalues, daily and monthly PCs
    eofs = np.zeros((len(lev), len(lat)), dtype='d')
    eigs = np.zeros(len(lev), dtype='d')
    pcs_da = np.zeros((n_tim, len(lev)), dtype='d')

    # Calendar-independent monthly mean
    sta_mon = []  # first day of the month
    mid_mon = []  # 15th of the month
    end_mon = []  # last day of the month (add +1 when slicing)

    mon = 999
    idate = 0

    while idate < len(date):

        # Save first day of the month
        if date[idate].month != mon:
            mon = date[idate].month
            sta_mon.append(idate)

        # Save month mid-day
        if date[idate].day == 15:
            mid_mon.append(idate)

        # Save last day of the month
        if idate == len(date) - 1 or (date[idate].month == mon and
                                      date[idate + 1].month != mon):
            end_mon.append(idate)

        idate += 1

    pcs_mo = np.zeros((len(date[mid_mon]), len(lev)), dtype='d')

    # Perform analysis by level
    for i_lev in np.arange(len(lev)):
        # Latitude weighting
        if lat_weighting is True:
            for j_lat in np.arange(len(lat)):
                zg_da_lp[:, i_lev,
                         j_lat] *= np.sqrt(abs(np.cos(lat[j_lat] * deg_to_r)))

        zg_da_lp_an = zg_da_lp[:, i_lev, :] - np.mean(zg_da_lp[:, i_lev, :],
                                                      axis=0)
        cov = np.dot(zg_da_lp_an.T, zg_da_lp_an) / (n_tim - 1)

        # Compute eigenvectors and eigenvalues
        eigenval, eigenvec = np.linalg.eig(cov)

        sum_eigenval = np.sum(eigenval)

        eigenval_norm = eigenval[:] / sum_eigenval

        # Largest eigenvalue
        max_eigenval = eigenval_norm.argmax()

        # PC calculation
        pc = np.dot(zg_da_lp_an[:, :], eigenvec)

        # Latitude de-weighting
        if lat_weighting is True:
            for i_lat in np.arange(len(lat)):
                eigenvec[i_lat, :] /= np.sqrt(
                    abs(np.cos(lat[i_lat] * deg_to_r)))

        # Retain leading standardized PC & EOF
        lead_pc_mean = np.mean(pc[:, max_eigenval], axis=0)
        lead_pc_std = np.std(pc[:, max_eigenval], ddof=1, axis=0)
        lead_pc = (pc[:, max_eigenval] - lead_pc_mean) / lead_pc_std
        lead_eof = eigenvec[:, max_eigenval]

        max_lat = max(range(len(lat)), key=lambda x: lat[x])
        min_lat = min(range(len(lat)), key=lambda x: lat[x])

        if lead_eof[max_lat] > lead_eof[min_lat]:
            lead_pc *= -1
            lead_eof *= -1

        lead_pc_mo = np.zeros(len(date[mid_mon]), dtype='d')
        time_mo = np.zeros(len(date[mid_mon]), dtype='d')

        # print(lead_pc)

        for k_mo in range(len(date[mid_mon])):
            lead_pc_mo[k_mo] = np.mean(lead_pc[sta_mon[k_mo]:end_mon[k_mo] +
                                               1])
            time_mo[k_mo] = time[mid_mon[k_mo]]

        # Store PC/EOF for this level (no time dependent)
        eigs[i_lev] = eigenval_norm[max_eigenval]
        eofs[i_lev, :] = lead_eof
        pcs_da[:, i_lev] = lead_pc
        pcs_mo[:, i_lev] = lead_pc_mo

    # Save output files

    # (1) daily PCs
    fname = outdir + '_'.join(src_props) + '_pc_da.nc'
    outfiles.append(fname)
    with netCDF4.Dataset(fname, mode='w') as file_out:
        file_out.title = 'Zonal mean annular mode (1)'
        file_out.contact = 'F. Serva (federico.serva@artov.ismar.cnr.it); \
                            C. Cagnazzo (chiara.cagnazzo@cnr.it)'

        file_out.createDimension('time', None)
        file_out.createDimension('plev', np.size(lev))
        file_out.createDimension('lat', np.size(lat))
        file_out.createDimension('lon', np.size(lon))

        time_var = file_out.createVariable('time', 'd', ('time', ))
        if time_lnam:
            time_var.setncattr('long_name', time_lnam)
        if time_snam:
            time_var.setncattr('standard_name', time_snam)
        time_var.setncattr('units', time_uni)
        time_var.setncattr('calendar', time_cal)
        time_var[:] = time_dim[:]

        lev_var = file_out.createVariable('plev', 'd', ('plev', ))
        if lev_lnam:
            lev_var.setncattr('long_name', lev_lnam)
        if lev_snam:
            lev_var.setncattr('standard_name', lev_snam)
        lev_var.setncattr('units', lev_uni)
        lev_var.setncattr('positive', lev_pos)
        lev_var.setncattr('axis', lev_axi)
        lev_var[:] = lev[:]

        pcs_da_var = file_out.createVariable('PC_da', 'd', (
            'time',
            'plev',
        ))
        pcs_da_var.setncattr('long_name', 'Daily zonal mean annular mode PC')
        pcs_da_var.setncattr(
            'comment',
            'Reference: Baldwin and Thompson (2009), doi:10.1002/qj.479')
        pcs_da_var[:] = pcs_da[:, :]

    # (2) monthly PCs
    fname = outdir + '_'.join(src_props) + '_pc_mo.nc'
    outfiles.append(fname)
    with netCDF4.Dataset(fname, mode='w') as file_out:
        file_out.title = 'Zonal mean annular mode (2)'
        file_out.contact = 'F. Serva (federico.serva@artov.ismar.cnr.it); \
        C. Cagnazzo (chiara.cagnazzo@cnr.it)'

        file_out.createDimension('time', None)
        file_out.createDimension('plev', np.size(lev))

        time_var = file_out.createVariable('time', 'd', ('time', ))
        if time_lnam:
            time_var.setncattr('long_name', time_lnam)
        if time_snam:
            time_var.setncattr('standard_name', time_snam)
        time_var.setncattr('units', time_uni)
        time_var.setncattr('calendar', time_cal)
        time_var[:] = time_mo

        lev_var = file_out.createVariable('plev', 'd', ('plev', ))
        if lev_lnam:
            lev_var.setncattr('long_name', lev_lnam)
        if lev_snam:
            lev_var.setncattr('standard_name', lev_snam)
        lev_var.setncattr('units', lev_uni)
        lev_var.setncattr('positive', lev_pos)
        lev_var.setncattr('axis', lev_axi)
        lev_var[:] = lev[:]

        pcs_mo_var = file_out.createVariable('PC_mo', 'd', (
            'time',
            'plev',
        ))
        pcs_mo_var.setncattr('long_name', 'Monthly zonal mean annular mode PC')
        pcs_mo_var.setncattr(
            'comment',
            'Reference: Baldwin and Thompson (2009), doi:10.1002/qj.479')
        pcs_mo_var[:] = pcs_mo[:, :]

    # (3) EOFs and explained varianceo
    fname = outdir + '_'.join(src_props) + '_eofs.nc'
    outfiles.append(fname)
    with netCDF4.Dataset(fname, mode='w') as file_out:
        file_out.title = 'Zonal mean annular mode (3)'
        file_out.contact = 'F. Serva (federico.serva@artov.ismar.cnr.it); \
        C. Cagnazzo (chiara.cagnazzo@cnr.it)'

        file_out.createDimension('time', None)
        file_out.createDimension('plev', np.size(lev))
        file_out.createDimension('lat', np.size(lat))
        file_out.createDimension('lon', np.size(lon))

        time_var = file_out.createVariable('time', 'd', ('time', ))
        if time_lnam:
            time_var.setncattr('long_name', time_lnam)
        if time_snam:
            time_var.setncattr('standard_name', time_snam)
        time_var.setncattr('units', time_uni)
        time_var.setncattr('calendar', time_cal)
        time_var[:] = 0
        #
        lev_var = file_out.createVariable('plev', 'd', ('plev', ))
        if lev_lnam:
            lev_var.setncattr('long_name', lev_lnam)
        if lev_snam:
            lev_var.setncattr('standard_name', lev_snam)
        lev_var.setncattr('units', lev_uni)
        lev_var.setncattr('positive', lev_pos)
        lev_var.setncattr('axis', lev_axi)
        lev_var[:] = lev[:]
        #
        lat_var = file_out.createVariable('lat', 'd', ('lat', ))
        lat_var.setncattr('units', lat_uni)
        lev_var.setncattr('axis', lat_axi)
        lat_var[:] = lat[:]
        #
        lon_var = file_out.createVariable('lon', 'd', ('lon', ))
        lon_var.setncattr('units', lon_uni)
        lon_var.setncattr('axis', lon_axi)
        lon_var[:] = lon[:]
        #
        eofs_var = file_out.createVariable('EOF', 'd', ('plev', 'lat'))
        eofs_var.setncattr('long_name', 'Zonal mean annular mode EOF')
        eofs_var.setncattr(
            'comment',
            'Reference: Baldwin and Thompson (2009), doi:10.1002/qj.479')
        eofs_var[:] = eofs[:, :]
        #
        eigs_var = file_out.createVariable('eigenvalues', 'd', ('plev'))
        eigs_var.setncattr('long_name',
                           'Zonal mean annular mode EOF explained variance')
        eigs_var.setncattr(
            'comment',
            'Reference: Baldwin and Thompson (2009), doi:10.1002/qj.479')
        eigs_var[:] = eigs[:]

    return outfiles
