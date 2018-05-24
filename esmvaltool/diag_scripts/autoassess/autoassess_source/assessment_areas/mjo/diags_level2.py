'''
Converting NCL diagnostics_cam.ncl
Functions include

1. decompose2SymAsym
2. rmvAnnualCycle
3. resolveWavesHayashi
4. genDispersionCurves
5. addHorVertLines
6. replaceChars
7. statAsymSym
8. wkSpaceTime
9. mjo_spectra_season
10. mjo_spectra
11. mjo_xcor_lag_season
12. mjo_wavenum_freq_season
'''
import numpy as np
import iris
import sys
import math
import os
import mjo_utils as mu
import mjo_plots as mjp

from scipy import signal


def _makecube(var, wave, freq):
    #===========================================================================
    # # Make a 2D cube of wavenumber & frequency dimensions
    #===========================================================================
    var_cube = iris.cube.Cube(var)
    var_cube.rename('spectra')
    wave_coord = iris.coords.DimCoord(wave, long_name='wavenumber')
    wave_coord.guess_bounds()
    freq_coord = iris.coords.DimCoord(freq, long_name='frequency')
    freq_coord.guess_bounds()
    var_cube.add_dim_coord(freq_coord, 0)
    var_cube.add_dim_coord(wave_coord, 1)
    return var_cube


def _makespeccube(var, lat, wave, freq):
    #===========================================================================
    # # Make a 3D cube of Latitude, wavenumber & frequency dimensions
    #===========================================================================
    var_cube = iris.cube.Cube(var)
    var_cube.rename('spectra')
    lat_coord = iris.coords.DimCoord(lat, long_name='latitude')
    lat_coord.guess_bounds()
    wave_coord = iris.coords.DimCoord(wave, long_name='wavenumber')
    wave_coord.guess_bounds()
    freq_coord = iris.coords.DimCoord(freq, long_name='frequency')
    freq_coord.guess_bounds()
    var_cube.add_dim_coord(lat_coord, 0)
    var_cube.add_dim_coord(freq_coord, 1)
    var_cube.add_dim_coord(wave_coord, 2)
    return var_cube


def _makecube_season_pow(var, wave, freq, name='spectra'):
    #===========================================================================
    # Make a cube of seasonal power
    #===========================================================================
    var_cube = iris.cube.Cube(var)
    var_cube.rename('spectra')
    #var_cube.long_name = long_name
    wave_coord = iris.coords.DimCoord(wave, long_name='wavenumber')
    wave_coord.guess_bounds()
    freq_coord = iris.coords.DimCoord(freq, long_name='frequency')
    freq_coord.guess_bounds()
    var_cube.add_dim_coord(wave_coord, 0)
    var_cube.add_dim_coord(freq_coord, 1)
    return var_cube


def _taper(ts, alpha=0.1, iopt=0):
    #===========================================================================
    #
    # this function applies split-cosine-bell tapering to the series x.
    # the series will be tapered to the mean of x.
    # See:
    #     Fourier Analysis of Time Series
    #     Peter Bloomfield
    #     Wiley-Interscience, 1976
    #  This is used prior performing an fft on non-cyclic data
    #  arguments:
    #     ts   series to be tapered (tapering done in place)
    #         **missing data not allowed**
    #     alpha   the proportion of the time series to be tapered
    #         [p=0.10 means 10 %]
    #     iopt  iopt=0 taper to series mean
    #     iopt  iopt=1 means *force* taper to 0.0
    #===========================================================================
    if all(x == ts[0] for x in ts):
        print 'all values equal'
        iopt = 1
    if iopt == 0:
        tsmean = np.mean(ts)
    else:
        tsmean = 0.
    n = len(ts)
    m = max(1, int(alpha * n + 0.5) / 2)
    pim = np.pi / m

    tst = ts.copy()
    for i in range(1, m + 1):
        weight = 0.5 - 0.5 * np.cos(pim * (i - 0.5))
        tst[i - 1] = (ts[i - 1] - tsmean) * weight + tsmean
        tst[n - i] = (ts[n - i] - tsmean) * weight + tsmean
    return tst


def _rmvAnnualCycle(var, nDayTot, fCrit, spd=1, rmvMeans=False):
    #===========================================================================
    # Prewhiten the data: eg remove the annual cycle.
    # Actually, this will remove all time periods less than
    #   those corresponding to 'fCrit'.
    # Note: The original fortran code provided by JET did not remove
    #   the grid point means so .... rmvMeans=False
    #       I assume that Matt Wheeler's code did that also.
    #===========================================================================
    print fCrit
    # Take a copy for metadata
    var_cut = var.copy()

    # mean for later
    varMean = var.collapsed('time', iris.analysis.MEAN)

    ntime, nlat, nlon = var.data.shape

    # compute FFT
    cf = np.fft.fft(var.data, axis=0)

    freq = np.fft.fftfreq(ntime)
    x = freq.copy()

    # cutting frequencies
    cf[np.abs(x) < fCrit] = 0.

    # inverse FFT
    var_cut.data = np.fft.ifft(cf, axis=0).astype(float)
    if not rmvMeans:
        var_cut.data += varMean.data
    return var_cut


def _decompose2SymAsym(var, axis=1):
    #===========================================================================
    # This code decomposes the data in to symmetric and anti-symmetric components
    # with respect to latitude axis. It is assumed that the latitude dimension is
    # along axis = 1 as default.
    #
    # antisymmetric part is stored in one hemisphere [eg: Northern Hemisphere]
    #          xOut( lat) = (x(lat)-x(-lat))/2
    # symmetric part is stored in other hemisphere [eg: Southern Hemisphere]
    #          xOut(-lat) = (x(lat)+x(-lat))/2
    #===========================================================================
    varSA = var.copy()  # copy to store the output
    nlat = var.shape[axis]  # get the number of latitude points

    N2 = nlat / 2
    N = N2
    if nlat % 2 == 1:
        N = N2 + 1  # offset to handle the Equator

    if axis == 1:
        for nl in np.arange(0, N2):
            print (nlat - 1 - nl), nl
            varSA.data[:, nl] = 0.5 * (var.data[:, nlat - 1 - nl] + var.data[:, nl])
            varSA.data[:, nlat - 1 - nl] = 0.5 * (var.data[:, nlat - 1 - nl] - var.data[:, nl])

        return varSA
    else:
        print 'Modify the code to accommodate other axes...'
        print 'Exiting...'
        sys.exit(0)


def resolveWavesHayashi(varfft, nDayWin, spd=1):
    #===========================================================================
    #
    #  Create array PEE(NL+1,NT+1) which contains the (real) power spectrum.
    #  all the following assume indexing starting with 0
    #  In this array, the negative wavenumbers will be from pn=0 to NL/2-1
    #  The positive wavenumbers will be for pn=NL/2+1 to NL.
    #  Negative frequencies will be from pt=0 to NT/2-1
    #  Positive frequencies will be from pt=NT/2+1 to NT  .
    #  Information about zonal mean will be for pn=NL/2  .
    #  Information about time mean will be for pt=NT/2  .
    #  Information about the Nyquist Frequency is at pt=0 and pt=NT
    #
    #  In PEE, define the
    #  WESTWARD waves to be either +ve frequency
    #           and -ve wavenumber or -ve freq and +ve wavenumber.
    #  EASTWARD waves are either +ve freq and +ve wavenumber
    #           OR -ve freq and -ve wavenumber.
    #
    #  Note that frequencies are returned from fftpack are ordered like so
    #     input_time_pos [ 0    1   2    3     4      5    6   7  ]
    #     ouput_fft_coef [mean 1/7 2/7  3/7 nyquist -3/7 -2/7 -1/7]
    #                     mean,pos freq to nyq,neg freq hi to lo
    #
    #  Rearrange the coef array to give you power array of freq and wave number east/west
    #  Note east/west wave number *NOT* eq to fft wavenumber see Hayashi '71
    #  Hence, NCL's 'cfftf_frq_reorder' can *not* be used.
    #
    #  For ffts that return the coefficients as described above, here is the algorithm
    #  coeff array varfft(2,n,t)   dimensioned (2,0:numlon-1,0:numtim-1)
    #  new space/time pee(2,pn,pt) dimensioned (2,0:numlon  ,0:numtim  )
    #
    #  NOTE: one larger in both freq/space dims
    #  the initial index of 2 is for the real (indx 0) and imag (indx 1) parts of the array
    #
    #
    #       if  |  0 <= pn <= numlon/2-1    then    | numlon/2 <= n <= 1
    #           |  0 <= pt < numtim/2-1             | numtim/2 <= t <= numtim-1
    #
    #       if  |  0         <= pn <= numlon/2-1    then    | numlon/2 <= n <= 1
    #           |  numtime/2 <= pt <= numtim                | 0        <= t <= numtim/2
    #
    #       if  |  numlon/2  <= pn <= numlon    then    | 0  <= n <= numlon/2
    #           |  0         <= pt <= numtim/2          | numtim/2 <= t <= 0
    #
    #       if  |  numlon/2   <= pn <= numlon    then    | 0        <= n <= numlon/2
    #           |  numtim/2+1 <= pt <= numtim            | numtim-1 <= t <= numtim/2
    #
    #===========================================================================
    N, mlon = varfft.shape
    pee = np.ones([N + 1, mlon + 1]) * -999.  # initialize
                                            # -999 scaling is for testing
                                            # putposes
    # Create the real power spectrum pee = sqrt(real^2+imag^2)^2
    varfft = np.abs(varfft) ** 2
    pee[:N / 2, :mlon / 2] = varfft[N / 2:N, mlon / 2:0:-1]
    pee[N / 2 :, :mlon / 2] = varfft[:N / 2 + 1, mlon / 2:0:-1]
    pee[:N / 2 + 1, mlon / 2:] = varfft[N / 2::-1, :mlon / 2 + 1]
    pee[N / 2 + 1:, mlon / 2:] = varfft[N - 1:N / 2 - 1:-1, :mlon / 2 + 1]

    return pee


def wk_smooth121(var):
    #===========================================================================
    # Special 1-2-1 smoother
    # Smooths vv by passing it through a 1-2-1 filter.
    # The first and last points are given 3-1 (1st) or 1-3 (last)
    # weightings (Note that this conserves the total sum).
    # The routine also skips-over missing data (np.nan)
    #===========================================================================
    varf = var.copy()
    nt = len(var)

    for i in range(0, nt, 1):
        if i == 0:
            if not np.isnan(var[i + 1]):
                varf[i] = (3.0 * var[i] + var[i + 1]) / 4.0
        elif np.isnan(var[i - 1]):
            if not np.isnan(var[i + 1]):
                varf[i] = (3.0 * var[i] + var[i + 1]) / 4.0
        elif (i == nt - 1) or (np.isnan(var[i + 1])):
            if not np.isnan(var[i - 1]):
                varf[i] = (var[i - 1] + 3.0 * var[i]) / 4.0
        else:
            varf[i] = (1.0 * var[i - 1] + 2.0 * var[i] + 1.0 * var[i + 1]) / 4.0
    return varf


def ClosestIndex(array, value):
    #===========================================================================
    # Find the closest index to a given value in an array
    #===========================================================================
    idx = (np.abs(array - value)).argmin()
    return idx


def ComputeBackground(peeAS, wave, freq, minwav4smth, maxwav4smth):
    #===========================================================================
    # Derive the background spectrum (red noise) ************
    # [1] Sum power over all latitude
    # [2] Put fill value in mean
    # [3] Apply smoothing to the spectrum. This smoothing DOES include wavenumber zero.
    #
    #===========================================================================
    psumb = np.sum(peeAS, axis=0)  #  sum over all latitudes
    N, mlon = psumb.shape
    smthlen = maxwav4smth - minwav4smth + 1

    for tt in range(N / 2 + 1, N):
        if freq[tt] < 0.1 :
            for i in range(1, 6):
                psumb[tt, minwav4smth:maxwav4smth + 1] = \
                wk_smooth121(psumb[tt, minwav4smth:maxwav4smth + 1])
        if freq[tt] >= 0.1 and freq[tt] < 0.2 :
            for i in range(1, 11):
                psumb[tt, minwav4smth:maxwav4smth + 1] = \
                wk_smooth121(psumb[tt, minwav4smth:maxwav4smth + 1])
        if freq[tt] >= 0.2 and freq[tt] < 0.3 :
            for i in range(1, 21):
                psumb[tt, minwav4smth:maxwav4smth + 1] = \
                wk_smooth121(psumb[tt, minwav4smth:maxwav4smth + 1])
        if freq[tt] >= 0.3 :
            for i in range(1, 41):
                psumb[tt, minwav4smth:maxwav4smth + 1] = \
                wk_smooth121(psumb[tt, minwav4smth:maxwav4smth + 1])

    pt8cpd = min([ClosestIndex(freq, 0.8), len(freq) - 1])

    # smth frequency up to .8 cycles per day
    for nw in range(minwav4smth, maxwav4smth + 1):
        smthlen = pt8cpd - (N / 2 + 1) + 1
        for i in range(1, 11):
            psumb[N / 2 + 1:pt8cpd + 1, nw] = \
            wk_smooth121(psumb[N / 2 + 1:pt8cpd + 1, nw])
    return psumb


def genDispersionCurves(nWaveType, nEquivDepth, nPlanetaryWave, rlat, Ahe, fillval):
    #---------------------------------------------------------------
    # Theoretical shallow water dispersion curves
    #--------------------------------------------------------------
    pi = 4.0 * math.atan(1.0)
    re = 6.37122e06     # [m]   average radius of earth
    g = 9.80665        # [m/s] gravity at 45 deg lat used by the WMO
    omega = 7.292e-05      # [1/s] earth's angular vel
    U = 0.0
    Un = 0.0   # since Un = U*T/L
    ll = 2.*pi * re * math.cos(abs(rlat))
    Beta = 2.*omega * math.cos(abs(rlat)) / re
    maxwn = nPlanetaryWave

    Apzwn = np.zeros([nWaveType, nEquivDepth, nPlanetaryWave], dtype=np.double)
    Afreq = np.zeros([nWaveType, nEquivDepth, nPlanetaryWave], dtype=np.double)

    for ww in range(1, nWaveType + 1):    #  wave type
        for ed in range(1, nEquivDepth + 1):   #  equivalent depth
            he = Ahe[ed - 1]
            T = 1. / math.sqrt(Beta) * (g * he) ** (0.25)
            L = (g * he) ** (0.25) / math.sqrt(Beta)

            for wn in range(1, nPlanetaryWave + 1):     # planetary wave number
                s = -20.*(wn - 1) * 2. / (nPlanetaryWave - 1) + 20.
                k = 2.*pi * s / ll
                kn = k * L

                # Anti-symmetric curves
                if ww == 1 :  # MRG wave
                    if k <= 0:
                        delx = math.sqrt(1. + (4.*Beta) / (k ** 2 * math.sqrt(g * he)))
                        deif = k * math.sqrt(g * he) * (0.5 - 0.5 * delx)
                    if k == 0:
                        deif = math.sqrt(math.sqrt(g * he) * Beta)
                    if k > 0 :
                        deif = fillval

                if ww == 2:     # n=0 IG wave
                    if k < 0:
                        deif = fillval
                    if k == 0:
                        deif = math.sqrt(math.sqrt(g * he) * Beta)
                    if k > 0:
                        delx = math.sqrt(1. + (4.0 * Beta) / (k ** 2 * math.sqrt(g * he)))
                        deif = k * math.sqrt(g * he) * (0.5 + 0.5 * delx)

                if ww == 3:      # n=2 IG wave
                    n = 2.
                    delx = (Beta * math.sqrt(g * he))
                    deif = math.sqrt((2.*n + 1.) * delx + (g * he) * k ** 2)
                    # do some corrections to the above calculated frequency.......
                    for i in range(1, 6):
                        deif = math.sqrt((2.*n + 1.) * delx + (g * he) * k ** 2 + g * he * Beta * k / deif)

                # symmetric curves
                if ww == 4:              # n=1 ER wave
                    n = 1.
                    if k < 0:
                        delx = (Beta / math.sqrt(g * he)) * (2.*n + 1.)
                        deif = -Beta * k / (k ** 2 + delx)
                    else:
                        deif = fillval
                if ww == 5:              # Kelvin wave
                    deif = k * math.sqrt(g * he)
                if ww == 6:             #  n=1 IG wave
                    n = 1.
                    delx = (Beta * math.sqrt(g * he))
                    deif = math.sqrt((2.*n + 1.) * delx + (g * he) * k ** 2)
                    # do some corrections to the above calculated frequency.......
                    for i in range(1, 6):
                        deif = math.sqrt((2.*n + 1.) * delx + (g * he) * k ** 2 + g * he * Beta * k / deif)

                eif = deif  # + k*U since  U=0.0
                P = 2.*pi / (eif * 24.*60.*60.)
                dps = deif / k
                R = L
                Rdeg = (180.*R) / (pi * 6.37e6)
                Apzwn[ww - 1, ed - 1, wn - 1] = s
                if deif != fillval:
                    P = 2.*pi / (eif * 24.*60.*60.)
                    Afreq[ww - 1, ed - 1, wn - 1] = 1. / P
                else:
                    Afreq[ww - 1, ed - 1, wn - 1] = fillval

    return Afreq, Apzwn


def wkSpaceTime(var, outdir, runid, latBound=15, spd=1, nDayWin=96, nDaySkip=15):
    """Create Wheeler-Kiladis Space-Time  plots.

     Note_1: The full logitudinal domain is used.
             This means that every planetary
             wavenumber will be represented.
     Note_2: Tapering in time is done to make the variable periodic.

     The calculations are also only made for the latitudes
     between '-latBound' and 'latBound'.

    ********************   REFERENCES  *******************************
     Wheeler, M., G.N. Kiladis Convectively Coupled Equatorial Waves:
        Analysis of Clouds and Temperature in the Wavenumber-Frequency
        Domain J. Atmos. Sci., 1999,  56: 374-399.
    ---
     Hayashi, Y. A Generalized Method of Resolving Disturbances into
        Progressive and Retrogressive Waves by Space and Fourier and
        TimeCross Spectral Analysis J. Meteor. Soc. Japan, 1971, 49: 125-128.
    """

    with iris.FUTURE.context(cell_datetime_objects=True):
        start_date = var.coord('time').cell(0).point.strftime('%Y-%m-%d')
        end_date = var.coord('time').cell(-1).point.strftime('%Y-%m-%d')

    varname = var.name()
    if varname == 'x_wind':
        assert len(var.coord('pressure').points) == 1
        pressure_level = var.coord('pressure').points[0]
        if pressure_level == 850:
            varname = 'x_wind_850hPa'
        if pressure_level == 200:
            varname = 'x_wind_200hPa'

    ntim, nlat, mlon = var.shape
    latN = latBound
    latS = -latBound  # make symmetric about the equator

    lonL = 0  # -180
    lonR = 360  #  180
    fCrit = 1. / nDayWin  # remove all contributions 'longer'

    tim_taper = 0.1  # time taper      [0.1   => 10%]
    lon_taper = 0.0  # longitude taper [0.0 for globe  only global supported]

    if lon_taper > 0.0 or lonR - lonL != 360.:
        print 'Code does currently allow lon_taper>0 or (lonR-lonL)<360'
        sys.exit(0)

    nDayTot = ntim / spd  # of days (total) for input variable
    nSampTot = nDayTot * spd  # of samples (total)
    nSampWin = nDayWin * spd  # of samples per temporal window
    nSampSkip = nDaySkip * spd  # of samples to skip between window segments
                                # neg means overlap
    nWindow = (nSampTot - nSampWin) / (nSampWin + nSampSkip) + 1
    N = nSampWin  # convenience [historical]

    if nDayTot < nDayWin :
        print "nDayTot=" + nDayTot + " is less the nDayWin=" + nDayWin
        print "        This is not allowed !!       "
        sys.exit(0)
    #-------------------------------------------------------------------
    #  Remove dominant signals
    # (a) Explicitly remove *long term* linear trend
    #      For consistency with JET code keep the grid point means.
    #      This necessitates that 'dtrend_msg' be used because 'dtrend'
    #      always removes the mean(s).
    #  (b) All variations >= approx 'nDayWin' days if full year available
    # -------------------------------------------------------------------

    # subset the data for 15S-15N
    constraint = iris.Constraint(latitude=lambda cell: latS <= cell <= latN)
    var = var.extract(constraint)

    ntim, nlat, mlon = var.shape
    peeAS = np.zeros([nlat, nSampWin + 1, mlon + 1])  # initialize

    # Wave numbers
    wave = np.arange(-mlon / 2, mlon / 2 + 1, 1)
    # Frequencies
    freq = np.linspace(-1 * nDayWin * spd / 2, nDayWin * spd / 2, nDayWin * spd + 1) / nDayWin

    wave = wave.astype(float)
    freq = freq.astype(float)
    lats = var.coord('latitude').points

    # Time mean (later to be added to the trend)
    varmean = var.collapsed('time', iris.analysis.MEAN)

    # remove linear trend
    var.data = signal.detrend(var.data, axis=0) + varmean.data # Mean added
    # print 'Passed detrend test!'

    print 'nDayTot = ' + str(nDayTot)

    if nDayTot >= 365:  #  remove dominant signals
        rmvMeans = True  #  original code did not remove
        var = _rmvAnnualCycle(var, nDayTot, fCrit, spd=1, rmvMeans=False)
    else:
        print 'Length of the variable is shorter than 365. Can not continue!'
        sys.exit(1)

    #-------------------------------------------------------------------
    #  Decompose to Symmetric and Asymmetric parts
    # -------------------------------------------------------------------
    xAS = _decompose2SymAsym(var)  #   create Asym and Sym parts
    # print 'Passed decomposition test.'

    # -------------------------------------------------------------------
    #  Because there is the possibility of overlapping *temporal* segments,
    #  we must use a less efficient approach and detrend/taper
    #  each window segment as it arises.
    #           t0   t1   t2   t3   t4  .................. t(N)
    #  lon(0):  x00  x01  x02  x03  x04 .................. x0(N)
    #      :    :   :   :   :   :                     :
    #  lon(M):  xM0  xM1  xM2  xM3  xM4 .................. xM(N)
    # -------------------------------------------------------------------
    #  q     - temporary array to hold the 2D complex results
    #          for each longitude/time (lon,time) window that is fft'd.
    #          This is one instance [realization] of space-time decomposition.
    #
    #  peeAS - symmetric and asymmetric power values in each latitude hemisphere.
    #          Add extra lon/time to match JET
    # -------------------------------------------------------------------
    print 'nSampWin = ' + str(nSampWin)

    for nl in range(nlat):
        nw = 0
        print 'Latitude: nl = ' + str(nl)
        ntStrt = 0
        ntLast = nSampWin
        while ntLast < nDayTot:
            if nl == 0:
                print 'nw = %s, ntStrt = %s, ntLast =%s ' % (nw, ntStrt, ntLast)
            work = xAS[ntStrt:ntLast, nl].copy()

            # detrend the window
            work.data = signal.detrend(xAS.data[ntStrt:ntLast, nl], axis=0)

            # taper the window along time axis
            # equivalent to NCL taper function described as
            # split-cosine-bell tapering.
            for lo in range(mlon):
                work.data[:, lo] = _taper(work.data[:, lo], \
                                         alpha=tim_taper, \
                                         iopt=0)
            # print 'Passed Tapering test'

            # Do actual FFT work
            ft = work.copy()
            ft.data = np.fft.fft2(work.data) / mlon / nSampWin

            # Shifting FFTs
            pee = resolveWavesHayashi(ft.data, nDayWin, spd)

            # Average
            peeAS[nl, :, :] = peeAS[nl, :, :] + (pee / nWindow)
            nw += 1
            ntStrt = ntLast + nSampSkip  #     set index for next temporal window
            ntLast = ntStrt + nSampWin

    peeAS_cube = _makespeccube(peeAS, lats, wave, freq)

    # -------------------------------------------------------------------
    #  now that we have the power array for sym and asym: use to
    #     1) plot raw power spectrum (some smoothing)
    #     2) derive and plot the background spectrum (lots of smoothing)
    #     3) derive a denoised spectrum that is raw power/background power
    # -------------------------------------------------------------------
    #  psumanti and psumsym will contain the symmetric and asymmetric power
    #  summed over latitude
    # -------------------------------------------------------------------
    if nlat % 2 == 0:
        psumanti = np.sum(peeAS[nlat / 2:nlat], axis=0)
        psumsym = np.sum(peeAS[:nlat / 2], axis=0)
    else:
        psumanti = np.sum(peeAS[nlat / 2 + 1:nlat], axis=0)
        psumsym = np.sum(peeAS[:nlat / 2 + 1], axis=0)
    # -------------------------------------------------------------------
    #  since summing over half the array (symmetric,asymmetric) the
    #  total variance is 2x the half sum
    # -------------------------------------------------------------------
    psumanti = 2.0 * psumanti
    psumsym = 2.0 * psumsym
    #-------------------------------------------------------------------
    # set the mean to missing to match original code
    #------------------------------------------------------------------
    zeroind = np.where(freq == 0.)[0][0]

    psumanti[zeroind, :] = np.nan
    psumsym[zeroind, :] = np.nan
    psumanti = np.ma.masked_invalid(psumanti)
    psumsym = np.ma.masked_invalid(psumsym)

    # -------------------------------------------------------------------
    #  Apply smoothing to the spectrum. smooth over limited wave numbers
    #  Smoothing in frequency only (check if mean should be smoothed
    #  not smoothing now)
    # --
    #  Smoothing parameters set these larger than the plotting
    #  wavenumbers to avoid smoothing artifacts
    # -------------------------------------------------------------------
    minwav4smth = -27
    maxwav4smth = 27

    indStrt = np.where(minwav4smth == wave)[0][0]
    indLast = np.where(maxwav4smth == wave)[0][0]

    for wv in np.arange(indStrt, indLast + 1):
        psumanti[N / 2 + 1:N, wv] = wk_smooth121(psumanti[N / 2 + 1:N, wv])
        psumsym[N / 2 + 1:N, wv] = wk_smooth121(psumsym[N / 2 + 1:N, wv])
    # -------------------------------------------------------------------
    #  Log10 scaling
    # -------------------------------------------------------------------
    psumanti_nolog = np.ma.masked_array(psumanti)
    psumsym_nolog = np.ma.masked_array(psumsym)

    psumanti = np.ma.log10(psumanti)
    psumsym = np.ma.log10(psumsym)

    # Creating Iris cube, assigning metadata
    psumanti_cube = _makecube(psumanti, wave, freq)
    psumsym_cube = _makecube(psumsym, wave, freq)

    # -----------------------------------------------------------------------------
    #  ******  now derive and plot the background spectrum (red noise) ************
    #  [1] Sum power over all latitude
    #  [2] Put fill value in mean
    #  [3] Apply smoothing to the spectrum. This smoothing DOES include
    #      wavenumber zero.
    # -----------------------------------------------------------------------------
    # print("======> BACKGROUND <=====")

    psumb = ComputeBackground(peeAS, wave, freq, indStrt, indLast)
    psumb_nolog = np.ma.masked_array(psumb)
    psumb = np.ma.log10(psumb)


    # -------------------------------------------------------------------------------
    #  Plot section
    # Fig.1 - Raw spectra Symmetric and Anti-symmetric
    #

    # Define contour levels for plots
    if varname == 'toa_outgoing_longwave_flux':
        levels = np.array([-1.3, -1.2, -1.1, -1, -0.8, -0.6, -0.4, -0.2, \
                           0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.2, 1.3])
    elif varname == 'precipitation_flux':
        levels = np.array([-0.5, -0.4, -0.3, -0.2, -0.1, 0, \
                           0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    elif varname == 'x_wind_850hPa':
        levels = np.arange(-3.25, 0.5, 0.25)
    elif varname == 'x_wind_200hPa':
        levels = np.arange(-3.3, 1.2, 0.3)

    title = 'Anti-symmetric ' + varname + ' log(power)' + '\n'
    title += ' '.join([runid, start_date, '-', end_date, '[15S-15N]'])
    forename = runid + '_' + varname + '_Raw_Spec_Asym'
    figname = os.path.join(outdir, "%s.png" % forename)
    # mjp.PlotSpecRaw(psumanti, freq, wave, levels=levels,
    #         title=title, figname=figname)
    ncname = os.path.join(outdir, "%s.nc" % forename)
    iris.save(psumanti_cube, ncname, netcdf_format="NETCDF3_CLASSIC")

    title = 'Symmetric ' + varname + ' log(power)' + '\n'
    title += ' '.join([runid, start_date, '-', end_date, '[15S-15N]'])
    forename = runid + '_' + varname + '_Raw_Spec_Sym'
    figname = os.path.join(outdir, "%s.png" % forename)
    # mjp.PlotSpecRaw(psumsym, freq, wave, levels=levels,
    #         title=title, figname=figname)
    ncname = os.path.join(outdir, "%s.nc" % forename)
    iris.save(psumsym_cube, ncname, netcdf_format="NETCDF3_CLASSIC")

    # Fig.2 - Background spectra
    #
    title = 'Background power ' + varname + ' log(power)' + '\n'
    title += ' '.join([runid, start_date, '-', end_date, '[15S-15N]'])
    forename = runid + '_' + varname + '_BG_Spec'
    figname = os.path.join(outdir, "%s.png" % forename)
    # mjp.PlotSpecRaw(psumb, freq, wave, levels=levels,
    #         title=title, figname=figname)
    psumb_cube = _makecube(psumb, wave, freq)
    ncname = os.path.join(outdir, "%s.nc" % forename)
    iris.save(psumb_cube, ncname, netcdf_format="NETCDF3_CLASSIC")

    # *************************************************************
    #  Fig 3a, 3b:  psum_nolog/psumb_nolog  [ratio]
    # ***************************************************************
    psumanti_nolog = np.ma.masked_array(psumanti_nolog / psumb_nolog)
    psumsym_nolog = np.ma.masked_array(psumsym_nolog / psumb_nolog)  #  (wave,freq)

    # Theoretical dispersion curves
    rlat = 0.0
    Ahe = np.array([50., 25., 12.])
    nWaveType = 6
    nPlanetaryWave = 50
    nEquivDepth = Ahe.size
    fillval = 1e20
    Afreq, Apzwn = genDispersionCurves(nWaveType, \
                   nEquivDepth, nPlanetaryWave, rlat, \
                   Ahe, fillval)
    # removing missing values
    Afreq = np.ma.masked_values(Afreq, fillval)
    Apzwn = np.ma.masked_values(Apzwn, fillval)

    # Define contour levels for plots
    if varname == 'toa_outgoing_longwave_flux':
        levels = np.array([0.2, .3, .4, .5, .6, .7, .8, .9, 1., 1.1, 1.2, 1.3,
            1.4, 1.5, 1.6, 1.7, 1.8])
    elif varname == 'precipitation_flux':
        levels = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.15, 1.2, 1.25,
            1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7])
    elif varname == 'x_wind_850hPa':
        levels = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2,
            1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9])
    elif varname == 'x_wind_850hPa':
        levels = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2,
            1.3, 1.4, 1.5, 1.6,  1.7, 1.8, 2])

    title = 'Anti-symmetric/Background ' + varname + ' log(power)' + '\n'
    title += ' '.join([runid, start_date, '-', end_date, '[15S-15N]'])
    forename = runid + '_' + varname + '_Ratio_Spec_Asym'
    figname = os.path.join(outdir, "%s.png" % forename)
    mjp.plotAntiSymmetric(psumanti_nolog, freq, wave, Apzwn, Afreq, \
        levels=levels, title=title, figname=figname)

    psumanti_nolog_cube = _makecube(psumanti_nolog, wave, freq)
    ncname = os.path.join(outdir, "%s.nc" % forename)
    iris.save(psumanti_nolog_cube, ncname, netcdf_format="NETCDF3_CLASSIC")

    # Define contour levels for plots
    if varname == 'toa_outgoing_longwave_flux':
        levels = np.array([0.2, .3, .4, .5, .6, .7, .8, .9, 1., 1.1, 1.2, 1.4,
            1.7, 2., 2.4, 2.8, 3.2])
    elif varname == 'precipitation_flux':
        levels = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.15, 1.2, 1.25,
            1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7])
    elif varname == 'x_wind_850hPa':
        levels = np.array([0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.3, 1.4, 1.5, 1.6,
            1.7, 1.8, 2, 2.2, 2.4, 2.6, 2.8])
    elif varname == 'x_wind_850hPa':
        levels = np.array([0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.3, 1.4, 1.5, 1.6,
            1.7, 1.8, 2, 2.2, 2.4, 2.6, 2.8])

    title = 'Symmetric/Background ' + varname + ' log(power)' + '\n'
    title += ' '.join([runid, start_date, '-', end_date, '[15S-15N]'])
    forename = runid + '_' + varname + '_Ratio_Spec_Sym'
    figname = os.path.join(outdir, "%s.png" % forename)
    mjp.plotSymmetric(psumsym_nolog, freq, wave, Apzwn, Afreq, \
        levels=levels, title=title, figname=figname)

    psumsym_nolog_cube = _makecube(psumsym_nolog, wave, freq)
    ncname = os.path.join(outdir, "%s.nc" % forename)
    iris.save(psumsym_nolog_cube, ncname, netcdf_format="NETCDF3_CLASSIC")


def mjo_wavenum_freq_season(var, seaName):
    #
    #  For each 'seaName' (say, 'winter') over a number
    #  of different years, calculate the
    #  pooled space-time spectra
    #  MJO CLIVAR: wavenumber-frequency spectra
    #              winter starts Nov 1  [180 days]
    #              summer starts May 1  [180 days]
    #              annual starts Jan 1  [365 days]
    ntim, mlon = var.shape
    time = var.coord('time')
    years, months, days = mu.getDates(time)
    mmdd = months * 100 + days
    yyyymmdd = years * 10000 + mmdd
    nDay = 180
    if seaName == "winter":
        mmddStrt = 1101
    if seaName == "summer":
        mmddStrt = 501
    if seaName == "annual":
        nDay = 365
        mmddStrt = 101

    iSea = [i for i in range(len(mmdd)) if  mmdd[i] == mmddStrt]
    nYear = len(iSea)

    '''
    # *****************************************************************
    #  For a specific season, calculate spectra via averaging
    #  over each seasonal segment.
    #  MJO Clivar says "no" to detrending/tapering.
    #  Hence, the following are just 'place holders'
    # *****************************************************************
    '''
    # detrend overall series in time
    # Time mean (later to be added to the trend)
    varmean = var.collapsed('time', iris.analysis.MEAN)

    # remove linear trend
    var.data = signal.detrend(var.data, axis=0) #+ varmean.data # Mean added

    # Initialise
    power = np.zeros([mlon + 1, nDay + 1])  # initialize
    work = var.data
    xAvgSea = 0.0
    xVarSea = 0.0                            #  variance (raw)
    xVarTap = 0.0                            #  variance after tapering
    kSea = 0                                #  count of seasons used
    N = nDay                            #  convenience

    for ny in range(nYear):
        iStrt = iSea[ny]                  # start index for current season
        iLast = iSea[ny] + nDay      # last
        if iLast < ntim - 1:
            print ny, kSea, iStrt, iLast, yyyymmdd[iStrt], yyyymmdd[iLast]
            xSeason = work[iStrt:iLast, :]
            xAvg = np.average(xSeason)              #  season average all time/lon
            xSeason = xSeason - xAvg              #  remove season time-lon mean
            xVarSea = xVarSea + np.var(xSeason)    #  overall variance
            kSea = kSea + 1

            for lo in range(mlon):
                xSeason[:, lo] = _taper(xSeason[:, lo], \
                                             alpha=0.1, \
                                             iopt=0)
            xVarTap = xVarTap + np.var(xSeason)    # variance after tapering
            # do 2d fft
            ft = xSeason.copy()
            ft = np.fft.fft2(xSeason.T) / mlon / nDay

            # Shifting FFTs
            power = power + resolveWavesHayashi(ft, nDay, spd=1)     #  (wave,freq)

    xVarSea = xVarSea / kSea            #   pooled seasonal variance
    xVarTap = xVarTap / kSea            #    pooled seasonal variance
    power = np.ma.masked_array(power / kSea)                   #   pooled spectra

    wave = np.arange(-mlon / 2, mlon / 2 + 1, 1)
    freq = np.linspace(-1 * nDay / 2, nDay / 2, nDay + 1) / nDay
    wave = wave.astype(float)
    freq = freq.astype(float)
    return power, wave, freq


def diagnos_level2(var, outdir, runid):

    with iris.FUTURE.context(cell_datetime_objects=True):
        start_date = var.coord('time').cell(0).point.strftime('%Y-%m-%d')
        end_date = var.coord('time').cell(-1).point.strftime('%Y-%m-%d')

    varname = var.name()
    if varname == 'x_wind':
        assert len(var.coord('pressure').points) == 1
        pressure_level = var.coord('pressure').points[0]
        if pressure_level == 850:
            varname = 'x_wind_850hPa'
        elif pressure_level == 200:
            varname = 'x_wind_200hPa'

    ### Wheeler-Kiladis spectral calculation
    latBound = 15
    coords = [0, -latBound, 360, latBound]
    xvar = var.extract(mu.region(coords))

    spd = 1
    nDayWin = 96  #  Wheeler-Kiladis [WK] temporal window length (days)
    nDaySkip = -65
    xs = wkSpaceTime(xvar, outdir, runid, latBound, spd, nDayWin, nDaySkip)


    ### Frequency spectra for each season
    # The wavenumber - frequency spectra for each season are calculated
    # via the mjo_wavenum_freq_season function.

    latBound = 10
    coords = [0, -latBound, 360, latBound]
    avar = var.extract(mu.region(coords))
    avar = var.collapsed('latitude', iris.analysis.MEAN)

    metrics = {}
    for season in ['winter', 'summer']:
        # This is the NCL method of computing the spectra which uses anomalies
        power, wave, freq = mjo_wavenum_freq_season(avar, season)
        pow_cube = _makecube_season_pow(power, wave, freq)

        # Define contour levels for plots
        if varname == 'toa_outgoing_longwave_flux':
            levels = np.arange(0.0, 2.4, 0.2)
        elif varname == 'precipitation_flux':
            levels = np.arange(0.0, 0.055, 0.005)
        elif varname == 'x_wind_850hPa':
            levels = np.arange(0.007, 0.07, 0.007)
        elif varname == 'x_wind_200hPa':
            levels = np.arange(0.05, 0.5, 0.05)

        title = runid + ' ' + season + ' daily ' + varname + '[10S-10N]'
        title += '\n' + start_date + ' - ' + end_date
        forename = runid + '_' + varname + '_wavenum_freq_season_' + season
        figname = os.path.join(outdir, "%s.png" % forename)
        mjp.mjo_wavenum_freq_season_plot(pow_cube, levels=levels,
                                         title=title, figname=figname)

        ncname = os.path.join(outdir, "%s.nc" % forename)
        iris.save(pow_cube, ncname, netcdf_format="NETCDF3_CLASSIC")  # so xconv can read it

        # Calculate metrics
        _metrics = east_west_power_ratio(pow_cube, season, varname)
        metrics.update(_metrics)

    return metrics


def east_west_power_ratio(power_spectra_cube, season, var_name):
    """Compute Power ratio of eastward and westward propagating waves.

    Args:
        power_spectra_cube: Cube containing power spectrum.
        season: Season name.
        var_name: Name of the variable from which the power spectrum was
            calculated.

    Returns: metrics dict

    Reference:
        D. Kim, et al (2014). Process-Oriented MJO Simulation Diagnostic:
        Moisture Sensitivity of Simulated Convection. J. Climate, 27, 5379-5395.
        DOI: 10.1175/JCLI-D-13-00497.1
    """
    metrics = {}
    freq_constraint = iris.Constraint(frequency=lambda cell: 0.0167 <= cell <= 0.0333)
    east_constraint = iris.Constraint(wavenumber=lambda cell: 1 <= cell <= 3)
    west_constraint = iris.Constraint(wavenumber=lambda cell: -3 <= cell <= -1)

    east = power_spectra_cube.extract(freq_constraint & east_constraint)
    west = power_spectra_cube.extract(freq_constraint & west_constraint)
    east = east.data.sum().sum()
    west = west.data.sum().sum()
    metrics[' '.join(['East power', season, var_name])] = round(east, 2)
    metrics[' '.join(['East/west power', season, var_name])] = round(east / west, 2)

    return metrics
