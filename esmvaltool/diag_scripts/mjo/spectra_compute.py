import math
import os
import sys

import cf_units as unit
import iris
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import scipy

from esmvaltool.diag_scripts.shared import save_figure


class WKSpectra:
    def __init__(self, cfg, attributes, check_missing=True):
        self.cfg = cfg
        self.spd = 1
        self.nDayWin = 96  # Wheeler-Kiladis [WK] temporal window length (days)
        self.nDaySkip = -65  # Negative means overlap
        self.latBound = 15
        self.filename = attributes["filename"]
        self.out_dir = os.path.dirname(self.filename)
        self.runid = attributes["dataset"]
        self.label = f"{attributes['dataset']}_{attributes['ensemble']}"
        self.plot_dir = attributes["plot_dir"]
        self.work_dir = attributes["work_dir"]

        self.cube = attributes["cube"]
        self.varname = attributes["varname"]

        if check_missing:
            # Checking for any missing data in the cube.
            # If found interpolating along longitude to fill gaps.
            if np.any(self.cube.data.mask):
                self.cube = self.interpolate_along_axis(self.cube, "longitude")
        print(self.cube)

    def interpolate_along_axis(self, cube, coord_name):
        """
        Interpolate masked values in a masked 3D array along a specified axis.

        Parameters:
        arr (np.ma.MaskedArray): The masked array with missing values.
        axis (int): The axis along which to perform the interpolation.

        Returns:
        np.ndarray: The array with interpolated values along the specified axis.
        """
        interpolated_cube = cube.copy()

        # Get the shape of the array
        shape = cube.shape

        # Initialize an empty array to store interpolated data
        interpolated_data = np.copy(cube.data)  # Copy the unmasked data

        # Find the longitude dimension index
        axis = cube.coord_dims(coord_name)[0]

        # Loop over all the axes except the one we're interpolating along
        # Create an iterator over the other dimensions
        it = np.nditer(np.ones(np.delete(shape, axis)), flags=["multi_index"])

        while not it.finished:
            # Create slices to fix all dimensions except the one we're interpolating along
            slices = list(it.multi_index)
            slices.insert(axis, slice(None))
            slices = tuple(slices)

            # Extract the 1D slice along the given axis
            slice_data = cube.data[slices]

            # Check if there are any masked values in this slice
            if np.any(slice_data.mask):
                # Get the indices of the non-masked values
                valid = ~slice_data.mask
                valid_indices = np.where(valid)[0]

                # Interpolate only if there are enough valid points to perform interpolation
                if valid_indices.size > 1:
                    # Perform interpolation along the axis
                    interp_func = scipy.interpolate.interp1d(
                        valid_indices,
                        slice_data[valid],
                        kind="linear",
                        fill_value="extrapolate",
                    )

                    # Fill the masked values with the interpolated values
                    interp_indices = np.arange(slice_data.shape[0])
                    interpolated_data[slices] = interp_func(interp_indices)

            it.iternext()  # Move to the next slice

        interpolated_cube.data = interpolated_data
        return interpolated_cube

    def split_time(self, date):
        """Split date string into yyyy, mm, dd integers"""
        d = str(date)
        year = int(d[0:4])
        month = int(d[5:7])
        day = int(d[8:10])
        return year, month, day

    def get_dates(self, time):
        """splits the iris time coordinate into yyyy, mm, dd"""
        if time.units.calendar == "gregorian":
            dates = unit.num2date(
                time.points, str(time.units), unit.CALENDAR_GREGORIAN
            )
        if time.units.calendar == "standard":
            dates = unit.num2date(
                time.points, str(time.units), unit.CALENDAR_STANDARD
            )
        else:
            dates = unit.num2date(
                time.points, str(time.units), time.units.calendar
            )
        try:
            dates
        except NameError:
            print(
                dates + " WASN'T defined after all!"
            )  # todo Effective but probs frowned upon
        else:
            year = np.zeros(len(dates), dtype=int)
            month = np.zeros(len(dates), dtype=int)
            day = np.zeros(len(dates), dtype=int)
            for i, date in enumerate(dates):
                year[i], month[i], day[i] = self.split_time(date)
        return year, month, day

    def makecube_season_pow(self, var, wave, freq, name="spectra"):
        # ===========================================================================
        # Make a cube of seasonal power
        # ===========================================================================
        var_cube = iris.cube.Cube(var)
        var_cube.rename("spectra")
        # var_cube.long_name = long_name
        wave_coord = iris.coords.DimCoord(wave, long_name="wavenumber")
        wave_coord.guess_bounds()
        freq_coord = iris.coords.DimCoord(freq, long_name="frequency")
        freq_coord.guess_bounds()
        var_cube.add_dim_coord(wave_coord, 0)
        var_cube.add_dim_coord(freq_coord, 1)
        return var_cube

    def remove_annual_cycle(self, var, nDayTot, fCrit, spd=1, rmvMeans=False):
        # ===========================================================================
        # Prewhiten the data: eg remove the annual cycle.
        # Actually, this will remove all time periods less than
        #   those corresponding to 'fCrit'.
        # Note: The original fortran code provided by JET did not remove
        #   the grid point means so .... rmvMeans=False
        #       I assume that Matt Wheeler's code did that also.
        # ===========================================================================
        print(fCrit)
        # Take a copy for metadata
        var_cut = var.copy()

        # mean for later
        varMean = var.collapsed("time", iris.analysis.MEAN)

        ntime, nlat, nlon = var.data.shape

        # compute FFT
        cf = np.fft.fft(var.data, axis=0)

        freq = np.fft.fftfreq(ntime)
        x = freq.copy()

        # cutting frequencies
        cf[np.abs(x) < fCrit] = 0.0

        # inverse FFT
        var_cut.data = np.fft.ifft(cf, axis=0).astype(float)
        if not rmvMeans:
            var_cut.data += varMean.data
        return var_cut

    def decompose_sym_asym(self, var, axis=1):
        # ===========================================================================
        # This code decomposes the data in to symmetric and anti-symmetric components
        # with respect to latitude axis. It is assumed that the latitude dimension is
        # along axis = 1 as default.
        #
        # antisymmetric part is stored in one hemisphere [eg: Northern Hemisphere]
        #          xOut( lat) = (x(lat)-x(-lat))/2
        # symmetric part is stored in other hemisphere [eg: Southern Hemisphere]
        #          xOut(-lat) = (x(lat)+x(-lat))/2
        # ===========================================================================
        varSA = var.copy()  # copy to store the output
        nlat = var.shape[axis]  # get the number of latitude points

        N2 = nlat // 2
        N = N2
        if nlat % 2 == 1:
            N = N2 + 1  # offset to handle the Equator

        if axis == 1:
            for nl in np.arange(0, N2):
                print((nlat - 1 - nl), nl)
                varSA.data[:, nl] = 0.5 * (
                    var.data[:, nlat - 1 - nl] + var.data[:, nl]
                )
                varSA.data[:, nlat - 1 - nl] = 0.5 * (
                    var.data[:, nlat - 1 - nl] - var.data[:, nl]
                )

            return varSA
        else:
            print("Modify the code to accommodate other axes...")
            print("Exiting...")
            sys.exit(0)

    def taper(self, ts, alpha=0.1, iopt=0):
        # ===========================================================================
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
        # ===========================================================================
        if all(x == ts[0] for x in ts):
            print("all values equal")
            iopt = 1
        if iopt == 0:
            tsmean = np.mean(ts)
        else:
            tsmean = 0.0
        n = len(ts)
        m = max(1, int(alpha * n + 0.5) // 2)
        pim = np.pi / m

        tst = ts.copy()
        for i in range(1, m + 1):
            weight = 0.5 - 0.5 * np.cos(pim * (i - 0.5))
            tst[i - 1] = (ts[i - 1] - tsmean) * weight + tsmean
            tst[n - i] = (ts[n - i] - tsmean) * weight + tsmean
        return tst

    def resolve_waves_hayashi(self, varfft, nDayWin, spd):
        # ===========================================================================
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
        # ===========================================================================
        N, mlon = varfft.shape
        pee = np.ones([N + 1, mlon + 1]) * -999.0  # initialize
        # -999 scaling is for testing
        # purpose
        # Create the real power spectrum pee = sqrt(real^2+imag^2)^2
        varfft = np.abs(varfft) ** 2
        pee[: N // 2, : mlon // 2] = varfft[N // 2 : N, mlon // 2 : 0 : -1]
        pee[N // 2 :, : mlon // 2] = varfft[: N // 2 + 1, mlon // 2 : 0 : -1]
        pee[: N // 2 + 1, mlon // 2 :] = varfft[N // 2 :: -1, : mlon // 2 + 1]
        pee[N // 2 + 1 :, mlon // 2 :] = varfft[
            N - 1 : N // 2 - 1 : -1, : mlon // 2 + 1
        ]

        return pee

    def wk_smooth121(self, var):
        # ===========================================================================
        # Special 1-2-1 smoother
        # Smooths vv by passing it through a 1-2-1 filter.
        # The first and last points are given 3-1 (1st) or 1-3 (last)
        # weightings (Note that this conserves the total sum).
        # The routine also skips-over missing data (np.nan)
        # ===========================================================================
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
                varf[i] = (
                    1.0 * var[i - 1] + 2.0 * var[i] + 1.0 * var[i + 1]
                ) / 4.0
        return varf

    def make_spec_cube(self, var, lat, wave, freq):
        # ===========================================================================
        # # Make a 3D cube of Latitude, wavenumber & frequency dimensions
        # ===========================================================================
        var_cube = iris.cube.Cube(var)
        var_cube.rename("spectra")
        lat_coord = iris.coords.DimCoord(lat, long_name="latitude")
        lat_coord.guess_bounds()
        wave_coord = iris.coords.DimCoord(wave, long_name="wavenumber")
        wave_coord.guess_bounds()
        freq_coord = iris.coords.DimCoord(freq, long_name="frequency")
        freq_coord.guess_bounds()
        var_cube.add_dim_coord(lat_coord, 0)
        var_cube.add_dim_coord(freq_coord, 1)
        var_cube.add_dim_coord(wave_coord, 2)
        return var_cube

    def make_cube(self, var, wave, freq):
        # ===========================================================================
        # # Make a 2D cube of wavenumber & frequency dimensions
        # ===========================================================================
        var_cube = iris.cube.Cube(var)
        var_cube.rename("spectra")
        wave_coord = iris.coords.DimCoord(wave, long_name="wavenumber")
        wave_coord.guess_bounds()
        freq_coord = iris.coords.DimCoord(freq, long_name="frequency")
        freq_coord.guess_bounds()
        var_cube.add_dim_coord(freq_coord, 0)
        var_cube.add_dim_coord(wave_coord, 1)
        return var_cube

    def closest_index(self, array, value):
        # ===========================================================================
        # Find the closest index to a given value in an array
        # ===========================================================================
        return (np.abs(array - value)).argmin()

    def compute_background(self, peeAS, wave, freq, minwav4smth, maxwav4smth):
        # ===========================================================================
        # Derive the background spectrum (red noise) ************
        # [1] Sum power over all latitude
        # [2] Put fill value in mean
        # [3] Apply smoothing to the spectrum. This smoothing DOES include wavenumber zero.
        #
        # ===========================================================================
        psumb = np.sum(peeAS, axis=0)  # sum over all latitudes
        N, mlon = psumb.shape
        smthlen = maxwav4smth - minwav4smth + 1

        for tt in range(N // 2 + 1, N):
            if freq[tt] < 0.1:
                for i in range(1, 6):
                    psumb[tt, minwav4smth : maxwav4smth + 1] = (
                        self.wk_smooth121(
                            psumb[tt, minwav4smth : maxwav4smth + 1]
                        )
                    )
            if freq[tt] >= 0.1 and freq[tt] < 0.2:
                for i in range(1, 11):
                    psumb[tt, minwav4smth : maxwav4smth + 1] = (
                        self.wk_smooth121(
                            psumb[tt, minwav4smth : maxwav4smth + 1]
                        )
                    )
            if freq[tt] >= 0.2 and freq[tt] < 0.3:
                for i in range(1, 21):
                    psumb[tt, minwav4smth : maxwav4smth + 1] = (
                        self.wk_smooth121(
                            psumb[tt, minwav4smth : maxwav4smth + 1]
                        )
                    )
            if freq[tt] >= 0.3:
                for i in range(1, 41):
                    psumb[tt, minwav4smth : maxwav4smth + 1] = (
                        self.wk_smooth121(
                            psumb[tt, minwav4smth : maxwav4smth + 1]
                        )
                    )

        pt8cpd = min([self.closest_index(freq, 0.8), len(freq) - 1])

        # smth frequency up to .8 cycles per day
        for nw in range(minwav4smth, maxwav4smth + 1):
            smthlen = pt8cpd - (N // 2 + 1) + 1
            for i in range(1, 11):
                psumb[N // 2 + 1 : pt8cpd + 1, nw] = self.wk_smooth121(
                    psumb[N // 2 + 1 : pt8cpd + 1, nw]
                )
        return psumb

    def generate_dispersion_curves(self):
        # Theoretical dispersion curves
        rlat = 0.0
        Ahe = np.array([50.0, 25.0, 12.0])
        nWaveType = 6
        nPlanetaryWave = 50
        nEquivDepth = Ahe.size
        fillval = 1e20
        # ---------------------------------------------------------------
        # Theoretical shallow water dispersion curves
        # --------------------------------------------------------------
        pi = 4.0 * math.atan(1.0)
        re = 6.37122e06  # [m]   average radius of earth
        g = 9.80665  # [m/s] gravity at 45 deg lat used by the WMO
        omega = 7.292e-05  # [1/s] earth's angular vel
        U = 0.0
        Un = 0.0  # since Un = U*T/L
        ll = 2.0 * pi * re * math.cos(abs(rlat))
        Beta = 2.0 * omega * math.cos(abs(rlat)) / re
        maxwn = nPlanetaryWave

        Apzwn = np.zeros(
            [nWaveType, nEquivDepth, nPlanetaryWave], dtype=np.double
        )
        Afreq = np.zeros(
            [nWaveType, nEquivDepth, nPlanetaryWave], dtype=np.double
        )

        for ww in range(1, nWaveType + 1):  # wave type
            for ed in range(1, nEquivDepth + 1):  # equivalent depth
                he = Ahe[ed - 1]
                T = 1.0 / math.sqrt(Beta) * (g * he) ** (0.25)
                L = (g * he) ** (0.25) / math.sqrt(Beta)

                for wn in range(
                    1, nPlanetaryWave + 1
                ):  # planetary wave number
                    s = -20.0 * (wn - 1) * 2.0 / (nPlanetaryWave - 1) + 20.0
                    k = 2.0 * pi * s / ll
                    kn = k * L

                    # Anti-symmetric curves
                    if ww == 1:  # MRG wave
                        if k <= 0:
                            delx = math.sqrt(
                                1.0 + (4.0 * Beta) / (k**2 * math.sqrt(g * he))
                            )
                            deif = k * math.sqrt(g * he) * (0.5 - 0.5 * delx)
                        if k == 0:
                            deif = math.sqrt(math.sqrt(g * he) * Beta)
                        if k > 0:
                            deif = fillval

                    if ww == 2:  # n=0 IG wave
                        if k < 0:
                            deif = fillval
                        if k == 0:
                            deif = math.sqrt(math.sqrt(g * he) * Beta)
                        if k > 0:
                            delx = math.sqrt(
                                1.0 + (4.0 * Beta) / (k**2 * math.sqrt(g * he))
                            )
                            deif = k * math.sqrt(g * he) * (0.5 + 0.5 * delx)

                    if ww == 3:  # n=2 IG wave
                        n = 2.0
                        delx = Beta * math.sqrt(g * he)
                        deif = math.sqrt(
                            (2.0 * n + 1.0) * delx + (g * he) * k**2
                        )
                        # do some corrections to the above calculated frequency.......
                        for i in range(1, 6):
                            deif = math.sqrt(
                                (2.0 * n + 1.0) * delx
                                + (g * he) * k**2
                                + g * he * Beta * k / deif
                            )

                    # symmetric curves
                    if ww == 4:  # n=1 ER wave
                        n = 1.0
                        if k < 0:
                            delx = (Beta / math.sqrt(g * he)) * (2.0 * n + 1.0)
                            deif = -Beta * k / (k**2 + delx)
                        else:
                            deif = fillval
                    if ww == 5:  # Kelvin wave
                        deif = k * math.sqrt(g * he)
                    if ww == 6:  # n=1 IG wave
                        n = 1.0
                        delx = Beta * math.sqrt(g * he)
                        deif = math.sqrt(
                            (2.0 * n + 1.0) * delx + (g * he) * k**2
                        )
                        # do some corrections to the above calculated frequency.......
                        for i in range(1, 6):
                            deif = math.sqrt(
                                (2.0 * n + 1.0) * delx
                                + (g * he) * k**2
                                + g * he * Beta * k / deif
                            )

                    eif = deif  # + k*U since  U=0.0
                    P = 2.0 * pi / (eif * 24.0 * 60.0 * 60.0)
                    dps = deif / k
                    R = L
                    Rdeg = (180.0 * R) / (pi * 6.37e6)
                    Apzwn[ww - 1, ed - 1, wn - 1] = s
                    if deif != fillval:
                        P = 2.0 * pi / (eif * 24.0 * 60.0 * 60.0)
                        Afreq[ww - 1, ed - 1, wn - 1] = 1.0 / P
                    else:
                        Afreq[ww - 1, ed - 1, wn - 1] = fillval

        # removing missing values
        Afreq = np.ma.masked_values(Afreq, fillval)
        Apzwn = np.ma.masked_values(Apzwn, fillval)
        return Afreq, Apzwn

    def spread_colorbar(self, C):
        x = np.linspace(0, 256, len(C))
        xnew = np.arange(256)
        fr = scipy.interpolate.interp1d(x, C[:, 0])
        fg = scipy.interpolate.interp1d(x, C[:, 1])
        fb = scipy.interpolate.interp1d(x, C[:, 2])
        C = np.array(
            [fr(xnew).astype(int), fg(xnew).astype(int), fb(xnew).astype(int)]
        ).T
        return C

    def get_colors(self, reverse=False):
        # Provided RGB values
        rgb_values_str = (
            ";R   G   B\n"
            "130 32  240\n"
            "0   0   150\n"
            "0   0   205\n"
            "65  105 225\n"
            "30  144 255\n"
            "0   191 255\n"
            "160 210 255\n"
            "210 245 255\n"
            "255 255 200\n"
            "255 225 50\n"
            "255 170 0\n"
            "255 110 0\n"
            "255 0   0\n"
            "200 0   0\n"
            "160 35  35\n"
            "255 105 180"
        )

        # Split the string into lines
        lines = rgb_values_str.split("\n")

        # Remove the header line
        lines.pop(0)

        # Initialize an empty list to store RGB values
        rgb_values = []

        # Parse each line and extract RGB values
        for line in lines:
            # Split the line into individual RGB values
            r, g, b = map(int, line.split())
            # Append the RGB values to the list
            rgb_values.append([r, g, b])

        # Convert the list to a NumPy array
        C = np.array(rgb_values)

        if len(C) < 256:
            C = self.spread_colorbar(C)
            cm = mpl.colors.ListedColormap(C / 256.0)
            if reverse:
                cm = mpl.colors.ListedColormap(C[::-1, :] / 256.0)
        return cm

    def get_provenance_record(
        self, caption
    ):  # Credit to AutoAssess _plot_mo_metrics.py
        """Create a provenance record describing the diagnostic data and plot."""
        # Get the list of input filenames
        filenames = [
            item["filename"] for item in self.cfg["input_data"].values()
        ]

        # Write the provenance dictionary using the provided caption
        record = {
            "caption": caption,
            "statistics": ["mean"],
            "domains": ["global"],
            "plot_types": ["zonal"],
            "authors": [
                "xavier_prince",
            ],
            "ancestors": filenames,
        }

        return record

    def plot_anti_symmetric(
        self,
        spec,
        freq,
        wave,
        Apzwn,
        Afreq,
        levels=np.arange(0, 2, 0.2),
        title="",
        figname="specAntiSym_test.ps",
    ):
        plt.clf()  # Not sure this is still needed
        minfrq4plt = 0.0
        maxfrq4plt = 0.8
        minwav4plt = -15
        maxwav4plt = 15

        minfrq = minfrq4plt
        maxfrq = min([maxfrq4plt, max(freq)])
        F, W = np.meshgrid(wave, freq)

        cmap = self.get_colors()

        norm = colors.BoundaryNorm(levels, len(cmap.colors))

        # Initialize the plot
        fig, ax = plt.subplots()

        CS = ax.contourf(
            F, W, spec, levels=levels, cmap=cmap, norm=norm, extend="both"
        )
        bar = fig.colorbar(CS, ax=ax)

        tick_locs = levels
        tick_labels = levels
        bar.locator = ticker.FixedLocator(tick_locs)
        bar.formatter = ticker.FixedFormatter(tick_labels)
        bar.update_ticks()

        # set axes range
        ax.set_xlim(minwav4plt, maxwav4plt)
        ax.set_ylim(minfrq, maxfrq)

        # Lines
        ax.plot([0, 0], [0, 0.5], "k--", lw=0.5)

        # Line markers of periods
        frqs = [80, 30, 6, 3]
        for frq in frqs:
            ax.plot([-15, 15], [1.0 / frq, 1.0 / frq], "k--", lw=0.5)
            ax.text(-14.7, 1.0 / frq, str(frq) + " days", {"color": "k"})

        ax.set_title(title)
        ax.set_xlabel("Westward     Zonal Wave Number     Eastward")
        ax.set_ylabel("Frequency (CPD)")

        # Symmetric waves
        # Equatorial Rossby
        for i in range(3):
            for j in range(3):
                ax.plot(Apzwn[i, j, :], Afreq[i, j, :], "k", lw=0.5)

        ax.text(-10.0, 0.15, "MRG", {"color": "k", "backgroundcolor": "w"})
        ax.text(-3.0, 0.58, "n=2 IG", {"color": "k", "backgroundcolor": "w"})
        ax.text(6.0, 0.4, "n=0 EIG", {"color": "k", "backgroundcolor": "w"})
        ax.text(-3.0, 0.475, "h=12", {"color": "k", "backgroundcolor": "w"})

        # Add provenance information
        caption = f"{figname}, [or other caption for antisymmetric]"  # TODO
        provenance_dict = self.get_provenance_record(caption)

        # Save the figure (also closes it)
        save_figure(
            figname,
            provenance_dict,
            self.cfg,
            figure=fig,
            close=True,
        )
        print(f"Plotted {figname}")

    def plot_symmetric(
        self,
        spec,
        freq,
        wave,
        Apzwn,
        Afreq,
        levels=np.arange(0, 5, 0.5),
        title="",
        figname="specSym_test.ps",
    ):
        plt.clf()  # Again, maybe not needed
        minfrq4plt = 0.0
        maxfrq4plt = 0.8
        minwav4plt = -15
        maxwav4plt = 15

        minfrq = minfrq4plt
        maxfrq = min([maxfrq4plt, max(freq)])
        F, W = np.meshgrid(wave, freq)

        cmap = self.get_colors()
        norm = colors.BoundaryNorm(levels, len(cmap.colors))

        # Initialize the plot
        fig, ax = plt.subplots()

        CS = ax.contourf(
            F, W, spec, levels=levels, cmap=cmap, norm=norm, extend="both"
        )
        bar = fig.colorbar(CS, ax=ax)

        tick_locs = levels
        tick_labels = levels
        bar.locator = ticker.FixedLocator(tick_locs)
        bar.formatter = ticker.FixedFormatter(tick_labels)
        bar.update_ticks()

        # set axes range
        ax.set_xlim(minwav4plt, maxwav4plt)
        ax.set_ylim(minfrq, maxfrq)

        # Lines
        ax.plot([0, 0], [0, 0.5], "k--", lw=0.5)

        # Line markers of periods
        frqs = [80, 30, 6, 3]
        for frq in frqs:
            ax.plot([-15, 15], [1.0 / frq, 1.0 / frq], "k--", lw=0.5)
            ax.text(-14.7, 1.0 / frq, str(frq) + " days", {"color": "k"})

        ax.set_title(title)  #
        ax.set_xlabel("Westward     Zonal Wave Number     Eastward")
        ax.set_ylabel("Frequency (CPD)")

        # Symmetric waves
        # Equatorial Rossby
        for i in range(3, 6):
            for j in range(3):
                ax.plot(Apzwn[i, j, :], Afreq[i, j, :], "k", lw=0.5)

        ax.text(11.5, 0.4, "Kelvin", {"color": "k", "backgroundcolor": "w"})
        ax.text(-10.7, 0.07, "n=1 ER", {"color": "k", "backgroundcolor": "w"})
        ax.text(-3.0, 0.45, "n=1 IG", {"color": "k", "backgroundcolor": "w"})
        ax.text(-14.0, 0.46, "h=12", {"color": "k", "backgroundcolor": "w"})

        # Add provenance information
        caption = f"{figname}, [or other caption for symmetric]"  # TODO
        provenance_dict = self.get_provenance_record(caption)

        # Save the figure (also closes it)
        save_figure(
            figname,
            provenance_dict,
            self.cfg,
            figure=fig,
            close=True,
        )
        print(f"Plotted {figname}")

    def wkSpaceTime(self):
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

        # if self.varname == 'x_wind':
        #    assert len(self.cube.coord('pressure').points) == 1
        #    pressure_level = self.cube.coord('pressure').points[0]
        #    if pressure_level == 850:
        #        varname = 'x_wind_850hPa'
        #    if pressure_level == 200:
        #        varname = 'x_wind_200hPa'

        ntim, nlat, mlon = self.cube.shape
        latN = self.latBound
        latS = -1 * self.latBound  # make symmetric about the equator

        lonL = 0  # -180
        lonR = 360  # 180
        fCrit = 1.0 / self.nDayWin  # remove all contributions 'longer'

        tim_taper = 0.1  # time taper      [0.1   => 10%]
        lon_taper = (
            0.0  # longitude taper [0.0 for globe  only global supported]
        )

        if lon_taper > 0.0 or lonR - lonL != 360.0:
            print("Code does currently allow lon_taper>0 or (lonR-lonL)<360")
            sys.exit(0)

        nDayTot = ntim / self.spd  # of days (total) for input variable
        nSampTot = nDayTot * self.spd  # of samples (total)
        nSampWin = self.nDayWin * self.spd  # of samples per temporal window
        nSampSkip = (
            self.nDaySkip * self.spd
        )  # of samples to skip between window segments
        # neg means overlap
        nWindow = (nSampTot - nSampWin) / (nSampWin + nSampSkip) + 1
        N = nSampWin  # convenience [historical]

        if nDayTot < self.nDayWin:
            print(
                "nDayTot=" + nDayTot + " is less the nDayWin=" + self.nDayWin
            )
            print("        This is not allowed !!       ")
            sys.exit(0)
        # -------------------------------------------------------------------
        #  Remove dominant signals
        # (a) Explicitly remove *long term* linear trend
        #      For consistency with JET code keep the grid point means.
        #      This necessitates that 'dtrend_msg' be used because 'dtrend'
        #      always removes the mean(s).
        #  (b) All variations >= approx 'nDayWin' days if full year available
        # -------------------------------------------------------------------

        # subset the data for 15S-15N
        constraint = iris.Constraint(
            latitude=lambda cell: latS <= cell <= latN
        )
        self.cube = self.cube.extract(constraint)
        ntim, nlat, mlon = self.cube.shape

        peeAS = np.zeros([nlat, nSampWin + 1, mlon + 1])  # initialize

        # Wave numbers
        wave = np.arange(-mlon / 2, mlon / 2 + 1, 1)
        # Frequencies
        freq = (
            np.linspace(
                -1 * self.nDayWin * self.spd / 2,
                self.nDayWin * self.spd / 2,
                self.nDayWin * self.spd + 1,
            )
            / self.nDayWin
        )

        wave = wave.astype(float)
        freq = freq.astype(float)
        lats = self.cube.coord("latitude").points

        # Time mean (later to be added to the trend)
        varmean = self.cube.collapsed("time", iris.analysis.MEAN)

        # remove linear trend
        self.cube.data = (
            scipy.signal.detrend(self.cube.data, axis=0) + varmean.data
        )  # Mean added

        print("nDayTot = " + str(nDayTot))

        if nDayTot >= 365:  # remove dominant signals
            self.cube = self.remove_annual_cycle(
                self.cube, nDayTot, fCrit, spd=1, rmvMeans=False
            )
        else:
            print(
                "Length of the variable is shorter than 365. Can not continue!"
            )
            sys.exit(1)

        # -------------------------------------------------------------------
        #  Decompose to Symmetric and Asymmetric parts
        # -------------------------------------------------------------------
        xAS = self.decompose_sym_asym(self.cube)  # create Asym and Sym parts

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
        print("nSampWin = " + str(nSampWin))

        for nl in range(nlat):
            nw = 0
            print("Latitude: nl = " + str(nl))
            ntStrt = 0
            ntLast = nSampWin
            while ntLast < nDayTot:
                if nl == 0:
                    print(
                        "nw = %s, ntStrt = %s, ntLast =%s "
                        % (nw, ntStrt, ntLast)
                    )
                work = xAS[ntStrt:ntLast, nl].copy()

                # Check for missing data
                # masked_inds = np.where(work.data.mask)
                # if not len(masked_inds[0]) > 0:

                # detrend the window
                work.data = scipy.signal.detrend(
                    xAS.data[ntStrt:ntLast, nl], axis=0
                )

                # taper the window along time axis
                # equivalent to NCL taper function described as
                # split-cosine-bell tapering.
                for lo in range(mlon):
                    work.data[:, lo] = self.taper(
                        work.data[:, lo], alpha=tim_taper, iopt=0
                    )
                # print 'Passed Tapering test'

                # Do actual FFT work
                ft = work.copy()
                ft.data = np.fft.fft2(work.data) / mlon / nSampWin

                # Shifting FFTs
                pee = self.resolve_waves_hayashi(
                    ft.data, self.nDayWin, self.spd
                )

                # Average
                peeAS[nl, :, :] = peeAS[nl, :, :] + (pee / nWindow)

                nw += 1

                # else:
                #    print('Missing data detected. Skipping to the next window...')

                ntStrt = (
                    ntLast + nSampSkip
                )  # set index for next temporal window
                ntLast = ntStrt + nSampWin

        peeAS_cube = self.make_spec_cube(peeAS, lats, wave, freq)

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
            psumanti = np.sum(
                peeAS[nlat // 2 : nlat], axis=0
            )  # // for integer result
            psumsym = np.sum(peeAS[: nlat // 2], axis=0)
        else:
            psumanti = np.sum(peeAS[nlat // 2 + 1 : nlat], axis=0)
            psumsym = np.sum(peeAS[: nlat // 2 + 1], axis=0)
        # -------------------------------------------------------------------
        #  since summing over half the array (symmetric,asymmetric) the
        #  total variance is 2 x the half sum
        # -------------------------------------------------------------------
        psumanti = 2.0 * psumanti
        psumsym = 2.0 * psumsym
        # -------------------------------------------------------------------
        # set the mean to missing to match original code
        # ------------------------------------------------------------------
        zeroind = np.where(freq == 0.0)[0][0]

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
            psumanti[N // 2 + 1 : N, wv] = self.wk_smooth121(
                psumanti[N // 2 + 1 : N, wv]
            )
            psumsym[N // 2 + 1 : N, wv] = self.wk_smooth121(
                psumsym[N // 2 + 1 : N, wv]
            )
        # -------------------------------------------------------------------
        #  Log10 scaling
        # -------------------------------------------------------------------
        psumanti_nolog = np.ma.masked_array(psumanti)
        psumsym_nolog = np.ma.masked_array(psumsym)

        psumanti = np.ma.log10(psumanti)
        psumsym = np.ma.log10(psumsym)

        # Creating Iris cube, assigning metadata
        psumanti_cube = self.make_cube(psumanti, wave, freq)
        psumsym_cube = self.make_cube(psumsym, wave, freq)

        # -----------------------------------------------------------------------------
        #  ******  now derive and plot the background spectrum (red noise) ************
        #  [1] Sum power over all latitude
        #  [2] Put fill value in mean
        #  [3] Apply smoothing to the spectrum. This smoothing DOES include
        #      wavenumber zero.
        # -----------------------------------------------------------------------------
        # print("======> BACKGROUND <=====")

        psumb = self.compute_background(peeAS, wave, freq, indStrt, indLast)
        psumb_nolog = np.ma.masked_array(psumb)
        psumb = np.ma.log10(psumb)
        psumb_cube = self.make_cube(psumb, wave, freq)

        # -------------------------------------------------------------------------------
        #  Plot section

        # Generate dispersion cuves
        Afreq, Apzwn = self.generate_dispersion_curves()

        # Fig.1 - Raw spectra Symmetric and Anti-symmetric
        #
        # Define contour levels for plots
        levels_dict = {
            "toa_outgoing_longwave_flux": np.array(
                [
                    -1.3,
                    -1.2,
                    -1.1,
                    -1,
                    -0.8,
                    -0.6,
                    -0.4,
                    -0.2,
                    0.0,
                    0.2,
                    0.4,
                    0.6,
                    0.8,
                    1.0,
                    1.1,
                    1.2,
                    1.3,
                ]
            ),
            "Precipitation": np.array(
                [-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
            ),
            "x_wind_850hPa": np.arange(-3.25, 0.5, 0.25),
            "x_wind_200hPa": np.arange(-3.3, 1.2, 0.3),
        }

        # Anti-symmetric
        title = f"{self.label}_{self.varname} \n  Anti-symmetric log(power) [15S-15N]"
        forename = f"{self.runid}_{self.varname}_Raw_Spec_Asym"
        figname = os.path.join(self.plot_dir, f"{forename}.png")
        self.plot_anti_symmetric(
            psumanti,
            freq,
            wave,
            Apzwn,
            Afreq,
            levels=levels_dict[self.varname],
            title=title,
            figname=figname,
        )
        ncname = os.path.join(self.work_dir, f"{forename}.nc")
        iris.save(psumanti_cube, ncname)

        # Symmetric
        title = (
            f"{self.label}_{self.varname} \n Symmetric log(power) [15S-15N]"
        )
        forename = f"{self.runid}_{self.varname}_Raw_Spec_Sym"
        figname = os.path.join(self.plot_dir, f"{forename}.png")
        self.plot_symmetric(
            psumsym,
            freq,
            wave,
            Apzwn,
            Afreq,
            levels=levels_dict[self.varname],
            title=title,
            figname=figname,
        )
        ncname = os.path.join(self.work_dir, f"{forename}.nc")
        iris.save(psumsym_cube, ncname)

        # Background spectra
        title = f"{self.label} {self.varname} \n Background power log(power) [15S-15N]"
        forename = f"{self.runid}_{self.varname}_BG_Spec"
        figname = f"{forename}.png"
        ncname = os.path.join(self.work_dir, f"{forename}.nc")
        iris.save(psumb_cube, ncname)

        # *************************************************************
        #  Fig 3a, 3b:  psum_nolog/psumb_nolog  [ratio]
        # ***************************************************************
        psumanti_nolog = np.ma.masked_array(psumanti_nolog / psumb_nolog)
        psumsym_nolog = np.ma.masked_array(
            psumsym_nolog / psumb_nolog
        )  # (wave,freq)
        # Make cubes
        psumanti_nolog_cube = self.make_cube(psumanti_nolog, wave, freq)
        psumsym_nolog_cube = self.make_cube(psumsym_nolog, wave, freq)

        # Anti-symmetric
        # Define contour levels for plots
        levels_dict = {
            "toa_outgoing_longwave_flux": np.array(
                [
                    0.2,
                    0.3,
                    0.4,
                    0.5,
                    0.6,
                    0.7,
                    0.8,
                    0.9,
                    1.0,
                    1.1,
                    1.2,
                    1.3,
                    1.4,
                    1.5,
                    1.6,
                    1.7,
                    1.8,
                ]
            ),
            "Precipitation": np.array(
                [
                    0.5,
                    0.6,
                    0.7,
                    0.8,
                    0.9,
                    1.0,
                    1.1,
                    1.15,
                    1.2,
                    1.25,
                    1.3,
                    1.35,
                    1.4,
                    1.45,
                    1.5,
                    1.6,
                    1.7,
                ]
            ),
            "x_wind_850hPa": np.array(
                [
                    0.3,
                    0.4,
                    0.5,
                    0.6,
                    0.7,
                    0.8,
                    0.9,
                    1.0,
                    1.1,
                    1.2,
                    1.3,
                    1.4,
                    1.5,
                    1.6,
                    1.7,
                    1.8,
                    1.9,
                ]
            ),
            "x_wind_200hPa": np.array(
                [
                    0.3,
                    0.4,
                    0.5,
                    0.6,
                    0.7,
                    0.8,
                    0.9,
                    1.0,
                    1.1,
                    1.2,
                    1.3,
                    1.4,
                    1.5,
                    1.6,
                    1.7,
                    1.8,
                    2,
                ]
            ),
        }

        title = f"{self.label} {self.varname} \n Anti-symmetric/Background log(power) [15S-15N]"
        forename = f"{self.runid}_{self.varname}_Ratio_Spec_Asym"
        figname = os.path.join(self.plot_dir, f"{forename}.png")
        self.plot_anti_symmetric(
            psumanti_nolog,
            freq,
            wave,
            Apzwn,
            Afreq,
            levels=levels_dict[self.varname],
            title=title,
            figname=figname,
        )
        ncname = os.path.join(self.work_dir, f"{forename}.nc")
        iris.save(psumanti_nolog_cube, ncname)

        # Symmetric
        # Define contour levels for plots
        levels_dict = {
            "toa_outgoing_longwave_flux": np.array(
                [
                    0.2,
                    0.3,
                    0.4,
                    0.5,
                    0.6,
                    0.7,
                    0.8,
                    0.9,
                    1.0,
                    1.1,
                    1.2,
                    1.4,
                    1.7,
                    2.0,
                    2.4,
                    2.8,
                    3.2,
                ]
            ),
            "Precipitation": np.array(
                [
                    0.5,
                    0.6,
                    0.7,
                    0.8,
                    0.9,
                    1.0,
                    1.1,
                    1.15,
                    1.2,
                    1.25,
                    1.3,
                    1.35,
                    1.4,
                    1.45,
                    1.5,
                    1.6,
                    1.7,
                ]
            ),
            "x_wind_850hPa": np.array(
                [
                    0.2,
                    0.4,
                    0.6,
                    0.8,
                    1.0,
                    1.2,
                    1.3,
                    1.4,
                    1.5,
                    1.6,
                    1.7,
                    1.8,
                    2,
                    2.2,
                    2.4,
                    2.6,
                    2.8,
                ]
            ),
            "x_wind_200hPa": np.array(
                [
                    0.2,
                    0.4,
                    0.6,
                    0.8,
                    1.0,
                    1.2,
                    1.3,
                    1.4,
                    1.5,
                    1.6,
                    1.7,
                    1.8,
                    2,
                    2.2,
                    2.4,
                    2.6,
                    2.8,
                ]
            ),
        }

        title = f"{self.label} {self.varname} \n Symmetric/Background log(power) [15S-15N]"
        forename = f"{self.runid}_{self.varname}_Ratio_Spec_Sym"
        figname = os.path.join(self.plot_dir, f"{forename}.png")
        self.plot_symmetric(
            psumsym_nolog,
            freq,
            wave,
            Apzwn,
            Afreq,
            levels=levels_dict[self.varname],
            title=title,
            figname=figname,
        )
        ncname = os.path.join(self.work_dir, f"{forename}.nc")
        iris.save(psumsym_nolog_cube, ncname)

    def mjo_wavenum_freq_season(self, seaName):
        #
        #  For each 'seaName' (say, 'winter') over a number
        #  of different years, calculate the
        #  pooled space-time spectra
        #  MJO CLIVAR: wavenumber-frequency spectra
        #              winter starts Nov 1  [180 days]
        #              summer starts May 1  [180 days]
        #              annual starts Jan 1  [365 days]
        latBound = 10
        var = self.cube.intersection(
            latitude=(-latBound, latBound), longitude=(0, 360)
        )
        var = var.collapsed("latitude", iris.analysis.MEAN)

        ntim, mlon = var.shape
        time = var.coord("time")
        years, months, days = self.get_dates(time)
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

        iSea = [i for i in range(len(mmdd)) if mmdd[i] == mmddStrt]
        nYear = len(iSea)

        """
        # *****************************************************************
        #  For a specific season, calculate spectra via averaging
        #  over each seasonal segment.
        #  MJO Clivar says "no" to detrending/tapering.
        #  Hence, the following are just 'place holders'
        # *****************************************************************
        """
        # detrend overall series in time
        # Time mean (later to be added to the trend)
        varmean = var.collapsed("time", iris.analysis.MEAN)

        # remove linear trend
        var.data = scipy.signal.detrend(
            var.data, axis=0
        )  # + varmean.data # Mean added

        # Initialise
        power = np.zeros([mlon + 1, nDay + 1])  # initialize
        work = var.data
        xAvgSea = 0.0
        xVarSea = 0.0  # variance (raw)
        xVarTap = 0.0  # variance after tapering
        kSea = 0  # count of seasons used
        N = nDay  # convenience

        for ny in range(nYear):
            iStrt = iSea[ny]  # start index for current season
            iLast = iSea[ny] + nDay  # last
            if iLast < ntim - 1:
                print(ny, kSea, iStrt, iLast, yyyymmdd[iStrt], yyyymmdd[iLast])
                xSeason = work[iStrt:iLast, :]
                xAvg = np.average(xSeason)  # season average all time/lon
                xSeason = xSeason - xAvg  # remove season time-lon mean
                xVarSea = xVarSea + np.var(xSeason)  # overall variance
                kSea = kSea + 1

                for lo in range(mlon):
                    xSeason[:, lo] = self.taper(
                        xSeason[:, lo], alpha=0.1, iopt=0
                    )
                xVarTap = xVarTap + np.var(xSeason)  # variance after tapering
                # do 2d fft
                ft = xSeason.copy()
                ft = np.fft.fft2(xSeason.T) / mlon / nDay

                # Shifting FFTs
                power = power + self.resolve_waves_hayashi(
                    ft, nDay, spd=1
                )  # (wave,freq)

        xVarSea = xVarSea / kSea  # pooled seasonal variance
        xVarTap = xVarTap / kSea  # pooled seasonal variance
        power = np.ma.masked_array(power / kSea)  # pooled spectra

        wave = np.arange(-mlon / 2, mlon / 2 + 1, 1)
        freq = np.linspace(-1 * nDay / 2, nDay / 2, nDay + 1) / nDay
        wave = wave.astype(float)
        freq = freq.astype(float)
        pow_cube = self.makecube_season_pow(power, wave, freq)
        return pow_cube

    def mjo_wavenum_freq_season_plot(
        self,
        pow_cube,
        levels=np.arange(0, 3, 0.2),
        title="",
        figname="wavenum_freq_season_plot.ps",
    ):
        NW = 6
        fMin = -0.05
        fMax = 0.05

        constraint = iris.Constraint(
            frequency=lambda cell: fMin <= cell <= fMax,
            wavenumber=lambda cell: 0 <= cell <= NW,
        )
        pow_cube = pow_cube.extract(constraint)

        freq = pow_cube.coord("frequency").points
        wave = pow_cube.coord("wavenumber").points

        W, F = np.meshgrid(freq, wave)

        # set zeroth frequency to minimum value
        izero = np.where(freq == 0)[0][0]
        pow_cube.data[:, izero] = np.min(pow_cube.data)  # 0th freq

        # Initialize the plot
        fig, ax = plt.subplots()

        CS = ax.contourf(
            W, F, pow_cube.data, levels=levels, cmap="YlGnBu", extend="both"
        )
        fig.colorbar(CS, ax=ax)

        # set axes range
        ax.set_ylim(0, NW)
        ax.set_xlim(fMin, fMax)

        # Line markers of periods
        frqs = [80, 30]
        for frq in frqs:
            ax.plot([1.0 / frq, 1.0 / frq], [-0, 15], "k--", lw=0.5)
            ax.text(1.0 / frq, 5.5, str(frq) + "d", {"color": "k"})

        ax.plot([0, 0], [-0, 15], "k:", lw=0.25)
        ax.set_title(title)
        ax.set_xlabel("Westward     Frequency     Eastward")
        ax.set_ylabel("Zonal wavenumber")

        # Add provenance information
        caption = f"{figname}, [or other caption for wavenum freq season plot]"  # TODO
        provenance_dict = self.get_provenance_record(caption)

        # Save the figure (also closes it)
        save_figure(
            figname,
            provenance_dict,
            self.cfg,
            figure=fig,
            close=True,
        )
        print(f"Plotted {figname}")

    def SpectraSeason(self):
        for season in ["winter", "summer"]:
            print(season)
            # This is the NCL method of computing the spectra which uses anomalies
            pow_cube = self.mjo_wavenum_freq_season(season)

            forename = (
                f"{self.runid}_{self.varname}_wavenum_freq_season_{season}"
            )
            ncname = os.path.join(self.work_dir, f"{forename}.nc")
            iris.save(pow_cube, ncname)

            # Define contour levels for plots
            levels_dict = {
                "toa_outgoing_longwave_flux": np.arange(0.0, 2.4, 0.2),
                "Precipitation": np.arange(0.0, 0.055, 0.005),
                "x_wind_850hPa": np.arange(0.007, 0.07, 0.007),
                "x_wind_200hPa": np.arange(0.05, 0.5, 0.05),
            }

            title = f"{self.label} \n {season}  daily {self.varname} [10S-10N]"
            figname = os.path.join(self.plot_dir, f"{forename}.png")
            self.mjo_wavenum_freq_season_plot(
                pow_cube,
                levels=levels_dict[self.varname],
                title=title,
                figname=figname,
            )
