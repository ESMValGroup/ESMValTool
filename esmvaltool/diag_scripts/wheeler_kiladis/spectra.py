"""
Core Wheeler-Kiladis spectral utilities.

This module provides functions for preprocessing equatorial longitude-time
fields and computing Wheeler-Kiladis (WK) wavenumber-frequency spectra
following Wheeler and Kiladis (1999).

The workflow consists of:

1. Spatial preprocessing
   - Regularize longitude and latitude coordinates.
   - Compute a cosine-weighted equatorial mean.

2. Temporal preprocessing
   - Convert data to daily means.
   - Remove leap days.
   - Remove linear trends.
   - Remove the annual cycle using harmonic regression.
   - Optionally remove the zonal mean.

3. Spectral analysis
   - Split the data into overlapping segments.
   - Compute 2D FFT-based wavenumber-frequency spectra.
   - Average spectra across segments.
   - Estimate a smoothed background spectrum.
   - Compute normalized spectra.

Diagnostic options
------------------

remove_zonal_mean : bool, optional
    Remove the zonal mean at each time step prior to spectral analysis.
    This suppresses purely zonally symmetric variability and emphasizes
    propagating disturbances. The default is True because Wheeler-Kiladis
    spectra are typically interpreted in terms of eastward- and westward-
    propagating equatorial wave modes.

nharmonics : int, optional
    Number of annual harmonics removed during seasonal-cycle correction.
    The default value of 3 removes the annual cycle and its first two
    higher harmonics, which generally captures the dominant seasonal
    variability without removing intraseasonal signals of interest.

max_missing_fraction : float, optional
    Maximum fraction of missing data allowed before the diagnostic aborts.
    The default threshold is intentionally conservative because Fourier-
    based spectral analyses are sensitive to data gaps. Small isolated
    gaps can be interpolated with minimal impact, whereas extensive
    missing coverage can distort the spectral estimates.

segment_length : int, optional
    Length of the time segments used to compute individual spectra.
    The default value of 180 days follows common Wheeler-Kiladis
    implementations and provides a balance between frequency resolution
    and the number of independent realizations available for averaging.

segment_overlap : int, optional
    Overlap between adjacent time segments. The default value of
    90 days corresponds to 50% overlap, increasing the number of
    segments used in the averaging while maintaining substantial
    independence between consecutive estimates.

samples_per_day : float, optional
    Temporal sampling frequency of the preprocessed data. The current
    implementation assumes daily data and therefore uses a default
    value of 1.0 samples day^-1.

sigma_freq : float, optional
    Standard deviation of the Gaussian filter applied in the frequency
    direction when estimating the background spectrum. The default
    value of 4.0 provides sufficient smoothing to represent the broad
    spectral background while retaining the large-scale structure of
    the spectrum.

sigma_wn : float, optional
    Standard deviation of the Gaussian filter applied in the zonal
    wavenumber direction when estimating the background spectrum.
    The default value of 4.0 was chosen to produce a smoothly varying
    background spectrum suitable for normalization while preserving
    the large-scale spectral envelope.

Notes
-----

The default configuration follows the methodology commonly used for
Wheeler-Kiladis diagnostics, with daily sampling, overlapping
180-day segments, seasonal-cycle removal through harmonic regression,
and normalization by a smoothed background spectrum. These choices
are intended to provide robust identification of convectively coupled
equatorial wave signals while remaining computationally efficient.

Outputs
-------

raw_power
    Unsmoothed Wheeler-Kiladis power spectrum.

background_power
    Smoothed background spectrum.

normalized_power
    Power spectrum normalized by the background spectrum.

References
----------

Wheeler, M. C. and Kiladis, G. N. (1999):
Convectively Coupled Equatorial Waves: Analysis of Clouds and
Temperature in the Wavenumber-Frequency Domain.
Journal of the Atmospheric Sciences, 56, 374-399.
"""

import logging
from pathlib import Path

import numpy as np
import xarray as xr
from scipy import signal
from scipy.ndimage import gaussian_filter

logger = logging.getLogger(__name__)


def prepare_regular_lon_lat(data):
    """Sort regular lat-lon data and remove duplicate cyclic longitudes."""
    rename_dict = {}

    if "latitude" in data.dims:
        rename_dict["latitude"] = "lat"
    if "longitude" in data.dims:
        rename_dict["longitude"] = "lon"

    if rename_dict:
        data = data.rename(rename_dict)

    data = data.transpose("time", "lat", "lon")

    if data.lat[0] > data.lat[-1]:
        data = data.sortby("lat")

    if data.lon[0] > data.lon[-1]:
        data = data.sortby("lon")

    lon = data.lon.values

    if len(lon) > 1 and np.isclose(abs(lon[-1] - lon[0]), 360.0):
        data = data.isel(lon=slice(0, -1))

    lon_mod = np.mod(data.lon.values, 360.0)
    _, unique_idx = np.unique(np.round(lon_mod, 8), return_index=True)
    unique_idx = np.sort(unique_idx)

    if len(unique_idx) != data.sizes["lon"]:
        data = data.isel(lon=unique_idx)

    return data


def equatorial_mean_regular(data):
    """Create cosine-weighted tropical mean time-lon field."""
    data = prepare_regular_lon_lat(data)

    weights = np.cos(np.deg2rad(data.lat))
    weights = weights / weights.mean()

    out = data.weighted(weights).mean("lat")
    out = out.transpose("time", "lon")

    return out


def remove_feb29(data):
    """Remove leap days."""
    is_feb29 = (data.time.dt.month == 2) & (data.time.dt.day == 29)
    return data.sel(time=~is_feb29)


def detrend_xarray(data, dim="time"):
    """Remove linear trend along a dimension."""
    axis = data.get_axis_num(dim)
    values = signal.detrend(data.values, axis=axis, type="linear")

    return xr.DataArray(
        values,
        dims=data.dims,
        coords=data.coords,
        attrs=data.attrs,
        name=data.name,
    )


def remove_annual_harmonics(data, nharmonics=3, samples_per_year=365.0):
    """Remove annual cycle using harmonic regression."""
    data = data.transpose("time", "lon")

    y = data.values
    nt, _ = y.shape
    t = np.arange(nt, dtype=float)

    columns = [np.ones(nt)]

    for n in range(1, nharmonics + 1):
        columns.append(np.sin(2.0 * np.pi * n * t / samples_per_year))
        columns.append(np.cos(2.0 * np.pi * n * t / samples_per_year))

    design = np.stack(columns, axis=1)
    beta, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    seasonal = design @ beta
    anomalies = y - seasonal

    return xr.DataArray(
        anomalies,
        dims=data.dims,
        coords=data.coords,
        attrs=data.attrs,
        name=data.name,
    )


def preprocess_equatorial_field(
    data,
    remove_zonal_mean=True,
    nharmonics=3,
    max_missing_fraction=0.01,
):
    """Create anomaly field for Wheeler-Kiladis analysis.

    The Wheeler-Kiladis spectral decomposition requires a complete,
    regularly sampled longitude-time field. Some observational datasets
    (e.g. PERSIANN-CDR) may contain isolated missing values that cause
    scipy.signal.detrend() to fail and may subsequently contaminate the
    spectral analysis.

    To preserve the full record while minimizing spectral distortion:

    1. Remove leap days.
    2. Check the fraction of missing values and abort if data coverage
       is poor.
    3. Fill isolated gaps using linear interpolation along the time axis
       only. This maintains a regular sampling interval and avoids the
       artificial discontinuities that would be introduced by filling
       with zeros or nearest-neighbour values.
    4. Remove the temporal mean.
    5. Remove the linear trend.
    6. Remove the annual cycle represented by the first few harmonics.
    7. Optionally remove the zonal mean.

    The interpolation step is intended only for small data gaps.
    Datasets with substantial missing coverage should be rejected
    rather than reconstructed.

    Parameters
    ----------
    data : xarray.DataArray
        Daily longitude-time precipitation or OLR field.
    remove_zonal_mean : bool, optional
        Remove zonal mean before spectral analysis.
    nharmonics : int, optional
        Number of annual harmonics to remove.
    max_missing_fraction : float, optional
        Maximum allowed fraction of missing values before raising an
        error. Default is 1%.

    Returns
    -------
    xarray.DataArray
        Preprocessed anomaly field suitable for WK analysis.
    """

    data = remove_feb29(data)

    # Reject datasets with substantial missing coverage. WK spectra
    # are sensitive to large gaps and these should be treated as a
    # data-quality issue rather than filled automatically.
    missing_fraction = float(data.isnull().sum()) / data.size

    if missing_fraction > max_missing_fraction:
        raise ValueError(
            f"Missing-value fraction ({100 * missing_fraction:.2f}%) "
            f"exceeds allowed threshold "
            f"({100 * max_missing_fraction:.2f}%)."
        )

    # Fill isolated gaps using temporal interpolation. This preserves
    # regular time sampling required by detrending and FFT operations.
    if data.isnull().any():
        n_missing = int(data.isnull().sum())

        logger.warning(
            "Interpolating %d missing values prior to detrending.",
            n_missing,
        )

        data = data.interpolate_na(
            dim="time",
            method="linear",
            fill_value="extrapolate",
        )

    # Remove long-term mean
    data = data - data.mean("time")

    # Remove linear trend
    data = detrend_xarray(data, dim="time")

    # Remove seasonal cycle
    data = remove_annual_harmonics(
        data,
        nharmonics=nharmonics,
    )

    # Remove zonal mean if desired
    if remove_zonal_mean:
        data = data - data.mean("lon")

    return data


def segment_time(data, nperseg=180, noverlap=90, time_dim="time"):
    """Split time-lon data into overlapping segments."""
    step = nperseg - noverlap
    nt = data.sizes[time_dim]

    if step <= 0:
        raise ValueError("noverlap must be smaller than nperseg")

    if nt < nperseg:
        raise ValueError(
            f"Not enough time steps for nperseg={nperseg}. Available nt={nt}."
        )

    segments = []

    for start in range(0, nt - nperseg + 1, step):
        segment = data.isel({time_dim: slice(start, start + nperseg)})
        segments.append(segment)

    return segments


def wk_spectrum_one_segment(
    segment,
    samples_per_day=1.0,
    apply_time_window=True,
):
    """Compute wavenumber-frequency power spectrum for one time-lon segment."""
    segment = segment.transpose("time", "lon")

    arr = segment.values.astype(float)
    nt, nlon = arr.shape

    arr = arr - np.nanmean(arr, axis=0, keepdims=True)
    arr = arr - np.nanmean(arr, axis=1, keepdims=True)
    arr = np.nan_to_num(arr)

    if apply_time_window:
        window = np.hanning(nt)[:, None]
        arr = arr * window
        window_norm = np.mean(window[:, 0] ** 2)
    else:
        window_norm = 1.0

    fft = np.fft.fft2(arr, axes=(0, 1))
    power = np.abs(fft) ** 2
    power = power / (nt * nlon * window_norm)

    freqs = np.fft.fftfreq(nt, d=1.0 / samples_per_day)
    wavenumbers_raw = np.fft.fftfreq(nlon, d=1.0 / nlon)

    power = np.fft.fftshift(power, axes=1)
    wavenumbers = np.fft.fftshift(wavenumbers_raw)

    # Positive wavenumber means eastward propagation.
    wavenumbers = -wavenumbers

    order = np.argsort(wavenumbers)
    wavenumbers = wavenumbers[order]
    power = power[:, order]

    positive_freq = freqs > 0
    freqs = freqs[positive_freq]
    power = power[positive_freq, :]

    return xr.DataArray(
        power,
        dims=("frequency", "wavenumber"),
        coords={
            "frequency": freqs,
            "wavenumber": wavenumbers,
        },
        name="power",
    )


def average_wk_spectra(segments, samples_per_day=1.0):
    """Compute average WK spectrum over segments."""
    spectra = [
        wk_spectrum_one_segment(segment, samples_per_day=samples_per_day)
        for segment in segments
    ]

    return xr.concat(spectra, dim="segment").mean("segment")


def smooth_background_spectrum(
    power,
    sigma_freq=4.0,
    sigma_wn=4.0,
    floor_quantile=0.05,
):
    """Estimate smooth WK background spectrum in log space."""
    arr = power.values.astype(float)
    valid = np.isfinite(arr) & (arr > 0)

    if not np.any(valid):
        raise ValueError("No positive finite values found in spectrum.")

    floor = np.quantile(arr[valid], floor_quantile)

    arr_safe = np.where(valid, arr, floor)
    arr_safe = np.where(arr_safe > floor, arr_safe, floor)

    smooth_log = gaussian_filter(
        np.log(arr_safe),
        sigma=(sigma_freq, sigma_wn),
        mode="nearest",
    )

    background = np.exp(smooth_log)
    background = np.where(background > floor, background, floor)

    return xr.DataArray(
        background,
        dims=power.dims,
        coords=power.coords,
        name="background_power",
    )


def compute_wk_spectra(
    data,
    segment_length=180,
    segment_overlap=90,
    samples_per_day=1.0,
    sigma_freq=4.0,
    sigma_wn=4.0,
):
    """Compute raw, background, and normalized WK spectra."""
    segments = segment_time(
        data,
        nperseg=segment_length,
        noverlap=segment_overlap,
    )

    raw = average_wk_spectra(
        segments,
        samples_per_day=samples_per_day,
    )

    background = smooth_background_spectrum(
        raw,
        sigma_freq=sigma_freq,
        sigma_wn=sigma_wn,
    )

    normalized = raw / background
    normalized.name = "normalized_power"

    return raw, background, normalized


def mask_zero_wavenumber(spectrum):
    """Mask exact k=0 column."""
    k0 = spectrum.wavenumber.sel(wavenumber=0, method="nearest")
    return spectrum.where(spectrum.wavenumber != k0)


def save_spectra(raw, background, normalized, output_dir, basename):
    """Save spectra to NetCDF files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_file = output_dir / f"{basename}_raw_power.nc"
    background_file = output_dir / f"{basename}_background_power.nc"
    normalized_file = output_dir / f"{basename}_normalized_power.nc"

    raw.to_netcdf(raw_file)
    background.to_netcdf(background_file)
    normalized.to_netcdf(normalized_file)

    return raw_file, background_file, normalized_file
