"""Seasonal MJO/BSISO extended-EOF index from daily OLR anomalies.

This diagnostic ports the ACCESS-NRI MJO/BSISO notebook to the ESMValTool
recipe interface. The recipe supplies daily, tropical, day-of-year OLR
anomalies. A 139-tap Hamming-window FIR filter isolates 25–90-day variability.
Lagged fields are built on the continuous daily series before selecting DJF
and JJA training dates. The first two area-weighted EEOFs define separate
MJO and BSISO indices, which are projected onto every available day.

The EOF/PC sign is paired, as in the source notebook. Geographic phase names
are deliberately omitted pending independent sign-convention validation.
"""

import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from scipy.signal import firwin, lfilter
from sklearn.decomposition import PCA

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    group_metadata,
    run_diagnostic,
    save_figure,
)


def _as_daily_array(cube):
    """Convert a preprocessed Iris cube to a continuous time/lat/lon array."""
    data = xr.DataArray.from_iris(cube)
    rename = {
        long_name: short_name
        for long_name, short_name in (
            ("latitude", "lat"),
            ("longitude", "lon"),
        )
        if long_name in data.dims
    }
    data = data.rename(rename).transpose("time", "lat", "lon")
    times = data.time.values
    if len(times) < 200:
        raise ValueError("At least 200 daily values are required")
    if any(
        (right - left) != np.timedelta64(1, "D")
        for left, right in zip(times[:-1], times[1:])
    ):
        raise ValueError("Input time coordinate must be continuous and daily")
    values = np.asarray(data.values, dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError(
            "Input OLR anomalies contain missing or nonfinite data"
        )
    return values, times, data.lat.values, data.lon.values


def _bandpass(values, times, low_period, high_period, window):
    """Reproduce the notebook's centered Hamming FIR and trim both edges."""
    if window < 3 or window % 2 != 1:
        raise ValueError("filter_window must be an odd integer >= 3")
    if not 0 < low_period < high_period:
        raise ValueError("Expected 0 < low_period < high_period")
    if len(times) <= window:
        raise ValueError("Time series is shorter than the FIR window")
    weights = firwin(
        window,
        [1.0 / high_period, 1.0 / low_period],
        pass_zero=False,
        fs=1.0,
    )
    half = window // 2
    causal = lfilter(weights, 1.0, values, axis=0)
    return causal[2 * half :], times[half:-half]


def _lagged_matrix(values, times, lags, latitudes):
    """Build weighted lag fields on continuous days for later seasonal fits."""
    lags = np.asarray(lags, dtype=int)
    if not len(lags):
        raise ValueError("lags must contain at least one day offset")
    if len(np.unique(lags)) != len(lags):
        raise ValueError("lags must be distinct")
    first = max(0, -int(lags.min()))
    stop = len(times) - max(0, int(lags.max()))
    if stop <= first:
        raise ValueError("Not enough dates for the requested lags")
    centers = np.arange(first, stop)
    fields = np.stack([values[centers + lag] for lag in lags], axis=1)
    latitude_weights = np.sqrt(np.cos(np.deg2rad(latitudes)))
    fields *= latitude_weights[np.newaxis, np.newaxis, :, np.newaxis]
    return fields.reshape(len(centers), -1), times[centers]


def _fit_season(matrix, times, months, spatial_shape, flip_pc2=False):
    """Fit two normalized seasonal PCs and project all available days."""
    season = np.isin([time.month for time in times], months)
    if season.sum() < 3:
        raise ValueError(f"Too few training days for months {months}")
    pca = PCA(n_components=2)
    raw = pca.fit_transform(matrix[season])
    standard_deviation = raw.std(axis=0)
    if np.any(standard_deviation == 0):
        raise ValueError("A seasonal PC has zero standard deviation")
    training = raw / standard_deviation
    projected = pca.transform(matrix) / standard_deviation
    eeofs = (
        pca.components_ * standard_deviation[:, np.newaxis]
    ).reshape((2, *spatial_shape))
    if flip_pc2:
        training[:, 1] *= -1
        projected[:, 1] *= -1
        eeofs[1] *= -1
    np.testing.assert_allclose(
        projected[season], training, rtol=1e-5, atol=1e-5
    )
    return {
        "training": training,
        "projected": projected,
        "eeofs": eeofs,
        "variance": pca.explained_variance_ratio_ * 100,
        "season": season,
    }


def calculate_indices(cube, cfg):
    """Calculate the EEOF patterns, PCs, amplitudes, and monthly frequency."""
    values, times, latitudes, longitudes = _as_daily_array(cube)
    low_period = cfg.get("low_period", 25)
    high_period = cfg.get("high_period", 90)
    window = cfg.get("filter_window", 139)
    lags = cfg.get("lags", [-10, -5, 0])
    filtered, filtered_times = _bandpass(
        values, times, low_period, high_period, window
    )
    del values
    matrix, dates = _lagged_matrix(filtered, filtered_times, lags, latitudes)
    del filtered
    shape = (len(lags), len(latitudes), len(longitudes))
    mjo = _fit_season(matrix, dates, [12, 1, 2], shape)
    bsiso = _fit_season(matrix, dates, [6, 7, 8], shape, flip_pc2=True)
    del matrix

    amplitude_mjo = np.linalg.norm(mjo["projected"], axis=1)
    amplitude_bsiso = np.linalg.norm(bsiso["projected"], axis=1)
    months = np.asarray([date.month for date in dates])
    dominant_mjo = (amplitude_mjo >= 1) & (amplitude_mjo > amplitude_bsiso)
    dominant_bsiso = (amplitude_bsiso >= 1) & (
        amplitude_bsiso >= amplitude_mjo
    )
    days_in_month = np.asarray(
        [31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    )
    mjo_days = np.asarray(
        [dominant_mjo[months == month].mean() for month in range(1, 13)]
    ) * days_in_month
    bsiso_days = np.asarray(
        [dominant_bsiso[months == month].mean() for month in range(1, 13)]
    ) * days_in_month

    result = xr.Dataset(
        data_vars={
            "eeof_mjo": (("mode", "lag", "lat", "lon"), mjo["eeofs"]),
            "eeof_bsiso": (("mode", "lag", "lat", "lon"), bsiso["eeofs"]),
            "pc_mjo_training": (
                ("time_mjo", "mode"),
                mjo["training"],
            ),
            "pc_bsiso_training": (
                ("time_bsiso", "mode"),
                bsiso["training"],
            ),
            "pc_mjo": (("time", "mode"), mjo["projected"]),
            "pc_bsiso": (("time", "mode"), bsiso["projected"]),
            "amplitude_mjo": ("time", amplitude_mjo),
            "amplitude_bsiso": ("time", amplitude_bsiso),
            "dominant_mjo": ("time", dominant_mjo.astype(np.int8)),
            "dominant_bsiso": ("time", dominant_bsiso.astype(np.int8)),
            "mjo_days": ("month", mjo_days),
            "bsiso_days": ("month", bsiso_days),
            "variance_mjo": ("mode", mjo["variance"]),
            "variance_bsiso": ("mode", bsiso["variance"]),
        },
        coords={
            "mode": [1, 2],
            "lag": lags,
            "lat": latitudes,
            "lon": longitudes,
            "time": dates,
            "time_mjo": dates[mjo["season"]],
            "time_bsiso": dates[bsiso["season"]],
            "month": np.arange(1, 13),
        },
        attrs={
            "description": (
                "Seasonal MJO/BSISO EEOF index from daily OLR anomalies"
            ),
            "filter": (
                f"Hamming FIR, {window} taps, {low_period}-{high_period} days"
            ),
            "training_seasons": "MJO: DJF; BSISO: JJA",
            "pc_normalization": "unit variance over each training season",
            "bsiso_sign": "PC2 and EEOF2 multiplied by -1 together",
        },
    )
    result["eeof_mjo"].attrs["units"] = str(cube.units)
    result["eeof_bsiso"].attrs["units"] = str(cube.units)
    result["variance_mjo"].attrs["units"] = "%"
    result["variance_bsiso"].attrs["units"] = "%"
    return result


def _season_year(time):
    return np.asarray([date.year + (date.month == 12) for date in time])


def _year_position(time):
    """Return numeric model years, also valid for pre-1678 calendars."""
    return np.asarray(
        [
            date.year + (date.month - 1) / 12 + (date.day - 1) / 365
            for date in time
        ]
    )


def _plot_segmented(ax, time, values, label, **kwargs):
    years = _season_year(time)
    positions = _year_position(time)
    for index, year in enumerate(np.unique(years)):
        selected = years == year
        ax.plot(
            positions[selected],
            values[selected],
            label=label if index == 0 else None,
            **kwargs,
        )


def plot_eeofs(result, label):
    """Plot lagged spatial patterns in EEOF1, EEOF2 order."""
    n_lags = result.sizes["lag"]
    fig, axes = plt.subplots(
        2 * n_lags,
        2,
        figsize=(14, 5 * n_lags),
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)},
    )
    levels = np.linspace(-15, 15, 11)
    for section, prefix in enumerate(("mjo", "bsiso")):
        patterns = result[f"eeof_{prefix}"]
        variance = result[f"variance_{prefix}"].values
        for lag_index, lag in enumerate(result.lag.values):
            for mode_index in range(2):
                ax = axes[section * n_lags + lag_index, mode_index]
                field = patterns.isel(mode=mode_index, lag=lag_index)
                contour = ax.contourf(
                    result.lon.values,
                    result.lat.values,
                    field.values,
                    levels=levels,
                    cmap="RdBu_r",
                    extend="both",
                    transform=ccrs.PlateCarree(),
                )
                ax.coastlines(color="dimgray", linewidth=0.7)
                ax.set_extent([30, 240, -25, 30], crs=ccrs.PlateCarree())
                ax.text(
                    0.98, 0.82, f"day {lag}",
                    ha="right", transform=ax.transAxes,
                )
                if lag_index == 0:
                    ax.set_title(
                        f"{prefix.upper()} EEOF{mode_index + 1} "
                        f"({variance[mode_index]:.1f}%)"
                    )
    fig.colorbar(
        contour,
        ax=axes.ravel().tolist(),
        orientation="horizontal",
        fraction=0.025,
        pad=0.04,
        label=f"OLR anomaly ({result.eeof_mjo.attrs['units']})",
    )
    fig.suptitle(f"{label}: seasonal MJO and BSISO EEOFs", y=0.995)
    return fig


def plot_pc_timeseries(result, label):
    """Plot seasonal PCs and amplitudes without connecting seasonal gaps."""
    fig, axes = plt.subplots(3, 1, figsize=(14, 12))
    for ax, prefix, color in zip(
        axes[:2], ("mjo", "bsiso"), ("royalblue", "indianred")
    ):
        time = result[f"time_{prefix}"].values
        pcs = result[f"pc_{prefix}_training"].values
        _plot_segmented(ax, time, pcs[:, 0], "PC1", color=color, lw=1)
        _plot_segmented(ax, time, pcs[:, 1], "PC2", color="darkorange", lw=1)
        ax.axhline(0, color="black", ls="--", lw=0.8)
        ax.set_title(f"{prefix.upper()} seasonal PCs")
        ax.set_ylabel("Normalized PC")
        ax.legend()
    for prefix, color in (("mjo", "navy"), ("bsiso", "darkorange")):
        time = result[f"time_{prefix}"].values
        amplitude = np.linalg.norm(
            result[f"pc_{prefix}_training"].values, axis=1
        )
        _plot_segmented(
            axes[2], time, amplitude,
            f"{prefix.upper()} amplitude", color=color,
        )
    axes[2].axhline(1, color="red", ls=":", label="Amplitude threshold")
    axes[2].set_ylabel("Amplitude")
    axes[2].set_xlabel("Model year")
    axes[2].legend()
    fig.suptitle(f"{label}: seasonal principal components")
    fig.tight_layout()
    return fig


def plot_occurrence(result, label):
    """Plot monthly frequency of the predominant projected mode."""
    fig, ax = plt.subplots(figsize=(11, 6))
    month = result.month.values
    ax.bar(
        month, result.mjo_days.values,
        color="royalblue", label="MJO pattern",
    )
    ax.bar(
        month,
        -result.bsiso_days.values,
        color="crimson",
        label="BSISO pattern",
    )
    ax.axhline(0, color="black")
    ax.set_xticks(
        month,
        [
            "Jan", "Feb", "Mar", "Apr", "May", "Jun",
            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
        ],
    )
    ax.set_ylabel("Days per month; BSISO shown below zero")
    ax.set_title(f"{label}: predominant projected ISO mode")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_amplitudes(result, label):
    """Compare year-round MJO and BSISO projected amplitudes."""
    fig, ax = plt.subplots(figsize=(8, 8))
    months = np.asarray([date.month for date in result.time.values])
    for selected_months, color, season in (
        ([6, 7, 8, 9, 10], "red", "Jun–Oct"),
        ([12, 1, 2, 3, 4], "blue", "Dec–Apr"),
        ([5, 11], "gray", "May, Nov"),
    ):
        selected = np.isin(months, selected_months)
        ax.scatter(
            result.amplitude_bsiso.values[selected],
            result.amplitude_mjo.values[selected],
            s=4,
            alpha=0.4,
            color=color,
            label=season,
            rasterized=True,
        )
    ax.axhline(1, color="black", lw=0.8)
    ax.axvline(1, color="black", lw=0.8)
    ax.axline((0, 0), slope=1, color="black", ls="--", lw=0.8)
    ax.set_xlabel("BSISO amplitude")
    ax.set_ylabel("MJO amplitude")
    ax.set_title(f"{label}: projected daily ISO amplitudes")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_phase_space(result, label):
    """Plot seasonal PC phase space without joining separate seasons."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, prefix, cmap in zip(axes, ("mjo", "bsiso"), ("Blues", "Reds")):
        pcs = result[f"pc_{prefix}_training"].values
        times = result[f"time_{prefix}"].values
        years = _season_year(times)
        amplitude = np.linalg.norm(pcs, axis=1)
        for year in np.unique(years):
            selected = years == year
            ax.plot(
                pcs[selected, 0], pcs[selected, 1],
                color="gray", lw=0.7, alpha=0.3,
            )
        scatter = ax.scatter(pcs[:, 0], pcs[:, 1], c=amplitude, cmap=cmap, s=7)
        ax.add_artist(plt.Circle((0, 0), 1, fill=False, color="black"))
        ax.axhline(0, color="black", lw=0.7)
        ax.axvline(0, color="black", lw=0.7)
        ax.axline((0, 0), slope=1, color="black", lw=0.7)
        ax.axline((0, 0), slope=-1, color="black", lw=0.7)
        ax.set(
            xlabel="PC1", ylabel="PC2", title=prefix.upper(),
            xlim=(-5, 5), ylim=(-5, 5),
        )
        ax.set_aspect("equal")
        fig.colorbar(scatter, ax=ax, label="Amplitude")
    fig.suptitle(f"{label}: seasonal phase space (orientation unvalidated)")
    fig.tight_layout()
    return fig


def _provenance(caption, ancestor):
    return {
        "caption": caption,
        "ancestors": [ancestor],
        "authors": ["sullivan_arnold", "chun_felicity", "beucher_romain"],
        "references": ["kikuchi12climdyn"],
        "domains": ["trop"],
        "statistics": ["other"],
        "plot_types": ["geo"],
    }


def main(cfg):
    """Run the index for each preprocessed OLR dataset."""
    groups = group_metadata(
        cfg["input_data"].values(), "alias", sort="dataset"
    )
    for label, entries in groups.items():
        if len(entries) != 1:
            raise ValueError(
                f"Expected one OLR file for {label}, got {len(entries)}"
            )
        source = entries[0]["filename"]
        result = calculate_indices(iris.load_cube(source), cfg)
        basename = f"{label.replace(' ', '_').replace('/', '_')}_mjo_bsiso"
        data_file = get_diagnostic_filename(basename, cfg)
        result.to_netcdf(data_file)
        data_record = _provenance(
            f"{label}: MJO and BSISO EEOF patterns, PCs, amplitudes, "
            "and occurrence",
            source,
        )
        with ProvenanceLogger(cfg) as logger:
            logger.log(data_file, data_record)

        figures = {
            "eeofs": plot_eeofs,
            "pc_timeseries": plot_pc_timeseries,
            "occurrence": plot_occurrence,
            "amplitudes": plot_amplitudes,
            "phase_space": plot_phase_space,
        }
        for name, plot_function in figures.items():
            figure = plot_function(result, label)
            record = _provenance(
                f"{label}: {name.replace('_', ' ')} of the seasonal "
                "MJO/BSISO index",
                source,
            )
            save_figure(
                f"{basename}_{name}",
                record,
                cfg,
                figure=figure,
                dpi=150,
                bbox_inches="tight",
            )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
