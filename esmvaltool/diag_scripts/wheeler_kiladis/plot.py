"""Plot utilities for Wheeler-Kiladis diagnostic."""

import matplotlib.pyplot as plt
import numpy as np


def kelvin_curve(wavenumber, depth_m):
    """
    Kelvin wave dispersion curve in cycles per day.

    Parameters
    ----------
    wavenumber : array-like
        Dimensionless zonal wavenumber.
    depth_m : float
        Equivalent depth [m].
    """
    gravity = 9.80665  # m s^-2
    earth_radius = 6.371e6  # m

    phase_speed = np.sqrt(gravity * depth_m)  # m s^-1
    zonal_wavenumber = np.asarray(wavenumber) / earth_radius  # m^-1

    angular_frequency = phase_speed * zonal_wavenumber  # rad s^-1
    frequency = angular_frequency / (2.0 * np.pi) * 86400.0  # cycles day^-1

    return frequency


def equatorial_rossby_curve(wavenumber, depth_m, meridional_mode=1):
    """
    Compute the dispersion relation for equatorial Rossby waves.

    Parameters
    ----------
    wavenumber : array-like
        Dimensionless zonal wavenumber.
    depth_m : float
        Equivalent depth [m].
    meridional_mode : int, optional
        Meridional mode number of the equatorial Rossby wave.
        Default is 1.

    Returns
    -------
    numpy.ndarray
        Wave frequency [cycles day^-1] corresponding to the input
        zonal wavenumbers.
    """
    earth_rotation = 7.2921159e-5  # rad s^-1
    earth_radius = 6.371e6  # m
    gravity = 9.80665  # m s^-2

    beta = 2.0 * earth_rotation / earth_radius  # m^-1 s^-1
    phase_speed = np.sqrt(gravity * depth_m)  # m s^-1

    wavenumber = np.asarray(wavenumber)  # dimensionless zonal wavenumber
    zonal_wavenumber = wavenumber / earth_radius  # m^-1

    denominator = (
        zonal_wavenumber**2
        + (2 * meridional_mode + 1) * beta / phase_speed  # m^-2
    )

    angular_frequency = -beta * zonal_wavenumber / denominator  # rad s^-1
    frequency = angular_frequency / (2.0 * np.pi) * 86400.0  # cycles day^-1

    return frequency


def plot_wk_spectrum(
    spectrum,
    title,
    output_file,
    max_wn=15,
    max_freq=0.5,
    levels=None,
    period_ticks=(2, 3, 5, 10, 20, 30, 60, 100),
    depths=(8, 12, 25, 50),
    show_dispersion=True,
    caption=None,
):
    """Plot normalized Wheeler-Kiladis spectrum."""
    if levels is None:
        levels = np.linspace(0.8, 4.0, 12)

    plot_data = spectrum.sel(
        wavenumber=slice(-max_wn, max_wn),
        frequency=slice(0, max_freq),
    )

    fig, ax = plt.subplots(figsize=(12, 11.0))

    plot_data.plot.contourf(
        ax=ax,
        x="wavenumber",
        y="frequency",
        levels=levels,
        cmap="Spectral_r",
        extend="both",
        cbar_kwargs={"label": "normalized power"},
    )

    ax.axvline(0, color="k", linewidth=0.8)

    period_ticks = np.asarray(period_ticks)
    freq_ticks = 1.0 / period_ticks
    valid_ticks = (freq_ticks > 0) & (freq_ticks <= max_freq)

    ax.set_yticks(freq_ticks[valid_ticks])
    ax.set_yticklabels([f"{period:g}" for period in period_ticks[valid_ticks]])

    ax.set_ylabel("Frequency [cycles day$^{-1}$]")
    ax.set_xlabel("Zonal wavenumber")

    ax.set_xlim(-max_wn, max_wn)
    ax.set_ylim(0, max_freq)

    ax.set_title(title)

    for period in [30, 100]:
        frequency = 1.0 / period

        if frequency <= max_freq:
            ax.axhline(
                frequency,
                color="k",
                linestyle="--",
                linewidth=0.8,
                alpha=0.7,
            )

    if show_dispersion:
        kelvin_wn = np.linspace(0.05, max_wn, 800)
        rossby_wn = -np.linspace(0.05, max_wn, 800)

        for depth in depths:
            kelvin_freq = kelvin_curve(kelvin_wn, depth)
            rossby_freq = np.abs(
                equatorial_rossby_curve(
                    rossby_wn,
                    depth,
                    meridional_mode=1,
                )
            )

            kelvin_valid = (
                np.isfinite(kelvin_freq)
                & (kelvin_freq >= 0)
                & (kelvin_freq <= max_freq)
            )

            rossby_valid = (
                np.isfinite(rossby_freq)
                & (rossby_freq >= 0)
                & (rossby_freq <= max_freq)
            )

            ax.plot(
                kelvin_wn[kelvin_valid],
                kelvin_freq[kelvin_valid],
                color="k",
                linewidth=1.0,
            )

            ax.plot(
                rossby_wn[rossby_valid],
                rossby_freq[rossby_valid],
                color="k",
                linewidth=1.0,
            )
    if caption is not None:
        fig.text(
            0.5,
            0.01,
            caption,
            ha="center",
            va="bottom",
            fontsize=8,
            wrap=True,
        )
    output_file = str(output_file)
    fig.savefig(output_file, bbox_inches="tight", dpi=150)
    plt.close(fig)

    return output_file
