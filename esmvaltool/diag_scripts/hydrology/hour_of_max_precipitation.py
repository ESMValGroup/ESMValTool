"""Calculate and plot hour of daily maximum precipitation.

Description
-----------
This diagnostics calculates the hour of daily maximum precipitation and plots
it. The input data needs to be subdaily and of shape (hour, latitude,
longitude). The hour dimension should be in local solar time (see
:func:`~esmvalcore.preprocessor.local_solar_time`).

Author
------
Manuel Schlund (DLR, Germany)

Configuration options in recipe
-------------------------------
caption: str, optional
    Figure caption used for provenance tracking. By default, uses "Global map
    of the hour of the daily maximum precipitation for {alias}.".
cbar_label: str, optional (default: "Hour of daily maximum precipitation")
    Colorbar label.
cbar_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.colorbar`. By
    default, uses ``{orientation: "horizontal"}``.
figure_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.figure`. By
    default, uses ``constrained_layout: true, ticks: [0, 3, 6, 9, 12, 15, 18,
    21]``.
matplotlib_rc_params: dict, optional
    Optional :class:`matplotlib.RcParams` used to customize matplotlib plots.
    Options given here will be passed to :func:`matplotlib.rc_context` and used
    for all plots produced with this diagnostic.
method: str, optional (default: `"dft"`)
    Method to determine the hour of daily maximum precipitation. Possible
    options:

    - If `"dft"`, apply a discrete Fourier transform (DFT) to the data and use
      the peak of the first component (see
      https://mdtf-diagnostics.readthedocs.io/en/latest/sphinx_pods/precip_diurnal_cycle.html).
    - If `"harmonic_fit"`, fit a 12hr- plus 24hr-harmonic to the input data and
      use the peak of the 24hr-harmonic (see Dai, 2024;
      https://doi.org/10.1007/s00382-024-07182-6). This should give
      identical/very similar results to `"dft"`, but is much slower.
    - If `"max"`, simply use the maximum in the input data.
plot_kwargs: dict, optional
    Optional keyword arguments for :func:`iris.plot.pcolormesh`. By default,
    uses ``cmap: twilight``.
projection: str, optional (default: None)
    Projection used for the plot. Needs to be a valid projection class of
    :mod:`cartopy.crs`. Keyword arguments can be specified using the option
    ``projection_kwargs``.
projection_kwargs: dict, optional
    Optional keyword arguments for the projection given by ``projection``. For
    map plots, the default keyword arguments ``{central_longitude: 10}`` are
    used.
pyplot_kwargs: dict, optional
    Optional calls to functions of :mod:`matplotlib.pyplot`. Dictionary keys
    are functions of :mod:`matplotlib.pyplot`. Dictionary values are used as
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can contain format strings (e.g., ``"{dataset}
    ({project})"``).
savefig_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.savefig`. By
    default, uses ``bbox_inches: tight, dpi: 300, orientation: landscape``.
seaborn_settings: dict, optional
    Options for :func:`seaborn.set_theme` (affects all plots). By default, uses
    ``style: ticks``.
threshold: float, optional (default: 0.0)
    Mask grid points where precipitation is lower than the given threshold.

"""

import logging
import warnings
from copy import deepcopy
from pathlib import Path
from typing import Any

import cartopy.crs as ccrs
import iris
import iris.plot
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from iris.cube import Cube
from iris.warnings import IrisVagueMetadataWarning
from matplotlib.figure import Figure
from scipy.optimize import curve_fit

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared._base import save_data, save_figure

logger = logging.getLogger(Path(__file__).stem)


def _harmonic_function(
    x_val: float,
    constant: float,
    amplitude_24: float,
    phase_24: float,
    amplitude_12: float,
    phase_12: float,
) -> float:
    """Harmonic function which can be fitted to diurnal cycle."""
    diurnal_component = amplitude_24 * np.sin(
        2 * np.pi * (x_val - phase_24) / 24.0,
    )
    semidiurnal_component = amplitude_12 * np.sin(
        2 * np.pi * (x_val - phase_12) / 12.0,
    )
    return constant + diurnal_component + semidiurnal_component


def _get_max_of_24hr_harmonic(
    x_data: np.ma.MaskedArray,
    y_data: np.ma.MaskedArray,
) -> float:
    """Get maximum of 24hr-harmonic fit."""
    # curve_fit cannot handle masked data
    if np.ma.is_masked(y_data):
        return np.nan
    params = curve_fit(_harmonic_function, x_data, y_data, nan_policy="omit")[
        0
    ]
    amplitude_24hr = params[1]
    phase_24hr = params[2]

    # Since we are using sin() as fit function here, the maximum of the fit
    # function is located 1/4 of a periodic length (i.e., π/2; here: 6hrs) to
    # the right/left (depending on the sign of the amplitude) of the phase (=
    # first zero of sin()). Since the phase can be outside [0, 24] due to the
    # periodicity of sin(), make sure to return the maximum within [0, 24].
    maximum = phase_24hr + 6.0 if amplitude_24hr > 0 else phase_24hr - 6.0
    return maximum % 24.0


_v_get_max_of_24hr_harmonic = np.vectorize(
    _get_max_of_24hr_harmonic,
    signature="(t),(t)->()",
)


def _calculate_hour_of_max_precipitation(
    cube: Cube,
    cfg: dict[str, Any],
) -> Cube:
    """Calculate hour of daily maximum daily precipitation."""
    hour_dim = cube.coord_dims("hour")[0]

    # Mask values mean diurnal cycle is smaller than threshold
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=IrisVagueMetadataWarning,
            module="iris",
        )
        mean_diurnal_cycle = cube.collapsed("hour", iris.analysis.MEAN)
    mask = mean_diurnal_cycle.data < cfg["threshold"]
    broadcasted_mask = np.broadcast_to(
        np.expand_dims(mask, axis=hour_dim),
        cube.shape,
    )
    cube = cube.copy(np.ma.masked_array(cube.data, mask=broadcasted_mask))

    # Calculate hour of daily maximum precipitation
    if cfg["method"] == "dft":
        dft = np.fft.rfft(cube.data, axis=hour_dim)

        # The phase here is the phase of a cosine, so this is identical to the
        # maximum.
        max_component_1 = np.arctan2(-dft[1].imag, dft[1].real)

        # Transform from radians to hours and make sure we end up in [0, 24]
        max_24hrs = (max_component_1 * 24.0 / 2.0 / np.pi) % 24.0
        hour_of_max_precipitation_data = np.ma.array(max_24hrs, mask=mask)
    elif cfg["method"] == "harmonic_fit":
        # np.vectorized assumes that the core dimension is the rightmost
        # dimension
        new_order = [d for d in range(cube.ndim) if d != hour_dim] + [hour_dim]
        hour_of_max_precipitation_data = _v_get_max_of_24hr_harmonic(
            cube.coord("hour").points,
            cube.data.transpose(new_order),
        )
        hour_of_max_precipitation_data = np.ma.masked_invalid(
            hour_of_max_precipitation_data,
        )
    elif cfg["method"] == "max":
        hour_of_max_precipitation_data = np.ma.masked_array(
            cube.coord("hour").points[np.argmax(cube.data, axis=hour_dim)],
            mask=mask,
        )
    else:
        supported_methods = (
            "dft",
            "harmonic_fit",
            "max",
        )
        msg = f"Expected one of {supported_methods} for 'method', got '{cfg['method']}'"
        raise ValueError(msg)

    return mean_diurnal_cycle.copy(hour_of_max_precipitation_data)


def _get_default_cfg(cfg: dict) -> dict:
    """Get default options for configuration dictionary."""
    cfg = deepcopy(cfg)

    cfg.setdefault(
        "caption",
        "Global map of the hour of the daily maximum precipitation for {alias}.",
    )
    cfg.setdefault("cbar_label", "Hour of daily maximum precipitation")
    cfg.setdefault(
        "cbar_kwargs",
        {"orientation": "horizontal", "ticks": [0, 3, 6, 9, 12, 15, 18, 21]},
    )
    cfg.setdefault("figure_kwargs", {"constrained_layout": True})
    cfg.setdefault("matplotlib_rc_params", {})
    cfg.setdefault("method", "dft")
    cfg.setdefault("plot_kwargs", {"cmap": "twilight"})
    cfg.setdefault("projection", "Robinson")
    cfg.setdefault("projection_kwargs", {"central_longitude": 10})
    cfg.setdefault("pyplot_kwargs", {})
    cfg.setdefault(
        "savefig_kwargs",
        {
            "bbox_inches": "tight",
            "dpi": 300,
            "orientation": "landscape",
        },
    )
    cfg.setdefault("seaborn_settings", {"style": "ticks"})
    cfg.setdefault("threshold", 0.0)

    supported_methods = (
        "dft",
        "harmonic_fit",
        "max",
    )

    if cfg["method"] not in supported_methods:
        msg = f"Expected one of {supported_methods} for 'method', got '{cfg['method']}'"
        raise ValueError(msg)

    return cfg


def _get_projection(cfg: dict[str, Any]) -> Any:
    """Get plot projection."""
    projection = cfg["projection"]
    projection_kwargs = cfg["projection_kwargs"]

    if not hasattr(ccrs, projection):
        msg = f"Got invalid projection '{projection}', expected class of cartopy.crs"
        raise AttributeError(msg)

    return getattr(ccrs, projection)(**projection_kwargs)


def _get_provenance_record(
    ancestors: list[str],
    caption: str,
) -> dict[str, Any]:
    """Get provenance record."""
    return {
        "ancestors": ancestors,
        "authors": ["schlund_manuel"],
        "caption": caption,
        "plot_types": ["map"],
        "realms": ["atmos"],
        "references": ["dai24climdyn"],
        "themes": ["phys"],
    }


def _create_plot(cube: Cube, dataset: dict, cfg: dict[str, Any]) -> Figure:
    """Plot map."""
    fig = plt.figure(**cfg["figure_kwargs"])
    axes = fig.add_subplot(projection=_get_projection(cfg))
    plot_kwargs = cfg["plot_kwargs"]
    plot_kwargs["axes"] = axes
    map_plot = iris.plot.pcolormesh(cube, **plot_kwargs)
    axes.set_title(dataset["alias"])
    _process_pyplot_kwargs(cfg["pyplot_kwargs"], dataset)
    cbar = plt.colorbar(map_plot, ax=axes, **cfg["cbar_kwargs"])
    cbar.set_label(cfg["cbar_label"])
    return fig


def _process_pyplot_kwargs(pyplot_kwargs: Any, dataset: dict) -> None:
    """Process functions for :mod:`matplotlib.pyplot`."""
    for func, arg in pyplot_kwargs.items():
        if arg is None:
            getattr(plt, func)()
        elif isinstance(arg, dict):
            getattr(plt, func)(**arg)
        elif isinstance(arg, str):
            getattr(plt, func)(arg.format(**dataset))
        else:
            getattr(plt, func)(arg)


def main(cfg: dict) -> None:
    """Run diagnostic."""
    cfg = _get_default_cfg(cfg)
    sns.set_theme(**cfg["seaborn_settings"])

    for dataset in cfg["input_data"].values():
        filename = dataset["filename"]
        basename = Path(filename).stem
        caption = cfg["caption"].format(**dataset)
        provenance_record = _get_provenance_record([filename], caption)

        # Calculation
        logger.info("Loading %s", filename)
        cube = iris.load_cube(filename)
        cube = _calculate_hour_of_max_precipitation(cube, cfg)
        save_data(basename, provenance_record, cfg, cube)

        # Plot
        logger.info("Plotting map for %s", dataset["alias"])
        figure = _create_plot(cube, dataset, cfg)
        save_figure(
            basename,
            provenance_record,
            cfg,
            figure=figure,
            **cfg["savefig_kwargs"],
        )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
