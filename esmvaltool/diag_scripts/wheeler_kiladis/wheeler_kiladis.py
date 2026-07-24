"""
Wheeler-Kiladis tropical wave spectrum diagnostic.

This diagnostic computes Wheeler-Kiladis (WK) wavenumber-frequency
spectra from preprocessed tropical fields following Wheeler and
Kiladis (1999). The diagnostic is intended for the analysis of
convectively coupled equatorial waves in precipitation and outgoing
longwave radiation (OLR).

The workflow consists of:

1. Read preprocessed ESMValTool input data.
2. Compute a cosine-weighted equatorial mean.
3. Remove the temporal mean and linear trend.
4. Remove the annual cycle using harmonic regression.
5. Optionally remove the zonal mean.
6. Compute wavenumber-frequency spectra using overlapping time segments.
7. Estimate a smooth background spectrum.
8. Compute normalized spectra.
9. Save NetCDF outputs and diagnostic plots.

Diagnostic options
------------------

annual_harmonics : int, optional
    Number of annual harmonics removed during seasonal-cycle
    correction. The default value of 3 removes the annual cycle
    and its first two harmonics while retaining most intraseasonal
    variability.

remove_zonal_mean : bool, optional
    Remove the zonal mean prior to spectral analysis. This is
    commonly done in Wheeler-Kiladis diagnostics to emphasize
    propagating equatorial wave signals. Default: True.

segment_length : int, optional
    Length of each spectral segment [days]. The default value of
    180 days is commonly used in Wheeler-Kiladis analyses and
    provides a balance between frequency resolution and statistical
    sampling.

segment_overlap : int, optional
    Overlap between consecutive segments [days]. The default value
    of 90 days corresponds to 50% overlap and increases the number
    of spectra contributing to the ensemble average.

sampling_frequency_per_day : float, optional
    Sampling frequency of the input data [samples day^-1].
    The default value is 1.0, corresponding to daily data.

sigma_freq : float, optional
    Gaussian smoothing width in the frequency direction used
    for background-spectrum estimation. The default value of
    4.0 provides sufficient smoothing to represent the broad
    spectral background while retaining its large-scale structure.

sigma_wn : float, optional
    Gaussian smoothing width in the zonal-wavenumber direction
    used for background-spectrum estimation. The default value
    of 4.0 provides a smoothly varying background spectrum
    suitable for normalization while preserving the large-scale
    spectral envelope.

max_wavenumber : int, optional
    Maximum zonal wavenumber displayed in diagnostic plots.
    Default: 15.

max_frequency : float, optional
    Maximum frequency displayed in diagnostic plots
    [cycles day^-1]. Default: 0.5.

period_ticks : list of int, optional
    Periods [days] shown as labels on the frequency axis.
    Default: [2, 3, 5, 10, 20, 30, 60, 100].

equivalent_depths : list of float, optional
    Equivalent depths [m] used to draw theoretical shallow-water
    dispersion curves. The default values [8, 12, 25, 50] span
    the range commonly associated with convectively coupled
    equatorial waves.

show_dispersion : bool, optional
    If True, overlay theoretical Kelvin and equatorial Rossby
    wave dispersion curves on the normalized Wheeler-Kiladis
    spectrum. Default: True.

mask_zero_wavenumber : bool, optional
    If True, mask the zonal-wavenumber-zero column in the
    normalized spectrum plot. This is often useful when the
    zonal mean has been removed. Default: True.

Outputs
-------

The diagnostic produces:

* Raw Wheeler-Kiladis power spectrum.
* Smoothed background power spectrum.
* Normalized Wheeler-Kiladis power spectrum.
* PNG figures showing normalized spectra with optional
  shallow-water dispersion curves.

References
----------

Wheeler, M. C. and Kiladis, G. N. (1999):
Convectively Coupled Equatorial Waves: Analysis of Clouds and
Temperature in the Wavenumber-Frequency Domain.
Journal of the Atmospheric Sciences, 56, 374-399.
"""

import logging
from pathlib import Path

import xarray as xr

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    run_diagnostic,
)
from esmvaltool.diag_scripts.wheeler_kiladis.plot import (
    plot_wk_spectrum,
)
from esmvaltool.diag_scripts.wheeler_kiladis.spectra import (
    compute_wk_spectra,
    equatorial_mean_regular,
    mask_zero_wavenumber,
    preprocess_equatorial_field,
    save_spectra,
)

logger = logging.getLogger(Path(__file__).stem)


def _load_regular_input(filename, variable):
    """Load regular gridded ESMValTool preprocessed input as xarray."""
    logger.info("Loading ESMValTool-preprocessed input: %s", filename)

    dataset = xr.open_dataset(filename)

    if variable in dataset:
        data = dataset[variable]
    else:
        # Fallback for cases where variable name differs from short_name.
        data_vars = list(dataset.data_vars)

        if len(data_vars) != 1:
            raise KeyError(
                f"Could not find variable {variable!r} in {filename}. "
                f"Available variables are: {data_vars}"
            )

        data = dataset[data_vars[0]]

    return data


def _process_regular_dataset(cfg, attrs):
    """Process one ESMValTool regular-grid dataset."""
    variable = attrs["short_name"]
    label = attrs.get("alias", attrs.get("dataset", "observations"))
    filename = attrs["filename"]

    data = _load_regular_input(filename, variable)

    data_eq = equatorial_mean_regular(data)

    data_anom = preprocess_equatorial_field(
        data_eq,
        remove_zonal_mean=cfg.get("remove_zonal_mean", True),
        nharmonics=cfg.get("annual_harmonics", 3),
    )

    raw, background, normalized = compute_wk_spectra(
        data_anom,
        segment_length=cfg.get("segment_length", 180),
        segment_overlap=cfg.get("segment_overlap", 90),
        samples_per_day=cfg.get("sampling_frequency_per_day", 1.0),
        sigma_freq=cfg.get("sigma_freq", 4.0),
        sigma_wn=cfg.get("sigma_wn", 4.0),
    )

    if cfg.get("remove_zonal_mean", True) and cfg.get(
        "mask_zero_wavenumber", True
    ):
        normalized_plot = mask_zero_wavenumber(normalized)
    else:
        normalized_plot = normalized

    return (
        label,
        variable,
        filename,
        raw,
        background,
        normalized,
        normalized_plot,
    )


def _make_method_caption(cfg, attrs, label):
    """Create a method caption describing WK settings."""
    dataset = attrs.get("dataset", label)

    start_year = attrs.get("start_year")
    end_year = attrs.get("end_year")

    caption = f"Wheeler-Kiladis spectrum for {dataset}. "

    if start_year is not None and end_year is not None:
        caption += f"Period: {start_year}-{end_year}. "

    caption += (
        "Region: 15S-15N, all longitudes. "
        "Data are converted to daily means and equatorially averaged. "
        f"Anomalies are computed after removing the time mean, linear trend, "
        f"and the first {cfg.get('annual_harmonics', 3)} annual harmonics. "
    )

    if cfg.get("remove_zonal_mean", True):
        caption += (
            "The zonal mean is removed before the "
            "wavenumber-frequency transform. "
        )

    caption += (
        f"Spectra are computed using "
        f"{cfg.get('segment_length', 180)}-day windows "
        f"with {cfg.get('segment_overlap', 90)}-day overlap. "
        "Background normalization uses Gaussian smoothing with "
        f"sigma_freq={cfg.get('sigma_freq', 4.0)} and "
        f"sigma_wn={cfg.get('sigma_wn', 4.0)}. "
        "Positive zonal wavenumber denotes eastward propagation."
    )
    return caption


def _plot_and_save(
    cfg,
    attrs,
    label,
    variable,
    ancestor,
    raw,
    background,
    normalized,
    normalized_plot,
):
    """Save spectra and plot for one dataset."""
    work_dir = Path(cfg["work_dir"])
    plot_dir = Path(cfg["plot_dir"])

    safe_label = label.replace(" ", "_").replace("/", "_")
    basename = f"{safe_label}_{variable}_wk"

    raw_file, background_file, normalized_file = save_spectra(
        raw,
        background,
        normalized,
        output_dir=work_dir,
        basename=basename,
    )

    plot_file = plot_dir / f"{basename}_normalized.png"

    method_caption = _make_method_caption(
        cfg,
        attrs,
        label,
    )

    plot_wk_spectrum(
        normalized_plot,
        title=f"{label}: normalized WK spectrum of {variable}",
        output_file=plot_file,
        max_wn=cfg["max_wavenumber"],
        max_freq=cfg["max_frequency"],
        levels=None,
        period_ticks=tuple(
            cfg.get("period_ticks", [2, 3, 5, 10, 20, 30, 60, 100])
        ),
        depths=tuple(cfg.get("equivalent_depths", [8, 12, 25, 50])),
        show_dispersion=cfg.get("show_dispersion", True),
        caption=method_caption,
    )

    record = {
        "caption": method_caption,
        "statistics": ["other"],
        "domains": ["trop"],
        "plot_types": ["zonal"],
        "authors": ["londono-castillo_santiago"],
        "references": ["acknow_project"],
        "ancestors": [str(ancestor)],
    }

    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(str(plot_file), record)
        provenance_logger.log(str(raw_file), record)
        provenance_logger.log(str(background_file), record)
        provenance_logger.log(str(normalized_file), record)

    return {
        "label": label,
        "ancestor": ancestor,
        "raw": raw,
        "background": background,
        "normalized": normalized,
        "normalized_plot": normalized_plot,
        "plot_file": plot_file,
    }


def main(cfg):
    """Run diagnostic."""
    logger.info("Starting Wheeler-Kiladis diagnostic")
    logger.info("Configuration: %s", cfg)

    metadata = list(cfg["input_data"].values())

    if not metadata:
        logger.warning("No ESMValTool input metadata found.")

    results = []

    for attrs in metadata:
        (
            label,
            variable,
            ancestor,
            raw,
            background,
            normalized,
            normalized_plot,
        ) = _process_regular_dataset(cfg, attrs)

        result = _plot_and_save(
            cfg,
            attrs,
            label,
            variable,
            ancestor,
            raw,
            background,
            normalized,
            normalized_plot,
        )

        results.append(result)
    logger.info("Wheeler-Kiladis diagnostic finished")


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
