# Copyright (C) 2026 ESMValTool development team
"""May-October precipitation and 850 hPa wind mean-state comparison.

The recipe supplies preprocessed reference precipitation and winds, plus
precipitation from historical and control experiments. Dataset names and
periods are taken from the input metadata so the diagnostic can be reused
with other datasets.
"""

import cartopy.crs as ccrs
import iris
import iris.analysis
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import regrid
from matplotlib import patches

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    group_metadata,
    run_diagnostic,
    save_figure,
)

REGIONS = {
    "BOB": (85, 10, 5, 5),
    "NEQ": (85, 5, 5, 5),
    "EEIO": (75, -5, 10, 10),
    "WP": (120, 10, 30, 10),
}


def _load_inputs(cfg):
    """Return one cube and its dataset alias for each variable group."""
    grouped = group_metadata(cfg["input_data"].values(), "variable_group")
    required = ("obs_pr", "obs_ua", "obs_va", "historical_pr", "control_pr")
    if set(grouped) != set(required):
        msg = f"Expected variable groups {required}, got {tuple(grouped)}"
        raise ValueError(
            msg,
        )
    cubes = {}
    aliases = {}
    ancestors = []
    for name in required:
        entries = grouped[name]
        if len(entries) != 1:
            msg = f"Expected one input for {name}, got {len(entries)}"
            raise ValueError(
                msg,
            )
        metadata = entries[0]
        filename = metadata["filename"]
        cubes[name] = iris.load_cube(filename)
        aliases[name] = metadata.get("alias") or metadata["dataset"]
        ancestors.append(filename)
    return cubes, aliases, ancestors


def calculate_biases(cubes, aliases):
    """Compute seasonal means and model-minus-reference precipitation bias."""
    means = {
        name: iris.util.squeeze(cube.collapsed("time", iris.analysis.MEAN))
        for name, cube in cubes.items()
    }
    mean_obs_pr = means["obs_pr"]
    historical_bias = means["historical_pr"] - regrid(
        mean_obs_pr,
        means["historical_pr"],
        scheme="linear",
    )
    control_bias = means["control_pr"] - regrid(
        mean_obs_pr,
        means["control_pr"],
        scheme="linear",
    )
    historical_bias.rename(
        f"{aliases['historical_pr']} minus {aliases['obs_pr']} precipitation",
    )
    control_bias.rename(
        f"{aliases['control_pr']} minus {aliases['obs_pr']} precipitation",
    )
    historical_bias.var_name = "historical_bias"
    control_bias.var_name = "control_bias"
    mean_obs_pr.var_name = "obs_pr"
    means["obs_ua"].var_name = "obs_ua"
    means["obs_va"].var_name = "obs_va"
    return means, historical_bias, control_bias


def _period_label(cube):
    """Return the years actually present after recipe preprocessing."""
    time = cube.coord("time")
    first = time.units.num2date(time.points[0]).year
    last = time.units.num2date(time.points[-1]).year
    return f"{first}-{last}"


def plot_biases(means, historical_bias, control_bias, periods, aliases):
    """Render the notebook's three-panel mean state and bias figure."""
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(10, 18),
        subplot_kw={"projection": ccrs.PlateCarree()},
    )
    observed_levels = np.linspace(0, 16, 17)
    bias_levels = np.unique(
        np.round(
            np.concatenate(
                (
                    np.arange(-6, -1, 2),
                    np.arange(-0.4, 0.5, 0.2),
                    np.arange(2, 7, 2),
                ),
            ),
            1,
        ),
    )
    observed = iplt.contourf(
        means["obs_pr"],
        axes=axes[0],
        levels=observed_levels,
        cmap="YlGnBu",
        extend="max",
    )
    coarse_u = regrid(means["obs_ua"], target_grid="10x10", scheme="nearest")
    coarse_v = regrid(means["obs_va"], target_grid="10x10", scheme="nearest")
    arrows = iplt.quiver(
        coarse_u,
        coarse_v,
        axes=axes[0],
        color="white",
        scale=150,
        width=0.005,
    )
    axes[0].quiverkey(
        arrows,
        0.9,
        1.05,
        10,
        "10 m/s",
        labelpos="E",
        coordinates="axes",
        color="black",
    )
    for name, (left, bottom, width, height) in REGIONS.items():
        rectangle = patches.Rectangle(
            (left, bottom),
            width,
            height,
            linewidth=2.5,
            edgecolor="red",
            facecolor="none",
            transform=ccrs.PlateCarree(),
            zorder=5,
        )
        axes[0].add_patch(rectangle)
        axes[0].text(
            left + width / 2,
            bottom + height / 2,
            name,
            color="red",
            fontweight="bold",
            fontsize=9,
            ha="center",
            va="center",
            bbox={
                "facecolor": "white",
                "alpha": 0.85,
                "edgecolor": "none",
                "pad": 1,
            },
            transform=ccrs.PlateCarree(),
            zorder=6,
        )
    axes[0].set_title(
        f"{aliases['obs_pr']} ({periods['obs_pr']})\n"
        "May-October precipitation and 850 hPa winds",
    )

    biases = (
        (
            historical_bias,
            (
                f"{aliases['historical_pr']} "
                f"({periods['historical_pr']}) - {aliases['obs_pr']}"
            ),
        ),
        (
            control_bias,
            (
                f"{aliases['control_pr']} "
                f"({periods['control_pr']}) - {aliases['obs_pr']}"
            ),
        ),
    )
    for axis, (cube, title) in zip(axes[1:], biases, strict=False):
        image = iplt.contourf(
            cube,
            axes=axis,
            levels=bias_levels,
            cmap="RdBu",
            extend="both",
        )
        axis.set_title(title)
    for axis in axes:
        axis.coastlines()
        axis.set_extent([50, 180, -45, 30], crs=ccrs.PlateCarree())
        grid = axis.gridlines(draw_labels=True, linestyle="--", alpha=0.5)
        grid.top_labels = grid.right_labels = False
    fig.tight_layout(rect=[0, 0, 0.85, 0.96])
    fig.colorbar(
        observed,
        cax=fig.add_axes([0.88, 0.70, 0.03, 0.22]),
        label="mm/day",
    )
    fig.colorbar(
        image,
        cax=fig.add_axes([0.88, 0.15, 0.03, 0.45]),
        label="Model minus reference (mm/day); blue = wetter",
    )
    fig.suptitle(
        "May-October mean state and precipitation biases",
        fontsize=18,
    )
    return fig


def main(cfg):
    """Save the plotted mean fields, biases, and figure with provenance."""
    cubes, aliases, ancestors = _load_inputs(cfg)
    periods = {name: _period_label(cube) for name, cube in cubes.items()}
    means, historical_bias, control_bias = calculate_biases(cubes, aliases)
    output = iris.cube.CubeList(
        [
            means["obs_pr"],
            means["obs_ua"],
            means["obs_va"],
            historical_bias,
            control_bias,
        ],
    )
    data_file = get_diagnostic_filename("mjo_bsiso_mean_state_bias", cfg)
    iris.save(output, data_file)
    provenance = {
        "caption": (
            f"May-October {aliases['obs_pr']} precipitation and 850 hPa "
            f"winds, with {aliases['historical_pr']} and "
            f"{aliases['control_pr']} precipitation biases relative to "
            "that reference"
        ),
        "ancestors": ancestors,
        "authors": ["sullivan_arnold", "chun_felicity", "beucher_romain"],
        "domains": ["trop"],
        "statistics": ["mean", "diff"],
        "plot_types": ["map"],
    }
    with ProvenanceLogger(cfg) as logger:
        logger.log(data_file, provenance)
    figure = plot_biases(
        means,
        historical_bias,
        control_bias,
        periods,
        aliases,
    )
    save_figure(
        "mjo_bsiso_mean_state_bias",
        provenance,
        cfg,
        figure=figure,
        dpi=150,
        bbox_inches="tight",
    )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
