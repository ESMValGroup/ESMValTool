# Copyright (C) 2026 ESMValTool development team
"""May-October precipitation and 850 hPa wind mean-state comparison.

The recipe supplies five already subset and season-selected inputs: daily
NCEP-NCAR-R1 precipitation and winds for 2000-2020, and daily precipitation
from the final 30 years of ACCESS-ESM1-5 historical and ACCESS-CM2 piControl.
The two model periods and experiments are intentionally labelled separately;
this is a descriptive comparison rather than a matched-experiment test.
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
    """Return exactly one cube for each named recipe variable group."""
    grouped = group_metadata(cfg["input_data"].values(), "variable_group")
    required = ("obs_pr", "obs_ua", "obs_va", "historical_pr", "control_pr")
    if set(grouped) != set(required):
        msg = f"Expected variable groups {required}, got {tuple(grouped)}"
        raise ValueError(
            msg,
        )
    cubes = {}
    ancestors = []
    for name in required:
        entries = grouped[name]
        if len(entries) != 1:
            msg = f"Expected one input for {name}, got {len(entries)}"
            raise ValueError(
                msg,
            )
        filename = entries[0]["filename"]
        cubes[name] = iris.load_cube(filename)
        ancestors.append(filename)
    return cubes, ancestors


def calculate_biases(cubes):
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
    historical_bias.rename("ACCESS-ESM1-5 historical precipitation bias")
    control_bias.rename("ACCESS-CM2 piControl precipitation bias")
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


def plot_biases(means, historical_bias, control_bias, periods):
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
        f"NCEP-NCAR-R1 ({periods['obs_pr']})\n"
        "May-October precipitation and 850 hPa winds",
    )

    biases = (
        (
            historical_bias,
            (
                "ACCESS-ESM1-5 historical "
                f"({periods['historical_pr']}) - NCEP-NCAR-R1"
            ),
        ),
        (
            control_bias,
            (f"ACCESS-CM2 piControl ({periods['control_pr']}) - NCEP-NCAR-R1"),
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
    cubes, ancestors = _load_inputs(cfg)
    periods = {name: _period_label(cube) for name, cube in cubes.items()}
    means, historical_bias, control_bias = calculate_biases(cubes)
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
            "May-October NCEP-NCAR-R1 precipitation and 850 hPa winds, "
            "with ACCESS-ESM1-5 historical and ACCESS-CM2 piControl "
            "precipitation biases relative to that reference"
        ),
        "ancestors": ancestors,
        "authors": ["sullivan_arnold", "chun_felicity", "beucher_romain"],
        "domains": ["trop"],
        "statistics": ["mean", "diff"],
        "plot_types": ["map"],
    }
    with ProvenanceLogger(cfg) as logger:
        logger.log(data_file, provenance)
    figure = plot_biases(means, historical_bias, control_bias, periods)
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
