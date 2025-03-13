"""Python example diagnostic."""

import logging
from pathlib import Path

import iris
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.shared import run_diagnostic, save_figure

logger = logging.getLogger(Path(__file__).stem)


def main(cfg):
    """Plot part of figure_9.3a from IPCC AR6."""
    colors = {
        "historical-ssp126": "#2a3652",
        "historical-ssp585": "#78333a",
    }
    fill_colors = {
        "historical-ssp126": "#d2d5dc",
        "historical-ssp585": "#ddced2",
    }
    labels = {
        "historical-ssp126": "Historical and SSP1-2.6",
        "historical-ssp585": "Historical and SSP5-8.5",
    }

    # Group input data by experiment
    groups = {}
    for filename, attributes in cfg["input_data"].items():
        exp = attributes["exp"]
        if exp not in groups:
            groups[exp] = {}
        groups[exp][attributes["dataset"]] = filename

    # Loop over experiments to populate plot
    for exp, group in groups.items():
        mean = iris.load_cube(group["MultiModelMean"])
        iris.quickplot.plot(
            mean,
            color=colors.get(exp),
            label=labels.get(exp, exp),
        )

        p17 = iris.load_cube(group["MultiModelPercentile17"])
        p83 = iris.load_cube(group["MultiModelPercentile83"])
        time_coord = mean.coord("time")
        time_axis = time_coord.units.num2date(time_coord.core_points())
        plt.fill_between(
            time_axis,
            p17.core_data(),
            p83.core_data(),
            color=fill_colors.get(exp),
            label="Likely (17% - 83%) ranges",
        )

    plt.title("Sea surface temperature anomaly")
    plt.legend(loc="upper left")

    filename = "IPCC_AR6_figure_9.3a_1850-2100"
    provenance_record = {
        "caption": "Part of figure 9.3a from IPCC AR6.",
        "authors": [
            "kalverla_peter",
            "andela_bouwe",
        ],
        "references": ["fox-kemper21ipcc"],
        "ancestors": list(cfg["input_data"].keys()),
    }
    save_figure(filename, provenance_record, cfg, dpi=300)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
