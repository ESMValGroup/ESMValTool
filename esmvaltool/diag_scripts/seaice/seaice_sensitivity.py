import logging
from pathlib import Path

import iris
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from scipy import stats

from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    save_figure,
    select_metadata,
)

logger = logging.getLogger(Path(__file__).stem)


# This is stolen from the example in AutoAssess _plot_mo_metrics.py
# TODO: Work out what should be written here
def get_provenance_record(cfg):
    """Create a provenance record describing the diagnostic data and plot."""
    filenames = [item["filename"] for item in cfg["input_data"].values()]

    region = [item["diagnostic"] for item in cfg["input_data"].values()][0]

    record = {
        "caption": f"Plots of {region.title()} sea ice sensitivity",
        "plot_type": "metrics",
        "authors": [
            "sellar_alistair",
        ],
        "ancestors": filenames,
    }

    return record


def list_datasets(data):
    """Actually returns a set of datatsets, to avoid duplication."""
    logger.debug("listing datasets")
    datasets = set()
    for element in data:
        datasets.add(element["dataset"])
    return datasets


def extract_cube(data, variable_group):
    """
    Select records from data based on the variable_group, assert that there is
    only one matching record, then return an iris cube for that record.
    """
    logger.debug("extracting %s cube from %s", variable_group, data)

    # Load any data in the variable_group
    selection = select_metadata(data, variable_group=variable_group)

    # Ensure there is only one file in the list
    assert len(selection) == 1, (
        f"None or too many matching files found for {variable_group}"
    )

    # Load the cube, [0] is because selection returns a list
    cube = iris.load_cube(selection[0]["filename"])

    return cube


def calculate_regression(independent, dependent, slope_only):
    """
    Use SciPy stats to calculate the least-squares regression
    """
    logger.debug(
        "Calculating linear relationship between %s and %s",
        dependent,
        independent,
    )

    # Use SciPy stats to calculate the regression
    # result = slope, intercept, rvalue, pvalue, stderr
    result = stats.linregress(independent, dependent)

    # Return only the slope
    if slope_only:
        return result.slope

    # Or else return everything
    return result


def calculate_annual_trends(data):
    """
    Calculate annual trends for surface air temperature (tas) and sea ice area (siconc),
    and the r and p values from the regression of siconc as a function of tas
    """
    logger.debug("calculating annual trends")

    # Load the preprocessed cubes
    si_cube = extract_cube(data, "siconc")
    tas_cube = extract_cube(data, "tas")

    # Calculate the individual trends over time
    years = tas_cube.coord("year").points
    si_trend = calculate_regression(years, si_cube.data, slope_only=True)
    tas_trend = calculate_regression(years, tas_cube.data, slope_only=True)

    # Calculate the direct regression for r and p values
    direct_regression = calculate_regression(
        tas_cube.data, si_cube.data, slope_only=False
    )

    dictionary = {
        "si_ann_trend": si_trend,
        "tas_ann_trend": tas_trend,
        "direct_r_val": direct_regression.rvalue,
        "direct_p_val": direct_regression.pvalue,
    }

    return dictionary


def calculate_direct_sensitivity(data):
    """
    Calculate slope of sea ice area over global mean temperature
    """
    logger.debug("calculating sensitivity")

    # Load the preprocessed cubes
    si_cube = extract_cube(data, "siconc")
    tas_cube = extract_cube(data, "tas")

    # Calculate the slope of the direct regression, NOT via time
    sensitivity = calculate_regression(
        tas_cube.data, si_cube.data, slope_only=True
    )

    return sensitivity


def create_titles_dict(data):
    """
    Create a dictionary of appropriate observations and titles
    depending on whether the plot is for the Arctic or Antarctic
    and assuming the recipe used September Arctic sea ice data
    and annually mean averaged Antarctic sea ice data
    """
    dictionary = {}

    first_variable = next(iter(data))

    if first_variable["diagnostic"] == "arctic":
        # Ed Blockley's dSIA/dT values in million sq km per Kelvin
        dictionary["obs"] = {
            "1979-2014": {"mean": -4.01, "std_dev": 0.32, "plausible": 1.28}
        }
        # Setting titles for plots
        dictionary["titles"] = {
            "notz_fig_title": "September Arctic Sea Ice Sensitivity",
            "notz_ax_title": f"dSIA/dGMST obs from {next(iter(dictionary['obs'].keys()))}",
            "notz_plot_filename": "September Arctic sea ice sensitivity",
            "roach_fig_title": "Trends in Annual Mean Temperature And September Arctic Sea Ice",
            "roach_plot_filename": "September Arctic sea ice trends",
        }

    elif first_variable["diagnostic"] == "antarctic":
        # We don't have equivalent observations  # TODO: there are observations on Roach 3e
        dictionary["obs"] = {
            "no years": {"mean": 0.0, "std_dev": 0.0, "plausible": 0.0}
        }
        # Setting titles for plots
        dictionary["titles"] = {
            "notz_fig_title": "Annually Meaned Antarctic Sea Ice Sensitivity",
            "notz_ax_title": "dSIA/dGMST",
            "notz_plot_filename": "Annual Antarctic sea ice sensitivity",
            "roach_fig_title": "Trends in Annual Mean Temperature And Annual Antarctic Sea Ice",
            "roach_plot_filename": "Annual Antarctic sea ice trends",
        }

    return dictionary


def notz_style_plot_from_dict(data_dictionary, titles_dictionary, cfg):
    """
    Save a plot of sensitivities and observations
    """
    # Read from (Ed Blockley's) observations dictionary
    obs_dict = titles_dictionary["obs"]
    obs_years = list(obs_dict)[0]
    obs_mean = obs_dict[obs_years]["mean"]
    obs_std_dev = obs_dict[obs_years]["std_dev"]
    obs_plausible = obs_dict[obs_years]["plausible"]

    # Set up the figure
    fig, ax = plt.subplots(figsize=(4, 6), layout="constrained")
    fig.suptitle(titles_dictionary["titles"]["notz_fig_title"])
    ax.set_title(titles_dictionary["titles"]["notz_ax_title"])  # Ed's title

    # Iterate over the dictionary
    for dataset, inner_dict in data_dictionary.items():
        ax.plot(
            0.5,
            inner_dict["direct_sensitivity"],
            label=dataset,
            marker="_",
            markersize=20,
        )

    # Add observations (style taken from Ed's code)
    ax.hlines(obs_mean, 0, 1, linestyle="--", color="black", linewidth=2)
    ax.fill_between(
        [0, 1],
        obs_mean - obs_std_dev,
        obs_mean + obs_std_dev,
        facecolor="k",
        alpha=0.15,
    )
    ax.hlines(
        obs_mean + obs_plausible, 0, 1, linestyle=":", color="0.5", linewidth=1
    )
    ax.hlines(
        obs_mean - obs_plausible, 0, 1, linestyle=":", color="0.5", linewidth=1
    )

    # Tidy the figure
    ax.set_xticks([])
    ax.set_ylabel("dSIA/dGMST (million km$^2$ K$^{-1}$)")
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))

    # Save the figure (also closes it)
    provenance_record = get_provenance_record(cfg)
    save_figure(
        titles_dictionary["titles"]["notz_plot_filename"],
        provenance_record,
        cfg,
        figure=fig,
        close=True,
    )


def roach_style_plot_from_dict(data_dictionary, titles_dictionary, cfg):
    """
    Save a plot of trend in SIA against trend in GMST to the given filename
    """
    # Set up the figure
    fig, ax = plt.subplots(figsize=(10, 6), layout="constrained")
    fig.suptitle(titles_dictionary["titles"]["roach_fig_title"])
    ax.set_title("Decades are scaled up from annual data")

    # Set up for colouring the points
    norm = Normalize(vmin=-1, vmax=1)
    cmap = plt.get_cmap("PiYG_r")

    # Set up the axes
    ax.set_xlim(-0.02, 0.505)
    ax.set_ylim(-1.1, 0.21)
    ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_yticks([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2])
    ax.hlines(0, -0.02, 0.505, color="black", alpha=0.5)
    ax.vlines(0, -1.1, 0.21, color="black", alpha=0.5)
    ax.set_xlabel("Trend in GMST (K decade$^{-1}$)")
    ax.set_ylabel("Trend in SIA (million km$^2$ decade$^{-1}$)")

    # Iterate over the dictionary
    for dataset, inner_dict in data_dictionary.items():
        # Determine the position of the point
        x = 10 * inner_dict["tas_trend"]  # for equivalence to decades
        y = 10 * inner_dict["siconc_trend"]  # for equivalence to decades

        # Determine the colour of the point
        r2 = inner_dict["direct_r_val"] ** 2 * np.sign(
            inner_dict["direct_r_val"]
        )  # TODO: check this is right

        # Decide if the point should be hatched
        if inner_dict["direct_p_val"] > 0.05:  # TODO: also check this
            h = 5 * "/"  # This is a hatch pattern
        else:
            h = None

        # Plot the point
        plt.scatter(
            x, y, marker="o", s=150, c=[r2], hatch=h, cmap=cmap, norm=norm
        )

        # Label with the dataset
        plt.annotate(dataset, xy=(x, y), xytext=(x + 0.02, y + 0.02))

    # Add a colour bar
    plt.colorbar(label="R^2 * sign(R)")

    # Save the figure (also closes it)
    provenance_record = get_provenance_record(cfg)
    save_figure(
        titles_dictionary["titles"]["roach_plot_filename"],
        provenance_record,
        cfg,
        figure=fig,
        close=True,
    )


def main(cfg):
    """
    Create two plots per diagnostic from preprocessed data
    """
    # Get the data from the cfg
    logger.info("Getting data from the config")
    input_data = cfg["input_data"].values()

    # Titles and observations depend on the diagnostic being plotted
    titles_and_obs_dict = create_titles_dict(input_data)

    # Initialize blank data dictionary to send to plotting codes later
    data_dict = {}

    # Get list of datasets from cfg
    logger.info("Listing datasets in the data")
    datasets = list_datasets(input_data)

    # Iterate over each dataset
    for dataset in datasets:
        # Select only data from that dataset
        logger.info("Selecting data from %s", dataset)
        selection = select_metadata(input_data, dataset=dataset)

        # Add the dataset to the dictionary with a blank inner dictionary
        data_dict[dataset] = {}

        # Calculations for the Notz-style plot
        logger.info("Calculating data for Notz-style plot")
        sensitivity = calculate_direct_sensitivity(selection)
        # Add to dictionary
        data_dict[dataset]["direct_sensitivity"] = sensitivity

        # Calculations for the Roach-style plot
        logger.info("Calculating data for Roach-style plot")
        trends = calculate_annual_trends(input_data)
        # Add to dictionary
        data_dict[dataset]["siconc_trend"] = trends["si_ann_trend"]
        data_dict[dataset]["tas_trend"] = trends["tas_ann_trend"]
        data_dict[dataset]["direct_r_val"] = trends["direct_r_val"]
        data_dict[dataset]["direct_p_val"] = trends["direct_p_val"]

    # Plot the sensitivities (and save and close the plots)
    logger.info("Creating Notz-style plot")
    notz_style_plot_from_dict(data_dict, titles_and_obs_dict, cfg)
    logger.info("Creating Roach-style plot")
    roach_style_plot_from_dict(data_dict, titles_and_obs_dict, cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
