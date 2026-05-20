"""Diagnostic that shows the sensitivity of sea ice area to global warming."""

import logging
from pathlib import Path

import iris
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import Normalize
from scipy import stats

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    run_diagnostic,
    save_figure,
    select_metadata,
)

logger = logging.getLogger(Path(__file__).stem)


# This is stolen from the example in AutoAssess _plot_mo_metrics.py
def get_provenance_record(cfg, caption):
    """Create a provenance record describing the diagnostic data and plot."""
    filenames = [item["filename"] for item in cfg["input_data"].values()]

    record = {
        "caption": caption,
        "plot_type": "other",
        "authors": [
            "sellar_alistair",
            "parsons_naomi",
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
    """Return the variable group's single iris cube from the data."""
    logger.debug("extracting %s cube from %s", variable_group, data)

    # Load any data in the variable_group
    selection = select_metadata(data, variable_group=variable_group)

    # Ensure there is only one file in the list
    if len(selection) != 1:
        raise ValueError(
            f"None or too many matching files found for {variable_group}",
        )

    # Load the cube, [0] is because selection returns a list
    cube = iris.load_cube(selection[0]["filename"])

    return cube


def calculate_regression(independent, dependent):
    """Use SciPy stats to calculate the least-squares regression."""
    logger.debug(
        "Calculating linear relationship between %s and %s",
        dependent,
        independent,
    )

    # Use SciPy stats to calculate the regression
    # result = slope, intercept, rvalue, pvalue, stderr
    result = stats.linregress(independent, dependent)

    # Return everything
    return result


def calculate_annual_trends(data):
    """
    Calculate annual trends for surface air temperature (tas) and sea ice area (siconc).

    Also used for the r and p values from the regression of siconc as a function of tas.
    """
    logger.debug("calculating annual trends")

    # Load the preprocessed cubes
    si_cube = extract_cube(data, "siconc")
    tas_cube = extract_cube(data, "tas")

    # Calculate the individual trends over time
    years = tas_cube.coord("year").points
    si_trend = calculate_regression(years, si_cube.data)
    tas_trend = calculate_regression(years, tas_cube.data)

    # Calculate the direct regression for r and p values
    direct_regression = calculate_regression(tas_cube.data, si_cube.data)

    dictionary = {
        "si_ann_trend": si_trend.slope,
        "tas_ann_trend": tas_trend.slope,
        "direct_r_val": direct_regression.rvalue,
        "direct_p_val": direct_regression.pvalue,
    }

    return dictionary


def calculate_direct_sensitivity(data):
    """Calculate slope of sea ice area over global mean temperature."""
    logger.debug("calculating sensitivity")

    # Load the preprocessed cubes
    si_cube = extract_cube(data, "siconc")
    tas_cube = extract_cube(data, "tas")

    # Calculate the slope of the direct regression, NOT via time
    sensitivity = calculate_regression(tas_cube.data, si_cube.data)

    return sensitivity.slope


def write_obs_from_cfg(cfg):
    """Write the observations from the recipe to a dictionary."""
    # Initialize the dictionary to hold observations
    obs_dict = {}

    # Add the observation period to the dictionary
    obs_dict["obs_period"] = cfg["observations"]["observation period"]

    # Add a blank dictionary for the Notz-style plot
    obs_dict["notz_style"] = {}

    # Add the observations values to the dictionary
    notz_values = cfg["observations"]["sea ice sensitivity (Notz-style plot)"]
    obs_dict["notz_style"]["mean"] = notz_values["mean"]
    obs_dict["notz_style"]["std_dev"] = notz_values["standard deviation"]
    obs_dict["notz_style"]["plausible"] = notz_values["plausible range"]

    # Add a blank dictionary for the Roach-style plot
    obs_dict["roach_style"] = {}

    # Add each observation point to the dictionary
    roach_values = cfg["observations"]["annual trends (Roach-style plot)"]
    for point in roach_values.keys():
        obs_dict["roach_style"][point] = {}

        # Add the individual values for the observation point
        obs_dict["roach_style"][point]["annual_tas_trend"] = roach_values[
            point
        ]["GMST trend"]
        obs_dict["roach_style"][point]["annual_siconc_trend"] = roach_values[
            point
        ]["SIA trend"]
        obs_dict["roach_style"][point]["r_value"] = roach_values[point][
            "Pearson CC of SIA over GMST"
        ]
        obs_dict["roach_style"][point]["p_value"] = roach_values[point][
            "significance of SIA over GMST"
        ]

    return obs_dict


def create_titles_dict(data, cfg):
    """
    Create a dictionary of appropriate observations and titles.

    Values depend on whether the plot is for the Arctic or Antarctic
    and assume the recipe used September Arctic sea ice data or
    annually mean averaged Antarctic sea ice data
    """
    dictionary = {}

    first_variable = next(iter(data))

    if first_variable["diagnostic"] == "arctic":
        dictionary["titles"] = {
            "notz_fig_title": "September Arctic sea-ice area sensitivity to global mean surface temperature",
            "notz_ax_title": "dSIA/dGMST",
            "notz_plot_filename": "September Arctic sea ice sensitivity",
            "roach_fig_title": "Trends in Annual Mean Temperature And September Arctic Sea Ice",
            "roach_plot_filename": "September Arctic sea ice trends",
        }

    elif first_variable["diagnostic"] == "antarctic":
        dictionary["titles"] = {
            "notz_fig_title": "Annually Meaned Antarctic Sea Ice Sensitivity",
            "notz_ax_title": "dSIA/dGMST",
            "notz_plot_filename": "Annual Antarctic sea ice sensitivity",
            "roach_fig_title": "Trends in Annual Mean Temperature And Annual Antarctic Sea Ice",
            "roach_plot_filename": "Annual Antarctic sea ice trends",
        }

    dictionary["obs"] = write_obs_from_cfg(cfg)

    return dictionary


def write_dictionary_to_csv(cfg, model_dict, filename):
    """Output the model dictionary to a csv file using Pandas."""
    # Read the work directory from the config and create the csv filepath
    csv_filepath = f"{cfg['work_dir']}/{filename}.csv"

    # Write the data to a csv file (via a Pandas DataFrame)
    pd.DataFrame.from_dict(model_dict, orient="index").to_csv(csv_filepath)
    logger.info("Wrote data to %s", csv_filepath)

    # Create a provenance record for the csv file
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(
            csv_filepath,
            get_provenance_record(cfg, "Annual (not decadal) figures"),
        )


def notz_style_plot_from_dict(data_dictionary, titles_dictionary, cfg):
    """Save a plot of sensitivities and observations."""
    # Read from observations dictionary
    obs_years = titles_dictionary["obs"]["obs_period"]
    obs_dict = titles_dictionary["obs"]["notz_style"]
    obs_mean = obs_dict["mean"]
    obs_std_dev = obs_dict["std_dev"]
    obs_plausible = obs_dict["plausible"]

    # Set up the figure
    fig, ax = plt.subplots(figsize=(3.5, 6), layout="constrained")
    fig.suptitle(titles_dictionary["titles"]["notz_fig_title"], wrap=True)
    ax.set_title(
        titles_dictionary["titles"]["notz_ax_title"],
        wrap=True,
        fontsize=10,
    )

    # Iterate over the dictionary
    for dataset, inner_dict in data_dictionary.items():
        ax.plot(
            0.25,
            inner_dict["direct_sensitivity_(notz-style)"],
            color="blue",
            marker="_",
            markersize=20,
        )

        # Label with the dataset if specified
        if inner_dict["label"] == "to_label":
            plt.annotate(
                dataset,
                xy=(0.25, inner_dict["direct_sensitivity_(notz-style)"]),
                xytext=(
                    0.35,
                    inner_dict["direct_sensitivity_(notz-style)"] - 0.05,
                ),
            )

    # Add observations only if obs_mean is a number
    if isinstance(obs_mean, int | float):
        ax.hlines(obs_mean, 0, 1, linestyle="--", color="black", linewidth=2)
        ax.fill_between(
            [0, 1],
            obs_mean - obs_std_dev,
            obs_mean + obs_std_dev,
            facecolor="k",
            alpha=0.15,
        )
        ax.hlines(
            obs_mean + obs_plausible,
            0,
            1,
            linestyle=":",
            color="0.5",
            linewidth=1,
        )
        ax.hlines(
            obs_mean - obs_plausible,
            0,
            1,
            linestyle=":",
            color="0.5",
            linewidth=1,
        )

    # Tidy the figure
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_ylabel(r"dSIA/dGMST ($million \ km^2 \ K^{-1}$)")

    # Create caption based on whether observation mean is presnt
    if isinstance(obs_mean, int | float):
        caption = (
            "Sensitivity of sea ice area to annual mean global warming."
            f"Mean (dashed), standard deviation (shaded) and plausible values from {obs_years}."
        )
    else:
        caption = "Sensitivity of sea ice area to annual mean global warming."

    # Save the figure (also closes it)
    save_figure(
        titles_dictionary["titles"]["notz_plot_filename"],
        get_provenance_record(cfg, caption),
        cfg,
        figure=fig,
        close=True,
    )


def roach_style_plot_from_dict(data_dictionary, titles_dictionary, cfg):
    """Save a plot of trend in SIA against trend in GMST to the given filename."""
    # Set up the figure
    fig, ax = plt.subplots(figsize=(10, 6), layout="constrained")
    fig.suptitle(titles_dictionary["titles"]["roach_fig_title"], wrap=True)

    # Set up for colouring the points
    norm = Normalize(vmin=-1, vmax=1)
    cmap = plt.get_cmap("PiYG_r")

    # Set up the axes
    ax.axhline(color="black", alpha=0.5)
    ax.axvline(color="black", alpha=0.5)
    ax.set_xlabel(r"Trend in GMST ($K \ decade^{-1}$)")
    ax.set_ylabel(r"Trend in SIA ($million \ km^2 \ decade^{-1}$)")

    # Iterate over the dictionary
    for dataset, inner_dict in data_dictionary.items():
        # Determine the position of the point
        x = 10 * inner_dict["annual_tas_trend"]  # for equivalence to decades
        y = (
            10 * inner_dict["annual_siconc_trend"]
        )  # for equivalence to decades

        # Determine the colour of the point
        r_corr = inner_dict["direct_r_val"]

        # Decide if the point should be hatched
        if inner_dict["direct_p_val"] >= 0.05:
            h = 5 * "/"  # This is a hatch pattern
        else:
            h = None

        # Plot the point
        plt.scatter(
            x,
            y,
            marker="o",
            s=150,
            c=[r_corr],
            hatch=h,
            cmap=cmap,
            norm=norm,
        )

        # Label with the dataset if specified
        if inner_dict["label"] == "to_label":
            plt.annotate(dataset, xy=(x, y), xytext=(x + 0.01, y - 0.005))

    # Read from observations dictionary
    obs_years = titles_dictionary["obs"]["obs_period"]
    obs_dict = titles_dictionary["obs"]["roach_style"]

    # Add the observations
    for point in obs_dict.keys():
        # Get the values for the point
        x = obs_dict[point]["annual_tas_trend"]
        y = obs_dict[point]["annual_siconc_trend"]
        r_corr = obs_dict[point]["r_value"]
        p_val = obs_dict[point]["p_value"]

        # Provide a default colour for the point if Pearson coefficient is missing
        if r_corr is None:
            r_corr = 0

        # Provide a pattern for the point if the p-value is present and sufficiently large
        if p_val is not None and p_val >= 0.05:
            h = 5 * "/"  # This is a hatch pattern
        else:
            h = None

        # Plot the point only if both x and y values are provided
        if x is not None and y is not None:
            plt.scatter(
                x,
                y,
                marker="s",
                s=150,
                c=[r_corr],
                hatch=h,
                cmap=cmap,
                norm=norm,
                zorder=0,
                edgecolors="black",
            )

    # Add a colour bar
    plt.colorbar(label="Pearson correlation coefficient")

    # Create caption based on whether observational temp trend is present
    if obs_dict["first point"]["annual_tas_trend"] is not None:
        caption = (
            "Decadal trends of sea ice area and global mean temperature."
            f"Observations from {obs_years} are plotted as squares."
        )
    else:
        caption = "Decadal trends of sea ice area and global mean temperature."

    # Save the figure (also closes it)
    save_figure(
        titles_dictionary["titles"]["roach_plot_filename"],
        get_provenance_record(cfg, caption),
        cfg,
        figure=fig,
        close=True,
    )


def main(cfg):
    """Create two plots per diagnostic from preprocessed data."""
    # Get the data from the cfg
    logger.info("Getting data from the config")
    input_data = cfg["input_data"].values()

    # Titles and observations depend on the diagnostic being plotted
    logger.info("Creating titles and observations dictionary")
    titles_and_obs_dict = create_titles_dict(input_data, cfg)
    logger.debug("Titles and observations dictionary: %s", titles_and_obs_dict)

    # Initialize blank data dictionary to send to plotting codes later
    data_dict = {}

    # Get list of datasets from cfg
    logger.info("Listing datasets in the data")
    datasets = list_datasets(input_data)

    # Iterate over each dataset
    for dataset in datasets:
        # Select only data from that dataset
        logger.debug("Selecting data from %s", dataset)
        selection = select_metadata(input_data, dataset=dataset)

        # Add the dataset to the dictionary with a blank inner dictionary
        data_dict[dataset] = {}

        # Add an entry to determine labelling in plots
        if "label_dataset" in selection[0]:
            data_dict[dataset]["label"] = "to_label"
            logger.info("Dataset %s will be labelled", dataset)
        else:
            data_dict[dataset]["label"] = "unlabelled"
            logger.info("Not labelling dataset %s in plots", dataset)

        # Calculations for the Notz-style plot
        logger.info("Calculating data for Notz-style plot")
        sensitivity = calculate_direct_sensitivity(selection)
        # Add to dictionary
        data_dict[dataset]["direct_sensitivity_(notz-style)"] = sensitivity

        # Calculations for the Roach-style plot
        logger.info("Calculating data for Roach-style plot")
        trends = calculate_annual_trends(selection)
        # Add to dictionary
        data_dict[dataset]["annual_siconc_trend"] = trends["si_ann_trend"]
        data_dict[dataset]["annual_tas_trend"] = trends["tas_ann_trend"]
        data_dict[dataset]["direct_r_val"] = trends["direct_r_val"]
        data_dict[dataset]["direct_p_val"] = trends["direct_p_val"]

    # Add the values to plot to a csv file
    logger.info("Writing values to csv")
    write_dictionary_to_csv(cfg, data_dict, "plotted_values")

    # Plot the sensitivities (and save and close the plots)
    logger.info("Creating Notz-style plot")
    notz_style_plot_from_dict(data_dict, titles_and_obs_dict, cfg)
    logger.info("Creating Roach-style plot")
    roach_style_plot_from_dict(data_dict, titles_and_obs_dict, cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
