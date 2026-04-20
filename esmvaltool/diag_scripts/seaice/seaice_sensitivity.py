"""Diagnostic that shows the sensitivity of sea ice area to global warming."""

import logging
from pathlib import Path

import iris
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import Normalize
from scipy.stats import linregress

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    run_diagnostic,
    save_figure,
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


def create_category_dict(cfg):
    """Create a structured dictionary for adding values to later on."""
    logger.debug("Creating blank dictionary.")
    # Create blank dictionary of correct structure
    category_dict = {
        "models": {},
        "tasa_obs": {},
        "siconc_obs": {},
        "cross-dataset-obs": {},
    }

    # Read the data from the config object
    input_data = cfg["input_data"].values()

    # Iterate over the datasets to add to the dictionary
    for section in input_data:
        # Check for tasa observations
        if section["variable_group"] == "tasa_obs":
            category_dict["tasa_obs"][section["dataset"]] = {}

        # Check for siconc observations
        elif section["variable_group"] == "siconc_obs":
            category_dict["siconc_obs"][section["dataset"]] = {}

        # Everything else should be a model
        else:
            # Add the model dataset if not already present (appears twice, for tas and siconc)
            if section["dataset"] not in category_dict["models"]:
                category_dict["models"][section["dataset"]] = {}

            # Add labelling info
            if section.get("label_dataset"):
                category_dict["models"][section["dataset"]]["label"] = (
                    "to_label"
                )
                logger.info("Dataset %s will be labelled", section["dataset"])
            else:
                category_dict["models"][section["dataset"]]["label"] = (
                    "unlabelled"
                )
                logger.info(
                    "Not labelling dataset %s in plots",
                    section["dataset"],
                )

    return category_dict


def fetch_cube(dataset, variable, cfg):
    """Fetch a data cube for a dataset and variable using info from the config."""
    logger.debug(
        "Fetching cube for dataset: %s, variable: %s",
        dataset,
        variable,
    )

    # Read the data from the config object
    input_data = cfg["input_data"].values()

    # Find the correct filepath for the dataset
    for section in input_data:
        # Check the dataset AND variable matches as models have two entries
        # Only matching the first three letters to avoid issues with sic vs siconc
        if (
            section["dataset"] == dataset
            and section["short_name"][:3] == variable[:3]
        ):
            filepath = section["filename"]
            break

    # Load the cube using iris
    cube = iris.load_cube(filepath, variable)
    return cube


def calculate_annual_trend(cube):
    """Calculate the linear trend of a cube over time using scipy.stats.linregress."""
    logger.debug("Calculating annual trend for cube %s.", cube.name())

    # Depending on preprocessor, coord may be 'year' or 'time'
    if "year" in cube.coords():
        no_years = list(range(len(cube.coord("years").points)))
    else:
        no_years = list(range(len(cube.coord("time").points)))

    # Return all of slope, intercept, rvalue, pvalue, stderr as hatching needs p
    return linregress(no_years, cube.data)


def calculate_direct_stats(dataset, cfg):
    """Calculate the direct sensitivity of siconc to tas for a given dataset."""
    logger.debug("Calculating direct sensitivity for dataset %s.", dataset)

    # Fetch the required cubes
    siconc_cube = fetch_cube(dataset, "siconc", cfg)
    tas_cube = fetch_cube(dataset, "tas", cfg)

    # regression (tas as independent) gives slope, intercept, rvalue, pvalue, stderr
    return linregress(tas_cube.data, siconc_cube.data)


def calculate_cross_dataset_stats(tasa_dataset, siconc_dataset, cfg):
    """Calculate the sensitivity of siconc to tasa across (obs) datasets."""
    logger.debug(
        "Calculating cross sensitivity for datasets %s, %s.",
        tasa_dataset,
        siconc_dataset,
    )

    # Fetch the required cubes
    siconc_cube = fetch_cube(siconc_dataset, "siconc", cfg)
    tasa_cube = fetch_cube(tasa_dataset, "tasa", cfg)

    # regression (tasa as independent) gives slope, intercept, rvalue, pvalue, stderr
    return linregress(tasa_cube.data, siconc_cube.data)


def write_values_to_dict(data_dict, cfg):
    """Calculate and write values to the structured dictionary."""
    logger.debug("Writing values to dictionary.")

    # Calculate all the values for the models
    for model_dataset in data_dict["models"]:
        # Calculate annual tas trend
        tas_cube = fetch_cube(model_dataset, "tas", cfg)
        ann_tas_trend = calculate_annual_trend(tas_cube)
        data_dict["models"][model_dataset]["annual_tas_trend"] = (
            ann_tas_trend.slope
        )

        # Calculate annual siconc trend
        siconc_cube = fetch_cube(model_dataset, "siconc", cfg)
        ann_siconc_trend = calculate_annual_trend(siconc_cube)
        # Add the slope for 2D positioning
        data_dict["models"][model_dataset]["annual_siconc_trend"] = (
            ann_siconc_trend.slope
        )
        # Add the p-value for hatching in 2D plot
        data_dict["models"][model_dataset]["annual_siconc_p-value"] = (
            ann_siconc_trend.pvalue
        )

        # Calculate direct sensitivity of siconc to tas
        direct_sensitivity = calculate_direct_stats(model_dataset, cfg)
        # Add the slope for sensitivity plot
        data_dict["models"][model_dataset]["direct_sensitivity"] = (
            direct_sensitivity.slope
        )
        # Add the r-value for colouring in 2D plot
        data_dict["models"][model_dataset]["direct_r-value"] = (
            direct_sensitivity.rvalue
        )

    # Calculate just the tasa trend for the tasa observations
    for obs_dataset in data_dict["tasa_obs"]:
        # Calculate annual tas trend
        tasa_cube = fetch_cube(obs_dataset, "tasa", cfg)
        ann_tasa_trend = calculate_annual_trend(tasa_cube)
        # Add the slope for 2D positioning
        data_dict["tasa_obs"][obs_dataset]["annual_tas_trend"] = (
            ann_tasa_trend.slope
        )

    # Calculate the siconc slope and p value  for the siconc observations
    for obs_dataset in data_dict["siconc_obs"]:
        # Calculate annual siconc trend
        siconc_cube = fetch_cube(obs_dataset, "siconc", cfg)
        ann_siconc_trend = calculate_annual_trend(siconc_cube)
        # Add the slope for 2D positioning
        data_dict["siconc_obs"][obs_dataset]["annual_siconc_trend"] = (
            ann_siconc_trend.slope
        )
        # Add the p-value for hatching in 2D plot
        data_dict["siconc_obs"][obs_dataset]["annual_siconc_p-value"] = (
            ann_siconc_trend.pvalue
        )

    # Calculate cross-dataset statistics between tasa and siconc observations
    for tasa_dataset in data_dict["tasa_obs"]:
        for siconc_dataset in data_dict["siconc_obs"]:
            # Determine structure of dictionary to store values
            key_name = f"{siconc_dataset}_to_{tasa_dataset}"
            data_dict["cross-dataset-obs"][key_name] = {}
            inner_dict = data_dict["cross-dataset-obs"][key_name]

            # Calculate cross-dataset sensitivity of siconc to tasa
            cross_sensitivity = calculate_cross_dataset_stats(
                tasa_dataset,
                siconc_dataset,
                cfg,
            )
            # Add the r-value for colouring in 2D plot
            inner_dict["direct_r-value"] = cross_sensitivity.rvalue
            # Store the direct sensitivity as the calculation was run anyway
            inner_dict["direct_sensitivity"] = cross_sensitivity.slope

    return data_dict


def write_dictionary_to_csv(sub_dict, filename, cfg):
    """
    Output a section of data dictionary to a csv file using Pandas.

    Only sections of the dictionary should be written at a time as otherwise
    the structure is too complex to easily convert to a DataFrame.
    """
    logger.debug("Writing dictionary to csv file.")

    # Create the csv filepath using info from the config
    csv_filepath = f"{cfg['work_dir']}/{filename}.csv"

    # Write the data to a csv file (via a Pandas DataFrame)
    dataframe = pd.DataFrame.from_dict(sub_dict, orient="index")
    dataframe.to_csv(csv_filepath)
    logger.info("Wrote data to %s", csv_filepath)


def write_obs_from_cfg(cfg):
    """Write the Notz-style observations from the recipe to a dictionary."""
    logger.debug("Writing observations from config file.")

    # Initialize the dictionary with observation period
    obs_dict = {"obs_period": cfg["observations"]["observation period"]}

    # Add the observations values to the dictionary
    notz_values = cfg["observations"]["sea ice sensitivity (Notz-style plot)"]
    obs_dict["mean"] = notz_values["mean"]
    obs_dict["std_dev"] = notz_values["standard deviation"]
    obs_dict["plausible"] = notz_values["plausible range"]

    return obs_dict


def create_titles_dict(cfg):
    """
    Create a dictionary of appropriate titles and hardcoded observations.
    Values depend on whether the plot is for the Arctic or Antarctic
    and assume the recipe used September Arctic sea ice data or
    annually mean averaged Antarctic sea ice data
    """
    logger.debug("Creating titles dictionary.")
    dictionary = {}

    data = cfg["input_data"].values()
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


def notz_style_plot_from_dict(data_dictionary, titles_dictionary, cfg):
    """Save a plot of sensitivities and observations for model datasets."""
    # Read from observations dictionary
    obs_years = titles_dictionary["obs"]["obs_period"]
    obs_dict = titles_dictionary["obs"]
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
            inner_dict["direct_sensitivity"],
            color="blue",
            marker="_",
            markersize=20,
        )

        # Label with the dataset if specified, offset correct by eye
        if inner_dict["label"] == "to_label":
            plt.annotate(
                dataset,
                xy=(0.25, inner_dict["direct_sensitivity"]),
                xytext=(
                    0.35,
                    inner_dict["direct_sensitivity"] - 0.05,
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

    # Create caption based on whether observation mean is present
    if isinstance(obs_mean, int | float):
        caption = (
            "Sensitivity of sea ice area to annual mean global warming."
            f"\nMean (dashed), standard deviation (shaded) and plausible values from {obs_years}."
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

    # Choose p-value to hatch
    min_p_to_hatch = 0.05

    # Set up the axes
    ax.axhline(color="black", alpha=0.5)
    ax.axvline(color="black", alpha=0.5)
    ax.set_xlabel(r"Trend in GMST ($K \ decade^{-1}$)")
    ax.set_ylabel(r"Trend in SIA ($million \ km^2 \ decade^{-1}$)")

    # Iterate over the models sub-dictionary
    for dataset, inner_dict in data_dictionary["models"].items():
        # Determine the position of the point
        x = 10 * inner_dict["annual_tas_trend"]  # for equivalence to decades
        y = (
            10 * inner_dict["annual_siconc_trend"]
        )  # for equivalence to decades

        # Determine the colour of the point
        r_corr = inner_dict["direct_r-value"]

        # Decide if the point should be hatched
        if inner_dict["annual_siconc_p-value"] >= min_p_to_hatch:
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

    # Add the observations
    siconc_dict = data_dictionary["siconc_obs"]
    tasa_dict = data_dictionary["tasa_obs"]

    # Iterate over the pairs in cross-dataset-obs
    for pair, inner_dict in data_dictionary["cross-dataset-obs"].items():
        # Retrieve the names of the datasets from the pair string
        siconc_ds, tasa_ds = pair.split("_to_")

        # Determine the position of the point, from other dictionaries
        x = (
            10 * tasa_dict[tasa_ds]["annual_tas_trend"]
        )  # This was labelled as tas, not tasa
        y = 10 * siconc_dict[siconc_ds]["annual_siconc_trend"]

        # Determine the colour of the point from the inner dictionary
        r_corr = inner_dict["direct_r-value"]

        # Decide if the point should be hatched
        if siconc_dict[siconc_ds]["annual_siconc_p-value"] >= min_p_to_hatch:
            h = 5 * "/"
        else:
            h = None

        # Plot the point
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

    # Save the figure (also closes it)
    caption = "Decadal trends of sea ice area and global mean temperature."
    save_figure(
        titles_dictionary["titles"]["roach_plot_filename"],
        get_provenance_record(cfg, caption),
        cfg,
        figure=fig,
        close=True,
    )


def main(cfg):
    # Create the structured dictionary
    data_dict = create_category_dict(cfg)

    # Calculate and write values to the dictionary
    logger.info("Calculating and writing values to dictionary.")
    data_dict = write_values_to_dict(data_dict, cfg)

    # Write the model and obs dictionaries to csv files
    logger.info("Writing dictionaries to csv files.")
    write_dictionary_to_csv(data_dict["models"], "models_values", cfg)
    write_dictionary_to_csv(data_dict["tasa_obs"], "tasa_obs_values", cfg)
    write_dictionary_to_csv(data_dict["siconc_obs"], "siconc_obs_values", cfg)

    # Write the cross-dataset obs dictionary to csv files (separately for each pair)
    for pair in data_dict["cross-dataset-obs"]:
        data_dict["cross-dataset-obs"][pair] = data_dict["cross-dataset-obs"][
            pair
        ]
        write_dictionary_to_csv(
            data_dict["cross-dataset-obs"][pair],
            f"{pair}",
            cfg,
        )

    # Create a single provenance record for the csv files
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(
            f"{cfg['work_dir']}/figures_as_csv",
            get_provenance_record(cfg, "Annual (not decadal) figures"),
        )

    # Titles and observations depend on the diagnostic being plotted
    logger.info("Creating titles and observations dictionary")
    titles_and_obs_dict = create_titles_dict(cfg)

    # Plot the sensitivities, uses model data only (and obs from recipe)
    logger.info("Creating Notz-style plot")
    notz_style_plot_from_dict(data_dict["models"], titles_and_obs_dict, cfg)
    logger.info("Creating Roach-style plot")
    roach_style_plot_from_dict(data_dict, titles_and_obs_dict, cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
