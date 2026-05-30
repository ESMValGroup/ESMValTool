"""Diagnostic that shows the sensitivity of sea ice area to global warming."""

import logging
from pathlib import Path
from typing import NamedTuple

import iris
import matplotlib.pyplot as plt
import numpy as np
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


def create_dataset_dict(cfg):
    """Create a dictionary for feeding to dataframe."""
    logger.debug("Creating dataset dictionary.")
    # Initialise a dictionary
    dataset_dict = {}

    tasa_obs = []
    siconc_obs = []

    # Iterate over the data in the cfg
    for section in cfg["input_data"].values():
        # Read relevant facets
        dataset = section["dataset"]
        group = section["variable_group"]

        # Create the dictionary's section for the dataset if it doesn't exist
        dataset_dict.setdefault(dataset, {})

        # Add the models, they are in tas and siconc
        if group == "tas":
            dataset_dict[dataset]["type"] = "model"

            # Add whether to label the models
            if section.get("label_dataset"):
                dataset_dict[dataset]["label"] = "to_label"
                logger.info(
                    "Dataset %s will be labelled",
                    dataset,
                )
            else:
                dataset_dict[dataset]["label"] = "unlabelled"
                logger.info(
                    "Not labelling dataset %s in plots",
                    dataset,
                )

        # List the GMST obs
        elif group == "tasa_obs":
            tasa_obs.append(dataset)

        # List the SIA obs
        elif group == "siconc_obs":
            siconc_obs.append(dataset)

    # Create a string for each pair of observations
    pairs = [f"{t}_v_{s}" for t in tasa_obs for s in siconc_obs]

    # Add the pairs to the dictionary
    for pair in pairs:
        dataset_dict.setdefault(pair, {})
        dataset_dict[pair]["type"] = "multi-obs"

        # Don't label the obs now, but change here later if needed
        dataset_dict[pair]["label"] = "unlabelled"
        logger.info(
            "Not labelling observations pair %s in plots",
            pair,
        )

    return dataset_dict


# Setting this up to query data and obs periods by name later
class Periods(NamedTuple):
    """Class to hold data and observational time ranges from recipe."""

    periods: list[str]
    obs_period: str | None
    data_period: str


def get_start_and_end_years_for_datasets(cfg):
    """Check the recipe or the preprocessed cubes for start and end years."""
    # List as sets as hopefully there is only one of each
    start_years = set()
    end_years = set()

    # Iterate over all the datasets
    datasets = cfg["input_data"].values()
    for dataset in datasets:
        # Read the start year from the recipe if possible
        if "start_year" in dataset and "end_year" in dataset:
            # Add to list
            start_years.add(dataset["start_year"])
            end_years.add(dataset["end_year"])
            logger.info(
                "Found start and end years specified for %s",
                dataset["dataset"],
            )

        # Read from the cube if not
        else:
            cube = iris.load_cube(dataset["filename"])
            points = cube.coord("time").points
            units = cube.coord("time").units
            datetimes = units.num2date(points)
            start_year = min(datetimes).year
            logger.info(
                "Dataset %s starts at year %", dataset["dataset"], start_year
            )
            end_year = max(datetimes).year
            logger.info(
                "Dataset %s ends at year %", dataset["dataset"], end_year
            )

            # Add to list
            start_years.add(start_year)
            end_years.add(end_year)

    # Warn if they aren't all the same
    if len(start_years) > 1 or len(end_years) > 1:
        logger.warning("Not all datasets model the same data period.")

    # Take the largest possible range
    start_year = min(start_years)
    end_year = max(end_years)
    logger.info(
        "Using start year % and end year % as data period",
        start_year,
        end_year,
    )

    return start_year, end_year


def retrieve_periods(cfg):
    """Retrieve both dataset and observational time ranges."""
    data_start, data_end = get_start_and_end_years_for_datasets(cfg)
    data_period = f"{data_start}-{data_end}"
    logger.debug("Data period = %s", data_period)

    # Set obs period to none if not present
    obs_period = None
    periods = [data_period]

    # This is only present in the arctic diagnostic
    if "observations" in cfg:
        obs_start = cfg["observations"]["observation_period"]["start_year"]
        obs_end = cfg["observations"]["observation_period"]["end_year"]
        obs_period = f"{obs_start}-{obs_end}"
        logger.debug("Observation period = %s", obs_period)

        # Return both periods if present and different
        if obs_period != data_period:
            periods = [obs_period, data_period]

    return Periods(
        periods=periods, obs_period=obs_period, data_period=data_period
    )


def create_df_columns(periods):
    """Create column headers for Pandas DataFrame."""
    # In future the historical period will extend beyond 2014
    num_periods = len(periods)
    logger.debug("Dataframe will cover %s periods.", num_periods)

    # Sections for types of regression
    regressions = ["gmst_over_time", "sia_over_time", "sia_over_gmst"]
    num_regressions = len(regressions)

    # Things that will be caluclated
    statistics = ["slope", "r_value", "p_value", "std_err_slope"]
    num_stats = len(statistics)

    # Produces 1979-2014 * 12 then 1979-2021 * 12
    period = [p for p in periods for _ in range(num_regressions * num_stats)]

    # Produces gmst_over_time * 4 then sia_over_time * 4 then sia_over_gmst * 4, two times
    regression = [
        r
        for _ in range(num_periods)
        for r in regressions
        for _ in range(num_stats)
    ]

    # Produces slope then r_value then p_value then std_err_slope, six times
    statistic = [
        s for _ in range(num_periods * num_regressions) for s in statistics
    ]

    # Create main column headers
    data_columns = pd.MultiIndex.from_arrays(
        [period, regression, statistic],
        names=["period", "regression", "statistic"],
    )

    # Add equivalent of single-level columns for dataset info
    first_columns = pd.MultiIndex.from_arrays(
        [["", ""], ["", ""], ["label", "type"]],
        names=["period", "regression", "statistic"],
    )

    # Concatenate columns
    return first_columns.append(data_columns)


def create_blank_dataframe(dataset_dict, columns):
    """Create a Pandas DataFrame for adding values to later on."""
    logger.debug("Creating blank dataframe.")

    # Create DataFrame from dataset_dict
    df = pd.DataFrame.from_dict(dataset_dict, orient="index")

    # Adjustment to account for multi-row headers
    df.columns = pd.MultiIndex.from_tuples(
        [("", "", col) for col in df.columns]
    )

    # Reindex to ensure all required columns are present
    df = df.reindex(columns=columns, fill_value=np.nan)

    # Fill missing "label" and "type" with empty string
    df[("", "", "label")] = df[("", "", "label")].fillna("")
    df[("", "", "type")] = df[("", "", "type")].fillna("")

    return df


def fetch_cube(dataset, variable, time_range, cfg):
    """Fetch a data cube and constrain it using the time range."""
    logger.debug(
        "Fetching cube for dataset: %s, variable: %s, time: %s",
        dataset,
        variable,
        time_range,
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

    # Read the time constraint from the period
    start_year, end_year = time_range.split("-")

    # The constraint will only limit data if data outside the period is available
    time_constraint = iris.Constraint(
        time=lambda cell: int(start_year) <= cell.point.year <= int(end_year)
    )

    # Load the cube using iris with the time constraint
    cube = iris.load_cube(filepath, variable & time_constraint)
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


def calculate_direct_stats(dataset, time_range, cfg):
    """Calculate the direct sensitivity of siconc to tas for a given dataset."""
    logger.debug("Calculating direct sensitivity for dataset %s.", dataset)

    # Fetch the required cubes
    siconc_cube = fetch_cube(dataset, "siconc", time_range, cfg)
    tas_cube = fetch_cube(dataset, "tas", time_range, cfg)

    # regression (tas as independent) gives slope, intercept, rvalue, pvalue, stderr
    return linregress(tas_cube.data, siconc_cube.data)


def calculate_cross_dataset_stats(
    tasa_dataset, siconc_dataset, time_range, cfg
):
    """Calculate the sensitivity of siconc to tasa across (obs) datasets."""
    logger.debug(
        "Calculating cross sensitivity for datasets %s, %s.",
        tasa_dataset,
        siconc_dataset,
    )

    # Fetch the required cubes
    siconc_cube = fetch_cube(siconc_dataset, "siconc", time_range, cfg)
    tasa_cube = fetch_cube(tasa_dataset, "tasa", time_range, cfg)

    # regression (tasa as independent) gives slope, intercept, rvalue, pvalue, stderr
    return linregress(tasa_cube.data, siconc_cube.data)


def add_values_to_df(df, period, cfg):
    """Calculate and write values to the DataFrame."""
    logger.info("Writing values to dataframe.")

    # Calculate all the values for the models
    models = df[df.loc[:, ("", "", "type")] == "model"]
    for dataset, _ in models.iterrows():
        # Calculate annual tas trend
        tas_cube = fetch_cube(dataset, "tas", period, cfg)
        ann_tas_trend = calculate_annual_trend(tas_cube)
        logger.debug("Dataset %s annual tas trend: %s", dataset, ann_tas_trend)
        # Add values to dataframe
        df.loc[dataset, (period, "gmst_over_time", "slope")] = (
            ann_tas_trend.slope
        )
        df.loc[dataset, (period, "gmst_over_time", "r_value")] = (
            ann_tas_trend.rvalue
        )
        df.loc[dataset, (period, "gmst_over_time", "p_value")] = (
            ann_tas_trend.pvalue
        )
        df.loc[dataset, (period, "gmst_over_time", "std_err_slope")] = (
            ann_tas_trend.stderr
        )

        # Calculate annual siconc trend
        siconc_cube = fetch_cube(dataset, "siconc", period, cfg)
        ann_siconc_trend = calculate_annual_trend(siconc_cube)
        logger.debug(
            "Dataset %s annual siconc trend: %s",
            dataset,
            ann_siconc_trend,
        )
        # Add values to dataframe
        df.loc[dataset, (period, "sia_over_time", "slope")] = (
            ann_siconc_trend.slope
        )
        df.loc[dataset, (period, "sia_over_time", "r_value")] = (
            ann_siconc_trend.rvalue
        )
        df.loc[dataset, (period, "sia_over_time", "p_value")] = (
            ann_siconc_trend.pvalue
        )
        df.loc[dataset, (period, "sia_over_time", "std_err_slope")] = (
            ann_siconc_trend.stderr
        )

        # Calculate direct sensitivity of siconc to tas
        direct_sensitivity = calculate_direct_stats(dataset, period, cfg)
        logger.debug("Dataset %s sensitivity: %s", dataset, direct_sensitivity)
        # Add values to dataframe
        df.loc[dataset, (period, "sia_over_gmst", "slope")] = (
            direct_sensitivity.slope
        )
        df.loc[dataset, (period, "sia_over_gmst", "r_value")] = (
            direct_sensitivity.rvalue
        )
        df.loc[dataset, (period, "sia_over_gmst", "p_value")] = (
            direct_sensitivity.pvalue
        )
        df.loc[dataset, (period, "sia_over_gmst", "std_err_slope")] = (
            direct_sensitivity.stderr
        )

    # Calculate all the values for the observations
    obs = df[df.loc[:, ("", "", "type")] == "multi-obs"]
    for combined_name, _ in obs.iterrows():
        # Calculate annual tasa trend
        gmst_dataset = combined_name.split("_v_")[0]
        tasa_cube = fetch_cube(gmst_dataset, "tasa", period, cfg)
        ann_tasa_trend = calculate_annual_trend(tasa_cube)
        logger.debug(
            "Dataset %s tasa annual trend: %s", gmst_dataset, ann_tasa_trend
        )
        # Add values to dataframe
        df.loc[combined_name, (period, "gmst_over_time", "slope")] = (
            ann_tasa_trend.slope
        )
        df.loc[combined_name, (period, "gmst_over_time", "r_value")] = (
            ann_tasa_trend.rvalue
        )
        df.loc[combined_name, (period, "gmst_over_time", "p_value")] = (
            ann_tasa_trend.pvalue
        )
        df.loc[combined_name, (period, "gmst_over_time", "std_err_slope")] = (
            ann_tasa_trend.stderr
        )

        # Calculate annual siconc trend
        sia_dataset = combined_name.split("_v_")[1]
        siconc_cube = fetch_cube(sia_dataset, "siconc", period, cfg)
        ann_siconc_trend = calculate_annual_trend(siconc_cube)
        logger.debug(
            "Dataset %s siconc annual trend: %s", sia_dataset, ann_siconc_trend
        )
        # Add values to dataframe
        df.loc[combined_name, (period, "sia_over_time", "slope")] = (
            ann_siconc_trend.slope
        )
        df.loc[combined_name, (period, "sia_over_time", "r_value")] = (
            ann_siconc_trend.rvalue
        )
        df.loc[combined_name, (period, "sia_over_time", "p_value")] = (
            ann_siconc_trend.pvalue
        )
        df.loc[combined_name, (period, "sia_over_time", "std_err_slope")] = (
            ann_siconc_trend.stderr
        )

        # Calculate sensitivity of siconc to tasa
        cross_dataset_stats = calculate_cross_dataset_stats(
            gmst_dataset, sia_dataset, period, cfg
        )
        logger.debug(
            "Obs pair %s sensitivity: %s", combined_name, cross_dataset_stats
        )
        # Add values to dataframe
        df.loc[combined_name, (period, "sia_over_gmst", "slope")] = (
            cross_dataset_stats.slope
        )
        df.loc[combined_name, (period, "sia_over_gmst", "r_value")] = (
            cross_dataset_stats.rvalue
        )
        df.loc[combined_name, (period, "sia_over_gmst", "p_value")] = (
            cross_dataset_stats.pvalue
        )
        df.loc[combined_name, (period, "sia_over_gmst", "std_err_slope")] = (
            cross_dataset_stats.stderr
        )

    return df


def write_df_to_csv(df, filename, cfg):
    """Copy DataFrame to a csv file."""
    logger.debug("Dataframe to write:\n %s", df)

    # Create the csv filepath using info from the config
    csv_filepath = f"{cfg['work_dir']}/{filename}.csv"

    # Write the data to a csv file
    df.to_csv(csv_filepath)
    logger.info("Wrote data to %s", csv_filepath)


def roach_style_plot_from_df(df, cfg):
    """Save a plot of trend in SIA against trend in GMST to the given filename."""
    logger.debug("DataFrame for Roach-style plot:\n %s", df)

    # Look up variable to determine plot title
    data = cfg["input_data"].values()
    first_variable = next(iter(data))

    # Set plot title and filename
    if first_variable["diagnostic"] == "arctic":
        title = (
            "Trends in Annual Mean Temperature And September Arctic Sea Ice"
        )
        save_as = "September Arctic sea ice trends"
    elif first_variable["diagnostic"] == "antarctic":
        title = (
            "Trends in Annual Mean Temperature And Annual Antarctic Sea Ice"
        )
        save_as = "Annual Antarctic sea ice trends"

    # Set up the figure
    fig, ax = plt.subplots(figsize=(10, 6), layout="constrained")
    fig.suptitle(title, wrap=True)

    # Set up the axes
    ax.axhline(color="black", alpha=0.5)
    ax.axvline(color="black", alpha=0.5)
    ax.set_xlabel(r"Trend in GMST ($K \ decade^{-1}$)")
    ax.set_ylabel(r"Trend in SIA ($million \ km^2 \ decade^{-1}$)")

    # Set up for colouring the points
    norm = Normalize(vmin=-1, vmax=1)
    cmap = plt.get_cmap("PiYG_r")

    # Choose p-value to hatch
    min_p_to_hatch = 0.05

    # Only plot the entire data period in this style
    data_period = retrieve_periods(cfg).data_period
    ax.set_title(data_period, loc="center")

    # Iterate over the values in the dataframe
    for dataset, _ in df.iterrows():
        # Look up the relevant values
        gmst_trend = df.loc[dataset, (data_period, "gmst_over_time", "slope")]
        sia_trend = df.loc[dataset, (data_period, "sia_over_time", "slope")]
        r_cross = df.loc[dataset, (data_period, "sia_over_gmst", "r_value")]
        p_sia = df.loc[dataset, (data_period, "sia_over_time", "p_value")]

        # The decadal figures (ten times bigger) are plotted
        x = 10 * gmst_trend
        y = 10 * sia_trend

        # Decide if the point should be hatched
        hatch = 5 * "/" if p_sia >= min_p_to_hatch else None

        # Shape is different for obs
        if df.loc[dataset, ("", "", "type")] == "model":
            shape = "o"
            edgecolor = None
            order = None
        else:
            shape = "s"
            edgecolor = "black"
            order = 0

        # Plot the point
        plt.scatter(
            x,
            y,
            marker=shape,
            s=150,
            c=[r_cross],
            hatch=hatch,
            cmap=cmap,
            norm=norm,
            edgecolors=edgecolor,
            zorder=order,
        )

        # Label with the dataset if specified
        if df.loc[dataset, ("", "", "label")] == "to_label":
            plt.annotate(dataset, xy=(x, y), xytext=(x + 0.01, y - 0.005))

    # Add a colour bar
    plt.colorbar(label="Pearson correlation coefficient")

    # Save the figure (also closes it)
    caption = "Decadal trends of sea ice area and global mean temperature."
    save_figure(
        save_as,
        get_provenance_record(cfg, caption),
        cfg,
        figure=fig,
        close=True,
    )


def notz_style_plot_from_df(df, cfg):
    """Save a plot of sensitivities and observations for model datasets."""
    logger.debug("DataFrame for Notz-style plot:\n %s", df)

    # Look up variable to determine plot title
    data = cfg["input_data"].values()
    first_variable = next(iter(data))

    # Set plot title and filename
    if first_variable["diagnostic"] == "arctic":
        title = "September Arctic sea-ice area sensitivity\ndSIA/dGMST"
        save_as = "September Arctic sea ice sensitivity"
    elif first_variable["diagnostic"] == "antarctic":
        title = "Annual Antarctic sea-ice area sensitivity\ndSIA/dGMST"
        save_as = "Annual Antarctic sea ice sensitivity"

    # Caption changes if obs are present
    caption = "Sensitivity of sea ice area to annual mean global warming."

    # Retrieve obs from config if present
    obs_period = retrieve_periods(cfg).obs_period
    if obs_period is not None:
        obs_mean = cfg["observations"]["sea_ice_sensitivity"]["mean"]
        obs_std_dev = cfg["observations"]["sea_ice_sensitivity"][
            "standard_deviation"
        ]
        obs_plausible = cfg["observations"]["sea_ice_sensitivity"][
            "plausible_range"
        ]
        # Also changes caption
        caption = (
            "Sensitivity of sea ice area to annual mean global warming.\n"
            "Mean (dashed), standard deviation (shaded) and plausible (dotted)"
            f" values from {obs_period} are shown in grey."
        )

    # Set up the figure
    fig, ax = plt.subplots(figsize=(3.5, 6), layout="constrained")
    fig.suptitle(title, wrap=True)

    # Set up the axes
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_ylabel(r"dSIA/dGMST ($million \ km^2 \ K^{-1}$)")

    # Check how many periods there are, for titles and positioning
    cfg_periods = retrieve_periods(cfg)
    num_periods = len(cfg_periods.periods)
    obs_period = cfg_periods.obs_period
    data_period = cfg_periods.data_period

    # Set up the spacing and titles accordingly
    if num_periods > 1:
        # Put obs period on the left (as it should end earlier)
        ax.set_title(obs_period, loc="left", fontsize=10, x=0.1)
        ax.set_title(data_period, loc="right", fontsize=10, x=0.9)

        # Four columns used for data
        lhs_model_x = 0.15
        lhs_obs_x = 0.35
        rhs_model_x = 0.65
        rhs_obs_x = 0.85

        # Labelling was too big
        label_font_size = 7
    else:
        # Only one axes title
        ax.set_title(data_period, loc="center")

        # Two columns used for data
        lhs_model_x = 0.35
        lhs_obs_x = 0.65

        # Labelling need not shrink
        label_font_size = 10

    # Add pre-defined observation values to the first half / all, if present
    section_width = 1 / num_periods
    if obs_period is not None:
        ax.hlines(
            obs_mean,
            0,
            section_width,
            linestyle="--",
            color="black",
            linewidth=2,
        )
        ax.fill_between(
            [0, section_width],
            obs_mean - obs_std_dev,
            obs_mean + obs_std_dev,
            facecolor="k",
            alpha=0.15,
        )
        ax.hlines(
            obs_mean + obs_plausible,
            0,
            section_width,
            linestyle=":",
            color="0.5",
            linewidth=1,
        )
        ax.hlines(
            obs_mean - obs_plausible,
            0,
            section_width,
            linestyle=":",
            color="0.5",
            linewidth=1,
        )

    # Function that may be called once or twice
    def add_points(ax, period, data_x, obs_x, font_size):
        # Iterate over the values in the dataframe
        for dataset, _ in df.iterrows():
            # Look up the sensitivity of SIA to GMST for the dataset
            sensitivity = df.loc[dataset, (period, "sia_over_gmst", "slope")]

            # Plotting for models
            if df.loc[dataset, ("", "", "type")] == "model":
                ax.plot(
                    data_x,
                    sensitivity,
                    color="blue",
                    marker="_",
                    markersize=20,
                )

                # Label with the dataset if specified, offset correct by eye
                if df.loc[dataset, ("", "", "label")] == "to_label":
                    plt.annotate(
                        dataset,
                        xy=(data_x, sensitivity),
                        xytext=(
                            data_x + 0.1,
                            sensitivity - 0.05,
                        ),
                        fontsize=font_size,
                    )

            # Plotting for computed observations
            if df.loc[dataset, ("", "", "type")] == "multi-obs":
                ax.plot(
                    obs_x,
                    sensitivity,
                    color="orange",
                    marker="_",
                    markersize=20,
                )

                # Shade around computed observations
                std_err = df.loc[
                    dataset, (period, "sia_over_gmst", "std_err_slope")
                ]
                ax.fill_between(
                    [obs_x - 0.05, obs_x + 0.05],
                    sensitivity - std_err,
                    sensitivity + std_err,
                    facecolor="orange",
                    alpha=0.15,
                )

    # Add the obs data period, or the full data period if there are no obs
    if obs_period is not None:
        add_points(ax, obs_period, lhs_model_x, lhs_obs_x, label_font_size)
    else:
        add_points(ax, data_period, lhs_model_x, lhs_obs_x, label_font_size)

    # Add the full data period on the right-hand side if different
    if num_periods > 1:
        add_points(ax, data_period, rhs_model_x, rhs_obs_x, label_font_size)

    # Save the figure (also closes it)
    save_figure(
        save_as,
        get_provenance_record(cfg, caption),
        cfg,
        figure=fig,
        close=True,
    )


def main(cfg):
    """Create two plots and one csv file per diagnostic."""
    # Log config dictionary
    logger.debug("cfg:\n%s", cfg)

    # Look at the datasets in the config object
    logger.info("Reading datasets.")
    datasets = create_dataset_dict(cfg)

    # Retrieve the data periods (in case obs period is different)
    logger.info("Checking data periods.")
    data_periods = retrieve_periods(cfg).periods

    # Create a dataframe with the right columns
    logger.info("Creating dataframe.")
    columns = create_df_columns(data_periods)
    df = create_blank_dataframe(datasets, columns)

    # Add the data for each period
    for data_period in data_periods:
        logger.info(
            "Calculating and writing values for period %s.", data_period
        )
        filled = add_values_to_df(df, data_period, cfg)

    # Write the dataframe to file, with provenance
    logger.info("Writing dataframe to csv file.")
    filename = "data_values"
    write_df_to_csv(filled, filename, cfg)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(
            f"{cfg['work_dir']}/{filename}",
            get_provenance_record(cfg, "Annual (not decadal) figures"),
        )

    # Plot the 2D figure
    logger.info("Creating Roach-style plot")
    roach_style_plot_from_df(filled, cfg)

    # Plot the 1D figure
    logger.info("Creating Notz-style plot")
    notz_style_plot_from_df(filled, cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
