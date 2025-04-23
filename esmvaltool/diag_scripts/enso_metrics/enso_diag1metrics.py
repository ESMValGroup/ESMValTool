"""Compute and plot ENSO characteristics metrics."""

import logging
import os

import iris
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import (
    climate_statistics,
    extract_month,
    extract_season,
    mask_above_threshold,
    mask_below_threshold,
)
from scipy.stats import linregress, skew

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
    select_metadata,
)

logger = logging.getLogger(os.path.basename(__file__))


def plot_level1(input_data, metricval, y_label, title, dtls):
    """Create plots for output data."""
    figure = plt.figure(figsize=(10, 6), dpi=300)

    if title in ["ENSO pattern", "ENSO lifecycle"]:
        # model first
        plt.plot(*input_data[0], label=dtls[0])
        plt.plot(*input_data[1], label=f"ref: {dtls[1]}", color="black")
        plt.text(
            0.5,
            0.95,
            f"RMSE: {metricval:.2f}",
            fontsize=12,
            ha="center",
            transform=plt.gca().transAxes,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
        )
        plt.legend()
    else:
        plt.scatter(
            range(len(input_data)),
            input_data,
            c=["black", "blue"],
            marker="D",
        )
        # obs first
        plt.xlim(-0.5, 2)
        plt.xticks([])
        plt.text(
            0.75,
            0.95,
            f"* {dtls[0]}",
            color="blue",
            transform=plt.gca().transAxes,
        )
        plt.text(
            0.75,
            0.9,
            f"* ref: {dtls[1]}",
            color="black",
            transform=plt.gca().transAxes,
        )
        plt.text(
            0.75,
            0.8,
            f"metric(%): {metricval:.2f}",
            fontsize=12,
            transform=plt.gca().transAxes,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
        )

    plt.title(title)  # metric name
    plt.grid(linestyle="--")
    plt.ylabel(y_label)

    if title == "ENSO pattern":
        plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(format_lon))
    elif title == "ENSO lifecycle":
        plt.axhline(y=0, color="black", linewidth=2)
        xticks = np.arange(1, 73, 6) - 36  # Adjust for lead/lag months
        xtick_labels = ["Jan", "Jul"] * (len(xticks) // 2)
        plt.xticks(xticks, xtick_labels)
        plt.yticks(np.arange(-2, 2.5, step=1))

    logger.info("%s : metric: %f", dtls[0], metricval)

    return figure


def lin_regress(cube_ssta, cube_nino34):
    """Linear regression for 1D output."""
    slope_ls = []
    for lon_slice in cube_ssta.slices(["time"]):
        res = linregress(cube_nino34.data, lon_slice.data)
        slope_ls.append(res[0])

    return cube_ssta.coord("longitude").points, slope_ls


def sst_regressed(n34_cube):
    """Regression function for sst_time_series on sst_enso(dec)."""
    n34_dec = extract_month(n34_cube, 12)
    n34_dec_years = [
        n34_dec.coord("time").units.num2date(time).year
        for time in n34_dec.coord("time").points
    ]
    # Ensure that the selected years are not the last years
    n34_selected = []
    for year in n34_dec_years[3:-3]:
        enso_epoch = [year - 2, year - 1, year, year + 1, year + 2, year + 3]

        # Select the data for the current year and append it to n34_selected
        year_enso = iris.Constraint(
            time=lambda cell, enso_epoch=enso_epoch: cell.point.year
            in enso_epoch,
        )
        cube_2 = n34_cube.extract(year_enso)
        n34_selected.append(cube_2.data.data)

    event_constr = iris.Constraint(
        time=lambda cell: cell.point.year in n34_dec_years[3:-3],
    )
    n34_dec_ct = n34_dec.extract(event_constr)

    # 1) linear regression of sst_time_series on sst_enso
    a_data = np.array(n34_selected)
    b_data = n34_dec_ct.data
    b_with_intercept = np.vstack([b_data, np.ones_like(b_data)]).T
    coefs, _, _, _ = np.linalg.lstsq(b_with_intercept, a_data, rcond=None)

    return coefs[0]


def data_to_cube(line_point, in_cube, metric):
    """Translate computed data to cube to save."""
    if in_cube is None:
        cube = iris.cube.Cube(
            line_point,
            long_name=metric,
        )
    else:
        if isinstance(in_cube, iris.cube.Cube):
            coord = in_cube.coord("longitude")
        else:
            coord = iris.coords.DimCoord(
                in_cube,
                long_name="months 6 year ENSO epoch",
            )
        cube = iris.cube.Cube(
            line_point,
            long_name=metric,
            dim_coords_and_dims=[(coord, 0)],
        )
    return cube


def compute_enso_metrics(input_pair, dt_ls, var_group, metric):
    """Compute values for each of the ENSO metrics.

    Takes groupings of datasets required for each ENSO metric sorted
    and iterated through in the main function to compute the metric.

    Args:
        input_pair: List of dictionaries [obs_datasets, model_datasets]
            where dictionary key of dataset is variable group(preprocessor).
        dt_ls: List of dataset names.
            Used for labels in plots.
        var_group: List referring to preprocessed group used for the metric;
            list length is 1,
            or 2 for pattern and diversity metrics for linear regrssion
        metric: Name of metric to calculate.
            eg. 09pattern, 10lifecyle.
    """
    data_values = []
    data_to_save = []
    fig, val = None, None
    if metric == "09pattern":
        model_ssta = input_pair[1][var_group[1]]
        model_nino34 = input_pair[1][var_group[0]]
        reg_mod = lin_regress(model_ssta, model_nino34)
        obs_ssta = input_pair[0][var_group[1]]
        obs_nino34 = input_pair[0][var_group[0]]
        reg_obs = lin_regress(obs_ssta, obs_nino34)

        val = np.sqrt(
            np.mean((np.array(reg_obs[1]) - np.array(reg_mod[1])) ** 2),
        )

        fig = plot_level1(
            [reg_mod, reg_obs],
            val,
            "reg(ENSO SSTA, SSTA)",
            "ENSO pattern",
            dt_ls,
        )

        data_to_save.append(data_to_cube(reg_mod[1], model_ssta, metric))
        data_to_save.append(data_to_cube(reg_obs[1], obs_ssta, metric))

    elif metric == "10lifecycle":
        model = sst_regressed(input_pair[1][var_group[0]])
        obs = sst_regressed(input_pair[0][var_group[0]])
        val = np.sqrt(np.mean((obs - model) ** 2))
        months = np.arange(1, 73) - 36

        fig = plot_level1(
            [(months, model), (months, obs)],
            val,
            "Degree C / C",
            "ENSO lifecycle",
            dt_ls,
        )

        data_to_save.append(data_to_cube(model, months, metric))
        data_to_save.append(data_to_cube(obs, months, metric))

    elif metric == "11amplitude":
        data_values = [
            input_pair[1][var_group[0]].data.item(),
            input_pair[0][var_group[0]].data.item(),
        ]
        val = compute(data_values[1], data_values[0])

        fig = plot_level1(
            data_values,
            val,
            "SSTA std (°C)",
            "ENSO amplitude",
            dt_ls,
        )

    elif metric == "12seasonality":
        for datas in input_pair:  # obs 0, mod 1
            preproc = {}
            for season in ["NDJ", "MAM"]:
                cube = extract_season(datas[var_group[0]], season)
                cube = climate_statistics(
                    cube,
                    operator="std_dev",
                    period="full",
                )
                preproc[season] = cube.data

            data_values.append(preproc["NDJ"] / preproc["MAM"])

        val = compute(data_values[1], data_values[0])
        fig = plot_level1(
            data_values,
            val,
            "SSTA std (NDJ/MAM)(°C/°C)",
            "ENSO seasonality",
            dt_ls,
        )
        data_to_save.append(data_to_cube(data_values[1], None, metric))
        data_to_save.append(data_to_cube(data_values[0], None, metric))

    elif metric == "13asymmetry":
        model_skew = skew(input_pair[1][var_group[0]].data, axis=0)
        obs_skew = skew(input_pair[0][var_group[0]].data, axis=0)
        data_values = [model_skew, obs_skew]

        val = compute(data_values[1], data_values[0])
        fig = plot_level1(
            data_values,
            val,
            "SSTA skewness(°C)",
            "ENSO skewness",
            dt_ls,
        )
        data_to_save.append(data_to_cube(model_skew, None, metric))
        data_to_save.append(data_to_cube(obs_skew, None, metric))

    elif metric == "14duration":
        model = sst_regressed(input_pair[1][var_group[0]])
        obs = sst_regressed(input_pair[0][var_group[0]])

        months = np.arange(1, 73) - 36
        # Calculate the number of months where slope > 0.25
        within_range = (months >= -30) & (months <= 30)
        for slopes in [model, obs]:
            slope_above_025 = slopes[within_range] > 0.25
            data_values.append(np.sum(slope_above_025))
        val = compute(data_values[1], data_values[0])

        fig = plot_level1(
            data_values,
            val,
            "Duration (reg > 0.25) (months)",
            "ENSO duration",
            dt_ls,
        )
        data_to_save.append(data_to_cube(data_values[0], None, metric))
        data_to_save.append(data_to_cube(data_values[1], None, metric))
    elif metric == "15diversity":
        for datas in input_pair:  # obs 0, mod 1
            events = enso_events(datas[var_group[0]])
            results_lon = diversity(datas[var_group[1]], events)
            results_lon["enso"] = results_lon["nino"] + results_lon["nina"]
            data_values.append(iqr(results_lon["enso"]))

        val = compute(data_values[1], data_values[0])
        fig = plot_level1(
            data_values,
            val,
            "IQR of min/max SSTA(°long)",
            "ENSO diversity",
            dt_ls,
        )
        data_to_save.append(data_to_cube(data_values[1], None, metric))
        data_to_save.append(data_to_cube(data_values[0], None, metric))

    return val, fig, data_to_save


def mask_to_years(events):
    """Get years from mask."""
    maskedtime = np.ma.masked_array(
        events.coord("time").points,
        mask=events.data.mask,
    )
    # return years
    return [
        events.coord("time").units.num2date(time).year
        for time in maskedtime.compressed()
    ]


def enso_events(cube):
    """Get enso event years from dataset.

    La Nina events are <-0.75 * standard deviation
    El Nino events are > 0.75 * standard deviation
    """
    std = cube.data.std()
    a_events = mask_to_years(mask_above_threshold(cube.copy(), -0.75 * std))
    o_events = mask_to_years(mask_below_threshold(cube.copy(), 0.75 * std))
    return {"nina": a_events, "nino": o_events}


def diversity(ssta_cube, events_dict):
    """Compute diversity from event years."""
    res_lon = {}
    for enso, events in events_dict.items():
        year_enso = iris.Constraint(
            time=lambda cell, events=events: cell.point.year in events,
        )
        cube = ssta_cube.extract(year_enso)
        if enso == "nina":
            cube = cube * -1
        # iterate through cube, each time get max/min value and return lon
        loc_ls = []
        for yr_slice in cube.slices(["longitude"]):
            indx = np.argmax(yr_slice.data)
            loc_ls.append(cube.coord("longitude").points[indx])

        res_lon[enso] = loc_ls
    return res_lon


def iqr(data):
    """Compute interquartile range."""
    qrt3, qrt1 = np.percentile(data, [75, 25])

    return qrt3 - qrt1


def format_lon(x_val, _):
    """Format longitude in plot axis."""
    if x_val > 180:
        return f"{(360 - x_val):.0f}°W"
    if x_val == 180:
        return f"{x_val:.0f}°"

    return f"{x_val:.0f}°E"


def compute(obs, mod):
    """Compute percentage metric value."""
    return abs((mod - obs) / obs) * 100


def group_obs_models(obs, models, metric, var_preproc, cfg):
    """Group obs to models for metric computation."""
    metricfile = get_diagnostic_filename("matrix", cfg, extension="csv")
    prov_record = get_provenance_record(
        metric,
        list(cfg["input_data"].keys()),
    )
    # obs datasets for each model
    obs_datasets = {
        dataset["variable_group"]: iris.load_cube(dataset["filename"])
        for dataset in obs
    }
    # group models by dataset
    for dataset, attributes in group_metadata(
        models,
        "dataset",
        sort="project",
    ).items():
        logger.info(
            "%s, dataset:%s",
            metric,
            dataset,
        )
        data_labels = [dataset, obs[0]["dataset"]]
        output = compute_enso_metrics(
            [
                obs_datasets,
                {
                    attr["variable_group"]: iris.load_cube(attr["filename"])
                    for attr in attributes
                },
            ],
            data_labels,
            var_preproc,
            metric,
        )
        # save returned cubes
        for i, cube in enumerate(output[2]):
            save_data(f"{data_labels[i]}_{metric}", prov_record, cfg, cube)

        if output[0]:
            with open(metricfile, "a+", encoding="utf-8") as fileo:
                fileo.write(f"{dataset},{metric},{output[0]}\n")

            save_figure(
                f"{dataset}_{metric}",
                prov_record,
                cfg,
                figure=output[1],
                dpi=300,
            )
        # clear value,fig
        output = None


def get_provenance_record(metric, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = {
        "09pattern": (
            "Zonal structure of sea surface temperature anomalies in the "
            + "equatorial Pacific (averaged between 5°S and 5°N)."
        ),
        "10lifecycle": (
            "Temporal evolution of sea surface temperature anomalies in "
            + "the central equatorial Pacific (Niño 3.4 region average), "
            + "illustrating the ENSO-associated variability."
        ),
        "11amplitude": (
            "Standard deviation of sea surface temperature anomalies in "
            + "the central equatorial Pacific (Niño 3.4 region average), "
            + "representing the amplitude of variability."
        ),
        "12seasonality": (
            "Ratio of winter to spring standard deviation of sea surface "
            + "temperature anomalies in the central equatorial Pacific, "
            + "illustrating the seasonal timing of SSTA."
        ),
        "13asymmetry": (
            "Skewness of sea surface temperature anomalies in the central "
            + "equatorial Pacific, illustrating the expected asymmetry where "
            + "positive SSTA values should typically be larger than negative "
            + "SSTA values (usually close to 0). "
        ),
        "14duration": (
            "Duration of the ENSO life cycle where SSTA exceeds 0.25, "
            + "illustrating the 'duration' of the SSTA event."
        ),
        "15diversity": (
            "Width of the zonal location of maximum (minimum) SSTA during "
            + "all El Niño (La Niña) events, illustrating the 'diversity' "
            + "of ENSO events."
        ),
        "values": "List of metric values.",
    }
    record = {
        "caption": caption[metric],
        "statistics": ["anomaly"],
        "domains": ["eq"],
        "plot_types": ["line"],
        "authors": [
            "chun_felicity",
            "beucher_romain",
            "sullivan_arnold",
        ],
        "references": [
            "planton2021",
        ],
        "ancestors": ancestor_files,
    }
    return record


def main(cfg):
    """Run ENSO metrics."""
    input_data = cfg["input_data"].values()

    # iterate through each metric and get variable group, select_metadata
    metrics = {
        "09pattern": ["tos_patdiv1", "tos_pat2"],
        "10lifecycle": ["tos_lifdur1"],
        "11amplitude": ["tos_amp"],
        "12seasonality": ["tos_seas_asym"],
        "13asymmetry": ["tos_seas_asym"],
        "14duration": ["tos_lifdur1"],
        "15diversity": ["tos_patdiv1", "tos_lifdurdiv2"],
    }

    # select twice with project to get obs, iterate through model selection
    for metric, var_preproc in metrics.items():
        logger.info("%s,%s", metric, var_preproc)
        obs, models = [], []
        for var_prep in var_preproc:
            obs += select_metadata(
                input_data,
                variable_group=var_prep,
                project="OBS",
            )
            obs += select_metadata(
                input_data,
                variable_group=var_prep,
                project="OBS6",
            )
            models += select_metadata(
                input_data,
                variable_group=var_prep,
                project="CMIP6",
            )

        group_obs_models(obs, models, metric, var_preproc, cfg)

    # write provenance for csv metrics
    metricfile = get_diagnostic_filename("matrix", cfg, extension="csv")
    prov = get_provenance_record("values", list(cfg["input_data"].keys()))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(metricfile, prov)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
