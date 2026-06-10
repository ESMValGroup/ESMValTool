"""Diagnostic script to plot ENSO sst nhf feedback metrics level 3 and 4."""

import logging
import os

import iris
import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import (
    climate_statistics,
    extract_month,
    mask_above_threshold,
    mask_below_threshold,
    rolling_window_statistics,
)

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_figure,
    select_metadata,
)

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def annual_structure_reg(nhf_cube, ts_cube):
    """Annual structure regression for net heat flux and temperature."""
    months_reg = []
    for i in range(1, 13):
        # extract for both cubes
        nhf = extract_month(nhf_cube, month=i)
        ts = extract_month(ts_cube, month=i)
        coefs = lin_regress_matrix(nhf, ts, level4=True)
        # collect array for months
        months_reg.append(coefs)

    # get month_number coordinate from an existing cube
    m_cube = climate_statistics(ts_cube, operator="max", period="month")

    # create cube, x = months, y = lon
    result_cube = iris.cube.Cube(
        np.array(months_reg),
        long_name="regression",
        dim_coords_and_dims=[
            (nhf.coord("longitude"), 1),
            (m_cube.coord("month_number"), 0),
        ],
    )
    return result_cube


def lin_regress_matrix(cubea, cubebsst, level4=False):
    """Perform linear regression between two Cubes."""
    a_data = cubea.data.reshape(
        cubea.shape[0],
        -1,
    )  # Shape (time, spatial_points)
    if cubea.shape[0] == cubebsst.shape[0]:
        b_data = cubebsst.data.flatten()  # or all
    else:
        b_data = cubebsst.data.compressed()  # masked threshold cube (time,)

    # Add intercept term by stacking a column of ones with cubeB
    b_with_intercept = np.vstack([b_data, np.ones_like(b_data)]).T

    logger.info(
        "least squares data shapes %s, %s",
        b_with_intercept.shape,
        a_data.shape,
    )
    # Solve the linear equations using least squares method
    coefs, _, _, _ = np.linalg.lstsq(b_with_intercept, a_data, rcond=None)

    if level4:
        return coefs[0]
    # Create a new Iris Cube for the regression results
    result_cube = iris.cube.Cube(
        coefs[0],
        long_name="regression A",
        dim_coords_and_dims=[(cubea.coord("longitude"), 0)],
    )

    return result_cube


def feedback_nonlin(sst_cube, tauu_cube, level4=False):
    """Calculate non-linear feedback metrics."""
    tauu_aux = tauu_cube.copy()
    sst_coord = iris.coords.AuxCoord(
        sst_cube.data,
        sst_cube.standard_name,
        sst_cube.long_name,
        sst_cube.var_name,
        sst_cube.units,
    )
    tauu_aux.add_aux_coord(sst_coord, 0)
    logger.info("non linear shapes %s, %s", sst_cube.shape, tauu_cube.shape)
    logger.info(tauu_aux.summary())
    below0 = iris.Constraint(
        coord_values={sst_cube.standard_name: lambda cell: cell < 0},
    )
    above0 = iris.Constraint(
        coord_values={sst_cube.standard_name: lambda cell: cell > 0},
    )
    ssta_neg = mask_above_threshold(sst_cube.copy(), 0)  # x<0
    ssta_pos = mask_below_threshold(sst_cube.copy(), 0)  # x=>0
    xbelow0 = tauu_aux.extract(below0)
    xabove0 = tauu_aux.extract(above0)

    if level4:
        outreg_cube = annual_structure_reg(xbelow0, ssta_neg)
        posreg_cube = annual_structure_reg(xabove0, ssta_pos)
        all_cube = annual_structure_reg(tauu_cube, sst_cube)
        return all_cube, outreg_cube, posreg_cube

    outreg_cube = lin_regress_matrix(xbelow0, ssta_neg)
    posreg_cube = lin_regress_matrix(xabove0, ssta_pos)
    return outreg_cube, posreg_cube, None


def obs_extract_overlap(obs_1, obs_2):
    """Extract overlapping time range from two observation datasets."""
    start_1 = obs_1.coord("time").cell(0).point
    end_1 = obs_1.coord("time").cell(-1).point
    start_2 = obs_2.coord("time").cell(0).point
    end_2 = obs_2.coord("time").cell(-1).point

    start_overlap = max(start_1, start_2)
    end_overlap = min(end_1, end_2)
    logger.info(
        "%s, %s obs time overlap: %s to %s",
        obs_1.standard_name,
        obs_2.standard_name,
        start_overlap,
        end_overlap,
    )
    obs1 = obs_1.extract(
        iris.Constraint(
            time=lambda t: start_overlap <= t.point <= end_overlap,
        ),
    )
    obs2 = obs_2.extract(
        iris.Constraint(
            time=lambda t: start_overlap <= t.point <= end_overlap,
        ),
    )

    return obs1, obs2


def format_longitude(x, _pos):
    """Format longitude values for plotting."""
    if x > 180:
        return f"{int(360 - x)}°W"
    if x == 180:
        return f"{int(x)}°"
    return f"{int(x)}°E"


def plot_level3(
    obs_ds,
    model_ds,
    metric_varls,
    ds_labels,
    title,
):
    """Plot level 3 diagnostics for ENSO feedback metrics."""
    figure = plt.figure(figsize=(10, 6), dpi=300)
    # plot whole regression
    cb = lin_regress_matrix(
        model_ds[metric_varls[1]],
        model_ds[metric_varls[0]],
    )
    qplt.plot(cb, color="black", linestyle="solid", label=ds_labels[1])

    # obs datasets can have different time range..
    obs1, obs2 = obs_extract_overlap(
        obs_ds[metric_varls[2]],
        obs_ds[metric_varls[0]],
    )

    cb2 = lin_regress_matrix(obs1, obs2)
    qplt.plot(cb2, color="black", linestyle="--", label=ds_labels[0])
    # process model data split
    neg, pos, _ = feedback_nonlin(
        model_ds[metric_varls[0]],
        model_ds[metric_varls[1]],
    )

    xvar = metric_varls[0].split("_")[0]  # ts, tauu, ssh
    yvar = metric_varls[2].split("_")[0]

    qplt.plot(neg, color="blue", linestyle="solid", label=f"{xvar.upper()}A<0")
    qplt.plot(pos, color="red", linestyle="solid", label=f"{xvar.upper()}A>0")
    # process obs data split
    neg, pos, _ = feedback_nonlin(obs2, obs1)
    qplt.plot(neg, color="blue", linestyle="--")
    qplt.plot(pos, color="red", linestyle="--")

    plt.xlim(170, 250)
    plt.xlabel("longitude")
    plt.ylabel(f"reg({xvar.upper()}A, {yvar.upper()}A)")
    plt.grid(linestyle="--")
    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(format_longitude))
    plt.legend()
    plt.title(title)
    return figure


def plot_level_4(obs_ds, model_ds, metric_varls, ds_labels):
    """Plot level 4 diagnostics for ENSO feedback metrics."""
    obs1, obs2 = obs_extract_overlap(
        obs_ds[metric_varls[2]],
        obs_ds[metric_varls[0]],
    )

    xvar = metric_varls[0].split("_")[0]  # ts, tauu, ssh
    yvar = metric_varls[2].split("_")[0]

    modells = feedback_nonlin(
        model_ds[metric_varls[0]],
        model_ds[metric_varls[1]],
        level4=True,
    )
    obsls = feedback_nonlin(obs2, obs1, level4=True)
    sub = ["", ">0", "<0"]

    figure = plt.figure(figsize=(8, 12), dpi=300)
    i = 321
    ax1 = plt.subplot(i)

    for isub, cbls in enumerate(zip(obsls, modells, strict=False)):
        for cb in cbls:
            if i == 321:
                ax1.set_yticks(range(1, 13, 4))
                ax1.set_yticklabels(["Jan", "May", "Sep"])
                ax1.set_ylabel("Months")
                ax1.set_title(
                    f"ref: {ds_labels[0]} \nreg({xvar.upper()}A, {yvar.upper()}A)",
                    loc="left",
                )  # label
            else:
                ax2 = plt.subplot(i, sharex=ax1, sharey=ax1)
                ax2.set_title(
                    f"reg({xvar.upper()}A{sub[isub]}, {yvar.upper()}A)",
                    loc="left",
                )
                if i == 322:
                    ax2.set_title(
                        f"{ds_labels[1]} \nreg({xvar.upper()}A, {yvar.upper()}A)",
                        loc="left",
                    )

            # contour plt data
            cf1 = iplt.contourf(
                cb,
                coords=["longitude", "month_number"],
                levels=np.arange(-50, 50, 5),
                extend="both",
                cmap="RdBu_r",
            )

            if i < 325:
                plt.tick_params("x", labelbottom=False)
            else:  # apply to bottom left =>325
                plt.gca().xaxis.set_major_formatter(
                    plt.FuncFormatter(format_longitude),
                )
                ax2.set_xlabel("longitude")

            if i % 2 == 0:
                plt.tick_params(labelleft=False)

            i += 1

    plt.subplots_adjust(top=0.95)
    cax = plt.axes([0.15, 0.04, 0.7, 0.02])
    cbar = figure.colorbar(
        cf1,
        cax=cax,
        orientation="horizontal",
        extend="both",
        ticks=np.arange(-50, 51, 25),
    )
    cbar.set_label("regression ((W/m2/°C)")

    return figure


def get_provenance_record(caption, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        "caption": caption,
        "statistics": ["anomaly"],
        "domains": ["eq"],
        "plot_types": ["line"],
        "authors": [
            "chun_felicity",
            "beucher_romain",
            "sullivan_arnold",
        ],
        "references": [
            "access-nri",
        ],
        "ancestors": ancestor_files,
    }
    return record


def _fix_time_coord(cube, _, _filename):
    """Set time points to same day and hour for monthly data."""
    t_coord = cube.coord("time")
    _unit = t_coord.units
    new_time = [
        d.replace(day=15, hour=0, minute=0, second=0, microsecond=0)
        for d in _unit.num2date(t_coord.points)
    ]
    t_coord.points = _unit.date2num(new_time).astype("float64")


def main(cfg):
    """Run ENSO feedback metrics."""
    input_data = cfg["input_data"].values()

    sst_nhf = ["sst_eqp", "nhf_eqp_mod", "nhf_eqp_obs"]
    metric = "sst_nhf"
    obs, models = [], []

    obs += select_metadata(input_data, project="OBS")
    obs += select_metadata(input_data, project="OBS6")
    models += select_metadata(input_data, project="CMIP6")

    # log
    msg = f"{metric} : observation datasets {len(obs)}, models {len(models)}"
    logger.info(msg)

    # group models by dataset
    model_ds = group_metadata(models, "dataset", sort="project")

    obs_ds = {
        dataset["variable_group"]: iris.load_cube(
            dataset["filename"],
            callback=_fix_time_coord,
        )
        for dataset in obs
    }
    obs_ds["nhf_eqp_obs"] = rolling_window_statistics(
        obs_ds["nhf_eqp_obs"],  # check index
        coordinate="longitude",
        operator="mean",
        window_length=30,
    )
    # dataset name
    for dataset, mod_ds in model_ds.items():
        logger.info(
            "%s, preprocessed cubes:%d, dataset:%s",
            metric,
            len(model_ds),
            dataset,
        )
        dt_files = [ds["filename"] for ds in obs] + [
            ds["filename"] for ds in mod_ds
        ]

        model = {
            attributes["variable_group"]: iris.load_cube(
                attributes["filename"],
            )
            for attributes in mod_ds
        }
        model["nhf_eqp_mod"] = -model["nhf_eqp_mod"]  # make negative
        model["nhf_eqp_mod"] = rolling_window_statistics(
            model["nhf_eqp_mod"],
            coordinate="longitude",
            operator="mean",
            window_length=30,
        )

        ds_labels = [
            f"{obs[0]['dataset']}_{obs[1]['dataset']}",
            mod_ds[0]["dataset"],
        ]
        title = "net heat flux feedback"
        # plot level 3
        prov_record = get_provenance_record(
            f"ENSO metrics {metric} feedback level 3",
            dt_files,
        )
        fig = plot_level3(obs_ds, model, sst_nhf, ds_labels, title)
        save_figure(
            f"{dataset}_{metric}_lvl3",
            prov_record,
            cfg,
            figure=fig,
            dpi=300,
        )

        prov_record = get_provenance_record(
            f"ENSO metrics {metric} feedback level 4",
            dt_files,
        )
        fig = plot_level_4(obs_ds, model, sst_nhf, ds_labels)
        save_figure(
            f"{dataset}_{metric}_lvl4",
            prov_record,
            cfg,
            figure=fig,
            dpi=300,
        )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
