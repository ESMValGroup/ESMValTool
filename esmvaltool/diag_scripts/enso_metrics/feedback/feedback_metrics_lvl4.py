"""diagnostic script to plot ENSO feedback metrics level 4."""

import logging
import os

import iris
import iris.plot as iplt
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


def feedback_nonlin_lvl4(sst_cube, tauu_cube):
    """Calculate ENSO feedback metrics for level 4."""
    tauu_aux = tauu_cube.copy()
    sst_coord = iris.coords.AuxCoord(
        sst_cube.data,
        sst_cube.standard_name,
        sst_cube.long_name,
        sst_cube.var_name,
        sst_cube.units,
    )
    tauu_aux.add_aux_coord(sst_coord, 0)
    logger.info(
        "Add aux coord: %s, %s, %s, %s",
        sst_cube.standard_name,
        sst_cube.long_name,
        sst_cube.var_name,
        sst_cube.units,
    )

    below0 = iris.Constraint(
        coord_values={sst_cube.standard_name: lambda cell: cell < 0},
    )
    above0 = iris.Constraint(
        coord_values={sst_cube.standard_name: lambda cell: cell > 0},
    )
    ssta_neg = mask_above_threshold(sst_cube.copy(), 0)  # x<0
    ssta_pos = mask_below_threshold(sst_cube.copy(), 0)  # x>0
    xbelow0 = tauu_aux.extract(below0)
    xabove0 = tauu_aux.extract(above0)

    logger.info("%s, %s", sst_cube.standard_name, tauu_cube.standard_name)
    logger.info(
        "%s: %s, %s: %s, %s: %s, %s: %s",
        "ssta_neg",
        ssta_neg.shape,
        "ssta_pos",
        ssta_pos.shape,
        "xbelow0",
        xbelow0.shape,
        "xabove0",
        xabove0.shape,
    )

    outreg_cube = annual_structure_reg(xbelow0, ssta_neg)
    posreg_cube = annual_structure_reg(xabove0, ssta_pos)
    all_cube = annual_structure_reg(tauu_cube, sst_cube)

    return all_cube, posreg_cube, outreg_cube


def lin_regress_4(cubea, cubebsst):
    """Perform linear regression of cubes."""
    a_data = cubea.data.reshape(
        cubea.shape[0],
        -1,
    )  # Shape (time, spatial_points)
    if cubea.shape[0] == cubebsst.shape[0]:
        b_data = cubebsst.data.flatten()
    else:
        b_data = cubebsst.data.compressed()  # compress masked

    b_with_intercept = np.vstack([b_data, np.ones_like(b_data)]).T
    logger.info(
        "B_with_intercept shape: %s, A_data shape: %s",
        b_with_intercept.shape,
        a_data.shape,
    )
    coefs, _, _, _ = np.linalg.lstsq(b_with_intercept, a_data, rcond=None)
    return coefs[0]


def annual_structure_reg(nhf_cube, ts_cube):
    """Calculate annual structure regression for cubes."""
    months_reg = []
    for i in range(1, 13):
        # extract for both cubes
        nhf = extract_month(nhf_cube, month=i)
        ts = extract_month(ts_cube, month=i)
        logger.info("nhf: %s, ts: %s, month: %d", nhf.shape, ts.shape, i)
        coefs = lin_regress_4(nhf, ts)
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


def plot_level_4(obs_ds, model_ds, metric_varls, ds_labels):
    """Plot level 4 diagnostics for ENSO feedback metrics."""
    tau_modcube = rolling_window_statistics(
        model_ds[metric_varls[1]],
        coordinate="longitude",
        operator="mean",
        window_length=30,
    )
    tau_obcube = rolling_window_statistics(
        obs_ds[metric_varls[1]],
        coordinate="longitude",
        operator="mean",
        window_length=30,
    )

    # obs datasets can have different time range..
    obs1, obs2 = obs_extract_overlap(tau_obcube, obs_ds[metric_varls[0]])

    xvar = metric_varls[0].split("_")[0]  # ts, tauu, ssh
    yvar = metric_varls[1].split("_")[0]
    rangesx = {
        "sst": (-20, 20, 2),
        "taux": (-0.5, 0.5, 0.05),
        "ssh": (-0.3, 0.3, 0.02),
    }

    modells = feedback_nonlin_lvl4(model_ds[metric_varls[0]], tau_modcube)
    obsls = feedback_nonlin_lvl4(obs2, obs1)
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
                levels=np.arange(*rangesx[xvar]),
                extend="both",
                cmap="RdBu_r",
            )

            if i < 325:  # bottom row
                plt.tick_params("x", labelbottom=False)
            else:
                plt.gca().xaxis.set_major_formatter(
                    plt.FuncFormatter(format_longitude),
                )  # apply to bottom left =>325
                ax2.set_xlabel("longitude")

            if i % 2 == 0:
                plt.tick_params(labelleft=False)

            i += 1

    plt.subplots_adjust(top=0.95)
    cax = plt.axes([0.15, 0.04, 0.7, 0.02])
    tickranges = (
        rangesx[xvar][0],
        rangesx[xvar][1] + rangesx[xvar][2],
        rangesx[xvar][1] * 0.5,
    )
    cbar = figure.colorbar(
        cf1,
        cax=cax,
        orientation="horizontal",
        extend="both",
        ticks=np.arange(*tickranges),
    )
    unitlabels = {"sst": "°C", "taux": "N/m2", "ssh": "cm", 1: "1e-3 ", 0: ""}
    n = 1 if xvar == "taux" or yvar == "taux" else 0
    cbar.set_label(
        f"regression ({unitlabels[n]} {unitlabels[yvar]}/{unitlabels[xvar]})",
    )  # title (1e-3 cm/N/m2) (°C/cm)

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
        d.replace(day=15, hour=0) for d in _unit.num2date(t_coord.points)
    ]
    t_coord.points = _unit.date2num(new_time).astype("float64")


def main(cfg):
    """Run ENSO feedback metrics."""
    input_data = cfg["input_data"].values()

    # iterate through each metric and get variable group, select_metadata, map to function call
    metrics = {
        "SST_TAUX": ["sst_east", "taux_eqp"],
        "TAUX_SSH": ["taux_west", "ssh_eqp"],
        "SSH_SST": ["ssh_eqp_area", "sst_eqp"],
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

        # log
        msg = (
            f"{metric} : observation datasets {len(obs)}, models {len(models)}"
        )
        logger.info(msg)

        # group models by dataset
        model_ds = group_metadata(models, "dataset", sort="project")

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
            prov_record = get_provenance_record(
                f"ENSO metrics {metric} feedback level 4",
                dt_files,
            )

            obs_ds = {
                dataset["variable_group"]: iris.load_cube(
                    dataset["filename"],
                    callback=_fix_time_coord,
                )
                for dataset in obs
            }
            model = {
                attributes["variable_group"]: iris.load_cube(
                    attributes["filename"],
                )
                for attributes in mod_ds
            }
            ds_labels = [
                f"{obs[0]['dataset']}_{obs[1]['dataset']}",
                mod_ds[0]["dataset"],
            ]

            # plot level 4
            fig = plot_level_4(obs_ds, model, var_preproc, ds_labels)
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
