#!/usr/bin/env  python
"""Calculate and plot relative area of drought events.

Creates timeseries of the spatial extend of all drought events. Different types
of events can be configured for specific ranges of drought indices. A
multimodel mean is calculated by averaging over all datasets event area ratios.
The datasets are required to be preprocessed accordingly to the same shape and
time axis.

The default parameters for this diagnostic are defined in
event_area_timeseries.yml.


Configuration options in recipe
-------------------------------
interval: int, optional (default: 240)
    Number of months per plot, for regular intervals. Set to 0 to create
    one plot for the full period.
intervals: list of tuples, optional (default: None)
    List of tuples with start and end time indices for each plot. If not set,
    `interval` will be used to generate regular ranges.
events: list of dicts
    List of event types with min and max index values, colors and labels.
    See event_area_timeseries.yml for an example.
fig_kwargs: dict, optional
    Additional keyword arguments for the figure creation. This can be set for
    specific plottypes as ``fullperiod.fig_kwargs`` or ``overview.fig_kwargs``.
overview: dict, optional
    Setup for a figure with multiple plots for selected intervals. The first
    interval is expected to be plotted once, all others are plotted for each
    scenario. plot_kwargs and fig_kwargs can be set as for this figure.
    Set ``overview.skip: True`` to not create this figure.
fullperiod: dict, optional
    Setup for a figure with one plot for each pair of scenario and region.
    ``plot_kwargs`` and ``fig_kwargs`` can be set as for this figure.
    Set ``fullperiod.skip: True`` to not create this figure.
regions: list of str, optional
    List of regions (acronyms) for which the area ratios are plotted.
    If not given, global data is used.
combine_regions: bool, optional (default: False)
    If true the list of acronyms is combined into one region,
    rather than plotting them individually.
basename: str, optional (default: "SPEI_{dataset}_{interval}")
    Format string for plot file names. The string will be formatted with the
    current meta data. "group", "interval" and "region" are be available
    in certain cases. For multimodel mean plots "MMM" is used as dataset name.
reuse: bool, optional (default: False)
    If True, the diagnostic will try to load existing data from the output of
    previous runs. This should be only set to true in temporary settings during
    development or adjusting plot parameters. Should not be changed in recipes.
"""

import logging
import os

import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import xarray as xr
import yaml
from esmvalcore import preprocessor as pp
from iris.analysis.cartography import area_weights
from matplotlib import gridspec
from matplotlib.dates import DateFormatter, YearLocator  # MonthLocator

import esmvaltool.diag_scripts.droughtindex.utils as ut
from esmvaltool.diag_scripts.shared import (
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

log = logging.getLogger(__file__)


def load_and_prepare(cfg, fname):
    """apply mask and guess lat/lon bounds."""
    cube = iris.load_cube(fname)
    ut.guess_lat_lon_bounds(cube)
    old_mask = cube.data.mask
    new_mask = get_2d_mask(cube, tile=True)
    diff_mask = np.logical_xor(old_mask, new_mask)
    cube.data.mask = new_mask
    cube.data.data[diff_mask] = np.NaN
    return cube


def get_intervals(cube, interval):
    if interval > 0:
        months = len(cube.coord("time").points)
        steps = int(months / interval)
        intervals = [(i * interval, (i + 1) * interval) for i in range(steps)]
        if months % interval > 0:  # plot last (incomplete) interval
            intervals += [(intervals[-1][-1], None)]
    else:
        intervals = [(0, None)]
    return intervals


def calc_ratio(cube, event, weights):
    """Calculates a timeseries of area ratio for specific index range.

    The fraction of area (or cells) with index values between min and max is
    calculated for each timestep. NaN values are ignored for each event.
    If imin and imax are set to "nan" the ratio of cells with nan/inf values is
    returned.

    Parameters
    ----------
    cube : iris.cube
        3D drought index, with area
    event : dict
        must contain floats or "nan" for keys: "min", "max"
    weights : np.ndarray
        weights array with the same shape as the cube
    nan_area: bool
        return the ratio of nan values ignoring imin and imax.
    """
    imin = event["min"]
    imax = event["max"]
    if imin == "nan" and imax == "nan":
        mask = np.isfinite(cube.data)
    else:
        mask = ~np.logical_and(cube.data >= imin, cube.data < imax)
    event_areas = ma.masked_array(weights, mask=mask)  # still 3D
    return np.sum(event_areas, axis=(1, 2))  # collapse lat/lon


def get_2d_mask(cube, mask_any=False, tile=False):
    """return a 2d (lat/lon) mask where any or all entries are masked.

    Parameters
    ----------
    cube : iris.cube
        3d cube with masked data
    mask_any: bool
        return true for any masked entrie along time dim, instead of all
    tile : bool
        return a 3d mask (matching cube) repeated along the time dim
    """
    if mask_any:
        mask2d = np.any(cube.data.mask, axis=0)
    else:
        mask2d = np.all(cube.data.mask, axis=0)
    if tile:
        mask2d = np.tile(mask2d, (cube.shape[0], 1, 1))
    return mask2d


def plot_area_ratios(cfg, meta, cube):
    """plot area ratio of given event types for a cube of index values

    The area weights are normalized on the masked cube data, resulting in the
    ratio between area with index values in a given range and the area of all
    grid cells, that are not masked (usually not glaciated land).

    Parameters
    ----------
    cube : iris.Cube
        3D index values.
    interval : int
        number of months per figure. negative values disable the split.
        Optional, -1 by default
    weighted : boolean
        calculate area weights based on grid boundaries. Masked cells are
        ignored in normalization.
    """
    # if cfg["weighted"]:
    #     # NOTE: area_weights does not apply cubes mask, so normalize manually
    #     weights = area_weights(cube, normalize=True)
    # else:
    #     weights = np.ones(cube.data.shape) / (cube.shape[1] * cube.shape[2])
    # # mask and normalize weights to sum up to 1 for unmasked (land/region)
    # mask = get_2d_mask(cube, tile=True)
    # weights = ma.masked_array(weights, mask=mask)
    # unmasked_area = np.sum(weights) / cube.shape[0]
    # weights = weights / unmasked_area
    # # calculate areas for each event
    # for event in cfg["events"]:
    #     event["area"] = calc_ratio(cube, event, weights=weights)
    #     mm_key = f'mm_{meta["region"]}' if "region" in meta else "mm"
    #     if mm_key not in event:
    #         event[mm_key] = []
    #     event[mm_key].append(event["area"])
    y = np.vstack([e["area"] for e in cfg["events"]])
    if cfg.get("intervals", None) is None:
        cfg["intervals"] = get_intervals(cube, cfg["interval"])
    if cfg.get("plot_models", True):
        for i in cfg["intervals"]:
            ftemp = f"{cube.name()}_{meta['dataset']}_{i[0]}-{i[1]}"
            if "region" in meta:
                ftemp += f"_{meta['region']}"
            fname = get_plot_filename(ftemp, cfg)
            plot(cfg, fname, i, y)


def plot(cfg, i, y, fname=None, fig=None, ax=None, label=None, full=False):
    """plot area ratios for given interval and events.
    pass either fname to save single plots, or fig and ax to plot one axis into
    existing figure.
    """
    log.debug("TIMES: %s", cfg["times"])
    if i is None:
        i = (0, len(cfg["times"]))
    t = cfg["times"][i[0] : i[1]]
    # dt = t[1] - t[0]
    if not fig:
        fig, ax = plt.subplots(figsize=cfg["figsize"], dpi=300)
    ax.xaxis.set_major_locator(YearLocator(5))
    ax.xaxis.set_major_formatter(DateFormatter("%Y"))
    ax.xaxis.set_minor_locator(YearLocator(1))
    ax.set_ylabel(label)
    plot_kwargs = dict(step="mid", colors=cfg["colors"], labels=cfg["labels"])
    plot_kwargs.update(cfg.get("plot_kwargs", {}))
    ax.stackplot(t, y[:, i[0] : i[1]], **plot_kwargs)
    ax.set_ylim(*cfg["ylim"])
    ax.set(**cfg.get("axes_properties", {}))
    ax.tick_params(direction="in", which="both")
    # ax.set_xlim(t[0]-dt.timedelta(days=20), t[-1])  # show first tick
    ax.set_xlim(t[0], t[-1])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=00, ha="center")
    print(t[0], t[-1])
    # if cfg["interval"] > 0 and not full:
    #     log.info("setting xlim for interval")
    #     ax.set_xlim(t[0] - dt, t[-1])
    #     if len(t) < cfg["interval"]:
    #         ax.set_xlim(t[0], t[0] + cfg["interval"] * dt)
    # else:
    #     log.info("setting xlim for full period %s %s", t[0], t[-1])
    #     ax.set_xlim(t[0], t[-1])
    if cfg.get("strip_plot", False):
        # standalone legend in seperate figure
        fig2, ax2 = plt.subplots(figsize=(5, 1), dpi=300)
        ax2.stackplot(t, y[:, i[0] : i[1]], **plot_kwargs)
        ax2.axis("off")
        lfname = get_plot_filename(f"{cfg['index']}_MMM_legend", cfg)
        fig2.legend(mode="expand", ncol=2, fancybox=False, framealpha=1)
        fig2.savefig(lfname)
    # elif not cfg.get("legend_latest", True) or i[1] is None:
    # add legend last or every figure
    # if not cfg.get("subplots") and not cfg.get("subplots_full"):
    #     fig.legend(
    #         loc="center right",
    #         fancybox=False,
    #         framealpha=1,
    #         labelspacing=0.3,
    #     )  # 'center right'
    if fname:
        fig.savefig(fname)
        plt.close()


def plot_mmm(cfg, region=None, fig=None, axs=None, label=None, meta={}):
    """plot multi model mean of area ratios for given interval and events."""
    mmm_key = f"mmm_{region}" if region else "mmm"
    mm_key = f"mm_{region}" if region else "mm"
    for event in cfg["events"]:
        event[mmm_key] = np.mean(np.array(event[mm_key]), axis=(0))
    y = np.vstack([e[mmm_key] for e in cfg["events"]])
    if cfg.get("intervals", None) is not None:
        for n, i in enumerate(cfg["intervals"]):
            meta["interval"] = f"{i[0]}-{i[1]}"
            basename = cfg["basename"].format(**meta)
            if cfg["subplots"]:
                plot(cfg, i, y, fig=fig, ax=axs[n], label=label)
            # else:
            #     fname = get_plot_filename(basename, cfg)
            #     plot(cfg, i, y, fname=fname)
    if cfg.get("subplots_full", False):
        plot(cfg, i, y, fig=fig, ax=axs, label=label)


def set_defaults(cfg):
    """update cfg with default values from diffmap.yml
    TODO: This could be a shared function reused in other diagnostics.
    """
    config_file = os.path.realpath(__file__)[:-3] + ".yml"
    with open(config_file, encoding="utf-8") as f:
        defaults = yaml.safe_load(f)
    for key, val in defaults.items():
        cfg.setdefault(key, val)
    cfg["colors"] = [e["color"] for e in cfg["events"]]
    cfg["labels"] = [e["label"] for e in cfg["events"]]
    cfg["mirror"] = len(cfg.get("mirror_events", [])) > 0
    cfg.setdefault("combine_regions", False)
    cfg.setdefault("basename", "SPEI_MMM_{interval}")
    cfg["overview"].setdefault("skip", False)
    cfg["fullperiod"].setdefault("skip", False)


def set_timecoords(cfg):
    """Read time coordinate points from reference dataset and store it in cfg.
    This is required to ensure that all datasets have the same time axis.
    TODO: maybe new regrid_time with common calendar could replace this?
    """
    meta = ut.select_single_meta(
        cfg["input_data"].values(),
        short_name=cfg["index"],
        # experiment="historical-ssp585",  # TODO:
        dataset=cfg["reference_dataset"],
        strict=False,
    )
    tc = iris.load_cube(meta["filename"]).coord("time")
    cfg["times"] = iplt._fixup_dates(tc, tc.points)
    print(cfg["times"])


def plot_overview(cfg, data, group="unnamed"):
    """prepare a figure with 1 histogrical and 3 future scenario intervals."""
    fig = plt.figure(**ut.sub_cfg(cfg, "overview", "fig_kwargs"), dpi=300)
    gs = gridspec.GridSpec(3, 2)
    scenario_axs = [
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 1]),
        fig.add_subplot(gs[2, 1]),
    ]
    hist_ax = fig.add_subplot(gs[0, 0], sharey=scenario_axs[0])
    # NOTE: sharey hides second yticklabels, use twin to show it
    twin = scenario_axs[0].twinx()
    twin.tick_params(axis="y", which="both", right=True, direction="in")
    # TODO: can this be moved to ax properties?
    twin.set_yticks(cfg.get("yticks", None))
    leg_ax = fig.add_subplot(gs[1:, 0])
    hist_plotted = False
    for n, ((exp), e_data) in enumerate(data.groupby(["exp"])):
        dat = e_data.squeeze()
        # pick first and last interval:
        if not hist_plotted:
            # first = dat.isel(time=slice(0, cfg["intervals"][0][1]))
            # print(first["event_ratio"])
            plot(
                cfg,
                cfg["intervals"][0],
                dat["event_ratio"].data,
                fig=fig,
                ax=hist_ax,
            )
        # last = dat.isel(time=slice(cfg["intervals"][-1][0], None))
        plot(
            cfg,
            cfg["intervals"][-1],
            dat["event_ratio"].data,
            fig=fig,
            ax=scenario_axs[n],
        )
    for ax in [hist_ax, *scenario_axs]:  # all plots
        ax.set_yticks(cfg.get("yticks", None))
        ax.grid(True, which="major", linestyle="--", linewidth=0.5)
        ax.tick_params(axis="x", which="both", top=True, bottom=True)
        ax.yaxis.tick_right()
    for ax in scenario_axs:  # right plots
        ax.tick_params(axis="y", which="both", left=True, right=True)
    for ax in scenario_axs[:-1]:  # disable xicks
        ax.set_xticklabels([])
    hist_ax.set_ylabel("Historical")
    hist_ax.yaxis.tick_right()
    hist_ax.set_yticklabels([])
    # axs[1][0].tick_params(axis='x', which='both', pad=10)
    leg_ax.axis("off")
    hands, labs = scenario_axs[0].get_legend_handles_labels()
    legend = leg_ax.legend(
        hands[0:8],
        labs[0:8],
        ncol=2,
        loc="center",
        bbox_to_anchor=(0.5, 0.34),
        fancybox=False,
        labelspacing=0,
        framealpha=0,
    )
    for txt in legend.get_texts():
        txt.set_linespacing(1.3)
    leg_ax.axis("off")
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.2)
    fig.savefig(get_plot_filename(f"overview_{cfg['index']}_MMM_{group}", cfg))
    return


def plot_full_periods(cfg, data):
    """prepare a figure with full time series for each scenario/region."""
    # setup figure with 1 row for legend and 1 for each scenario/region pair
    fig = plt.figure(**ut.sub_cfg(cfg, "fullperiod", "fig_kwargs"), dpi=300)
    rows = len(data["exp"]) * len(data["region"])
    gs = gridspec.GridSpec(rows + 1, 1, height_ratios=[0.3] + [1] * (rows))
    leg_ax = fig.add_subplot(gs[0, 0])
    axs = []
    # loop through data slices and plot to axis
    for n, ((exp, reg), dat) in enumerate(data.groupby(["exp", "region"])):
        ax = fig.add_subplot(gs[n + 1, 0])
        dat = dat.squeeze()
        y = dat["event_ratio"].data
        # hardcode full interval for this plot
        i = [0, None]
        fname = get_plot_filename(f"event_area_{exp}_{reg}", cfg)
        plot(cfg, i, y, fname=fname, fig=fig, ax=ax)
        axs.append(ax)
        # if cfg.get("intervals", None) is None:
        #     cfg["intervals"] = get_intervals(cube, cfg["interval"])
        # if cfg.get("plot_models", True):
        #     for i in cfg["intervals"]:
        #         ftemp = f"{cube.name()}_{meta['dataset']}_{i[0]}-{i[1]}"
        #         if "region" in meta:
        #             ftemp += f"_{meta['region']}"
        #         fname = get_plot_filename(ftemp, cfg)
        #         plot(cfg, fname, i, y)
    # for i, (exp, emetas) in enumerate(exp_metas.items()):
    #     exp_axs = [scenario_axs[i]]
    #     process_datasets(cfg, emetas, fig=fig, axs=exp_axs)

    for ax in axs:  # all plots
        ax.set_yticks(cfg.get("yticks", None))
        ax.grid(True, which="both", linestyle="--", linewidth=0.5)
        ax.tick_params(axis="x", which="both", top=True, bottom=True)
        # ax.yaxis.tick_right()
        ax.tick_params(axis="y", which="both", left=True, right=True)
    for ax in axs[:-1]:  # disable xicks
        ax.set_xticklabels([])
    leg_ax.axis("off")
    hands, labs = axs[0].get_legend_handles_labels()
    labs = [lab.replace("\n", " ") for lab in labs]
    ncol = 4
    legend = leg_ax.legend(
        hands[0:8],
        labs[0:8],
        ncol=ncol,
        loc="center",
        bbox_to_anchor=(0.5, 0.34),
        fancybox=False,
        labelspacing=0,
        framealpha=0,
    )
    for txt in legend.get_texts():
        txt.set_linespacing(1.3)
    leg_ax.axis("off")
    fig.tight_layout()
    # fig.subplots_adjust(hspace=0.2)
    fig.subplots_adjust(
        left=0.05, right=0.95, top=0.95, bottom=0.05, hspace=0.2
    )
    fig.savefig(get_plot_filename(f"{cfg['index']}_MMM", cfg))
    return


def plot_each_interval(cfg, exp_metas):
    """create an individual figure for each interval and each scenario."""
    for i, (exp, emetas) in enumerate(exp_metas):
        process_datasets(cfg, emetas, fig=None, axs=None)


def process_datasets(cfg, metas, fig=None, axs=None):
    """load all models and call event area calculation for each."""
    last_meta = None
    for meta in metas:
        fname = meta["filename"]
        if not meta["short_name"].lower() == cfg["index"]:
            log.info("Not matching index (skipped): %s", cfg["index"])
            continue
        cube = load_and_prepare(cfg, fname)
        if cfg.get("global", True):
            plot_area_ratios(cfg, meta, cube)
        if cfg.get("regions", False) and not cfg["combine_regions"]:
            for region in cfg["regions"]:
                log.info("-- region %s", region)
                meta["region"] = region
                extracted = pp.extract_shape(
                    cube, shapefile="ar6", ids={"Acronym": [region]}
                )
                plot_area_ratios(cfg, meta, extracted)
        elif cfg.get("regions", False):
            log.info("-- combined region")
            extracted = pp.extract_shape(
                cube, shapefile="ar6", ids={"Acronym": cfg["regions"]}
            )
            if "region" in meta:
                del meta["region"]
            plot_area_ratios(cfg, meta, extracted)
        last_meta = meta
    # multi model mean
    if cfg.get("plot_mmm", True):
        ylabel = cfg.get("ylabels", {}).get(meta["exp"], meta["exp"])
        if cfg.get("global", True) or cfg["combine_regions"]:
            print("plotting mmm global/combined")
            plot_mmm(cfg, fig=fig, axs=axs, label=ylabel, meta=last_meta)
        else:
            for region in cfg.get("regions", [None]):
                print("plotting mmm each region")
                plot_mmm(
                    cfg,
                    region=region,
                    fig=fig,
                    axs=axs,
                    label=ylabel,
                    meta=last_meta,
                )


def extract_regions(cfg, cube):
    """extract regions and return a list of cubes."""
    extracted = {}
    params = {"shapefile": "ar6", "ids": {"Acronym": cfg["regions"]}}
    if cfg["regions"] and cfg["combine_regions"]:
        log.info("extracting combined region")
        extracted["combined"] = pp.extract_shape(cube, **params)
    elif cfg["regions"]:
        for region in cfg["regions"]:
            log.info("extracting region %s", region)
            params["ids"]["Acronym"] = [region]
            extracted[region] = pp.extract_shape(cube, **params)
    else:
        extracted["global"] = cube
    return extracted


def regional_weights(cfg, cube):
    """calculate area weights normalized to the total unmasked area."""
    # NOTE: area_weights does not apply cubes mask, normalize manually
    if cfg["weighted"]:
        weights = area_weights(cube, normalize=True)
    else:
        weights = np.ones(cube.data.shape) / (cube.shape[1] * cube.shape[2])
    # mask and normalize weights to sum up to 1 for unmasked (land/region)
    mask = get_2d_mask(cube, tile=True)
    weights = ma.masked_array(weights, mask=mask)
    unmasked_area = np.sum(weights) / cube.shape[0]
    weights = weights / unmasked_area
    return weights


def calculate_event_ratios(cfg, metas, output):
    """load data and save calculated event ratio timelines."""
    # data: dataset x exp x region x event
    # data_mmm: exp x region x event
    if cfg["regions"] and cfg["combine_regions"]:
        regions = ["combined"]
    elif cfg["regions"]:
        regions = cfg["regions"]
    else:
        regions = ["global"]
    coords = {
        "dataset": list(group_metadata(metas.values(), "dataset").keys()),
        "region": regions,
        "exp": list(group_metadata(metas.values(), "exp").keys()),
        "event": [e["label"] for e in cfg["events"]],
        "time": cfg["times"],
    }
    dims = list(coords.keys())
    nans = np.full([len(c) for c in coords.values()], np.NaN)
    data = xr.Dataset({"event_ratio": (dims, nans)}, coords=coords)
    for meta in metas.values():
        cube = load_and_prepare(cfg, meta["filename"])
        extracted = extract_regions(cfg, cube)
        loc = {"dataset": meta["dataset"], "exp": meta["exp"]}
        for region, r_cube in extracted.items():
            loc["region"] = region
            log.info("calculating %s", region)
            weights = regional_weights(cfg, r_cube)
            for event in cfg["events"]:
                loc["event"] = event["label"]
                ratios = calc_ratio(r_cube, event, weights=weights)
                data["event_ratio"].loc[loc] = ratios
    fname = get_diagnostic_filename("event_area", cfg)
    output[fname] = {"filename": fname, "plottype": "event_area"}
    data.to_netcdf(fname)
    data_mmm = data.mean("dataset")
    # data_mmm = data_mmm.drop_vars("dataset")
    fname = get_diagnostic_filename("event_area_mmm", cfg)
    output[fname] = {"filename": fname, "plottype": "event_area_mmm"}
    data_mmm.to_netcdf(fname)
    return data, data_mmm


def load_event_ratios(cfg):
    """load calculated event ratios for datasets and MMM from files."""
    names = ["event_area", "event_area_mmm"]
    fnames = [get_diagnostic_filename(f, cfg) for f in names]
    datas = {}
    for fname in fnames:
        try:
            datas[fname] = xr.open_dataset(fname)
        except FileNotFoundError:
            log.error("File not found: %s", fname)
            datas[fname] = None
    return datas[fnames[0]], datas[fnames[1]]


def main(cfg):
    """get common time coordinates, execute the diagnostic code.
    Loop over experiments, than datasets.
    """
    set_defaults(cfg)
    set_timecoords(cfg)
    metas = cfg["input_data"]
    output = {}
    if cfg["reuse"]:
        data, data_mmm = load_event_ratios(cfg)
    else:
        data, data_mmm = calculate_event_ratios(cfg, metas, output)
    if not cfg["overview"]["skip"]:
        for (reg), reg_data in data_mmm.groupby("region"):
            plot_overview(cfg, data_mmm, group=reg)
    if not cfg["fullperiod"]["skip"]:
        plot_full_periods(cfg, data_mmm)
    # if cfg.get("plot_intervals", False):
    #     exp_metas = group_metadata(metas.values(), "exp")
    #     plot_each_interval(cfg, exp_metas)
    # TODO: when data is loaded its not written to metadata rn.
    ut.save_metadata(cfg, output)


if __name__ == "__main__":
    with run_diagnostic() as cfg:
        main(cfg)
