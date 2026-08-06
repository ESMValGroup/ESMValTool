"""Utility functions for plotting."""

import logging
import warnings
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import iris
import iris.analysis.cartography
import iris.plot as iplt
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import seaborn as sns
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from matplotlib.colors import ListedColormap
from wt_utils import get_driver, get_provenance_record, log_provenance

iris.FUTURE.datum_support = True

logger = logging.getLogger(Path(__file__).name)

# Ignoring a warning that is produced when selecting timesteps of a weathertype
warnings.filterwarnings("ignore", ".*Collapsing a non-contiguous coordinate*")


def generate_grayscale_hex_values(x):
    """Generate grayscale values for plotting seasonal occurrences.

    Parameters
    ----------
    x
        Amount of weathertypes.

    Returns
    -------
    list
        List of grayscale hex values
    """
    grayscale_values = np.linspace(0, 1, x)

    return [
        f"#{int(value * 255):02x}{int(value * 255):02x}{int(value * 255):02x}"
        for value in grayscale_values
    ]


def plot_seasonal_occurrence(
    cfg: dict,
    wt_cubes: iris.cube.Cube,
    data_info: dict,
    cube_path: str,
):
    """Plot seasonal occurrences of weathertypes.

    Parameters
    ----------
    cfg
        Configuration dictionary from recipe
    wt_cubes
        List of cubes of lwt, slwt_ERA5 and slwt_EOBS
    data_info
        Dictionary with info to dataset
    cube_path
        Paths to weathertype cubes (ancestors)
    """
    driver = get_driver(data_info)

    output_path = f"{cfg['plot_dir']}/seasonal_occurrence"

    if not Path(output_path).exists():
        Path(output_path).mkdir(parents=True, exist_ok=True)

    month_list = [
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ]

    relative_occurrences = {}
    # {wt_string: {month: {wt1: rel_occurrence, wt2: rel_occurrence, ....}}}
    # first do absolute occurrence, then relative occurrence
    for cube in wt_cubes:
        month_dict = {}
        for month in range(1, 13):
            array = cube.extract(
                iris.Constraint(time=iris.time.PartialDateTime(month=month)),
            ).data
            unique, counts = np.unique(array, return_counts=True)
            month_dict[month] = dict(
                zip(unique, counts / sum(counts), strict=False),
            )
        relative_occurrences[cube.long_name] = month_dict

    for wt_string in relative_occurrences:
        wt_numbers = max(
            len(value)
            for value in relative_occurrences.get(wt_string).values()
        )
        wt_stack = np.zeros((wt_numbers, 12))
        for month, month_value in relative_occurrences.get(wt_string).items():
            for wt in month_value:
                wt_stack[np.int8(wt - 1), month - 1] = month_value.get(wt)

        y = np.vstack(list(wt_stack))

        # plot
        _, ax_ = plt.subplots(figsize=(10, 10))

        ax_.set_title(
            f"Seasonal occurrence of {wt_string}, \
            {data_info.get('dataset')}, {data_info.get('timerange')}",
        )

        ax_.stackplot(
            month_list,
            y,
            colors=generate_grayscale_hex_values(wt_numbers),
        )

        ax_.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.05),
            fancybox=True,
            shadow=True,
            ncol=9,
            labels=tuple(f"WT {i + 1}" for i in range(wt_numbers)),
        )

        ax_.set(
            xlim=(0, 11),
            xticks=np.arange(0, 12),
            ylim=(0, 1),
            yticks=np.arange(0, 1.1, 0.1),
        )
        ax_.set_xlabel("Month")
        ax_.set_ylabel("Cumulative Relative occurrence")

        plt.savefig(
            f"{output_path}/{driver}{data_info.get('dataset')}_{data_info.get('ensemble', '')}_"
            f"{wt_string}_seasonal_occurrence_{data_info.get('timerange')}.png",
        )
        plt.savefig(
            f"{output_path}/{driver}{data_info.get('dataset')}_{data_info.get('ensemble', '')}_"
            f"{wt_string}_seasonal_occurrence_{data_info.get('timerange')}.pdf",
        )
        plt.close()

        # ancestors here are just the wt_cubes
        ancestors = [
            f"{cube_path}",
        ]
        provenance_record = get_provenance_record(
            f"Seasonal occurrences for {wt_string}, ",
            ancestors,
            ["wt occurrences"],
            ["seas"],
        )

        log_provenance(
            f"{output_path}/{driver}{data_info.get('dataset')}_{data_info.get('ensemble', '')}_"
            f"{wt_string}_seasonal_occurrence_{data_info.get('timerange')}.png",
            cfg,
            provenance_record,
        )


def set_gridlines(ax: plt.Axes):
    """Set gridlines for plotting maps.

    Parameters
    ----------
    ax
        Axes object to draw gridlines on.

    Returns
    -------
    gl : (ax.gridlines)
        Gridlines object
    """
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.5,
        color="gray",
        alpha=0.5,
        linestyle="--",
    )
    gl.left_labels = True
    gl.bottom_labels = True
    gl.top_labels = False
    gl.right_labels = False
    gl.xlines = True
    gl.ylocator = mticker.FixedLocator(np.arange(20, 70, 5))
    gl.xlocator = mticker.FixedLocator([-10, -5, 0, 5, 10, 15])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {"size": 8, "color": "black"}
    gl.ylabel_style = {"color": "black", "size": 8}

    return gl


def prepare_plot_title(
    data_info: dict,
    wt: int,
    var_name: str,
    mode: str,
) -> str:
    """Return formatted plot title.

    Parameters
    ----------
    data_info
        dict containing info on dataset
    wt
        weathertype
    var_name
        name of variable
    mode
        statistic used

    Returns
    -------
    str
        title of plot
    """
    ensemble = data_info.get("ensemble", "")
    timerange = data_info.get("timerange")
    dataset = data_info.get("dataset")
    if var_name == "pr" and dataset == "ERA5":
        title = f"{dataset} {ensemble}, total {var_name} {mode}\n{timerange}, wt: {wt}"
    elif var_name == "pr":
        title = f"{dataset} {ensemble}, {var_name} flux {mode}\n{timerange}, wt: {wt}"
    elif var_name == "tas":
        title = f"{dataset} {ensemble}, 1000 hPa {var_name} {mode}\n{timerange}, wt: {wt}"
    else:  # psl or others
        title = (
            f"{dataset} {ensemble}, {var_name} {mode}\n{timerange}, wt: {wt}"
        )
    return title


def get_unit(var_name: str, dataset: str) -> str:
    """Get unit of variables.

    Parameters
    ----------
    var_name
        name of variable
    dataset
        name of dataset

    Returns
    -------
    str
        unit of variable
    """
    if var_name == "psl":
        return "[hPa]"
    if var_name == "pr" and dataset == "ERA5":
        return "[m]"
    if var_name == "pr":
        return "[kg m-2 s-1]"
    if var_name == "tas":
        return "[K]"
    return ""


def plot_maps(
    wt: np.array,
    cfg: dict,
    cube: iris.cube.Cube,
    data_info: dict,
    mode: str,
):
    """Plot maps of means, std and anomalies.

    Parameters
    ----------
    wt
        WT array
    cfg
        Configuration dictionary from recipe
    cube
        Data to be plotted
    data_info
        Dictionary with info to dataset
    mode
        Statistics that is used
    """
    var_name = data_info.get("var")

    logger.info(
        "Plotting %s %s %s for %s %s",
        data_info.get("dataset"),
        var_name,
        mode,
        data_info.get("wt_string"),
        wt,
    )

    ax = plt.axes(projection=ccrs.PlateCarree())

    cmap = get_colormap(var_name if var_name != "pr" else "prcp")

    title = prepare_plot_title(data_info, wt, var_name, mode)
    plt.title(title)
    unit = get_unit(var_name, data_info.get("dataset"))

    im = iplt.contourf(cube if var_name != "psl" else cube / 100, cmap=cmap)
    cb = plt.colorbar(im)
    cb.ax.tick_params(labelsize=8)
    cb.set_label(label=f"{var_name} {mode} {unit}")

    set_gridlines(ax)

    ax.set_extent([-15, 20, 27.5, 62.5])

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=":")

    if not Path(f"{cfg.get('plot_dir')}/{mode}").exists():
        Path(f"{cfg.get('plot_dir')}/{mode}").mkdir(
            parents=True,
            exist_ok=True,
        )

    plt.savefig(
        f"{cfg.get('plot_dir')}/{mode}/{data_info.get('wt_string')}_{wt}{get_driver(data_info)}_{data_info.get('dataset')}_{data_info.get('ensemble', '')}"
        f"_{var_name}_{mode}_{data_info.get('timerange')}.png",
    )
    plt.savefig(
        f"{cfg.get('plot_dir')}/{mode}/{data_info.get('wt_string')}_{wt}{get_driver(data_info)}_{data_info.get('dataset')}_{data_info.get('ensemble', '')}_"
        f"{var_name}_{mode}_{data_info.get('timerange')}.pdf",
    )
    plt.close()

    # log provenance
    ancestors = [
        f"{data_info.get('preproc_path')}",
        f"{cfg.get('work_dir')}/ERA5.nc",
    ]
    provenance_record = get_provenance_record(
        f"{var_name} {mode} for {data_info.get('wt_string')}",
        ancestors,
        [var_name],
        ["map"],
        [mode],
    )

    local_path = f"{cfg.get('plot_dir')}/{mode}"

    log_provenance(
        f"{local_path}/{data_info.get('wt_string')}_{wt}{get_driver(data_info)}_"
        f"{data_info.get('dataset')}_{data_info.get('ensemble', '')}"
        f"_{var_name}_{mode}_{data_info.get('timerange')}.png",
        cfg,
        provenance_record,
    )


def plot_corr_rmse_heatmaps(
    cfg: dict,
    pattern_correlation_matrix: np.array,
    rmse_matrix: np.array,
    dataset: str,
    timerange: str,
):
    """Plot correlation and rmse heatmaps.

    Parameters
    ----------
    cfg
        Configuration dictionary from recipe
    pattern_correlation_matrix
        Pattern correlation matrix
    rmse_matrix
        RMSE matrix
    dataset
        Name of dataset
    timerange
        Time range for the calculation
    """
    output_path = f"{cfg.get('plot_dir')}/heatmaps"

    if not Path(f"{output_path}").exists():
        Path(f"{output_path}").mkdir(parents=True, exist_ok=True)

    labels = np.arange(1, 28)

    mask = np.zeros_like(pattern_correlation_matrix)
    mask[np.triu_indices_from(mask)] = True
    with sns.axes_style("white"):
        plt.figure(figsize=(10, 10))
        plt.title(f"Correlation Matrix , {dataset}, {timerange}")
        levels = np.linspace(
            np.min(pattern_correlation_matrix),
            np.max(pattern_correlation_matrix),
            9,
        )
        ax = sns.heatmap(
            pattern_correlation_matrix,
            mask=mask,
            square=True,
            annot=True,
            annot_kws={"size": 6},
            cmap="seismic",
            xticklabels=labels,
            yticklabels=labels,
            cbar_kws={"ticks": levels, "shrink": 0.8, "format": "%.2f"},
        )
        ax.set_xlabel("lwt", fontsize=8)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        ax.set_ylabel("lwt", fontsize=8)
        plt.tight_layout()
        plt.savefig(
            f"{output_path}/correlation_matrix_{dataset}_{timerange}.png",
        )
        plt.savefig(
            f"{output_path}/correlation_matrix_{dataset}_{timerange}.pdf",
        )
        plt.close()

    mask = np.zeros_like(rmse_matrix)
    mask[np.triu_indices_from(mask)] = True
    with sns.axes_style("white"):
        plt.figure(figsize=(10, 10))
        plt.title(f"RMSE Matrix, {dataset}, {timerange}")
        levels = np.linspace(np.min(rmse_matrix), np.max(rmse_matrix), 9)
        ax = sns.heatmap(
            rmse_matrix,
            mask=mask,
            square=True,
            annot=True,
            annot_kws={"size": 5},
            cmap="seismic",
            xticklabels=labels,
            yticklabels=labels,
            cbar_kws={"ticks": levels, "shrink": 0.8, "format": "%.2f"},
        )
        ax.set_xlabel("lwt", fontsize=8)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        ax.set_ylabel("lwt", fontsize=8)
        plt.tight_layout()
        plt.savefig(f"{output_path}/rmse_matrix_{dataset}_{timerange}.png")
        plt.savefig(f"{output_path}/rmse_matrix_{dataset}_{timerange}.pdf")
        plt.close()
    # log provenance
    # rmse matrix
    ancestors = [
        f"{cfg.get('work_dir')}/ERA5.nc",
    ]
    provenance_record = get_provenance_record(
        "rmse matrix",
        ancestors,
        ["rmse"],
        ["other"],
    )

    local_path = f"{cfg.get('plot_dir')}/heatmaps"

    log_provenance(
        f"{local_path}/rmse_matrix_{dataset}_{timerange}.png",
        cfg,
        provenance_record,
    )
    # correlation matrix
    provenance_record = get_provenance_record(
        "correlation matrix",
        ancestors,
        ["correlation"],
        ["other"],
    )

    log_provenance(
        f"{local_path}/correlation_matrix_{dataset}_{timerange}.png",
        cfg,
        provenance_record,
    )


def get_colormap(colormap_string: str) -> ListedColormap:
    """Get colormaps for plottings.

    Parameters
    ----------
    colormap_string
        string to identify colormap

    Returns
    -------
    ListedColormap
        Colormap for the specified variable
    """
    misc_seq_2_disc = [
        (230 / 255, 240 / 255, 240 / 255),
        (182 / 255, 217 / 255, 228 / 255),
        (142 / 255, 192 / 255, 226 / 255),
        (118 / 255, 163 / 255, 228 / 255),
        (116 / 255, 130 / 255, 222 / 255),
        (121 / 255, 97 / 255, 199 / 255),
        (118 / 255, 66 / 255, 164 / 255),
        (107 / 255, 40 / 255, 121 / 255),
        (86 / 255, 22 / 255, 75 / 255),
        (54 / 255, 14 / 255, 36 / 255),
    ]

    temp_seq_disc = [
        (254 / 255, 254 / 255, 203 / 255),
        (251 / 255, 235 / 255, 153 / 255),
        (244 / 255, 204 / 255, 104 / 255),
        (235 / 255, 167 / 255, 84 / 255),
        (228 / 255, 134 / 255, 80 / 255),
        (209 / 255, 98 / 255, 76 / 255),
        (164 / 255, 70 / 255, 66 / 255),
        (114 / 255, 55 / 255, 46 / 255),
        (66 / 255, 40 / 255, 24 / 255),
        (25 / 255, 25 / 255, 0 / 255),
    ]

    prec_seq_disc = [
        (255 / 255, 255 / 255, 229 / 255),
        (217 / 255, 235 / 255, 213 / 255),
        (180 / 255, 216 / 255, 197 / 255),
        (142 / 255, 197 / 255, 181 / 255),
        (105 / 255, 177 / 255, 165 / 255),
        (67 / 255, 158 / 255, 149 / 255),
        (44 / 255, 135 / 255, 127 / 255),
        (29 / 255, 110 / 255, 100 / 255),
        (14 / 255, 85 / 255, 74 / 255),
        (0 / 255, 60 / 255, 48 / 255),
    ]

    if colormap_string == "psl":
        return ListedColormap(misc_seq_2_disc)
    if colormap_string == "prcp":
        return ListedColormap(prec_seq_disc)
    if colormap_string == "temp":
        return ListedColormap(temp_seq_disc)

    return None
