"""Utility functions for plotting."""

# operating system manipulations (e.g. path constructions)
import logging
import os
import warnings

# plotting imports
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# to manipulate iris cubes
import iris
import iris.analysis.cartography
import iris.plot as iplt
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# general imports
import numpy as np
import seaborn as sns

# local imports
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from matplotlib.colors import ListedColormap

iris.FUTURE.datum_support = True

logger = logging.getLogger(os.path.basename(__file__))

# Ignoring a warning that is produced when selecting timesteps of a weathertype
warnings.filterwarnings("ignore", ".*Collapsing a non-contiguous coordinate*")


def generate_grayscale_hex_values(x):
    """Generate grayscale values for plotting seasonal occurrence.

    Args:
    ----
        x : int
            number of weathertypes

    Returns
    -------
        np.list
            array with grayscale values as hex
    """
    grayscale_values = np.linspace(0, 1, x)
    grayscale_hex = [
        f"#{int(value * 255):02x}{int(value * 255):02x}{int(value * 255):02x}"
        for value in grayscale_values
    ]

    return grayscale_hex


def plot_seasonal_occurrence(
    cfg: dict, wt_cubes: iris.cube.Cube, data_info: dict
):
    """Plot relative monthly/seasonal occurrence of weathertypes.

    Args:
    ----
        cfg : dict
            Configuration dictionary from recipe
        wt_cubes : iris.cube.Cube
            Cube with weathertypes
        data_info : dict
            Dictionary with relevant info to dataset
    """
    dataset_name = data_info.get("dataset")
    timerange = data_info.get("timerange")
    ensemble = data_info.get("ensemble", "")
    driver = data_info.get("driver", "")

    output_path = f"{cfg['plot_dir']}/seasonal_occurrence"

    if not os.path.exists(f"{output_path}"):
        os.makedirs(f"{output_path}")

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
            month_constraint = iris.Constraint(
                time=iris.time.PartialDateTime(month=month)
            )
            array = cube.extract(month_constraint).data
            unique, counts = np.unique(array, return_counts=True)
            count_dict = dict(zip(unique, counts / sum(counts), strict=False))
            month_dict[month] = count_dict
        relative_occurrences[cube.long_name] = month_dict

    x = month_list

    for wt_string in relative_occurrences:
        wt_numbers = max(
            len(value)
            for value in relative_occurrences.get(wt_string).values()
        )
        colors = generate_grayscale_hex_values(wt_numbers)
        wt_stack = np.zeros((wt_numbers, 12))
        for month, month_value in relative_occurrences.get(wt_string).items():
            print(month_value)
            for wt in month_value.keys():
                print(month_value.get(wt))
                wt_stack[np.int8(wt - 1), month - 1] = month_value.get(wt)

        y = np.vstack(list(wt_stack))

        # plot
        _, ax_ = plt.subplots(figsize=(10, 10))

        ax_.set_title(
            f"Seasonal occurence of {wt_string}, \
            {dataset_name}, {timerange}"
        )

        ax_.stackplot(x, y, colors=colors)

        ax_.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.05),
            fancybox=True,
            shadow=True,
            ncol=9,
            labels=tuple(f"WT {i + 1}" for i in range(0, wt_numbers)),
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
            f"{output_path}/{driver}{dataset_name}_{ensemble}_"
            f"{wt_string}_rel_occurrence_{timerange}.png"
        )
        plt.savefig(
            f"{output_path}/{driver}{dataset_name}_{ensemble}_"
            f"{wt_string}_rel_occurrence_{timerange}.pdf"
        )
        plt.close()


def plot_maps(
    wt: np.array, cfg: dict, cube: iris.cube.Cube, data_info: dict, mode: str
):
    """Plot maps.

    Args:
    ----
        wt : np.array
            weathertype number
        cfg : dict
            Configuration dicitonary from recipe
        cube : iris.cube.Cube
            Data to be plotted
        data_info : dict
            Dictionary with info to dataset
        mode : str
            Statistics that is used
    """
    dataset = data_info.get("dataset")
    var_name = data_info.get("var")
    wt_string = data_info.get("wt_string")
    ensemble = data_info.get("ensemble", "")
    timerange = data_info.get("timerange")
    driver = data_info.get("driver", "")
    if driver != "":
        driver = f"_{driver}"

    logger.info(
        "Plotting %s %s %s for %s %s", dataset, var_name, mode, wt_string, wt
    )

    local_path = f"{cfg.get('plot_dir')}/{mode}"

    if not os.path.exists(f"{local_path}"):
        os.makedirs(f"{local_path}")

    ax = plt.axes(projection=ccrs.PlateCarree())

    if var_name == "psl":
        psl_cmap = get_colormap("psl")
        plt.title(
            f"{dataset} {ensemble}, {var_name} {mode}\n"
            + f"{timerange}, wt: {wt}"
        )
        unit = "[hPa]"
        im = iplt.contourf(cube / 100, cmap=psl_cmap)
        cb = plt.colorbar(im)
        cb.ax.tick_params(labelsize=8)
        cb.set_label(label=f"{var_name} {mode} {unit}")
    elif var_name == "pr":
        prcp_cmap = get_colormap("prcp")
        if dataset == "ERA5":
            unit = "[m]"
            plt.title(
                f"{dataset} {ensemble}, total {var_name} {mode}\n"
                + f"{timerange}, wt: {wt}"
            )
        else:
            unit = "[kg m-2 s-1]"
            plt.title(
                f"{dataset} {ensemble}, {var_name} flux {mode}\n"
                + f"{timerange}, wt: {wt}"
            )
        im = iplt.contourf(cube, cmap=prcp_cmap)
        cb = plt.colorbar(im)
        cb.ax.tick_params(labelsize=8)
        cb.set_label(label=f"{var_name} {mode} {unit}")
    elif var_name == "tas":
        temp_cmap = get_colormap("temp")
        unit = "[K]"
        plt.title(
            f"{dataset} {ensemble}, 1000 hPa {var_name} {mode}\n"
            + f"{timerange}, wt: {wt}"
        )
        im = iplt.contourf(cube, cmap=temp_cmap)
        cb = plt.colorbar(im)
        cb.ax.tick_params(labelsize=8)
        cb.set_label(label=f"{var_name} {mode} {unit}")

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

    ax.set_extent([-15, 20, 27.5, 62.5])

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=":")

    plt.savefig(
        f"{local_path}/{wt_string}_{wt}{driver}_{dataset}_{ensemble}"
        f"_{var_name}_{mode}_{timerange}.png"
    )
    plt.savefig(
        f"{local_path}/{wt_string}_{wt}{driver}_{dataset}_{ensemble}_"
        f"{var_name}_{mode}_{timerange}.pdf"
    )
    plt.close()


def plot_corr_rmse_heatmaps(
    cfg: dict,
    pattern_correlation_matrix: np.array,
    rmse_matrix: np.array,
    dataset: str,
    timerange: str,
):
    """Plot heatmaps for correlation and rmse matrices.

    Args:
    ----
        cfg : dict
            Configuration dictionary from recipe
        pattern_correlation_matrix : np.array
            pattern correlation matrix
        rmse_matrix : np.array
            rmse matrix
        dataset : str
            string of dataset
        timerange : str
            string of timerange
    """
    output_path = f"{cfg.get('plot_dir')}/heatmaps"

    if not os.path.exists(f"{output_path}"):
        os.makedirs(f"{output_path}")

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
            f"{output_path}/correlation_matrix_{dataset}_{timerange}.png"
        )
        plt.savefig(
            f"{output_path}/correlation_matrix_{dataset}_{timerange}.pdf"
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


def get_colormap(colormap_string: str) -> ListedColormap:
    """Get colormaps based on string.

    Args:
    ----
        colormap_string : str
            String to get Colormaps for either
            psl, tas or precipitation.

    Returns
    -------
        ListedColormap
            Choosen Colormap
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
