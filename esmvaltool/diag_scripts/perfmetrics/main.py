import matplotlib as mpl
from esmvalcore import preprocessor as pp
import logging
import matplotlib.pyplot as plt
import esmvaltool.diag_scripts.shared as e
from esmvaltool.diag_scripts.shared import (
    # ProvenanceLogger,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
    select_metadata,
)

import iris
import numpy as np
import matplotlib.pyplot as plt


def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw=None, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if ax is None:
        ax = plt.gca()

    if cbar_kw is None:
        cbar_kw = {} 

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
            va="top", rotation_mode="anchor")

    # Turn spines off and create white grid.
    # ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="black", linestyle='-', linewidth=1)
    ax.tick_params(which="both", bottom=False, left=False)

    return im, cbar


def remove_reference(metas):
    return [meta for meta in metas if not meta.get("reference_for_metric", False)]


def plot(cfg, metas):
    x = cfg.get("x", "dataset")
    y = cfg.get("y", "variable_group")
    y_labels = group_metadata(metas, y).keys()
    x_labels = group_metadata(metas, x).keys()
    data = np.zeros((len(y_labels), len(x_labels)))
    for i, y_label in enumerate(y_labels):
        #y_metas = select_metadata(metas, variable_group=y_label)
        for j, x_label in enumerate(x_labels):
            selection = {x: x_label, y: y_label}
            meta = select_metadata(metas, **selection)[0]
            cube = iris.load_cube(meta["filename"])
            data[i][j] = cube.data 
    ax = heatmap(data, y_labels, x_labels, **{"cmap": "RdYlBu_r"}, cbar_kw={"extend": "both"})
    plt.savefig("heatmap.png")


def main(cfg):
    metas = cfg["input_data"].values()
    metas = remove_reference(metas)
    plot(cfg, metas)
    # grouped = group_metadata(metas, "project")
    # for group, group_metas in grouped.items():
    #     plot_group()
    

if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)