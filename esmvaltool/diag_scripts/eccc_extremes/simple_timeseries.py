"""
Diagnostic to reproduce fig. 8 from Malinina&Gillett (2024)

This diagnostic plots simple timeseries for individual models and
Multi-Model mean, weighted by the number of realizations in each model.
If there is more than one realization the min/maxshading is showed as
as a shading around the line.

The reqired variable groups are `reference` and other data variable
group. The preprocessed data should be a 1D timeseries cube.

Optional keywords for diagnostic are:
    `units`: unit flag which will be used for plotting
    `color_style`: name of teh color style file
    `mpl_style`: name of the matplotlib style file
    `figure_caption`: figure caption to be used
    `region`: name of the region

Author: Elizaveta Malinina (CCCma)
        elizaveta.malinina-rieger@ec.gc.ca
"""

import logging
import os

import iris
import matplotlib.pyplot as plt
import numpy as np

import esmvaltool.diag_scripts.shared.plot as eplot

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_figure,
)

# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def obtain_reference(data_group: list):
    """This function cretes a dictionary with reference data.

    Parameters:
    -----------
    data_group:
        list with metadata for reference group

    Returns:
    --------
    reference_dic:
        dictionary with the reference data. Dataset name is a keyword.
    ref_max:
        maximum value for all reference datasets
    ref_min:
        minimum value for all reference datasets
    """

    reference_dic = {}

    ref_fnames = group_metadata(data_group, "filename")

    mins = list()
    maxs = list()

    for dataset_f in ref_fnames.keys():
        dataset_n = ref_fnames[dataset_f][0]["dataset"]
        if len(group_metadata(data_group, "dataset")[dataset_n]) > 1:
            key = ref_fnames[dataset_f][0]["alias"]
        else:
            key = dataset_n
        ref_cb = iris.load_cube(dataset_f)
        # make key dataset if there's only one filename
        # of this dataset otherwise alias
        reference_dic[key] = ref_cb
        mins.append(ref_cb.collapsed("time", iris.analysis.MIN).data)
        maxs.append(ref_cb.collapsed("time", iris.analysis.MAX).data)

    ref_max = np.asarray(maxs).max()
    ref_min = np.asarray(mins).min()

    return reference_dic, ref_max, ref_min


def create_provenance(caption: str):
    """Creates provenance dictionary."""
    provenance_dic = {
        "authors": ["malinina_elizaveta"],
        "caption": caption,
        "references": ["malinina24"],
    }

    return provenance_dic


def plot_timeseries(
    data_dic: dict,
    reference_dic: dict,
    cfg: dict,
    min_val: float,
    max_val: float,
):
    """
    This function plots timeseries and the max/min values.

    Parameters:
    -----------
    data_dic:
        dictionary with the data to plot
    reference_dic:
        dictionary with the refrence data
    cfg:
        config dictionary, comes from ESMValCore
    min_val:
        minimum value for the xlims
    max_val:
        maximum value for the xlims
    """

    b_max = np.ceil(max_val * 1.05)
    b_min = np.floor(min_val * 1.05)

    mpl_st_file = eplot.get_path_to_mpl_style(cfg.get("mpl_style"))
    plt.style.use(mpl_st_file)

    for dataset in data_dic.keys():
        fig_ts, ax_ts = plt.subplots(nrows=1, ncols=1)

        model_var_data = list()
        weights = list()

        for i in range(0, len(data_dic[dataset])):
            var_cb = data_dic[dataset][i]
            model_var_data.append(var_cb.data)
            weights.append(
                var_cb.attributes["ensemble_weight"]
                * var_cb.attributes["reverse_dtsts_n"]
            )
        weights = np.array(weights)
        model_var_data = np.array(model_var_data)
        mean_var_arr = np.average(model_var_data, axis=0)
        min_var_arr = np.min(model_var_data, axis=0)
        max_var_arr = np.max(model_var_data, axis=0)

        years = [
            var_cb.coord("time").cell(i).point.year
            for i in range(var_cb.coord("time").shape[0])
        ]

        color_st = eplot.get_dataset_style(dataset, cfg.get("color_style"))
        ax_ts.plot(
            years, mean_var_arr, c=color_st["color"], zorder=3, label=dataset
        )
        if len(data_dic[dataset]) > 1:
            ax_ts.fill_between(
                years,
                min_var_arr,
                max_var_arr,
                color=color_st["color"],
                lw=0,
                alpha=0.25,
            )

        for ref_data in reference_dic.keys():
            ref_cb = reference_dic.get(ref_data)
            ref_color_st = eplot.get_dataset_style(
                ref_data, cfg.get("color_style")
            )
            ref_years = [
                ref_cb.coord("time").cell(i).point.year
                for i in range(ref_cb.coord("time").shape[0])
            ]
            ax_ts.plot(
                ref_years,
                ref_cb.data,
                label=ref_data,
                c=ref_color_st["color"],
                zorder=3,
            )

        ax_ts.axhline(color="grey", zorder=1)
        ax_ts.set_ylim(b_min, b_max)
        ax_ts.legend(loc=0, ncols=4, fancybox=False, frameon=False)

        variable = (
            cfg["var_label"] if cfg.get("var_label") else var_cb.var_name
        )
        exp_variable = variable.replace("_", " ")
        units = cfg["units"] if cfg.get("units") else var_cb.units
        region = cfg["region"] if cfg.get("units") else "region"

        ax_ts.set_ylabel(f"{exp_variable}, {units}")
        ax_ts.set_xlabel("year")

        default_caption = f"Timeseries of {exp_variable} in {region}"

        caption = (
            cfg["figure_caption"]
            if cfg.get("figure_caption")
            else default_caption
        )

        caption += f" ({dataset})"

        fig_ts.suptitle(caption)

        prov_dic = create_provenance(caption)

        plt.tight_layout()

        fig_path = os.path.join(
            cfg["plot_dir"], f"timeseries_{variable}_{region}_{dataset}"
        )

        save_figure(fig_path, prov_dic, cfg, fig_ts, close=True)

    return


def main(cfg):
    """Initiates diagnostic"""
    input_data = cfg["input_data"]

    groups = group_metadata(input_data.values(), "variable_group", sort=True)

    maxs, mins = [], []
    data_dic = {}

    reference_dic, ref_max, ref_min = obtain_reference(groups.pop("reference"))

    maxs.append(ref_max)
    mins.append(ref_min)

    remaining_metadata = []
    for k in groups.keys():
        remaining_metadata.extend(groups[k])

    datasets = group_metadata(remaining_metadata, "dataset")
    ens_var_cubelist = iris.cube.CubeList()
    for dataset in datasets.keys():
        filepaths = list(group_metadata(datasets[dataset], "filename").keys())
        n_real = len(filepaths)
        mod_var_cubelist = iris.cube.CubeList()
        for filepath in filepaths:
            mod_cb = iris.load_cube(filepath)
            # adding weights to the data cubes
            mod_cb.attributes["ensemble_weight"] = 1 / n_real
            mod_cb.attributes["reverse_dtsts_n"] = 1 / len(datasets)
            ens_var_cubelist.append(mod_cb)
            mod_var_cubelist.append(mod_cb)
            mins.append(mod_cb.collapsed("time", iris.analysis.MIN).data)
            maxs.append(mod_cb.collapsed("time", iris.analysis.MAX).data)
        data_dic[dataset] = mod_var_cubelist
    data_dic["Multi-Model"] = ens_var_cubelist

    plot_timeseries(
        data_dic,
        reference_dic,
        cfg,
        min_val=np.asarray(mins).min(),
        max_val=np.asarray(maxs).max(),
    )

    logger.info("Success")


if __name__ == "__main__":
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)
