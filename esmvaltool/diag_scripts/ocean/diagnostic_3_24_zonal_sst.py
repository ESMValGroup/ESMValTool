"""
Diagnostic for fig 3.24 in Chapter 3 of IPCC AR6 WGI.

This diagnostic either calculates bias with respect to the reference
dataset or plots the data and the reference dataset.
The diagnostic accepts the groups of the 1D datasets (with exception
of the resolved mask).
The additional keywords for the diagnostic provided through the recipe:

bias: True # if bias should be calculated (optional), default: false
mask: # if the data should be masked (optional), default: false
    flag: true
    type: simple # accepted values: 'simple' and 'resolved'
statistics:
    best_guess: mean # best guess statistic
    borders: # borders for the shading
    - operator: percentile
        percent: 5
    - operator: percentile
        percent: 95
mpl_style: ipcc_ar6_fgd # name of the matplotlib style file (optional)
caption: Zonal mean SST bias # Figure caption (optional)
color_style: sst_bias # name of the color style, colors are defined for
                      # variable groups

The resolved mask can have as many dimensions as one needs.

Initial development by: Lee de Mora (PML)
                        ledm@pml.ac.uk
Revised and corrected (15.01.2021): Elizaveta Malinina (CCCma)
                        elizaveta.malinina-rieger@ec.gc.ca
Rewritten and adapted to main branch v2.12+ (08.03.2025):
                            Elizaveta Malinina (CCCma)
                            elizaveta.malinina-rieger@ec.gc.ca
"""

import logging
import os
import sys

import esmvalcore.preprocessor as eprep
import iris
import iris.plot
import matplotlib.pyplot as plt

import esmvaltool.diag_scripts.shared.plot as eplot
from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_figure,
    select_metadata,
)

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


class Data4Analyis:
    """Data class which later will be used for plotting.

    Attributes
    ----------
    name: str
        name of the data group for which class is created
    bias: bool
        flag, showing if bias should be calculated
    mask: bool
        flag if data should be masked
    mask_type: str
        type of the mask: 'simple' or 'resolved'
    reference: str
        name of the reference dataset
    ref_cube: iris.cube.Cube
        iris cube with the reference dataset data
    data: iris.cube.CubeList
        iris cubelist with the data
    best_guess: iris.cube.Cube
        iris cube with the best guess to be plotted
    border1: iris.cube.Cube
        iris cube with the bottom/top for the shading
    border2: iris.cube.Cube
        iris cube with the top/bottom for the shading"""

    def __init__(
        self,
        name: str,
        group: list,
        cfg: dict,
        mask_meta: list[dict] | None,
    ):
        """Initialization of the class.

        Parameters:
        -----------
        name:
            Name of the data class, comes from the name of the group
        group:
            List with the metadata of the variable group
        cfg:
            Config dictionary coming from ESMValCore
        mask_type:
            Type of the mask to be used, if no mask to be used set to False
        mask_meta:
            Metadata for the mask group if mask_type='resolved'"""
        self.name = name
        self.bias = bool(cfg.get("bias"))
        self.mask = cfg.get("mask").get("flag") if cfg.get("mask") else False
        self.mask_type = cfg["mask"].get("type") if self.mask else False
        self.determine_reference(group)
        self.obtain_data(group, mask_meta)
        stats = cfg.get("data_statistics")
        if not stats:
            raise ValueError(
                "statistics dictionary should be provided in the recipe. "
                "The keywords 'best_guess' and 'borders' should be provided."
            )
        self.calculate_statistics(stats)

    def determine_reference(self, group: list):
        """Determines the reference dataset from the data group.

        Parameters:
        -----------
        group:
            List with the metadata of the variable group

        Raises:
        -------
        ValueError
            if more than one or no reference datasets have been provided"""
        if self.bias or self.mask:
            reference = list(group_metadata(group, "reference_dataset").keys())
            if len(reference) > 1:
                raise ValueError("More then one reference dataset was given")
            self.reference = reference[0]
            if reference:
                ref_metadata = select_metadata(group, dataset=self.reference)[
                    0
                ]
                self.ref_cube = iris.load_cube(ref_metadata["filename"])
                group.remove(ref_metadata)
            else:
                raise ValueError("No reference dataset has been provided")

    def obtain_data(self, group: list, mask_meta: list[dict] | None):
        """Obtains data for further statistics calculation.

        Parameters:
        -----------
        group:
            List with the metadata of the variable group
        mask_meta:
            List with the mask_metadata in case mask_type='resolved'"""
        files = list(group_metadata(group, "filename", sort=True).keys())
        data_cblst = iris.load(files)
        self.data = data_cblst

        if self.mask:
            self.mask_data(mask_meta)

        if self.bias:
            self.data = eprep.bias(self.data, reference=self.ref_cube)

        return

    def mask_data(self, mask_meta: list[dict] | None):
        """Masks data.

        Parameters:
        -----------
        mask_meta:
            List with the mask_metadata in case mask_type='resolved'

        Raises:
        -------
        ValueError
            if more than one datasets for resolved mask are provided
            or if data cubes have more than one dimension or if an
            unsupported type of mask is provided"""
        if self.mask_type == "simple":
            mask = self.ref_cube.data.mask
        elif self.mask_type == "resolved":
            mask_files = list(
                group_metadata(mask_meta, "filename", sort=True).keys()
            )
            if len(mask_files) > 1:
                raise ValueError(
                    "More than one dataset for the resolved mask "
                    "has been provided. Only one is supported."
                )
            mask_cb = iris.load_cube(mask_files[0])
            data_coord = self.data[0].dim_coords
            if len(data_coord) > 1:
                raise ValueError(
                    "The data cubes have more than one "
                    "cordinates. Only flat cubes are supported."
                )
            data_coord = data_coord[0]
            # determine the number of dimension over which the data is calculated
            dim = [c.name() for c in mask_cb.dim_coords].index(
                data_coord.name()
            )
            mask = []
            # if there is less than half of data over the dimension the cell
            # is masked
            for i in range(mask_cb.shape[dim]):
                mask.append(
                    mask_cb.data[(slice(None),) * dim + (i,)].count()
                    <= 0.5 * len(mask_cb.data[(slice(None),) * dim + (i,)])
                )
            # masking the reference cube
            self.ref_cube.data.mask = self.ref_cube.data.mask | mask
        else:
            raise ValueError(
                f"Mask type {self.mask_type} is not supported. "
                "Only 'simple' and 'resolved' are supported."
            )

        for n_cb in enumerate(self.data):
            self.data[n_cb[0]].data.mask = self.data[n_cb[0]].data.mask | mask

        return

    def calculate_statistics(self, stats: dict):
        """Calculates statistics which later will be plotted.
        Parameters:
        -----------
        stats:
            dictionary with the statistics which will be calculated.
            Dictionary should have keywords 'best_guess' and 'borders'."""
        if len(self.data) > 1:
            bg_dic = eprep.multi_model_statistics(
                self.data,
                span="full",
                statistics=[stats["best_guess"]],
                ignore_scalar_coords=True,
            )
            self.best_guess = bg_dic[list(bg_dic.keys())[0]]
            bord_dic = eprep.multi_model_statistics(
                self.data,
                span="full",
                statistics=stats["borders"],
                ignore_scalar_coords=True,
            )
            self.border1 = bord_dic[list(bord_dic.keys())[0]]
            self.border2 = bord_dic[list(bord_dic.keys())[1]]
        else:
            # if just one dataset is in the group there is no need
            # to calculate the borders
            self.best_guess = self.data[0]
            self.border1 = None
            self.border2 = None


def create_provenance(caption: str):
    """Creates provenance dictionary."""
    provenance_dic = {
        "authors": ["malinina_elizaveta", "demora_lee"],
        "caption": caption,
        "references": ["eyring21ipcc"],
    }

    return provenance_dic


def plot_bias_plot(data_list: list[Data4Analyis], cfg: dict):
    """Plots the diagnostic figure.
    Parameters:
    -----------
    data_list:
        List with the data classes which will be plotted
    cfg:
        Config dictionary coming from ESMValCore"""
    caption = cfg.get("caption") if cfg.get("caption") else ""
    prov_dic = create_provenance(caption=caption)

    st_file = eplot.get_path_to_mpl_style(cfg.get("mpl_style"))
    plt.style.use(st_file)

    fig = plt.figure(figsize=(6, 2.5))

    for data in data_list:
        if data.best_guess.dim_coords[0].name() == "longitude":
            # to make the Pacific and Atlantic continuous
            data.best_guess = data.best_guess.intersection(
                longitude=(20.0, 380.0)
            )
            data.border1 = data.border1.intersection(longitude=(20.0, 380.0))
            data.border2 = data.border2.intersection(longitude=(20.0, 380.0))
            data.ref_cube = data.ref_cube.intersection(longitude=(20.0, 380.0))
        data_col = eplot.get_dataset_style(data.name, cfg.get("color_style"))
        iris.plot.plot(data.best_guess, color=data_col["color"])
        iris.plot.fill_between(
            data.best_guess.dim_coords[0],
            data.border1,
            data.border2,
            alpha=0.2,
            linewidth=0,
            color=data_col["color"],
        )

    if cfg.get("bias"):
        xlim = fig.axes[0].get_xlim()
        plt.hlines(0, *xlim, colors="silver", linestyles="dashed", zorder=1)
        plt.xlim(*xlim)  # to make the xlims the same as before
    else:
        ref_col = eplot.get_dataset_style(
            data_list[0].reference, cfg.get("color_style")
        )
        iris.plot.plot(data_list[0].ref_cube, c=ref_col["color"])

    plt.title(caption)
    plt.tight_layout()
    fig_path = os.path.join(cfg["plot_dir"], "figure")
    save_figure(fig_path, prov_dic, cfg, fig, close=True)


def main(cfg: dict):
    """Diagnostic itself."""
    input_data = cfg["input_data"]

    groups = group_metadata(input_data.values(), "variable_group", sort=True)

    mask_type = cfg.get("mask").get("type") if cfg.get("mask") else False
    mask_group = (
        cfg.get("mask").get("group") if mask_type == "resolved" else None
    )
    mask_meta = groups.pop(mask_group) if mask_group else None

    data_list = []
    for group, group_cont in groups.items():
        group_data = Data4Analyis(
            name=group,
            group=group_cont,
            cfg=cfg,
            mask_meta=mask_meta,
        )
        data_list.append(group_data)

    plot_bias_plot(data_list, cfg)

    logger.info("Success")


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
