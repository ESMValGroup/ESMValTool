"""Lifetime diagnostic to show multiple datasets in one plot.

Description
-----------
This diagnostic can be used to visualize lifetime of multiple datasets in
one plot.

Input data needs to be 4D (time,
air_pressure/atmosphere_hybrid_sigma_pressure_coordinate, latitude, longitude).

For some plot types, a reference dataset can be defined. For this, use the
facet ``reference_for_lifetime_diags: true`` in the definition of the dataset
in the recipe. Note that at most one reference dataset is supported.

Currently supported plot types (use the option ``plots`` to specify them):
    - Time series (plot type ``timeseries``): all datasets are plotted in one
      single figure.

Author
------
Franziska Winterstein (DLR, Germany)

Configuration options in recipe
-------------------------------
facet_used_for_labels: str, optional (default: 'dataset')
    Facet used to label different datasets in plot titles and legends. For
    example, ``facet_used_for_labels: dataset`` will use dataset names in plot
    titles and legends; ``facet_used_for_labels: exp`` will use experiments in
    plot titles and legends. In addition, ``facet_used_for_labels`` is used to
    select the correct ``plot_kwargs`` for the different datasets (see
    configuration options for the different plot types below).
figure_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.figure`. By
    default, uses ``constrained_layout: true``.
m_air: float, optional (default: 28.970)
    Molar mass of air in g/mol.
molarmass: float
    Molar mass of reactant in g/mol.
oxidant: dict
    Oxidant reaction rates for chemical species, e.g., {"oh": {"A": 1.85e-20,
    "ER": 987.0, "b": 2.82}}.
plots: dict, optional
    Plot types plotted by this diagnostic (see list above). Dictionary keys
    must be ``timeseries``.
    Dictionary values are dictionaries used as options for the corresponding
    plot. The allowed options for the different plot types are given below.
plot_filename: str, optional
    Filename pattern for the plots.
    Defaults to ``{plot_type}_{real_name}_{dataset}_{mip}_{exp}_{ensemble}``.
    All tags (i.e., the entries in curly brackets, e.g., ``{dataset}``, are
    replaced with the corresponding tags).
plot_folder: str, optional
    Path to the folder to store figures. Defaults to
    ``{plot_dir}/../../{dataset}/{exp}/{modeling_realm}/{real_name}``.  All
    tags (i.e., the entries in curly brackets, e.g., ``{dataset}``, are
    replaced with the corresponding tags). ``{plot_dir}`` is replaced with the
    default ESMValTool plot directory (i.e.,
    ``output_dir/plots/diagnostic_name/script_name/``, see
    :ref:`esmvalcore:outputdata`).
reactant: str
    Name of reactant, e.g., "CH4".
regions: list of str, optional (default: ["TROP"])
    Regions considered for line plots. Currently, only "TROP" (troposphere) and
    "STRA" (stratosphere) are supported.
savefig_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.savefig`. By
    default, uses ``bbox_inches: tight, dpi: 300, orientation: landscape``.
seaborn_settings: dict, optional
    Options for :func:`seaborn.set_theme` (affects all plots). By default, uses
    ``style: ticks``.
units: str, optional (default: "years")
    Target units of the lifetime.
weight_type: str
    Type of weighting. If "mass CH4", use CH4 mass as weights. Otherwise, use
    equal weights.

Configuration options for plot type ``timeseries``
--------------------------------------------------
annual_mean: str, optional (default: False)
    Optional switch to turn on annual means to be displayed  'only' or
    additional 'both' to the original timeseries. If not set or set to
    ``False`` only the original timeseries is shown.
annual_mean_kwargs: dict, optional
    Optional keyword arguments for :func:`iris.plot.plot` for plotting annual
    means. These keyword arguments update (and potentially overwrite) the
    ``plot_kwargs`` for the annual mean plots. Use ``annual_mean_kwargs`` to
    not show annual means.
by_timestep: bool, optional (default: False)
    Calculate lifetime for each time step individually. Slower, but less
    memory-intensive.
display_mean: bool, optional (default: False)
    Show mean over entire period in plot label.
gridline_kwargs: dict, optional
    Optional keyword arguments for grid lines. By default, ``color: lightgrey,
    alpha: 0.5`` are used. Use ``gridline_kwargs: false`` to not show grid
    lines.
legend_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.legend`. Use
    ``legend_kwargs: false`` to not show legends.
plot_kwargs: dict, optional
    Optional keyword arguments for :func:`iris.plot.plot`. Dictionary keys are
    elements identified by ``facet_used_for_labels`` or ``default``, e.g.,
    ``CMIP6`` if ``facet_used_for_labels: project`` or ``historical`` if
    ``facet_used_for_labels: exp``. Dictionary values are dictionaries used as
    keyword arguments for :func:`iris.plot.plot`. String arguments can include
    facets in curly brackets which will be derived from the corresponding
    dataset, e.g., ``{project}``, ``{short_name}``, ``{exp}``. Examples:
    ``default: {linestyle: '-', label: '{project}'}, CMIP6: {color: red,
    linestyle: '--'}, OBS: {color: black}``.
pyplot_kwargs: dict, optional
    Optional calls to functions of :mod:`matplotlib.pyplot`. Dictionary keys
    are functions of :mod:`matplotlib.pyplot`. Dictionary values are used as
    single argument for these functions. String arguments can include facets in
    curly brackets which will be derived from the datasets plotted in the
    corresponding plot, e.g., ``{short_name}``, ``{exp}``. Facets like
    ``{project}`` that vary between the different datasets will be transformed
    to something like  ``ambiguous_project``. Examples: ``title: 'Awesome Plot
    of {long_name}'``, ``xlabel: '{short_name}'``, ``xlim: [0, 5]``.

.. hint::

   Extra arguments given to the recipe are ignored, so it is safe to use yaml
   anchors to share the configuration of common arguments with other monitor
   diagnostic script.

"""

import logging
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import iris
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from iris.coord_categorisation import add_year

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.lifetime.lifetime_func import (
    calculate_gridmassdry,
    calculate_lifetime,
    calculate_reaction_rate,
    calculate_rho,
    climatological_tropopause,
    create_press,
)
from esmvaltool.diag_scripts.monitor.monitor_base import MonitorBase
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    io,
    run_diagnostic,
)

logger = logging.getLogger(Path(__file__).stem)


def add_model_level(cube):
    """Add a simplified level coordinate.

    Parameters
    ----------
    cube : iris.cube.Cube
        Input cube.

    """
    z_coord = cube.coord(axis="z", dim_coords=True)
    n_levels = z_coord.shape[0]
    levels = np.array(range(n_levels, 0, -1), ndmin=1)
    level_coords = iris.coords.AuxCoord(
        levels,
        bounds=None,
        standard_name="model_level_number",
        units="1",
    )
    cube.add_aux_coord(level_coords, cube.coord_dims(z_coord))


class CH4Lifetime(MonitorBase):
    """Diagnostic to plot ch4 lifetime."""

    def __init__(self, config):
        """Initialize class member."""
        super().__init__(config)

        # Get default stettings
        self.cfg = deepcopy(self.cfg)

        self.cfg.setdefault("facet_used_for_labels", "dataset")
        self.cfg.setdefault("figure_kwargs", {"constrained_layout": True})
        self.cfg.setdefault(
            "savefig_kwargs",
            {
                "bbox_inches": "tight",
                "dpi": 300,
                "orientation": "landscape",
            },
        )
        self.cfg.setdefault("seaborn_settings", {"style": "ticks"})
        logger.info(
            "Using facet '%s' to create labels",
            self.cfg["facet_used_for_labels"],
        )
        self.cfg.setdefault("regions", ["TROP"])
        self.cfg.setdefault("units", "years")

        # set default molarmasses
        self.cfg.setdefault("m_air", 28.970)  # [g_air/mol_air]

        # Load input data
        self.input_data_dataset = self._calculate_coefficients()

        # base info
        self.info = {"short_name": f"tau_{self._get_name('reactant').upper()}"}
        oxidants = [ox.upper() for ox in self._get_name("oxidant")]
        self.info["long_name"] = (
            "Lifetime of"
            f" {self._get_name('reactant').upper()}"
            " with respect to"
            f" {', '.join(oxidants)}"
        )
        self.units = self.cfg["units"]
        self.info["units"] = self.cfg["units"]

        # Check given plot types and set default settings for them
        self.supported_plot_types = [
            "timeseries",
        ]
        for plot_type, plot_options in self.plots.items():
            if plot_type not in self.supported_plot_types:
                raise ValueError(
                    f"Got unexpected plot type '{plot_type}' for option "
                    f"'plots', expected one of {self.supported_plot_types}",
                )
            if plot_options is None:
                self.plots[plot_type] = {}

            # Default options for the differ N MMNMJHJ JMHTHent plot types
            if plot_type == "timeseries":
                self.plots[plot_type].setdefault("annual_mean", False)
                self.plots[plot_type].setdefault("annual_mean_kwargs", {})
                self.plots[plot_type].setdefault("display_mean", False)
                self.plots[plot_type].setdefault("gridline_kwargs", {})
                self.plots[plot_type].setdefault("legend_kwargs", {})
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("pyplot_kwargs", {})
                self.plots[plot_type].setdefault("by_timestep", False)

        # Check that facet_used_for_labels is present for every dataset
        for dataset in self.input_data:
            if self.cfg["facet_used_for_labels"] not in dataset:
                raise ValueError(
                    f"facet_used_for_labels "
                    f"'{self.cfg['facet_used_for_labels']}' not present for "
                    f"the following dataset:\n{pformat(dataset)}",
                )

        # Load seaborn settings
        sns.set_theme(**self.cfg["seaborn_settings"])

    def _add_colorbar(
        self,
        plot_type,
        plot_left,
        plot_right,
        axes_left,
        axes_right,
        dataset_left,
        dataset_right,
    ):
        """Add colorbar(s) for plots."""
        fontsize = self.plots[plot_type]["fontsize"]
        cbar_kwargs = self._get_cbar_kwargs(plot_type)
        cbar_label_left = self._get_cbar_label(plot_type, dataset_left)
        cbar_label_right = self._get_cbar_label(plot_type, dataset_right)

        # Create one common colorbar for the top panels
        # Note: Increase aspect ratio for nicer looks
        if self.plots[plot_type]["common_cbar"]:
            if "aspect" in cbar_kwargs:
                cbar_kwargs["aspect"] += 20.0
            cbar = plt.colorbar(
                plot_left,
                ax=[axes_left, axes_right],
                **cbar_kwargs,
            )
            cbar.set_label(cbar_label_left, fontsize=fontsize)
            cbar.ax.tick_params(labelsize=fontsize)

        # Create two separate colorbars for the top panels
        else:
            cbar_left = plt.colorbar(plot_left, ax=axes_left, **cbar_kwargs)
            cbar_left.set_label(cbar_label_left, fontsize=fontsize)
            cbar_left.ax.tick_params(labelsize=fontsize)
            cbar_right = plt.colorbar(plot_right, ax=axes_right, **cbar_kwargs)
            cbar_right.set_label(cbar_label_right, fontsize=fontsize)
            cbar_right.ax.tick_params(labelsize=fontsize)

    def _get_custom_mpl_rc_params(self, plot_type):
        """Get custom matplotlib rcParams."""
        fontsize = self.plots[plot_type]["fontsize"]
        custom_rc_params = {
            "axes.titlesize": fontsize + 2.0,
            "axes.labelsize": fontsize,
            "xtick.labelsize": fontsize,
            "ytick.labelsize": fontsize,
        }
        return custom_rc_params

    def _get_label(self, dataset):
        """Get label of dataset."""
        return dataset[self.cfg["facet_used_for_labels"]]

    def _get_cbar_kwargs(self, plot_type, bias=False):
        """Get colorbar kwargs."""
        cbar_kwargs = deepcopy(self.plots[plot_type]["cbar_kwargs"])
        if bias:
            cbar_kwargs.update(self.plots[plot_type]["cbar_kwargs_bias"])
        return deepcopy(cbar_kwargs)

    def _get_cbar_label(self, plot_type, dataset, bias=False):
        """Get colorbar label."""
        if bias:
            cbar_label = self.plots[plot_type]["cbar_label_bias"]
            descr = f"cbar_label_bias of {plot_type} '{cbar_label}'"
        else:
            cbar_label = self.plots[plot_type]["cbar_label"]
            descr = f"cbar_label of {plot_type} '{cbar_label}'"
        cbar_label = self._fill_facet_placeholders(cbar_label, dataset, descr)
        return cbar_label

    def _get_gridline_kwargs(self, plot_type):
        """Get gridline kwargs."""
        gridline_kwargs = self.plots[plot_type]["gridline_kwargs"]
        return deepcopy(gridline_kwargs)

    def _get_plot_func(self, plot_type):
        """Get plot function."""
        plot_func = self.plots[plot_type]["plot_func"]
        if not hasattr(iris.plot, plot_func):
            raise AttributeError(
                f"Got invalid plot function '{plot_func}' for plotting "
                f"{plot_type}, expected function of iris.plot",
            )
        logger.info(
            "Creating %s plots using function '%s'",
            plot_type,
            plot_func,
        )
        return getattr(iris.plot, plot_func)

    def _get_plot_kwargs(self, plot_type, dataset, bias=False):
        """Get keyword arguments for plot functions."""
        all_plot_kwargs = self.plots[plot_type]["plot_kwargs"]
        all_plot_kwargs = deepcopy(all_plot_kwargs)

        # First get default kwargs, then overwrite them with dataset-specific
        # ones
        plot_kwargs = all_plot_kwargs.get("default", {})
        label = self._get_label(dataset)
        plot_kwargs.update(all_plot_kwargs.get(label, {}))

        # For bias plots, overwrite the kwargs with bias-specific option
        if bias:
            bias_kwargs = self.plots[plot_type]["plot_kwargs_bias"]
            plot_kwargs.update(bias_kwargs)

        # Replace facets with dataset entries for string arguments
        for key, val in plot_kwargs.items():
            if isinstance(val, str):
                val = self._fill_facet_placeholders(
                    val,
                    dataset,
                    f"plot_kwargs of {plot_type} '{key}: {val}'",
                )
                plot_kwargs[key] = val

        # Default settings for different plot types
        if plot_type == "timeseries":
            plot_kwargs.setdefault("label", label)

        return deepcopy(plot_kwargs)

    def _calculate_coefficients(self):
        """Load data and calculate coefficients."""
        self.input_data = list(self.cfg["input_data"].values())

        # loops over different variables ch4, oh, ta etc.
        input_data_dataset = {}
        for dataset in self._get_all_datasets(self.input_data):
            input_data_dataset[dataset] = self._get_dataset_data(dataset)[0]
            input_data_dataset[dataset]["dataset_data"] = (
                self._get_dataset_data(dataset)
            )

            variables = {}
            for variable in input_data_dataset[dataset]["dataset_data"]:
                filename = variable["filename"]

                logger.info("Loading %s", filename)
                cube = iris.load_cube(filename, variable["short_name"])

                # Fix time coordinate if present
                if cube.coords("time", dim_coords=True):
                    ih.unify_time_coord(cube)

                # Remove surface air pressure coordinate if necessary
                if cube.coords("Surface Air Pressure"):
                    cube.remove_coord("Surface Air Pressure")

                # Cube for each variable
                variables[variable["short_name"]] = cube

            rho = calculate_rho(variables)

            oxidant = {ox: variables[ox] for ox in self._get_name("oxidant")}
            self._set_oxidant_defaults(oxidant)

            reaction = self._calculate_reaction(
                oxidant,
                rho,
                variables["ta"],
                self._get_name("reactant"),
            )

            # Set Z-coordinate
            if reaction.coords("air_pressure", dim_coords=True):
                self.z_coord = "air_pressure"
                z_coord = reaction.coord("air_pressure", dim_coords=True)
                z_coord.attributes["positive"] = "down"
                z_coord.convert_units("Pa")
                use_z_coord = "air_pressure"
            elif reaction.coords(
                "atmosphere_hybrid_sigma_pressure_coordinate",
                dim_coords=True,
            ):
                self.z_coord = "atmosphere_hybrid_sigma_pressure_coordinate"
                z_coord = reaction.coord(
                    "atmosphere_hybrid_sigma_pressure_coordinate",
                    dim_coords=True,
                )
                z_coord.attributes["positive"] = "down"
                use_z_coord = "model_level_number"
            else:
                raise NotImplementedError(
                    "Lifetime calculation is not implemented for the present "
                    "type of vertical coordinate",
                )

            if not {"TROP", "STRA"}.isdisjoint(
                self.cfg["regions"],
            ) and not {"timeseries"}.isdisjoint(self.plots):
                # calculate climatological tropopause pressure (tp_clim)
                # but only if no tropopause is given by data
                if "ptp" not in variables and "tp_i" not in variables:
                    tropopause = climatological_tropopause(
                        variables["ta"][:, 0, :, :],
                    )

                # If z_coord is defined as:
                #     - air_pressure, use:
                #          - ptp and air_pressure
                #          - tp_clim and air_pressure
                #     - atmosphere_hybrid_sigma_pressure_coordinate, use:
                #          - tp_i and model_level_number
                #          - ptp and (derived) air_pressure
                #          - tp_clim and (derived) air_pressure
                if z_coord.name() == "air_pressure":
                    tropopause = variables.get("ptp", tropopause)
                elif (
                    z_coord.name()
                    == "atmosphere_hybrid_sigma_pressure_coordinate"
                ):
                    if "tp_i" in variables:
                        tropopause = variables["tp_i"]
                    else:
                        tropopause = variables.get("ptp", tropopause)
                        # fall back to air_pressure
                        use_z_coord = "air_pressure"

                input_data_dataset[dataset]["tropopause"] = tropopause

            weight = self._define_weight(variables)

            if reaction.coords(
                "atmosphere_hybrid_sigma_pressure_coordinate",
                dim_coords=True,
            ):
                add_model_level(weight)
                add_model_level(reaction)

            input_data_dataset[dataset]["z_coord"] = z_coord
            input_data_dataset[dataset]["use_z_coord"] = use_z_coord
            input_data_dataset[dataset]["variables"] = variables
            input_data_dataset[dataset]["reaction"] = reaction
            input_data_dataset[dataset]["weight"] = weight

        return input_data_dataset

    def _set_oxidant_defaults(self, oxidant):
        """Set the defaults for the oxidant reaction rates."""
        for name in oxidant:
            # set default reaction
            if name == "oh":
                self.cfg["oxidant"][name].setdefault("A", 1.85e-20)
                self.cfg["oxidant"][name].setdefault("ER", 987.0)
                self.cfg["oxidant"][name].setdefault("b", 2.82)
            else:
                if (
                    "A" not in self.cfg["oxidant"][name]
                    or "ER" not in self.cfg["oxidant"][name]
                ):
                    raise KeyError(
                        "No sufficient reaction coefficients are given",
                    )
                self.cfg["oxidant"][name].setdefault("b", None)

    def _get_all_datasets(self, input_data):
        """Return (sorted) list of datasets."""
        datasets = []
        for data in input_data:
            datasets.append(data[self.cfg["facet_used_for_labels"]])
        datasets = list(set(datasets))
        datasets.sort()
        return datasets

    def _get_dataset_data(self, dataset):
        """Return input data corresponding to the chosen dataset."""
        input_data = []
        for data in self.input_data:
            if data[self.cfg["facet_used_for_labels"]] == dataset:
                input_data.append(data)

        return input_data

    def _get_name(self, case="reactant"):
        """Return variable name.

        Return the name of the reactant or the oxidant respectively.
        """
        if isinstance(self.cfg[case], dict):
            name = list(self.cfg[case])
        else:
            name = self.cfg[case]

        return name

    def _calculate_reaction(self, oxidant, rho, temp, name_reactant):
        """Calculate product of reaction rate and oxidant."""
        reaction = 0.0
        for name, oxid in oxidant.items():
            oxid_in_molec = oxid * rho
            reaction_rate = calculate_reaction_rate(
                temp,
                f"{name_reactant.upper()}+{name.upper()}",
                self.cfg["oxidant"][name]["A"],
                self.cfg["oxidant"][name]["ER"],
                self.cfg["oxidant"][name]["b"],
            )
            reaction = reaction + reaction_rate * oxid_in_molec
        # test for correct units
        if not reaction.units == "s-1":
            raise ValueError(
                "The units of the reaction rate is"
                " not consistent. Check input variables.",
            )

        return reaction

    def _define_weight(self, variables):
        """Define used weights in the lifetime calculation.

        Currently only one weight type is implemented. Any other options result
        in weights = 1. (no weighting)

        weight type:
         mass CH4: mass of CH4 -> convert to kg per gridbox
        """
        if self.cfg["weight_type"] == "mass CH4":
            weight = self._convert_to_mass("ch4", variables)
        else:
            weight = 1.0

        return weight

    def _convert_to_mass(self, varname, variables):
        """Convert to kg per gridbox.

        Used constants:
        m_var Molarmass of current species
        """
        # molar mass constants
        m_var = self.cfg["molarmass"]

        if "grmassdry" in variables:
            grmassdry = variables["grmassdry"]
        else:
            hus = variables["hus"]
            if "press" in variables:
                press = variables["press"]
            else:
                # create 4D pressure variable from
                # pressure (auxilliary) coordinate
                press = create_press(hus)

            grmassdry = calculate_gridmassdry(press, hus, self.z_coord)

        var = variables[varname] * grmassdry * (m_var / self.cfg["m_air"])

        return var

    def get_regional_plot_path(self, plot_type, dataset, region):
        """Get plot path."""
        plot_path = Path(super().get_plot_path(plot_type, dataset))
        filename = f"{plot_path.stem}_{region}{plot_path.suffix}"
        return plot_path.parent / filename

    def _process_pyplot_kwargs(self, plot_type, dataset):
        """Process functions for :mod:`matplotlib.pyplot`."""
        pyplot_kwargs = self.plots[plot_type]["pyplot_kwargs"]
        for func, arg in pyplot_kwargs.items():
            if isinstance(arg, str):
                arg = self._fill_facet_placeholders(
                    arg,
                    dataset,
                    f"pyplot_kwargs of {plot_type} '{func}: {arg}'",
                )
            if arg is None:
                getattr(plt, func)()
            else:
                getattr(plt, func)(arg)

    @staticmethod
    def _fill_facet_placeholders(string, dataset, description):
        """Fill facet placeholders."""
        try:
            string = string.format(**dataset)
        except KeyError as exc:
            raise ValueError(
                f"Not all necessary facets in {description} available for "
                f"dataset\n{pformat(dataset)}",
            ) from exc
        return string

    @staticmethod
    def _get_multi_dataset_facets(datasets):
        """Derive common facets for multiple datasets."""
        all_keys = {key for dataset in datasets for key in dataset}
        multi_dataset_facets = {}
        for key in all_keys:
            if all(d.get(key) == datasets[0].get(key) for d in datasets):
                multi_dataset_facets[key] = datasets[0][key]
            else:
                multi_dataset_facets[key] = f"ambiguous_{key}"
        return multi_dataset_facets

    @staticmethod
    def _get_reference_dataset(datasets, short_name):
        """Extract reference dataset."""
        ref_datasets = [
            d
            for d in datasets.values()
            if d.get("reference_for_lifetime_diags", False)
        ]
        if len(ref_datasets) > 1:
            raise ValueError(
                f"Expected at most 1 reference dataset (with "
                f"'reference_for_lifetime_diags: true' for variable "
                f"'{short_name}', got {len(ref_datasets):d}",
            )
        if ref_datasets:
            return ref_datasets[0]
        return None

    def create_timeseries_plot(self, region, input_data, base_datasets):
        """Create time series plot."""
        plot_type = "timeseries"
        if plot_type not in self.plots:
            return

        if not base_datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        logger.info("Plotting %s", plot_type)
        fig = plt.figure(**self.cfg["figure_kwargs"])
        axes = fig.add_subplot()

        # Plot all datasets in one single figure
        ancestors = []
        cubes = {}

        for label, dataset in input_data.items():
            ancestors.extend(
                variable["filename"] for variable in dataset["dataset_data"]
            )

            # Call by timestep will take longer, but it is less memory
            # intensive
            if self.plots[plot_type]["by_timestep"]:
                slice_dataset = {}
                slice_dataset["z_coord"] = dataset["z_coord"]
                slice_dataset["use_z_coord"] = dataset["use_z_coord"]
                cube_slices = iris.cube.CubeList()
                for reaction_slice, weight_slice, tp_slice in zip(
                    dataset["reaction"].slices_over("time"),
                    dataset["weight"].slices_over("time"),
                    dataset["tropopause"].slices_over("time"),
                    strict=True,
                ):
                    slice_dataset["reaction"] = reaction_slice
                    slice_dataset["weight"] = weight_slice
                    slice_dataset["tropopause"] = tp_slice

                    cube_slice = calculate_lifetime(
                        slice_dataset,
                        plot_type,
                        region,
                    )
                    # make it real to reduce memory demand later
                    cube_slice.data  # noqa: B018
                    cube_slices.append(cube_slice)

                cube = cube_slices.merge_cube()
            else:
                cube = calculate_lifetime(dataset, plot_type, region)

            # convert units
            cube.convert_units(self.info["units"])

            cubes[label] = cube

            # Plot original time series
            plot_kwargs = self._get_plot_kwargs(
                plot_type,
                base_datasets[label],
            )
            plot_kwargs["axes"] = axes

            if self.plots[plot_type]["display_mean"] is not False:
                mean = cube.collapsed("time", iris.analysis.MEAN).data
                plot_kwargs["label"] = f"{plot_kwargs['label']} ({mean:.2f})"

            # Plot annual means if desired
            annual_mean = self.plots[plot_type]["annual_mean"]
            if annual_mean in [False, "both"]:
                iris.plot.plot(cube, **plot_kwargs)
                plot_kwargs.pop("label", None)
            if annual_mean in ["both", "only"]:
                logger.info("Plotting annual means")
                if not cube.coords("year"):
                    add_year(cube, "time")
                annual_mean_cube = cube.aggregated_by(
                    "year",
                    iris.analysis.MEAN,
                )

                plot_kwargs.update(self.plots[plot_type]["annual_mean_kwargs"])
                iris.plot.plot(annual_mean_cube, **plot_kwargs)
            else:
                raise ValueError(
                    "Unknown option for annual_mean; choose between False, "
                    "'both', or 'only'",
                )

        # Default plot appearance
        multi_dataset_facets = self._get_multi_dataset_facets(
            list(base_datasets.values()),
        )
        axes.set_title(f"{self.info['long_name']} in region {region}")
        axes.set_xlabel("Time")
        axes.set_ylabel(
            f"τ({self._get_name('reactant').upper()}) [{self.info['units']}]",
        )
        gridline_kwargs = self._get_gridline_kwargs(plot_type)
        if gridline_kwargs is not False:
            axes.grid(**gridline_kwargs)

        # Legend
        legend_kwargs = self.plots[plot_type]["legend_kwargs"]
        if legend_kwargs is not False:
            axes.legend(**legend_kwargs)

        # Customize plot appearance
        self._process_pyplot_kwargs(plot_type, multi_dataset_facets)

        # Save plot
        plot_path = self.get_regional_plot_path(
            plot_type,
            multi_dataset_facets,
            region,
        )
        fig.savefig(plot_path, **self.cfg["savefig_kwargs"])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save netCDF file
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)
        var_attrs = {
            "short_name": self.info["short_name"],
            "long_name": self.info["long_name"],
            "units": self.info["units"],
        }
        io.save_1d_data(cubes, netcdf_path, "time", var_attrs)

        # Provenance tracking
        caption = (
            f"Time series of {multi_dataset_facets['long_name']} for various "
            f"datasets in region {region}."
        )
        provenance_record = {
            "ancestors": ancestors,
            "authors": ["schlund_manuel", "winterstein_franziska"],
            "caption": caption,
            "plot_types": ["line"],
            "long_names": [var_attrs["long_name"]],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def compute(self):
        """Plot preprocessed data."""
        input_data = self.input_data_dataset
        base_datasets = {
            label: dataset["dataset_data"][0]
            for label, dataset in input_data.items()
        }

        # At the moment regions only apply to TROP and STRAT
        for region in self.cfg["regions"]:
            logger.info("Plotting lifetime for region %s", region)
            self.create_timeseries_plot(region, input_data, base_datasets)


def main():
    """Run diagnostic."""
    with run_diagnostic() as config:
        CH4Lifetime(config).compute()


if __name__ == "__main__":
    main()
