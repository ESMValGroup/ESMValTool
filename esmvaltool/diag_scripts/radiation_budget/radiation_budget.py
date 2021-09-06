"""PLACEHOLDER modified from AA docs: The aim of monitoring the energy budget
is to understand the (im)balance of energy flux between the atmosphere and the
surface of a model due to its link with the hydrological cycle and climate
change."""

# To run the doctests:
# % cd ESMValTool/esmvaltool/
# % python -m doctest diag_scripts/radiation_budget/radiation_budget.py
import logging
import os

import iris
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

CWD = os.path.abspath(os.path.dirname(__file__))
STEPHENS_FILENAME = "Stephens_et_al_2012_obs_Energy_Budget.txt"
DEMORY_FILENAME = "Demory_et_al_2014_obs_Energy_Budget.txt"


def load_data(filenames):
    """Return the loaded cubes.

    Parameters
    ----------
    filenames : list of strings
        The filenames to load.

    Returns
    -------
    :class:`iris.cube.Cube`
        The loaded cubes.
    """
    cubes = iris.load(filenames)
    return cubes


def get_filenames(group):
    """Return all the filenames for the group.

    Parameters
    ----------
    group : list(dict)
        The grouped metadata describing preprocessed data.

    Returns
    -------
    list
        All the filenames for the group.
    """
    filenames = [item["filename"] for item in group]
    return filenames


def derive_additional_variables(cubes):
    """Return input ``cubes`` with the additional cubes.

    ``cubes`` must contain the variables specified in the recipe.

    The additional cubes derived from the cubes in ``cubes`` are as follows:

    * total_sw_cloud_forcing
    * upward_sw_reflected_surface
    * sw_reflected_clouds
    * sw_absorbed_atm
    * upward_lw_emitted_surface
    * total_lw_cloud_forcing
    * net_surface_radiation
    * radiation_adsorbed_surface
    * radiation_net_toa

    Parameters
    ----------
    cubes : :class:`iris.cube.CubeList`
        The cubes corresponding with the variables in the recipe.

    Returns
    -------
    :class:`iris.cube.CubeList`
        The input ``cubes`` with the additional cubes.
    """
    def _constraint(var_name):
        return iris.Constraint(
            cube_func=lambda cube: cube.var_name == var_name)

    rss = cubes.extract_cube(_constraint("rss"))
    rsdt = cubes.extract_cube(_constraint("rsdt"))
    rsut = cubes.extract_cube(_constraint("rsut"))
    rsutcs = cubes.extract_cube(_constraint("rsutcs"))
    rsds = cubes.extract_cube(_constraint("rsds"))
    rls = cubes.extract_cube(_constraint("rls"))
    rlut = cubes.extract_cube(_constraint("rlut"))
    rlutcs = cubes.extract_cube(_constraint("rlutcs"))
    rlds = cubes.extract_cube(_constraint("rlds"))
    hfss = cubes.extract_cube(_constraint("hfss"))
    hfls = cubes.extract_cube(_constraint("hfls"))

    total_sw_cloud_forcing = rsut - rsutcs
    upward_sw_reflected_surface = rsds - rss
    sw_reflected_clouds = rsut - upward_sw_reflected_surface
    sw_absorbed_atm = rsdt - sw_reflected_clouds - rsds
    upward_lw_emitted_surface = rlds - rls
    total_lw_cloud_forcing = rlutcs - rlut
    net_surface_radiation = rss + rls
    radiation_adsorbed_surface = rss + rls - hfss - hfls
    radiation_net_toa = rsdt - rsut - rlut

    total_sw_cloud_forcing.standard_name = ""
    total_sw_cloud_forcing.long_name = "total_sw_cloud_forcing"

    upward_sw_reflected_surface.standard_name = ""
    upward_sw_reflected_surface.long_name = "upward_sw_reflected_surface"

    sw_reflected_clouds.standard_name = ""
    sw_reflected_clouds.long_name = "sw_reflected_clouds"

    sw_absorbed_atm.standard_name = ""
    sw_absorbed_atm.long_name = "sw_absorbed_atm"

    upward_lw_emitted_surface.standard_name = ""
    upward_lw_emitted_surface.long_name = "upward_lw_emitted_surface"

    total_lw_cloud_forcing.standard_name = ""
    total_lw_cloud_forcing.long_name = "total_lw_cloud_forcing"

    net_surface_radiation.standard_name = ""
    net_surface_radiation.long_name = "net_surface_radiation"

    radiation_adsorbed_surface.standard_name = ""
    radiation_adsorbed_surface.long_name = "radiation_adsorbed_surface"

    radiation_net_toa.standard_name = ""
    radiation_net_toa.long_name = "radiation_net_toa"

    additional_cubes = [
        total_sw_cloud_forcing,
        upward_sw_reflected_surface,
        sw_reflected_clouds,
        sw_absorbed_atm,
        upward_lw_emitted_surface,
        total_lw_cloud_forcing,
        net_surface_radiation,
        radiation_adsorbed_surface,
        radiation_net_toa,
    ]

    cubes.extend(additional_cubes)
    return cubes


def validate_variable_data(variable_data, name, unit):
    """Return the variable from ``variable_data`` that has the same name and
    units as provided by ``name`` and ``unit``.

    If ``name`` doesn't exist in ``variable_data``, the returned variable will
    have a name and unit equal to ``name`` and ``unit`` and data equal to
    'NaN'.

    Parameters
    ----------
    variable_data : list of dictionaries
        The data to check where each dictionary corresponds
        to a variable and the key of the dictionary is the
        metadata attribute name.
    name : string
        The name of the variable to validate.
    unit : string
        The unit of the variable to validate.

    Raises
    ------
    RuntimeError
        If multiple ``name`` exist in ``variable_data``.
        If ``unit`` doesn't match the unit in ``variable_data``.

    Returns
    -------
    dictionary
        The validated variable.

    Examples
    --------
    >>> var1 = {"name": "sw_reflected_clouds", "unit": "W m-2", "data": 79.0}
    >>> var2 = {"name": "toa_outgoing_longwave_flux", "unit": "W m-2",
    ...         "data": 239.0}
    >>> variable_data = [var1, var2]
    >>> name = "sw_reflected_clouds"
    >>> unit = "W m-2"
    >>> validated_variable = validate_variable_data(variable_data, name, unit)
    >>> assert validated_variable == var1
    """
    items = [item for item in variable_data if item["name"] == name]

    if not items:
        variable = {"name": name, "unit": unit, "data": np.nan}

    if len(items) == 1:
        variable = items[0]

    if len(items) > 1:
        raise RuntimeError(f"Multiple '{name}' exist in '{items}'.")

    if variable["unit"] != unit:
        raise RuntimeError(
            f"Unit {unit} does not match the unit {variable['unit']} "
            f"in {variable} for {name}.")

    return variable


def order_data(cubes, obs_names, obs_units):
    """Return the data from the cubes in the order defined by ``obs_names``.

    The units from the cubes are checked against ``obs_units``.

    Parameters
    ----------
    cubes : :class:`iris.cube.CubeList`
        The cubes in a random order.
    obs_names : list
        The ordered names from the observation files.
    obs_units : list
        The ordered units from the observation files.

    Returns
    -------
    list
        The ordered data from the model cubes.
    """
    variable_data = []
    for cube in cubes:
        variable = {}
        variable["name"] = cube.name()
        variable["unit"] = cube.units
        if np.ma.isMaskedArray(cube.data):
            variable["data"] = cube.data.data
        else:
            variable["data"] = cube.data
        variable_data.append(variable)

    ordered_model_data = []
    for index, obs_name in enumerate(obs_names):
        obs_unit = obs_units[index]
        validated_variable = validate_variable_data(variable_data, obs_name,
                                                    obs_unit)
        ordered_model_data.append(validated_variable["data"])

    return ordered_model_data


def read_text_file(filepath):
    """Return contents of text file.

    Parameters
    ----------
    filepath : string
        The full path to the text file.

    Returns
    -------
    list of dictionaries
        The contents of the text file where each dictionary corresponds
        to a line in the file and the key of the dictionary is the name
        of the column.
    """
    contents = []
    with open(filepath, "r") as file_handle:
        lines = file_handle.readlines()
    for line in lines:
        line_dict = {}
        items = line.split()
        line_dict["name"] = items[0]
        line_dict["unit"] = " ".join(items[1:3]).strip(",")
        line_dict["data"] = float(items[3])  # float('NaN') == np.nan
        try:
            line_dict["error"] = float(items[4])
        except IndexError:
            pass
        contents.append(line_dict)
    return contents


def load_obs_data():
    """Return the names, units, data and error from the Stephens and Demory
    observation files.

    The observation files should exist in the same directory as this
    module.

    Returns
    -------
    tuple of lists
        The names, units, stephens data, stephens error and demory data
        from the observation files.
    """
    # Stephens data contains name, units part 1, units part 2, data, error.
    stephens_filepath = os.path.join(CWD, STEPHENS_FILENAME)
    stephens_contents = read_text_file(stephens_filepath)

    # Demory data contains name, units part 1, units part 2, data.
    demory_filepath = os.path.join(CWD, DEMORY_FILENAME)
    demory_contents = read_text_file(demory_filepath)

    # Arbitrarily use the order as defined in the Stephens filename.
    names = []
    units = []
    stephens_data = []
    stephens_error = []
    demory_data = []
    for line in stephens_contents:
        name = line["name"]
        unit = line["unit"]
        names.append(name)
        units.append(unit)
        stephens_data.append(line["data"])
        stephens_error.append(line["error"])

        demory_line = validate_variable_data(demory_contents, name, unit)
        demory_data.append(demory_line["data"])

    return names, units, stephens_data, stephens_error, demory_data


def plot_data(obs_names, obs_units, stephens_data, stephens_error, demory_data,
              ceres_data, model_data, plot_filename):
    """Produce and save the radiation budget comparison plot.

    Parameters
    ----------
    obs_names : list
        The names of variables included in the observation data.
    obs_units : list
        The units of variables included in the observation data.
    stephens_data : list
        Stephens observation data values.
    stephens_error : list
        Stephens observation data error values.
    demory_data : list
        Demory observation data values.
    ceres_data : list
        Ceres observation data values.
    model_data : list
        Data values from the model for which this comparison plot is being
        generated.
    plot_filename : string
        The name of the image file to be saved.
    """
    model_minus_stephens = np.array(model_data) - np.array(stephens_data)
    model_minus_demory = np.array(model_data) - np.array(demory_data)
    model_minus_ceres = np.array(model_data) - np.array(ceres_data)

    n_groups = 20
    fig = plt.Figure(figsize=(12, 8))
    # plt.FigureCanvas(fig)
    axis = fig.add_subplot(111)
    index = np.arange(n_groups) * 2.0

    bar_width, opacity = 0.5, 0.4
    axis.set_ylim(-20, 20)
    axis.bar(
        index + 0.2,
        model_minus_stephens,
        bar_width,
        alpha=opacity + 0.2,
        color="cornflowerblue",
        label="MODEL" + "(" + "season" + ")" + " - Stephens(2012)",
        yerr=stephens_error,
    )
    axis.bar(
        index + 0.2 + bar_width,
        model_minus_ceres,
        bar_width,
        alpha=opacity + 0.2,
        color="orange",
        label="MODEL" + "(" + "season" + ")" + " - CERES" + "(" + "season" +
        ")",
    )
    axis.bar(
        index + 0.2 + bar_width * 2,
        model_minus_demory,
        bar_width,
        alpha=opacity - 0.2,
        color="black",
        label="MODEL" + "(" + "season" + ")" + " - Demory(2014)",
    )
    axis.spines["bottom"].set_position(("data", 0))
    axis.spines["top"].set_position(("data", 0))

    axis.set_xticks(index + bar_width + 0.5)
    axis.set_xticklabels(obs_names, ha="center", rotation=90, fontsize=10)
    axis.legend(frameon=False, fontsize=10, loc="upper left")

    # TODO: use static filename, and put additional information into figure
    # title
    fig.savefig(plot_filename)
    fig.clear()


def main(config):
    """Radiation budget for HadGEM3 vs UKESM1.

    Parameters
    ----------
    config : dict
        The ESMValTool configuration.
    """
    logger = logging.getLogger(__name__)

    input_data = config["input_data"]
    datasets = group_metadata(input_data.values(), "dataset")

    (
        obs_names,
        obs_units,
        stephens_data,
        stephens_error,
        demory_data,
    ) = load_obs_data()

    ceres_group = datasets.pop("CERES-EBAF")
    ceres_filenames = get_filenames(ceres_group)
    raw_ceres_data = load_data(ceres_filenames)
    ceres_data = order_data(raw_ceres_data, obs_names, obs_units)

    for model_dataset, group in datasets.items():
        # 'model_dataset' is the name of the model dataset.
        # 'group' is a list of dictionaries containing metadata.
        logger.info("Processing data for %s", model_dataset)
        filenames = get_filenames(group)
        unordered_model_data = load_data(filenames)
        all_model_data = derive_additional_variables(unordered_model_data)
        model_data = order_data(all_model_data, obs_names, obs_units)
        plot_filename = get_plot_filename(model_dataset, config)
        plot_data(
            obs_names,
            obs_units,
            stephens_data,
            stephens_error,
            demory_data,
            ceres_data,
            model_data,
            plot_filename,
        )


if __name__ == "__main__":
    with run_diagnostic() as CONFIG:
        main(CONFIG)
