"""Plot the global radiation budget."""

# To run the doctests:
# % cd ESMValTool/esmvaltool/
# % python -m doctest diag_scripts/radiation_budget/radiation_budget.py
import logging
import os

import iris
import matplotlib.pyplot as plt
import numpy as np
from common import get_filenames, load_data

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_figure,
)

CWD = os.path.abspath(os.path.dirname(__file__))
STEPHENS_FILENAME = "Stephens_et_al_2012_obs_Energy_Budget.txt"
DEMORY_FILENAME = "Demory_et_al_2014_obs_Energy_Budget.txt"


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

    # Derivations for the following two cloud_forcing variables are
    # performed this way so that they match with the observational data
    # (all positive), the convention used is to treat SW as positive
    # downward and LW as positive upward.
    total_sw_cloud_forcing = rsut - rsutcs
    total_lw_cloud_forcing = rlutcs - rlut
    upward_sw_reflected_surface = rsds - rss
    sw_reflected_clouds = rsut - upward_sw_reflected_surface
    sw_absorbed_atm = rsdt - sw_reflected_clouds - rsds
    upward_lw_emitted_surface = rlds - rls
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
    KeyError
        If multiple ``name`` exist in ``variable_data``.
    ValueError
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
        raise KeyError(f"Multiple '{name}' exist in '{items}'.")

    if variable["unit"] != unit:
        raise ValueError(
            f"Unit {unit} does not match the unit {variable['unit']} "
            f"in {variable} for {name}.")

    return variable


def order_data(cubes, obs_names, obs_unit):
    """Return the data from the cubes in the order defined by ``obs_names``.

    The units from the cubes are checked against ``obs_units``.

    Parameters
    ----------
    cubes : :class:`iris.cube.CubeList`
        The cubes in a random order.
    obs_names : list
        The ordered names from the observation files.
    obs_unit : string
        The unit of the observation variables.

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
    for obs_name in obs_names:
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

    if len(set(units)) == 1:
        unit = units[0]
    else:
        raise RuntimeError("Not all observations have the same unit.")

    return names, unit, stephens_data, stephens_error, demory_data


def plot_data(
    model_dataset,
    model_data,
    model_period,
    obs_names,
    obs_unit,
    stephens_data,
    stephens_error,
    demory_data,
    ceres_dataset,
    ceres_data,
    ceres_period,
):
    """Produce and save the radiation budget comparison plot.

    Parameters
    ----------
    model_dataset : string
        The name of the model.
    model_data : list
        Data values from the model for which this comparison plot is being
        generated.
    model_period : string
        The start and end years of the model dataset.
    obs_names : list
        The names of variables included in the observation data.
    obs_unit : list
        The unit of variables included in the observation data.
    stephens_data : list
        Stephens observation data values.
    stephens_error : list
        Stephens observation data error values.
    demory_data : list
        Demory observation data values.
    ceres_dataset : string
        The name of the CERES observation data.
    ceres_data : list
        CERES observation data values.
    ceres_period : string
        The start and end years of the CERES observation data.

    Returns
    -------
    :class:`matplotlib.figure.Figure`
        The figure containing the plot.
    """
    model_minus_stephens = np.array(model_data) - np.array(stephens_data)
    model_minus_demory = np.array(model_data) - np.array(demory_data)
    model_minus_ceres = np.array(model_data) - np.array(ceres_data)

    figure, axes = plt.subplots(figsize=(12, 8))
    title = f"Radiation budget for {model_dataset}"
    y_label = f"Difference between model output and observations [{obs_unit}]"
    y_lim = (-20, 20)
    axes.set(title=title, ylabel=y_label, ylim=y_lim)

    num_x_ticks = len(obs_names)
    x_ticks = np.arange(0, num_x_ticks * 2, 2)

    bar_width = 0.5
    opacity = 0.6
    axes.bar(
        x_ticks + 0.2,
        model_minus_stephens,
        bar_width,
        alpha=opacity,
        color="cornflowerblue",
        label=f"{model_dataset} ({model_period}) - Stephens et al. (2012)",
        yerr=stephens_error,
    )
    axes.bar(
        x_ticks + 0.2 + bar_width,
        model_minus_ceres,
        bar_width,
        alpha=opacity,
        color="orange",
        label=(f"{model_dataset} ({model_period}) - {ceres_dataset} "
               f"({ceres_period})"),
    )
    axes.bar(
        x_ticks + 0.2 + bar_width * 2,
        model_minus_demory,
        bar_width,
        alpha=opacity,
        color="darkgrey",
        label=f"{model_dataset} ({model_period}) - Demory et al. (2014)",
    )
    axes.spines["bottom"].set_position(("data", 0))
    axes.spines["top"].set_position(("data", 0))

    axes.set_xticks(x_ticks + bar_width + 0.5)
    axes.set_xticklabels(obs_names, ha="center", rotation=90, fontsize=10)

    axes.legend(frameon=False, fontsize=10, loc="upper left")

    return figure


def get_provenance_record(filenames):
    """Return a provenance record describing the plot.

    Parameters
    ----------
    filenames : list of strings
        The filenames containing the data used to create the plot.

    Returns
    -------
    dictionary
        The provenance record describing the plot.
    """
    record = {
        'ancestors': filenames,
    }
    return record


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
        obs_unit,
        stephens_data,
        stephens_error,
        demory_data,
    ) = load_obs_data()

    ceres_dataset = "CERES-EBAF"
    ceres_group = datasets.pop(ceres_dataset)
    ceres_filenames = get_filenames(ceres_group)
    raw_ceres_data = load_data(ceres_filenames)
    ceres_data = order_data(raw_ceres_data, obs_names, obs_unit)
    ceres_period = (f"{ceres_group[0]['start_year']} - "
                    f"{ceres_group[0]['end_year']}")

    for model_dataset, group in datasets.items():
        # 'model_dataset' is the name of the model dataset.
        # 'group' is a list of dictionaries containing metadata.
        logger.info("Processing data for %s", model_dataset)
        filenames = get_filenames(group)
        unordered_model_data = load_data(filenames)
        all_model_data = derive_additional_variables(unordered_model_data)
        model_data = order_data(all_model_data, obs_names, obs_unit)
        model_period = f"{group[0]['start_year']} - {group[0]['end_year']}"
        figure = plot_data(
            model_dataset,
            model_data,
            model_period,
            obs_names,
            obs_unit,
            stephens_data,
            stephens_error,
            demory_data,
            ceres_dataset,
            ceres_data,
            ceres_period,
        )
        provenance_record = get_provenance_record(filenames)
        save_figure(model_dataset,
                    provenance_record,
                    config,
                    figure,
                    close=True)


if __name__ == "__main__":
    with run_diagnostic() as CONFIG:
        main(CONFIG)
