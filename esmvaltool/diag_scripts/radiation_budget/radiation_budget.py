"""Plot the global radiation budget."""

# To run the doctests:
# % cd ESMValTool/esmvaltool/
# % python -m doctest diag_scripts/radiation_budget/radiation_budget.py
import logging
import os

import iris
import matplotlib.pyplot as plt
import numpy as np
import yaml
from esmvalcore.iris_helpers import var_name_constraint

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_figure,
)

CWD = os.path.abspath(os.path.dirname(__file__))
STEPHENS_FILENAME = "Stephens_et_al_2012_obs_Energy_Budget.yml"
DEMORY_FILENAME = "Demory_et_al_2014_obs_Energy_Budget.yml"


def derive_additional_variables(cubes, is_model=True):
    """Return input ``cubes`` with the additional cubes.

    ``cubes`` must contain the variables specified in the recipe.

    When ``cubes`` contains observational data, the additional cubes
    derived are as follows:

    * total_sw_cloud_forcing
    * total_lw_cloud_forcing

    When ``cubes`` contains model data, the following cubes are derived in
    addition to those above:

    * upward_sw_reflected_surface
    * sw_reflected_clouds
    * sw_absorbed_atm
    * upward_lw_emitted_surface
    * net_surface_radiation
    * radiation_adsorbed_surface
    * radiation_net_toa

    Parameters
    ----------
    cubes : :class:`iris.cube.CubeList`
        The cubes corresponding with the variables in the recipe.
    is_model : bool
        Flag that dictates whether additional model variables are derived,
        or just those common to models and observations.

    Returns
    -------
    :class:`iris.cube.CubeList`
        The input ``cubes`` with the additional cubes.
    """
    common_variables = ['rsut', 'rsutcs', 'rlut', 'rlutcs']
    model_variables = ['rss', 'rsdt', 'rsds', 'rls', 'rlds', 'hfss', 'hfls']
    var = {}

    for variable in common_variables:
        var[variable] = cubes.extract_cube(var_name_constraint(variable))

    if is_model:
        for variable in model_variables:
            var[variable] = cubes.extract_cube(var_name_constraint(variable))

    # Derivations for the following two cloud_forcing variables are
    # performed this way so that they match with the observational data
    # (all positive), the convention used is to treat SW as positive
    # downward and LW as positive upward.
    total_sw_cloud_forcing = var['rsut'] - var['rsutcs']
    total_lw_cloud_forcing = var['rlutcs'] - var['rlut']

    _update_long_name(total_sw_cloud_forcing, "total_sw_cloud_forcing")
    _update_long_name(total_lw_cloud_forcing, "total_lw_cloud_forcing")
    additional_cubes = [total_sw_cloud_forcing, total_lw_cloud_forcing]

    if is_model:
        upward_sw_reflected_surface = var['rsds'] - var['rss']
        sw_reflected_clouds = var['rsut'] - upward_sw_reflected_surface
        sw_absorbed_atm = var['rsdt'] - sw_reflected_clouds - var['rsds']
        upward_lw_emitted_surface = var['rlds'] - var['rls']
        net_surface_radiation = var['rss'] + var['rls']
        radiation_adsorbed_surface = (var['rss'] + var['rls'] - var['hfss'] -
                                      var['hfls'])
        radiation_net_toa = var['rsdt'] - var['rsut'] - var['rlut']

        _update_long_name(upward_sw_reflected_surface,
                          'upward_sw_reflected_surface')
        _update_long_name(sw_reflected_clouds, 'sw_reflected_clouds')
        _update_long_name(sw_absorbed_atm, 'sw_absorbed_atm')
        _update_long_name(upward_lw_emitted_surface,
                          'upward_lw_emitted_surface')
        _update_long_name(net_surface_radiation, 'net_surface_radiation')
        _update_long_name(radiation_adsorbed_surface,
                          'radiation_adsorbed_surface')
        _update_long_name(radiation_net_toa, 'radiation_net_toa')

        additional_cubes.extend([
            upward_sw_reflected_surface,
            sw_reflected_clouds,
            sw_absorbed_atm,
            upward_lw_emitted_surface,
            net_surface_radiation,
            radiation_adsorbed_surface,
            radiation_net_toa,
        ])

    cubes.extend(additional_cubes)
    return cubes


def _update_long_name(cube, long_name):
    cube.standard_name = ""
    cube.long_name = long_name


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
        If ``unit`` does not match the unit in ``variable_data``.

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


def read_yaml_file(filepath):
    """Return contents of a yaml file.

    Parameters
    ----------
    filepath : string
        The full path to the yaml file.

    Returns
    -------
    list of dictionaries
        The contents of the yaml file where each dictionary corresponds
        to a line in the file and the key of the dictionary is the name
        of the column.
    """
    with open(filepath, "r") as stream:
        contents = yaml.safe_load(stream)
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
    # Stephens data contains name, units, data, error.
    stephens_filepath = os.path.join(CWD, STEPHENS_FILENAME)
    stephens_contents = read_yaml_file(stephens_filepath)

    # Demory data contains name, units, data.
    demory_filepath = os.path.join(CWD, DEMORY_FILENAME)
    demory_contents = read_yaml_file(demory_filepath)

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
    """Radiation budget comparison for models defined in the radiation_budget
    recipe file.

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
    ceres_filenames = [item["filename"] for item in ceres_group]
    raw_ceres_data = iris.load(ceres_filenames)
    all_ceres_data = derive_additional_variables(raw_ceres_data,
                                                 is_model=False)
    ceres_data = order_data(all_ceres_data, obs_names, obs_unit)
    ceres_period = (f"{ceres_group[0]['start_year']} - "
                    f"{ceres_group[0]['end_year']}")

    for model_dataset, group in datasets.items():
        # 'model_dataset' is the name of the model dataset.
        # 'group' is a list of dictionaries containing metadata.
        logger.info("Processing data for %s", model_dataset)
        filenames = [item["filename"] for item in group]
        unordered_model_data = iris.load(filenames)
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
