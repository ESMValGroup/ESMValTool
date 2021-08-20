import logging
import os

import iris

from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic

CWD = os.path.abspath(os.path.dirname(__file__))
STEPHENS_FILENAME = "Stephens_et_al_2012_obs_Energy_Budget.txt"
DEMORY_FILENAME = "Demory_et_al_2014_obs_Energy_Budget.txt"


def load_model_data(filenames):
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
        If a single ``name`` doesn't exist in ``variable_data``.
        If ``unit`` doesn't match the unit in ``variable_data``.

    Returns
    -------
    dictionary
        The validated variable.
    """
    items = [item for item in variable_data if item["name"] == name]

    if len(items) != 1:
        raise RuntimeError(f"A single '{name}' does not exist in '{items}'.")

    variable = items[0]
    variable_unit = variable["unit"]

    if variable_unit != unit:
        raise RuntimeError(
            f"Unit {unit} does not match the unit {variable_unit} "
            f"in {variable} for {name}")

    return variable


def organise_variables(cubes, obs_names, obs_units):
    """Return variables in the order defined by ``obs_names``.

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
        The organised cubes.
    """
    logger = logging.getLogger(__name__)

    named_cubes = {cube.name(): cube for cube in cubes}

    organised_cubes = [named_cubes[name] for name in obs_names]

    for cube in organised_cubes:
        logger.info(f"{cube.name()}, {cube.units}, {cube.data}")

    return organised_cubes


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
        line_dict["data"] = items[3]
        try:
            line_dict["error"] = items[4]
        except IndexError:
            pass
        contents.append(line_dict)
    return contents


def load_ceres_data():
    pass


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
    logger = logging.getLogger(__name__)

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

    for index, name in enumerate(names):
        logger.info(f"{name}, {units[index]}, {stephens_data[index]}, "
                    f"{stephens_error[index]}")
        logger.info(f"{name}, {units[index]}, {demory_data[index]}")

    return names, units, stephens_data, stephens_error, demory_data


def plot_data(cubes):
    pass


def main(config):
    """Radiation budget for HadGEM3 vs UKESM1.

    Parameters
    ----------
    config : dict
        The ESMValTool configuration.
    """
    logger = logging.getLogger(__name__)

    input_data = config["input_data"]
    model_datasets = group_metadata(input_data.values(), "dataset")

    obs_names, obs_units, stephens_data, stephens_error, demory_data = (
        load_obs_data())

    for model_dataset, group in model_datasets.items():
        # 'model_dataset' is the name of the model dataset.
        # 'group' is a list of dictionaries containing metadata.
        logger.info(f"Processing data for {model_dataset}")
        filenames = get_filenames(group)
        model_data = load_model_data(filenames)
        all_model_data = derive_additional_variables(model_data)
        organised_model_data = organise_variables(all_model_data, obs_names,
                                                  obs_units)
        plot_data(organised_model_data)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
