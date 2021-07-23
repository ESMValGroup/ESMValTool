import logging

import iris

from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic


def load_cubes(filenames):
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


def organise_variables(cubes):
    """Return variables in an order that matches AutoAssess plots.

    Parameters
    ----------
    cubes : :class:`iris.cube.CubeList`
        The cubes to plot in a random order.

    Returns
    -------
    list
        The organised cubes to plot.
    """
    logger = logging.getLogger(__name__)

    organised_names = [
        "radiation_net_toa",
        "toa_incoming_shortwave_flux",
        "toa_outgoing_shortwave_flux",
        "toa_outgoing_shortwave_flux_assuming_clear_sky",
        "total_sw_cloud_forcing",
        "surface_downwelling_shortwave_flux_in_air",
        "surface_net_downward_shortwave_flux",
        "upward_sw_reflected_surface",
        "sw_reflected_clouds",
        "sw_absorbed_atm",
        "toa_outgoing_longwave_flux",
        "toa_outgoing_longwave_flux_assuming_clear_sky",
        "total_lw_cloud_forcing",
        "surface_downwelling_longwave_flux_in_air",
        "surface_net_downward_longwave_flux",
        "upward_lw_emitted_surface",
        "net_surface_radiation",
        "surface_upward_sensible_heat_flux",
        "surface_upward_latent_heat_flux",
        "radiation_adsorbed_surface",
    ]

    named_cubes = {cube.name(): cube for cube in cubes}

    organised_cubes = [named_cubes[name] for name in organised_names]

    for cube in organised_cubes:
        logger.info(f"{cube.name()}, {cube.units}, {cube.data}")

    return organised_cubes


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

    for dataset, group in datasets.items():
        # 'dataset' is the name of the dataset.
        # 'group' is a list of dictionaries containing metadata.
        logger.info(f"Processing data for {dataset}")
        filenames = get_filenames(group)
        cubes = load_cubes(filenames)
        all_cubes = derive_additional_variables(cubes)
        _ = organise_variables(all_cubes)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
