"""Base class for lifetime diagnostics."""

import logging
from copy import deepcopy

import dask.array as da
import iris
from iris.util import broadcast_to_shape
from scipy.constants import N_A, R, g

# set standard molar masses
M_AIR = 28.970  # [g_air/mol_air]
M_H2O = 18.02  # [g_air/mol_air]

logger = logging.getLogger(__name__)


def create_press(var):
    """Create a pressure variable."""
    resolver = iris.common.resolve.Resolve(var, var)
    if var.coord("air_pressure").shape == var.shape:
        press = resolver.cube(var.coord("air_pressure").lazy_points())
    else:
        press = resolver.cube(
            broadcast_to_shape(
                var.coord("air_pressure").lazy_points(),
                var.shape,
                var.coord_dims("air_pressure"),
            ),
        )
    press.long_name = "air_pressure"
    press.var_name = "air_pressure"
    press.units = "Pa"

    return press


def calculate_gridmassdry(press, hus, z_coord):
    """Calculate the dry gridmass according to pressure and humidity.

    Used constants:
    g  standard acceleration of gravity from scipy
    """
    # calculate minimum in cubes
    if isinstance(press, iris.cube.Cube):
        pmin = da.nanmin(press.core_data())
    else:
        pmin = da.nanmin(press)

    # surface air pressure as cube
    surface = press.coord(z_coord)[0].points
    pmax = press.extract(
        iris.Constraint(coord_values=dict.fromkeys([z_coord], surface)),
    ).copy()
    if "surface_air_pressure" in press.coords():
        pmax.data = press.coord("surface_air_pressure").lazy_points()
    else:
        pmax.data = da.full_like(
            press.extract(
                iris.Constraint(
                    coord_values=dict.fromkeys([z_coord], surface),
                ),
            ),
            101525.0,
        )

    # grid area -> m2
    lat_lon_dims = sorted(
        tuple(set(hus.coord_dims("latitude") + hus.coord_dims("longitude"))),
    )
    lat_lon_slice = next(hus.slices(["latitude", "longitude"], ordered=False))
    area_2d = iris.analysis.cartography.area_weights(lat_lon_slice)
    area = deepcopy(press)
    area.data = broadcast_to_shape(
        da.array(area_2d),
        hus.shape,
        lat_lon_dims,
    )

    # delta pressure levels
    delta_p = dpres_plevel(press, pmin, pmax, z_coord=z_coord)

    # gridmass of dry air
    gridmassdry = deepcopy(press)

    # Following steps are executed:
    # - calculation of delta pressure to
    #     grid mass per square meter > kg m-2
    # - multiplied with area per grid -> mass per grid
    # - remove humidity for dry weight
    #
    gridmassdry = area * (delta_p / g) * (1.0 - hus)

    # attach metadata
    gridmassdry.long_name = "gridmassdry"
    gridmassdry.var_name = "gridmassdry"
    gridmassdry.units = "kg"

    return gridmassdry


def calculate_rho(variables):
    """Calculate number density (rho).

    Calculate number density with respect to the given input
    variables.
    """
    logger.info("Calculate number density (rho)")

    # model levels
    if "grmassdry" in variables and "grvol" in variables:
        rho = _number_density_dryair_by_grid(
            variables["grmassdry"],
            variables["grvol"],
        )
    # pressure levels
    elif "ta" in variables and "hus" in variables:
        rho = _number_density_dryair_by_press(
            variables["ta"],
            variables["hus"],
        )
    else:
        raise NotImplementedError(
            "The necessary variables"
            " to calculate number"
            " density of dry air"
            " are not provided.\n"
            "Provide either:\n"
            " - grmassdry and grvol\n"
            " or\n"
            " - ta and hus",
        )

    return rho


def _number_density_dryair_by_press(temp, hus, press=None):
    """
    Calculate number density of dry air.

    Used to convert from mol / mol_dry into molec / cm3
    by using present temperature and humidity.

    ##qqq
    Should there be an option to provide the simulated pressure field
    rather than the derived pressure from interpolation?
    ###

    Used constants:
    - N_A    Avogrado constant from scipy
    - R      gas constant from scipy
    - M_AIR   Molarmass of Air
    - M_H2O   Molarmass of watervapor
    """
    logger.info("Calculate number density of dry air by pressure")

    if not press:
        logger.info("Pressure not given")
        press = create_press(temp)

    rho = (
        N_A
        / 10.0**6
        * press
        / (R * temp)
        * (1.0 - hus)
        / (1.0 + hus * (M_AIR / M_H2O - 1.0))
    )

    # correct metadata [ 1 / cm^3 ]
    rho.var_name = "rho"
    rho.units = "cm-3"

    return rho


def _number_density_dryair_by_grid(grmassdry, grvol):
    """
    Calculate number density of dry air.

    Used to convert from mol / mol_dry into molec / cm3
    by using present gridmass of dry air and gridvolume.

    Since gridvolume might not be appropriate after interpolation
    to pressure coordinates, this version should only be used
    on model levels.

    Used constants:
    - N_A    Avogrado constant from scipy
    - M_AIR   Molarmass of Air
    """
    logger.info("Calculate number density of dry air by grid information")
    rho = (grmassdry / grvol) * (N_A / M_AIR) * 10 ** (-3)  # [ 1 / cm^3 ]
    # correct metadata
    rho.var_name = "rho"
    rho.units = "cm-3"

    return rho


def dpres_plevel(plev, pmin, pmax, z_coord="air_pressure"):
    """Calculate pressure levels."""
    cube_coords = [coord.name() for coord in plev.dim_coords]

    if [z_coord, "latitude", "longitude"] == cube_coords:
        dplev = dpres_plevel_3d(plev, pmin, pmax, z_coord)
    elif ["time", z_coord, "latitude", "longitude"] == cube_coords:
        dplev = dpres_plevel_4d(plev, pmin, pmax, z_coord)
    else:
        raise NotImplementedError(
            "Pressure level calculation is not implemented for the present "
            "coordinates",
        )

    return dplev


def dpres_plevel_3d(plev, pmin, pmax, z_coord="air_pressure"):
    """Calculate delta pressure levels.

    The delta pressure levels are based
    on the given pressure as a
    four dimensional cube.
    """
    axis = plev.coord_dims(z_coord)

    # shifted vectors:
    # - p_m1: shifted by one index down
    #         (the value of the last index is the former first)
    # - p_p1: shifted by one index up
    #         (the value of the first index is the former last)
    # both vector values are divided by two
    p_p1 = da.roll(plev.lazy_data(), 1, axis=axis) / 2.0
    p_m1 = da.roll(plev.lazy_data(), -1, axis=axis) / 2.0

    # modify the first entry in p_p1
    # note: p_p1[1, :, :] .eq. plev[0, :, :] / 2.
    p_p1[0, :, :] = pmax.lazy_data() - p_p1[1, :, :]

    # and the last entry in p_m1
    # note: p_m1[-2, :, :] .eq. plev[-1, :, :] / 2.
    p_m1[-1, :, :] = pmin - p_m1[-2, :, :]

    # calculate difference
    dplev = deepcopy(plev)
    dplev.data = p_p1 - p_m1

    return dplev


def dpres_plevel_4d(plev, pmin, pmax, z_coord="air_pressure"):
    """Calculate delta pressure levels.

    The delta pressure levels are based
    on the given pressure as a
    four dimensional cube.
    """
    axis = plev.coord_dims(z_coord)

    # shifted vectors:
    # - p_m1: shifted by one index down
    #         (the value of the last index is the former first)
    # - p_p1: shifted by one index up
    #         (the value of the first index is the former last)
    # both vector values are divided by two
    p_p1 = da.roll(plev.lazy_data(), 1, axis=axis) / 2.0
    p_m1 = da.roll(plev.lazy_data(), -1, axis=axis) / 2.0

    # modify the first entry in p_p1
    # note: p_p1[:, 1, :, :] .eq. plev[:, 0, :, :] / 2.
    p_p1[:, 0, :, :] = pmax.lazy_data() - p_p1[:, 1, :, :]

    # and the last entry in p_m1
    # note: p_m1[:, -2, :, :] .eq. plev[:, -1, :, :] / 2.
    p_m1[:, -1, :, :] = pmin - p_m1[:, -2, :, :]

    # calculate difference
    dplev = deepcopy(plev)
    dplev.data = p_p1 - p_m1

    return dplev


def calculate_lifetime(dataset, plot_type, region=None):
    """Calculate the lifetime for the given plot_type and region."""
    # Extract region from weights and reaction
    if region is not None:
        reaction = extract_region(dataset, region, case="reaction")
        weight = extract_region(dataset, region, case="weight")
    else:
        reaction = dataset["reaction"]
        weight = dataset["weight"]

    # Calculate nominator and denominator and sum of nominator and denominator
    # via plot_type dimensions
    nominator = sum_up_to_plot_dimensions(weight, plot_type)
    denominator = sum_up_to_plot_dimensions(weight * reaction, plot_type)

    return nominator / denominator


def extract_region(dataset, region, case="reaction"):
    """Return cube with everything outside region set to zero.

    Current aware regions:
    - TROP: troposphere (excl. tropopause)
    - STRA: stratosphere (incl. tropopause)
    """
    var = dataset[case].copy()
    z_coord = dataset["use_z_coord"]

    # mask regions outside
    z_4d = broadcast_to_shape(
        var.coord(z_coord).lazy_points(),
        var.shape,
        var.coord_dims(z_coord),
    )

    tp_4d = broadcast_to_shape(
        dataset["tropopause"].core_data(),
        var.shape,
        tuple(
            var.coord_dims(item)[0]
            for item in var.dim_coords
            if item != dataset["z_coord"]
        ),
    )

    if region == "TROP":
        var.data = da.ma.masked_array(var.core_data(), mask=z_4d <= tp_4d)
    elif region == "STRA":
        var.data = da.ma.masked_array(var.core_data(), mask=z_4d > tp_4d)
    else:
        raise NotImplementedError(f"region '{region}' is not supported")

    return var


def climatological_tropopause(cube):
    """Return cube with climatological tropopause pressure."""
    if not cube.coords("latitude", dim_coords=True):
        raise NotImplementedError(
            "The provided cube must have a latitude cooridnate",
        )

    tpp = (
        300.0
        - 215.0
        * (da.cos(da.deg2rad(cube.coord("latitude").lazy_points())) ** 2)
    ) * 100.0

    tp_clim = cube.copy()
    tp_clim.data = broadcast_to_shape(
        tpp,
        cube.shape,
        cube.coord_dims("latitude"),
    )
    tp_clim.standard_name = "tropopause_air_pressure"
    tp_clim.var_name = "tp_clim"
    tp_clim.long_name = "climatological tropopause pressure"
    tp_clim.units = "Pa"

    return tp_clim


def sum_up_to_plot_dimensions(var, plot_type):
    """Return the cube summed over the appropriate dimensions."""
    if plot_type == "timeseries":
        if var.coords("air_pressure", dim_coords=True):
            z_coord = var.coords("air_pressure", dim_coords=True)[0]
        elif var.coords("lev", dim_coords=True):
            z_coord = var.coords("lev", dim_coords=True)[0]
        elif var.coords(
            "atmosphere_hybrid_sigma_pressure_coordinate",
            dim_coords=True,
        ):
            z_coord = var.coords(
                "atmosphere_hybrid_sigma_pressure_coordinate",
                dim_coords=True,
            )[0]
        else:
            msg = f"No valid Z-coord available for data\n{var}"
            raise ValueError(msg)
        cube = var.collapsed(
            ["longitude", "latitude", z_coord],
            iris.analysis.SUM,
            weights=None,
        )
    else:
        raise NotImplementedError(
            f"The sum to plot dimensions for plot_type '{plot_type}' is "
            f"currently not implemented",
        )

    return cube


def calculate_reaction_rate(
    temp,
    reaction_type,
    coeff_a,
    coeff_er,
    coeff_b=None,
):
    """Calculate the reaction rate.

    Calculated in Arrhenius form or in a given special form
    depending on the oxidation partner.
    """
    reaction_rate = deepcopy(temp)
    reaction_rate.units = "unknown"

    # special reaction rate
    if coeff_b is not None:
        reaction_rate = coeff_a * iris.analysis.maths.exp(
            coeff_b * iris.analysis.maths.log(reaction_rate)
            - coeff_er / reaction_rate,
        )
    else:
        # standard reaction rate (arrhenius)
        reaction_rate = coeff_a * iris.analysis.maths.exp(
            coeff_er / reaction_rate,
        )

    # set units
    reaction_rate.units = "cm3 s-1"
    reaction_rate.var_name = "reaction_rate"
    reaction_rate.long_name = f"Reaction rate of {reaction_type}"
    return reaction_rate
