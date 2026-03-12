"""Part of the ESMValTool trace gas surface concentration diagnostics.

This module contains utility functions commonly used by
trace gas surface concentration assessment routines.

Classes and functions largely adapted from the AOD-AERONET diagnostic
at /diag_scripts/aerosols/aero_utils.py:
- FlaskAnsError
- _add_bounds
- _extract_pt
"""

import logging
from pathlib import Path

import iris
import iris.exceptions
import numpy as np

logger = logging.getLogger(Path(__file__).stem)

TRACE_GASES_FACTOR = {
    "ch4": 1e9,
    "co2": 1e6,
    "n2o": 1e9,
}
TRACE_GASES_UNITS = {
    "ch4": "ppb",
    "co2": "ppm",
    "n2o": "ppb",
}


class FlaskAnsError(Exception):
    """Exception class for errors raised in extract_pt module."""


def _add_bounds(cube):
    """Add bounds to a cube's latitude and longitude coordinates.

    Parameters
    ----------
    cube : Iris cube
        Iris cube with latitude and longitude coordinates.

    Returns
    -------
    cube : Iris cube.
        Iris cube with bounds added to the latitude and longitude coordinates.
    """
    if not cube.coord("latitude").has_bounds():
        cube.coord("latitude").guess_bounds()
    if not cube.coord("longitude").has_bounds():
        cube.coord("longitude").guess_bounds()
    return cube


def _extract_pt(icube, pt_lat, pt_lon, height=None, level=None, nearest=False):
    """Extract given location(s) (3-D) from a cube.

    Method
    ------
    Uses Iris module Analysis.Interpolate to extract values,
    initially based on horizontal coordinates, and then based on
    height, if specified.

    If height ("altitude") is requested, checks if cube heights
    include orography, i.e. HybridHeights have been derived.

    Parameters
    ----------
    icube : iris.cube.Cube
        Input iris cube.
    pt_lat : float, list/array of floats
        Latitude coordinates of desired points.
    pt_lon : float, list/array of floats
        Longitude coordinates of desired points.
    height : float, list/array of floats
        Altitude (above geoid) of point. Initialized to None.
    level : int
        Model level or pseudo level or tile number.
        Initialized to None, meaning that all available levels in
        the cube are used.
    nearest : bool
        Specify whether to use 'nearest neighbour', instead
        of 'linear' method while extracting data. Default is False.

    Returns
    -------
    data_out : List
        List of single point values, corresponding to each point specified.

    Raises
    ------
    FlaskAnsError : If the number of latitude and longitude points are
        mismatched. OR if both level and height are passed as args.
        OR if the cube contains a time coordinate. OR if a pseudo level
        coordinate is requested, but not present in the cube. OR if the numbers
        of latitude/longitude and height points are mismatched. OR if height
        is requested but the cube does not contain an altitude coordinate.
    """
    # Check that input data is a (single) cube
    if not isinstance(icube, iris.cube.Cube):
        m = "Extract_pt:First argument must be a single cube."
        raise FlaskAnsError(m)

    # Check if the cube contains a time dimension, which is
    # currently unsupported.
    if icube.coords()[0].name() == "time":
        m = "Extract_pt:Cannot handle time dimension at present."
        raise FlaskAnsError(m)

    # Check that equal number of lat/lon pairs are passed in point coordinates.
    # Convert arguments to lists for easier processing if necessary.
    pt_lat1 = []
    pt_lon1 = []

    if not isinstance(pt_lat, list):
        pt_lat1.append(pt_lat)
        pt_lon1.append(pt_lon)

    else:
        for n_lat in np.arange(len(pt_lat)):
            pt_lat1.append(pt_lat[n_lat])
            pt_lon1.append(pt_lon[n_lat])

    if len(pt_lat1) != len(pt_lon1):
        m = "Extract_pt:Mismatch in number of lat/long values."
        raise FlaskAnsError(m)

    # Check that both level and height haven't been requested.
    if level is not None and height is not None:
        m = "Extract_pt: Both Level and Height requested."
        raise FlaskAnsError(m)

    # Check that the cube has a level coordinate if level has been requested.
    if (
        level is not None
        and not icube.coord("model_level_number")
        and not icube.coord("pseudo_level")
    ):
        m = "Extract_pt:Level requested, but not found in cube."
        raise FlaskAnsError(m)

    # Check that the number of height points is equal to the number of
    # lat/lon pairs. Convert the argument to a list for easier
    # processing if necessary.
    if height is not None:
        pt_hgt = []
        if isinstance(height, list):
            pt_hgt.extend(height)
        else:
            pt_hgt.append(height)

        if len(pt_lat1) != len(pt_hgt):
            m = "Extract_pt:Mismatch in number of points for lat/long/height."
            raise FlaskAnsError(m)

        # Check that heights have been merged with orography.
        if not icube.coords("altitude"):
            m = "Extract_pt:Height requested but input data does not contain altitude coordinate."
            raise FlaskAnsError(m)

        # Store the min and max altitudes from cube data so that user
        # cannot request points located below/ above that.
        # Will extract =min/max if beyond limits.
        hgt_min = icube.coord("altitude").points.min()
        hgt_max = icube.coord("altitude").points.max()

    # ---------- Finished checks -- begin processing -------------------------

    # If level specified, extract slice first
    if level is not None:
        try:
            icube = icube.extract(iris.Constraint(model_level_number=level))

        except iris.exceptions.ConstraintMismatchError:
            logger.debug("Model level number not available. Use pseudo level.")

        else:
            icube = icube.extract(iris.Constraint(pseudo_level=level))

    # Extract values for specified points lat/lon
    # NOTE: Does not seem to handle multiple points if 3-D
    data_out = []

    # Set lat/lon coordinates for model grid cell interpolation
    for n_lat1 in np.arange(len(pt_lat1)):
        latlon_coords = [
            ("latitude", pt_lat1[n_lat1]),
            ("longitude", pt_lon1[n_lat1]),
        ]

        if nearest:
            tcube = icube.interpolate(latlon_coords, iris.analysis.Nearest())
        else:
            tcube = icube.interpolate(latlon_coords, iris.analysis.Linear())

        # If height specified, interpolate to requested height
        if height is not None:
            # Set vertical coordinates for model grid cell interpolation
            point = max(pt_hgt[n_lat1], hgt_min)
            point = min(pt_hgt[n_lat1], hgt_max)
            hgt_coords = [("altitude", point)]

            if nearest:
                tcube = tcube.interpolate(hgt_coords, iris.analysis.Nearest())
            else:
                tcube = tcube.interpolate(hgt_coords, iris.analysis.Linear())

        data_out.append(tcube.data)

    return data_out


def _compute_taylor_statistics(ref, model):
    """Compute statistics necessary for Taylor diagram.

    Parameters
    ----------
    ref : numpy.array
        Array containing the values of the observations.
    model : numpy.array
        Array containing the values of the model.

    Returns
    -------
    std_model : float
        Standard deviation of the model data.
    corr_coeff : float
        Correlation coefficient between obs and model data.
    """
    mask = ~np.isnan(ref) & ~np.isnan(model)
    if np.sum(mask) < 3:
        return np.nan, np.nan, np.nan
    ref_clean = ref[mask]
    model_clean = model[mask]
    std_model = np.std(model_clean)
    corr_coeff = np.corrcoef(ref_clean, model_clean)[0, 1]
    return std_model, corr_coeff


def _latitude_weights(latitudes):
    """Compute weights based on latitude coordinates.

    Parameters
    ----------
    latitudes : numpy.array
        Array containing the latitude points for the observations.

    Returns
    -------
    weights : numpy.array
        Array containing the latitude weights for the observations.
    """
    weights = np.cos(np.radians(latitudes))
    weights[np.isnan(weights)] = 0
    return weights / np.nansum(weights)


def _fisher_z_transform(corrs, weights):
    """Compute Fisher transformation and weighting of correlation coefficients.

    Parameters
    ----------
    corrs : numpy.array
        Array containing the correlation coefficients for the stations.
    weights : numpy.array
        Array containing the latitude weights for the stations.

    Returns
    -------
    corr_mean : float
        Fisher-transformed and weighted mean correlation coefficient.
    """
    corrs = np.clip(corrs, -0.9999, 0.9999)
    z_vals = np.arctanh(corrs)
    mean_z = np.nansum(weights * z_vals)
    return np.tanh(mean_z)


def _aggregate_model_stats(obs, model, station_lats):
    """Aggregate model statistics across stations.

    Parameters
    ----------
    obs : numpy.array
        Array containing the observation data at stations.
    model : numpy.array
        Array containing the model data at stations.
    station_lats : numpy.array
        Array containing the latitudes of the stations.

    Returns
    -------
    std_model_mean : float
        Mean weighted standard deviation of the model at stations.
    corr_mean : float
        Fisher-transformed and weighted mean correlation coefficient.
    """
    n_stations = obs.shape[0]
    std_models = []
    corrs = []
    valid = []
    for i in range(n_stations):
        std_m, corr = _compute_taylor_statistics(obs[i, :], model[i, :])
        std_models.append(std_m)
        corrs.append(corr)
        valid.append(~np.isnan(corr))
    std_models = np.array(std_models)
    corrs = np.array(corrs)
    valid = np.array(valid)
    # Filter latitudes and weights for valid stations only
    weights = _latitude_weights(station_lats[valid])
    std_model_mean = np.nansum(weights * std_models[valid])
    corr_mean = _fisher_z_transform(corrs[valid], weights)
    return std_model_mean, corr_mean


def _quick_fix_cube(input_file, trace_gas):
    """Fix the model data cube (add bounds and unit conversion).

    Parameters
    ----------
    input_file : str
        String containing the model cube data filename.
    trace_gas : str
        Trace gas name.

    Returns
    -------
    cube : iris.cube.Cube
        Cube containing the model data after the small fixes.
    """
    cube = iris.load_cube(input_file)
    # Add bounds for lat and lon if not present
    cube = _add_bounds(cube)
    # Put observations and model data on same scale for ppm/ppb
    cube = TRACE_GASES_FACTOR[trace_gas] * cube
    # Change units accordingly
    cube.attributes["unit"] = TRACE_GASES_UNITS[trace_gas]
    # Add month and year coordinates
    iris.coord_categorisation.add_year(cube, "time")
    iris.coord_categorisation.add_month(cube, "time")
    # Extract years
    years = sorted(set(cube.coord("year").points))
    # Realize the data to avoid reading it from disk multiple times.
    cube.data  # noqa B018
    return {
        "cube": cube,
        "time": cube.coord("time", dim_coords=True),
        "years": years,
    }


def _colocate_obs_model(obs, model, w_id=False):
    """Colocate model data at observations locations.

    Parameters
    ----------
    obs : iris.cube.Cube
        Cube containing the observation data.
    model : iris.cube.Cube
        Cube containing the model data.
    w_id : bool
        Flag to indicate if the station ids from the observation data should be
        returned by the function.

    Returns
    -------
    v_obs : list
        List containing the observation data at the extracted locations.
    v_mod : list
        List containing the model data at the extracted locations.
    v_id : list or None
        List containing the station ids of the extracted locations.
    """
    # Latitude/longitude entries
    lats = obs.coord("latitude").points.tolist()
    lons = obs.coord("longitude").points.tolist()
    # Colocate
    colocated_cube = _extract_pt(
        model,
        lats,
        lons,
        nearest=True,
    )
    # Filter only valid points
    valid_indices = ~(obs.data.mask | np.isnan(colocated_cube))
    v_obs = obs.data[valid_indices].tolist()
    v_mod = [
        tg.item()
        for i, tg in enumerate(colocated_cube)
        if (tg.item() is not None) and valid_indices[i]
    ]
    # Extract station indices if w_id = True
    if w_id:
        v_id = (
            obs.coord("Station index (arbitrary)", dim_coords=True)
            .points[valid_indices]
            .tolist()
        )
    else:
        v_id = None
    return v_obs, v_mod, v_id


def _setup_growth_cube(cube_yearly, config, cube_og, type_cube):
    """Set yearly growth cube.

    Parameters
    ----------
    cube_yearly : iris.cube.Cube
        Cube containing yearly means.
    config : dict
        ESMValTool recipe configuration.
    cube_og : iris.cube.Cube
        Original cube before yearly mean processing.
    type_cube : str
        String to indicate if the input cubes are for model or obs data.

    Returns
    -------
    cube : iris.cube.Cube
        Cube containing the yearly absolute growth.
    """
    dim_coords = None
    if type_cube == "obs":
        dim_coords = [
            (
                iris.coords.DimCoord.from_coord(cube_yearly.coord("year")[1:]),
                0,
            ),
            (cube_yearly.coord("Station index (arbitrary)"), 1),
        ]
    elif type_cube == "model":
        dim_coords = [
            (
                iris.coords.DimCoord.from_coord(cube_yearly.coord("year")[1:]),
                0,
            ),
            (cube_yearly.coord("latitude"), 1),
            (cube_yearly.coord("longitude"), 2),
        ]
    # Create growth cube from yearly mean cube
    cube = iris.cube.Cube(
        np.diff(cube_yearly.data, axis=0),
        long_name=f"{config['trace_gas']}_growth",
        units=TRACE_GASES_UNITS[config["trace_gas"]],
        dim_coords_and_dims=dim_coords,
    )
    # Fill in attributes from original cube
    cube.attributes = cube_og.attributes
    # Fill in auxiliary coordinates from original cube for obs
    if type_cube == "obs":
        station_index_og = cube_og.coord("Station index (arbitrary)").points
        station_index_new = cube.coord("Station index (arbitrary)").points
        index_map = np.where(np.isin(station_index_og, station_index_new))[0]
        for aux_coord_name in [
            "altitude",
            "latitude",
            "longitude",
            "platform_name",
        ]:
            try:
                aux_coord = cube_og.coord(aux_coord_name)
                if cube_og.coord_dims(aux_coord):
                    coord_dims = cube_og.coord_dims(aux_coord)
                    if coord_dims == (1,):
                        # Subset the aux coord points using the index map
                        new_points = aux_coord.points[index_map]
                        # If the coordinate has bounds
                        if aux_coord.bounds is not None:
                            new_bounds = aux_coord.bounds[index_map]
                        else:
                            new_bounds = None
                        # Create a new AuxCoord and add it to cube
                        new_aux = iris.coords.AuxCoord(
                            new_points,
                            standard_name=aux_coord.standard_name,
                            long_name=aux_coord.long_name,
                            var_name=aux_coord.var_name,
                            units=aux_coord.units,
                            bounds=new_bounds,
                        )
                        cube.add_aux_coord(
                            new_aux,
                            (cube.coord_dims("Station index (arbitrary)")[0],),
                        )
                    else:
                        # For scalar coordinates or broadcasted ones, copy as is
                        cube.add_aux_coord(aux_coord.copy(), ())
            except iris.exceptions.CoordinateNotFoundError:
                msg = f"Auxiliary coordinate {aux_coord_name} not found in obs_cube."
                logger.debug(msg)
    return cube
