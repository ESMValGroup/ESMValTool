"""Utility functions to run the evaluation for the climate drivers of fire.

Functions are taken from:
https://github.com/douglask3/Bayesian_fire_models/tree/AR7_REF
maintainer: Douglas Kelley kelley_douglas
    https://orcid.org/0000-0003-1413-4969
    https://github.com/douglask3

Functions from:
    - /libs directory
    - /fire_models/ConFire.py
The directory and filename are indicated above each function definition.

"""

from __future__ import annotations

import ast
import copy
import logging
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING

import arviz as az
import cartopy.crs as ccrs
import cf_units
import iris
import iris.coords
import iris.exceptions
import iris.quickplot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

from esmvaltool.diag_scripts.shared import ProvenanceLogger

if TYPE_CHECKING:
    from types import ModuleType

    import pytensor

logger = logging.getLogger(Path(__file__).stem)


def get_provenance_record(
    ancestors: str | list,
    model_name: str,
    project: str | list,
    experiment: str | list,
    timerange: str,
    var: str | None = None,
) -> dict:
    """Create a provenance record describing the diagnostic data and plot.

    Parameters
    ----------
    ancestors : str, list, dict
        List of ancestor files.
    project : str
        Project facet.
    model_name : str
        Model facet.
    experiment : str
        Experiment facet.
    timerange : str
        Timerange facet as start_year/end_year.
    var : str, default to None
        Variable of the plot.

    Returns
    -------
    record : dict
        Dictionary containing the ancestors.
    """
    spec = (
        f"for the {model_name} ({project}-{experiment}) "
        f"for the time period {timerange} "
    )
    captions = {
        "burnt_fraction": "Burnt area fraction "
        + spec
        + "as computed with the ConFire model (Jones et al., 2024).",
        "fire_weather_control": "Fire weather control "
        + spec
        + "as computed with the ConFire model (Jones et al., 2024).",
        "fuel_load_continuity_control": "Fuel load continuity control "
        + spec
        + "as computed with the ConFire model (Jones et al., 2024).",
    }
    return {
        "caption": captions[var] if var is not None else "",
        "model": model_name,
        "project": project,
        "experiment": experiment,
        "timerange": timerange,
        "authors": ["lenhardt_julien", "kelley_douglas"],
        "ancestors": ancestors,
    }


# /libs/select_key_or_default.py
def _select_key_or_default(
    dirc: dict,
    key: str,
    default: None | list | int = None,
    numpck: None | ModuleType = None,
    stack: bool = True,
) -> ModuleType:
    """Select specific key from dictionary.

    Parameters
    ----------
    dirc: dict
        Input dictionary to select from.
    key: str
        Key to select in the dictionary.
    default:
        Default value to return if key not present.
    stack: bool
        Boolean flag if the extracted results should be stacked.
    numpck: imported package
        Package to use to stack the data if necessary.

    Returns
    -------
    out: numpck instance
        numpck instance from the extracted in the dictionary.
    """
    if numpck is None:
        numpck = __import__("numpy")
    dirc = dict(sorted(dirc.items()))
    out = [dirc[name] for name in dirc if key in name]

    if len(out) == 0:
        out = default
    elif len(out) == 1:
        out = out[0]
    elif stack:
        out = numpck.stack([i[0] for i in out])

    if isinstance(out, list):
        try:
            if stack:
                out = numpck.stack(out)[:, 0]
        except (ValueError, IndexError) as expt:
            logger.debug("select_key_or_default error %s", expt)

    return out


# /libs/iris_plus.py
def _sort_time(
    cube: iris.cube.Cube,
    field: str,
    filename: str,
) -> iris.cube.Cube:
    """Sort time dimension in the iris cube.

    Parameters
    ----------
    cube: iris cube
        Input cube.
    field: str
        Variable name in cube.
    filename: str
        Filename of cube.

    Returns
    -------
    cube: iris cube
        Cube with sorted and added time dimensions.
    """
    logger.debug("Sorting time for variable %s in cube %s", field, filename)

    cube.coord("time").bounds = None
    tcoord = cube.coord("time")
    tcoord.units = cf_units.Unit(tcoord.units.origin, calendar="gregorian")
    tcoord.convert_units("days since 1661-01-01 00:00:00")
    tcoord.units = cf_units.Unit(
        tcoord.units.origin,
        calendar="proleptic_gregorian",
    )
    cube.remove_coord("time")
    cube.add_dim_coord(tcoord, 0)

    if not cube.coords("year"):
        iris.coord_categorisation.add_year(cube, "time")

    if not cube.coords("month"):
        iris.coord_categorisation.add_month_number(cube, "time", name="month")

    return cube


def _insert_data_into_cube(
    data: np.array,
    eg_cube: iris.cube.Cube,
    mask: np.array | None = None,
) -> iris.cube.Cube:
    """Insert data into cube following mask.

    Parameters
    ----------
    data: np.array
        data that we want to insert into the cube.
        Should have same shape of eg_cube, or same length as eg_cube
        or length equal to Trues in mask.
    eg_cube: iris cube
        The cube we want to insert data into.
    mask: Boolean array
        Array of shape or length x where True, will inster data.
        Default of None which means True for all points in eg_cube.

    Returns
    -------
    pred_cube: iris cube
        cube with data replaced by x
    """
    pred_cube = eg_cube.copy()
    # Only keep parent's attributes
    attrs_parent = {
        key: item
        for key, item in eg_cube.attributes.items()
        if key[:6] == "parent"
    }
    pred_cube.attributes = attrs_parent
    # Copy data
    pred = pred_cube.core_data().copy().flatten()

    if mask is None:
        pred[:] = data
    else:
        pred[mask] = data

    pred_cube.data = pred.reshape(pred_cube.core_data().shape)
    return pred_cube


# /libs/namelist_functions.py
def _read_variables_from_namelist(file_name: str) -> dict:
    """Read variables from a file and create them with their original names.

    Parameters
    ----------
    file_name: str
        The name of the file containing the variables.

    Returns
    -------
    dict: dict
        A dictionary of variable names and their values.

    Example Usage:
    --------------
    file_name = 'variables.txt'
    read_variables = read_variables_from_file(file_name)
    # Create variables with their original names and assign the values
    for variable_name, variable_value in read_variables.items():
        exec(f"{variable_name} = {variable_value}")
    # Now you have the variables with their original names and values
    print(variable1)  # Output: Hello
    print(variable2)  # Output: 42

    Example input file:
    -------------------
    variables.txt:
        y_filen:: 'filename.nc'
        x_filen_list:: ['file1.nc', 'file2.nc', 'file3.nc']
        CA_filen:: None
        dir:: 'path/to/driving/data/'
        filename_out:: '_filename_file_output'
        dir_outputs:: 'path/to/outputs/'
        subset_function:: function_name
        subset_function_args:: {'function_arg': [0, 1, 2]}
        out_file:: '_output_file'
        data_file:: 'path/to/data_file.nc'
        trace_file:: 'path/to/trace_file.nc'
        other_params_file:: 'path/to/param_file.txt'
        scale_file:: 'path/to/scale_file.csv'
    """
    variables = {}

    def _define_function(fun: str):
        return ast.literal_eval(fun.split("function ")[1].split(" at ")[0])

    with Path.open(file_name, encoding="utf-8") as file:
        for line in file:
            parts = line.strip().split("::")
            if len(parts) == 2:
                variable_name = parts[0].strip()
                variable_value = parts[1].strip()
                if callable(variable_value):
                    # If the variable is a function, save its name
                    variable_value_set = variable_value
                elif variable_value.startswith(
                    '"',
                ) and variable_value.endswith('"'):
                    # If the variable is a string, remove the quotes
                    variable_value_set = variable_value[1:-1]
                elif variable_value.startswith(
                    "[",
                ) and variable_value.endswith("]"):
                    # If the variable is a list, parse it
                    try:
                        variable_value_set = ast.literal_eval(variable_value)
                    except Exception as expt:
                        logger.debug(
                            "_read_variables_from_namelist error %s",
                            expt,
                        )
                        functions = variable_value.split(", ")
                        variable_value_set = [
                            _define_function(fun) for fun in functions
                        ]
                else:
                    try:
                        # Try to parse the variable as a dictionary
                        variable_value_set = ast.literal_eval(variable_value)
                    except (SyntaxError, NameError):
                        # If parsing fails, assume it's a non-string,
                        # or non-list variable
                        variable_value_set = variable_value
                if variable_name in variables:
                    if isinstance(variables[variable_name], list):
                        variables[variable_name].append(variable_value_set)
                    else:
                        variables[variable_name] = [
                            variables[variable_name],
                            variable_value_set,
                        ]
                else:
                    variables[variable_name] = variable_value_set
    return variables


# /libes/pymc_extras.py
def _select_post_param(trace: str) -> dict:
    """Select parameters from a pymc nc trace file.

    Parameters
    ----------
    trace: str
        pymc netcdf trace file as filename or already opened

    Returns
    -------
    dict of paramater values with each item names after the parameter
    """

    def _select_post_param_name(name: str) -> np.array:
        out = trace.posterior[name].to_numpy()
        a_shape = out.shape[0]
        b_shape = out.shape[1]
        new_shape = ((a_shape * b_shape), *out.shape[2:])
        return np.reshape(out, new_shape)

    try:
        trace = az.from_netcdf(trace, engine="netcdf4")
    except (ValueError, OSError) as expt:
        logger.debug("_select_post_param error %s", expt)
    params = trace.to_dict()["posterior"]
    params_names = params.keys()
    params = [_select_post_param_name(var) for var in params_names]
    return params, list(params_names)


def _construct_param_comb(
    i: int,
    params: list,
    params_names: list,
    extra_params: dict,
) -> dict:
    """Construct a new dictionary containing parameters.

    Parameters
    ----------
    i: int
        index for the parameter to add.
    params: list
        List of input parameters.
    params_names: list
        List of input parameters' names.
    extra_params: dict
        Dictionary of extra parameters to be added to the dictionary.

    Returns
    -------
    param_input: dict
        Dictionary of paramater values.
    """
    param_input = [
        param[i] if param.ndim == 1 else param[i, :] for param in params
    ]
    param_input = dict(zip(params_names, param_input, strict=False))
    param_input.update(extra_params)
    return param_input


# /libs/read_variable_from_netcdf.py
def _read_variable_from_netcdf(
    filename: str,
    directory: str | None = None,
    subset_function: Callable | list[Callable] | None = None,
    units: str | None = None,
    subset_function_args: dict | list[dict] | None = None,
    time_series: list | None = None,
    time_points: np.array | float | None = None,
    extent: np.array | list | None = None,
    *,
    make_flat: bool = False,
    return_time_points: bool = False,
    return_extent: bool = False,
) -> iris.cube.Cube:
    """Read data from a netCDF file.

    Assumes that the variables in the netcdf file all have the name
    "variable". Assumes that values < -9E9, you dont want. This could
    be different in some circumstances.

    Parameters
    ----------
    filename: str
        a string with filename or two element python list
        containing the name of the file and the target variable name.
        If just the string of "filename" assumes variable name is "variable"
    directory: str
        The directory the file is in. Path can be in "filename" and None
        means no additional directory path needed.
    subset_function: func or list of funcs
        a function or list of functions to be applied to each data set.
    subset_function_args: dict or list of dicts
        If subset_function is a function, dict arguments for that function.
        If subset_function is a list or dict containing arguments
        for those functions in turn.
    make_flat: bool
        Should the output variable to flattened or remain cube
    time_series: list
        List comtaining range of years. If making flat and
        returned a time series, checks if that time series contains year.

    Returns
    -------
    dataset: iris cube
        if make_flat, a numpy vector of the target variable, otherwise
        returns iris cube.
    """
    logger.info("Opening:")
    logger.info(filename)

    if filename[0] == "~" or filename[0] == "/" or filename[0] == ".":
        directory = ""

    if isinstance(filename, str):
        dataset = iris.load_raw(
            Path(directory) / filename,
            callback=_sort_time,
        )
    else:
        dataset = iris.load_raw(
            Path(directory) / filename[0],
            filename[1],
            callback=_sort_time,
        )
    dataset = dataset[0]

    coord_names = [coord.name() for coord in dataset.coords()]

    if dataset is None:
        return None

    if time_points is not None:
        if "time" in coord_names:
            dataset = dataset.interpolate(
                [("time", time_points)],
                iris.analysis.Linear(),
            )
        else:
            time_coord_to_add = [
                iris.coords.DimCoord(
                    np.array([time_point]),
                    standard_name="time",
                    units="days since 1661-01-01 00:00:00",
                )
                for time_point in time_points
            ]
            dataset_time = [
                dataset.add_aux_coord(t) for t in time_coord_to_add
            ]
            dataset = iris.cube.CubeList(dataset_time).merge_cube()

    if extent is not None:
        dataset = dataset.regrid(extent, iris.analysis.Linear())

    if units is not None:
        dataset.units = units

    if subset_function is not None:
        if isinstance(subset_function, list):
            for func, _args in zip(
                subset_function,
                subset_function_args,
                strict=False,
            ):
                try:
                    dataset = func(dataset, **_args)
                except ValueError as expt:
                    logger.debug("_read_variable_from_netcdf error %s", expt)
                    logger.debug(
                        "Warning! function: %s not applied to file: %s",
                        func.__name__,
                        directory + filename,
                    )
        else:
            dataset = subset_function(dataset, **subset_function_args)

    if return_time_points:
        time_points = dataset.coord("time").points

    if return_extent:
        extent = dataset[0]

    if make_flat:
        if time_series is not None:
            years = dataset.coord("year").points
        try:
            dataset = dataset.data.flatten()
        except ValueError as expt:
            logger.debug("_read_variable_from_netcdf error %s", expt)
        if time_series is not None:
            if years[0] != time_series[0]:
                dataset = np.append(
                    np.repeat(np.nan, years[0] - time_series[0]),
                    dataset,
                )
            if years[-1] != time_series[1]:
                dataset = np.append(
                    dataset,
                    np.repeat(np.nan, time_series[1] - years[-1]),
                )

    if return_time_points:
        dataset = (dataset, time_points)

    if return_extent:
        dataset += (extent,)

    return dataset


def _read_all_data_from_netcdf(
    y_filename: str | list,
    x_filename_list: list,
    ca_filename: list | None = None,
    y_threshold: float | None = None,
    scalers: np.array | None = None,
    frac_random_sample: float = 1.0,
    min_data_points_for_sample: int | None = None,
    *args: tuple,
    add_1s_columne: bool = False,
    x_normalise01: bool = False,
    check_mask: bool = True,
    **kw: dict,
) -> tuple[np.array]:
    """Read data from netCDF files.

    Parameters
    ----------
    y_filename: list
        a two element python list containing the name of the file and
        the target variable name.
    x_filename_list: list
        a python list of filename containing the feature variables.
    ca_filename: list
        a python list of filename containing the area of the cover type.
    y_threshold: float
        if converting y into boolean, the threshold we use to split into
        0's and 1's.
    add_1s_columne: bool
        useful for if using for regressions. Adds a variable
        of just 1's t rperesent y = SUM(a_i * x_i) + c
    x_normalise01: bool
        If True, then x's are normalised between 0 and 1.
    scalers: None or np.array
        None or numpy array of shape 2 by n. columns of x.
        Defines what scalers (min and max) to apply to each x column.
        If None, doesn't apply anything.
    check_mask: bool
        If True, simple checks if there are any large
        negtaive numbers and makes them out. Assunes that values < -9E9,
        you dont want. This could be different in some circumstances.
    frac_random_sample: int
        fraction of data to be returned
    args: tuple
        See _read_variable_from_netcdf comments.
    kw: dict
        See _read_variable_from_netcdf comments.

    Returns
    -------
    y: np.array
        a numpy array of the target variable.
    x: np.array
        an n-D numpy array of the feature variables.
    """
    y_var, time_points, extent = _read_variable_from_netcdf(
        *args,
        y_filename,
        make_flat=True,
        return_time_points=True,
        return_extent=True,
        **kw,
    )

    if ca_filename is not None:
        ca_var = _read_variable_from_netcdf(
            *args,
            ca_filename,
            make_flat=True,
            time_points=time_points,
            extent=extent,
            **kw,
        )

    # Create a new categorical variable based on the threshold
    if y_threshold is not None:
        y_var = np.where(y_var >= y_threshold, 0, 1)

    x_var = np.zeros([len(y_var), len(x_filename_list)])

    for i, filename in enumerate(x_filename_list):
        x_var[:, i] = _read_variable_from_netcdf(
            *args,
            filename,
            make_flat=True,
            time_points=time_points,
            extent=extent,
            **kw,
        )

    if add_1s_columne:
        # add a column of ones to x
        x_var = np.column_stack((x_var, np.ones(len(x_var))))

    if check_mask:
        if ca_filename is not None:
            cells_we_want = np.array(
                [
                    np.all(rw > -9e9) and np.all(rw < 9e9)
                    for rw in np.column_stack((x_var, y_var, ca_var))
                ],
            )
            ca_var = ca_var[cells_we_want]
        else:
            cells_we_want = np.array(
                [
                    np.all(rw > -9e9) and np.all(rw < 9e9)
                    for rw in np.column_stack((x_var, y_var))
                ],
            )
        y_var = y_var[cells_we_want]
        x_var = x_var[cells_we_want, :]

    if x_normalise01 and scalers is None:
        scalers = np.array([np.min(x_var), np.max(x_var)])
        squidge = (scalers[1, :] - scalers[0, :]) / (x_var.shape[0])
        scalers[0, :] = scalers[0, :] - squidge
        scalers[1, :] = scalers[1, :] + squidge

        test = scalers[1, :] == scalers[0, :]
        scalers[0, test] = 0.0
        scalers[1, test] = 1.0

    if frac_random_sample is None:
        frac_random_sample = 1000
    elif min_data_points_for_sample is not None:
        min_data_frac = min_data_points_for_sample / len(y_var)
        frac_random_sample = max(frac_random_sample, min_data_frac)

    if frac_random_sample < 1:
        m_x = x_var.shape[0]
        selected_rows = np.random.choice(
            m_x,
            size=int(m_x * frac_random_sample),
            replace=False,
        )
        y_var = y_var[selected_rows]
        x_var = x_var[selected_rows, :]
        if ca_filename is not None:
            ca_var = ca_var[selected_rows]

    if scalers is not None:
        x_var = (x_var - scalers[0, :]) / (scalers[1, :] - scalers[0, :])
        if check_mask and ca_filename is not None:
            return y_var, x_var, ca_var, cells_we_want, scalers
        return y_var, x_var, cells_we_want, scalers, None

    if (check_mask or frac_random_sample) and ca_filename is not None:
        return y_var, x_var, cells_we_want, None, ca_var

    return y_var, x_var, cells_we_want, None, None


# /fire_models/ConFire.py
class ConFire:
    """Class for a ConFire model object.

    This class contains all the necessary parameters to initialize a ConFire
    model instance and functions to run the evaluation.
    To initialize an object, you need to provide the parameters of the model.
    """

    def __init__(
        self,
        params: dict | list,
        *,
        inference: bool = False,
    ) -> None:
        """Initalise parameters and calculates the key variables.

        Parameters
        ----------
        params: dict or list of dict
        inference: bool
            Flag indicating to run inference or not.
        """
        self.inference = inference
        if self.inference:
            self.numpck = __import__("pytensor").tensor
        else:
            self.numpck = __import__("numpy")

        self.params = params

        def _select_param_or_default(*args: tuple, **kw: dict) -> ModuleType:
            return _select_key_or_default(
                *args,
                dirc=self.params,
                numpck=self.numpck,
                **kw,
            )

        self.controlid = self.params["controlID"]
        self.control_direction = self.params["control_Direction"]
        self.x0_param = _select_param_or_default(key="x0", default=[0])
        self.log_control = _select_param_or_default(
            key="log_control",
            default=[False] * len(self.control_direction),
        )
        self.betas = _select_param_or_default(
            key="betas",
            default=[[0]],
            stack=False,
        )
        self.powers = _select_param_or_default(
            key="powers",
            default=None,
            stack=False,
        )
        self.driver_direction = self.params["driver_Direction"]
        self.fmax = _select_param_or_default(
            key="Fmax",
            default=None,
            stack=False,
        )

    def burnt_area(
        self,
        data: np.array | pytensor.tensor,
        *,
        return_controls: bool = False,
        return_limitations: bool = False,
    ) -> np.array | pytensor.tensor:
        """Compute burnt area.

        Parameters
        ----------
        data: numpck instance
            Input drivers.
        return_controls: bool
        return_limitation: bool

        Returns
        -------
        ba: numpck instance
            Burnt area data.
        """

        # finds controls
        def _cal_control(cid: int = 0) -> np.array | pytensor.tensor:
            ids = self.controlid[cid]
            betas = self.betas[cid] * self.driver_direction[cid]

            x_i = data[:, ids]
            if self.powers is not None:
                x_i = self.numpck.power(x_i, self.powers[cid])

            out = self.numpck.sum(x_i * betas[None, ...], axis=-1)
            if self.log_control[cid]:
                out = self.numpck.log(out)
            return out + self.x0_param[cid]

        controls = [_cal_control(i) for i in range(len(self.controlid))]
        if return_controls:
            return controls

        def _sigmoid(
            data: np.array | pytensor.tensor,
            factor: float,
        ) -> np.array | pytensor.tensor:
            """Compute sigmoid.

            Parameters
            ----------
            data: numpck instance
                Input data.
            factor: float
                Exponential factor.

            Returns
            -------
            Applied sigmoid function value.
            """
            if factor == 0:
                return None
            if self.inference:
                return self.numpck.math.sigmoid(-data * factor)
            return 1.0 / (1.0 + self.numpck.exp(-data * factor))

        limitations = [
            _sigmoid(y, k)
            for y, k in zip(controls, self.control_direction, strict=False)
        ]

        if return_limitations:
            return limitations

        limitations = [lim for lim in limitations if lim is not None]

        b_a = self.numpck.prod(limitations, axis=0)
        if self.fmax is not None:
            b_a = _sigmoid(self.fmax, 1.0) * b_a

        return b_a


def _get_parameters(config: dict) -> tuple:
    """Get parameters necessary for ConFire run.

    Parameters
    ----------
    config: dict
        Dictionary from ESMValTool recipe.

    Returns
    -------
    output_dir: str
        Path to output directory.
    params: list
        List of parameters.
    params_names: list
        List of parameters names.
    extra_params: dict
        Dictionary of additional parameters to add to the output cubes.
    driving_data: numpy.array or iris.cube
        Object containing the concatenated drivers.
    lmask: numpy.array
        Land-sea mask.
    eg_cube: iris.cube
        Example of cube to be used to insert output data.
    control_direction: list
        List containing different experiment setups.
    """
    work_dir = config["work_dir"]
    confire_param = config["confire_param_dir"]
    # **Define Paths for Parameter Files  and for outputs**
    output_dir = work_dir + "/ConFire_outputs/"
    # Parameter files (traces, scalers, and other model parameters)
    param_file_trace = list(Path(confire_param).glob("trace*.nc"))[0]
    param_file_none_trace = list(
        Path(confire_param).glob("none_trace-params*.txt"),
    )[0]
    scale_file = list(Path(confire_param).glob("scalers*.csv"))[0]
    # **Load Variable Information and NetCDF Files**
    # Replace these lines with user-specified NetCDF files, ensuring variable
    # order is the same.
    nc_files = config["files_input"]
    nc_dir = ""
    # **Load Driving Data and Land Mask**
    logger.info("Loading data for ConFire model...")
    scalers = pd.read_csv(scale_file).to_numpy()
    _, driving_data, lmask, _, _ = _read_all_data_from_netcdf(
        y_filename=nc_files[0],
        x_filename_list=nc_files,
        scalers=scalers,
        directory=nc_dir,
    )
    # Load a sample cube (used for inserting data)
    eg_cube = _read_variable_from_netcdf(nc_files[0], directory=nc_dir)
    # **Extract Model Parameters**
    logger.info("Loading ConFire model parameters...")
    params, params_names = _select_post_param(param_file_trace)
    extra_params = _read_variables_from_namelist(param_file_none_trace)
    control_direction = extra_params["control_Direction"].copy()
    return (
        output_dir,
        params,
        params_names,
        extra_params,
        driving_data,
        lmask,
        eg_cube,
        control_direction,
    )


def _setup_cube_output(
    cube: iris.cube.Cube,
    output: str,
    provenance: dict,
) -> iris.cube.Cube:
    """Set up the output cube.

    Apply:
        - replace variable name
        - remove standard name
        - replace long_name
        - replace units
        - apply scaling factor

    Parameters
    ----------
    cube: iris cube
        Input cube to modify.
    output: str
        Output variable contained in the cube.
    provenance: dict
        Dictionary w/ provenance record.

    Returns
    -------
    cube: iris cube
        Modified output cube.
    """
    # Default attributes to fill in
    parameter_dict = {
        "burnt_fraction": {
            "units": "%",
            "long_name": "Burnt Fraction",
            "factor": 100.0,
        },
        "fire_weather_control": {
            "units": "1",
            "long_name": "Fire Weather control",
            "factor": 1.0,
        },
        "fuel_load_continuity_control": {
            "units": "1",
            "long_name": "Fuel Load/Continuity control",
            "factor": 1.0,
        },
    }
    if output == "burnt_fraction":
        cube.var_name = "burnt_fraction"
        cube.standard_name = None
        cube.long_name = parameter_dict["burnt_fraction"]["long_name"]
        cube.units = parameter_dict["burnt_fraction"]["units"]
        cube.data *= parameter_dict["burnt_fraction"]["factor"]
    elif output == "fire_weather_control":
        cube.var_name = "fire_weather_control"
        cube.standard_name = None
        cube.long_name = parameter_dict["fire_weather_control"]["long_name"]
        cube.units = parameter_dict["fire_weather_control"]["units"]
        cube.data *= parameter_dict["fire_weather_control"]["factor"]
    elif output == "fuel_load_continuity_control":
        cube.var_name = "fuel_load_continuity_control"
        cube.standard_name = None
        cube.long_name = parameter_dict["fuel_load_continuity_control"][
            "long_name"
        ]
        cube.units = parameter_dict["fuel_load_continuity_control"]["units"]
        cube.data *= parameter_dict["fuel_load_continuity_control"]["factor"]
    else:
        logger.debug("Output %s cannot be processed. Saving output as is.")
    cube.attributes["ancestors"] = provenance["ancestors"]
    return cube


def diagnostic_run_confire(
    config: dict,
    model_name: str = "model",
    timerange: str = "none",
    project: str | list[str] = "project",
    experiment: str | list[str] = "exp",
) -> list:
    """Run ConFire as a diagnostic.

    The outputs from the model run are saved and the plots are returned.

    Parameters
    ----------
    config: dict
        Dictionary containing the ESMValTool recip configuration.
    model_name: str
        Model name to include in plots.
    timerange: str
        Time range of the input data to include in plots.
    project: str or list of str
        Project of the model data.
    experiment: str or list of str
        Experiment of the model data.

    Returns
    -------
    figures: list
        List of matplotlib figures produced for the burnt area results.
    """
    # --------------------------------------------------------
    # This script runs the ConFire model using pre-generated parameter files.
    # It calculates burnt area under different control conditions and saves
    # results.
    # --------------------------------------------------------
    (
        output_dir,
        params,
        params_names,
        extra_params,
        driving_data,
        lmask,
        eg_cube,
        control_direction,
    ) = _get_parameters(config)

    # Number of samples to run from the trace file
    nsample_for_running = 100
    nexp = len(control_direction)
    # **Sample Iterations from Trace**
    nits = len(params[0])
    idx = range(0, nits, int(np.floor(nits / nsample_for_running)))
    # **Storage for Model Outputs (Including full model, physical controls)**
    out_cubes = [[] for _ in range(nexp + 1)]

    def _run_model_into_cube(
        param_in: dict,
        coord: iris.coords.Coord,
    ) -> iris.cube.Cube:
        """
        Run the ConFire model with given parameters.

        Results are inserted into an iris cube.

        Parameters
        ----------
        param_in: dict
            Dictionary of input parameters for the run.
        coord: iris coord
            Coordinate to add to the output cube.

        Returns
        -------
        cube: iris cube
            ConFire model output cube.
        """
        out = ConFire(param_in).burnt_area(driving_data)
        cube = _insert_data_into_cube(out, eg_cube, lmask)
        cube.add_aux_coord(coord)
        return cube

    # **Run ConFire Model with Different Control Scenarios**
    logger.info("Running ConFire model...")
    for index, i in zip(idx, range(len(idx)), strict=False):
        coord = iris.coords.DimCoord(i, "realization")
        param_in_full = _construct_param_comb(
            index,
            params,
            params_names,
            copy.deepcopy(extra_params),
        )

        # **Run Full Model**
        out_cubes[0].append(_run_model_into_cube(param_in_full, coord))

        # **Run Model with Individual Controls Turned On**
        for exp in range(nexp):
            param_exp = copy.deepcopy(param_in_full)
            param_exp["control_Direction"][:] = [0] * nexp
            param_exp["control_Direction"][exp] = control_direction[exp]
            out_cubes[exp + 1].append(_run_model_into_cube(param_exp, coord))

    # **Save Output Cubes**
    # In the current configuration, outputs from the model runs are:
    #   - burnt area (model)
    #   - fire weather (control direction 0)
    #   - fuel loads (control direction 1)
    logger.info("Saving ConFire output cubes...")
    Path.mkdir(output_dir, parents=True, exist_ok=True)
    timerange = timerange.replace("/", "-")
    ancestors = [
        config["files_input"][i][0]
        for i in range(len(config["files_input"]))
        if config["files_input"][i][1] != "vpd"
    ]
    if "vpd" in config["var_order"]:
        ancestors += config["provenance_record_vpd"]["ancestors"]
    for i, o_c in enumerate(out_cubes):
        filename = (
            Path(output_dir)
            / f"{config['filenames_out'][i]}_{model_name}_{timerange}.nc"
        )
        provenance = get_provenance_record(
            ancestors=ancestors,
            model_name=model_name,
            project=project,
            experiment=experiment,
            timerange=timerange,
            var=config["filenames_out"][i],
        )
        cubes = iris.cube.CubeList(o_c).merge_cube()
        cubes = _setup_cube_output(
            cubes,
            config["filenames_out"][i],
            provenance,
        )
        iris.save(cubes, filename)
        if not config["remove_confire_files"]:
            with ProvenanceLogger(config) as provenance_logger:
                provenance_logger.log(filename, provenance)

    # --------------------------------------------------------
    # **Visualization: Plot Resultant Maps**
    # --------------------------------------------------------
    # **Load and Plot Output Data**
    logger.info("Plotting model diagnostic outputs...")
    figures = []
    parameter_plot = {
        "burnt_fraction": {
            "vmin": 0.0,
            "vmax": 100.0,
        },
        "fire_weather_control": {
            "vmin": 0.0,
            "vmax": 1.0,
        },
        "fuel_load_continuity_control": {
            "vmin": 0.0,
            "vmax": 1.0,
        },
    }
    # Colormap setup
    oranges = plt.cm.get_cmap("Oranges", 256)
    colors = np.vstack(
        [
            [1, 1, 1, 1],
            oranges(np.linspace(0, 1, 256)),
        ],
    )
    cmap = ListedColormap(colors)
    # Plots
    for filename in config["filenames_out"]:
        filepath = Path(output_dir) / f"{filename}_{model_name}_{timerange}.nc"
        fig, axes = plt.subplots(
            nrows=1,
            ncols=2,
            figsize=(12, 8),
            subplot_kw={"projection": ccrs.Robinson()},
        )
        axes = axes.flatten()
        plotn = 0
        try:
            # Load saved NetCDF file
            cube = iris.load_cube(filepath)
            for pct in [5, 95]:
                pc_cube = cube.collapsed(
                    "time",
                    iris.analysis.MEAN,
                ).collapsed(
                    "realization",
                    iris.analysis.PERCENTILE,
                    percent=pct,
                )
                img = iris.quickplot.pcolormesh(
                    pc_cube,
                    axes=axes[plotn],
                    colorbar=False,
                    vmin=parameter_plot[filename]["vmin"],
                    vmax=parameter_plot[filename]["vmax"],
                    cmap=cmap,
                )
                axes[plotn].set_title(
                    (
                        filename.replace("_", " ").capitalize()
                        + f" [{pct!s}th percentile]\n"
                        + f"{model_name} "
                        + f"({'/'.join(project)} - "
                        + f"{'/'.join(experiment)}, "
                        + f"{timerange})"
                    ),
                )
                axes[plotn].coastlines()
                # Add colorbar
                fig.colorbar(
                    img,
                    ax=axes[plotn],
                    orientation="vertical",
                    extend="neither",
                    label=filename.replace("_", " ").capitalize()
                    + f" [{cube.units}]",
                    fraction=0.046,
                    pad=0.04,
                    shrink=0.3,
                )
                # Iterate plot number
                plotn = plotn + 1
            figures.append(fig)

        except ValueError as expt:
            logger.info("Skipping %s due to error: %s", filename, expt)
            logger.debug("Skipping %s due to error: %s", filename, expt)

    return figures
