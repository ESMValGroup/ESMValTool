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

import os
import logging
import ast
import glob
import cartopy.crs as ccrs
import cf_units
import iris
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import arviz as az


logger = logging.getLogger(os.path.basename(__file__))


# /libs/select_key_or_default.py
def select_key_or_default(
        dirc, key, default=None, stack=True,
        numpck=__import__('numpy')):
    """Select specific key from dictionary.

    Arguments:
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
    Returns:
        out: numpck instance
            numpck instance from the extracted in the dictionary.
    """
    dirc = dict(sorted(dirc.items()))
    out = [dirc[name] for name in dirc if key in name]

    if len(out) == 0:
        out = default
    elif len(out) == 1:
        out = out[0]
    else:
        if stack:
            out = numpck.stack([i[0] for i in out])

    if isinstance(out, list):
        try:
            if stack:
                out = numpck.stack(out)[:, 0]
        except Exception as expt:
            logger.debug('select_key_or_default error %s', expt)

    return out


# /libs/iris_plus.py
def sort_time(cube, field, filename):
    """Sort time dimension in the iris cube.

    Arguments:
        cube: iris cube
            Input cube.
        field: str
            Variable name in cube.
        filename: str
            Filename of cube.
    Returns:
        cube: iris cube
            Cube with sorted and added time dimensions.
    """
    logger.debug('Sorting time for variable %s in cube %s', field, filename)

    cube.coord("time").bounds = None
    tcoord = cube.coord("time")
    tcoord.units = cf_units.Unit(tcoord.units.origin, calendar="gregorian")
    tcoord.convert_units("days since 1661-01-01 00:00:00")
    tcoord.units = cf_units.Unit(
        tcoord.units.origin, calendar="proleptic_gregorian"
    )
    cube.remove_coord("time")
    cube.add_dim_coord(tcoord, 0)

    try:
        iris.coord_categorisation.add_year(cube, 'time')
    except Exception as expt:
        logger.debug('sort_time could not add year %s', expt)

    try:
        try:
            cube.remove_coord("month")
        except Exception as expt:
            logger.debug('sort_time could not remove month %s', expt)
        iris.coord_categorisation.add_month_number(cube, 'time', name='month')
    except Exception as expt:
        logger.debug('sort_time could not add month number %s', expt)

    try:
        del cube.attributes["history"]
    except Exception as expt:
        logger.debug('sort_time could not remove cube history %s', expt)

    return cube


def insert_data_into_cube(data, eg_cube, mask=None):
    """ Insert data into cube following mask.

    Arguments:
        x: np.array
            data that we want to insert into the cube.
            Should have same shape of eg_cube, or same length as eg_cube
            or length equal to Trues in mask.
        eg_cube: iris cube
            The cube we want to insert data into.
        mask: Boolean array
            Array of shape or length x where True, will inster data.
            Default of None which means True for all points in eg_cube.
    Returns:
        eg_cube: iris cube
            cube with data replaced by x
    """

    pred_cube = eg_cube.copy()
    pred = pred_cube.data.copy().flatten()

    if mask is None:
        pred[:] = data
    else:
        pred[mask] = data

    pred_cube.data = pred.reshape(pred_cube.data.shape)
    return pred_cube


# /libs/namelist_functions.py
def read_variables_from_namelist(file_name):
    """Read variables from a file and create them with their original names.

    Arguments:
        file_name: str
            The name of the file containing the variables.
    Returns:
        dict: dict
            A dictionary of variable names and their values.
    Example Usage:
        file_name = 'variables.txt'
        read_variables = read_variables_from_file(file_name)
        # Create variables with their original names and assign the values
        for variable_name, variable_value in read_variables.items():
            exec(f"{variable_name} = {variable_value}")
        # Now you have the variables with their original names and values
        print(variable1)  # Output: Hello
        print(variable2)  # Output: 42
    """
    variables = {}

    def define_function(fun):
        return ast.literal_eval(fun.split('function ')[1].split(' at ')[0])

    with open(file_name, 'r', encoding="utf-8") as file:
        for line in file:
            parts = line.strip().split("::")
            if len(parts) == 2:
                variable_name = parts[0].strip()
                variable_value = parts[1].strip()
                if callable(variable_value):
                    # If the variable is a function, save its name
                    variable_value_set = variable_value
                elif (variable_value.startswith('"') and
                      variable_value.endswith('"')):
                    # If the variable is a string, remove the quotes
                    variable_value_set = variable_value[1:-1]
                elif (variable_value.startswith('[') and
                      variable_value.endswith(']')):
                    # If the variable is a list, parse it
                    try:
                        variable_value_set = ast.literal_eval(variable_value)
                    except Exception as expt:
                        logger.debug(
                            'read_variables_from_namelist error %s', expt)
                        functions = variable_value.split(', ')
                        variable_value_set = [
                            define_function(fun) for fun in functions
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
                            variables[variable_name], variable_value_set]
                else:
                    variables[variable_name] = variable_value_set
    return variables


# /libes/pymc_extras.py
def select_post_param(trace):
    """Selects paramaeters from a pymc nc trace file.

    Arguments:
        trace -- pymc netcdf trace file as filename or already opened
    Returns:
        dict of paramater values with each item names after the parameter
    """
    def select_post_param_name(name):
        out = trace.posterior[name].values
        a_shape = out.shape[0]
        b_shape = out.shape[1]
        new_shape = ((a_shape * b_shape), * out.shape[2:])
        return np.reshape(out, new_shape)

    try:
        trace = az.from_netcdf(trace, engine='netcdf4')
    except Exception as expt:
        logger.debug('select_post_params error %s', expt)
    params = trace.to_dict()['posterior']
    params_names = params.keys()
    params = [select_post_param_name(var) for var in params_names]
    return params, list(params_names)


def construct_param_comb(i, params, params_names, extra_params):
    """Construct a new dictionary containing parameters.

    Arguments:
        i: int
            index for the parameter to add.
        params: list
            List of input parameters.
        params_names: list
            List of input parameters' names.
        extra_params: dict
            Dictionary of extra parameters to be added to the dictionary.
    Returns:
        param_in: dict
            Dictionary of paramater values.
    """
    param_in = [
        param[i] if param.ndim == 1 else param[i, :] for param in params
    ]
    param_in = dict(zip(params_names, param_in))
    param_in.update(extra_params)
    return param_in


# /libs/read_variable_from_netcdf.py
def read_variable_from_netcdf(
        filename, *args, directory=None, subset_function=None, make_flat=False,
        units=None, subset_function_args=None, time_series=None,
        time_points=None, extent=None, return_time_points=False,
        return_extent=False):
    """Read data from a netCDF file.
    Assumes that the variables in the netcdf file all have the name
    "variable". Assumes that values < -9E9, you dont want. This could
    be different in some circumstances.

    Arguments:
        filename: str
            a string with filename or two element python list
            containing the name of the file and the target variable name.
            If just the sting of "filename" assumes variable name is "variable"
        directory: str
            The directory the file is in. Path can be in "filename" and None
            means no additional directory path needed.
        subset_function: func or list of funcs
            a function or list of functions to be applied to each data set.
        subset_function_args: dict or list of dicts
            If subset_function is a function, dict arguments for that function.
            If subset_function is a list, a list or dict constaining arguments
            for those functions in turn.
        make_flat: bool
            Should the output variable to flattened or remain cube
        time_series: list
            List comtaining range of years. If making flat and
            returned a time series, checks if that time series contains year.
    Returns:
        dataset: iris cube
            if make_flat, a numpy vector of the target variable, otherwise
            returns iris cube.
    """
    logger.info("Opening:")
    logger.info(filename)

    if filename[0] == '~' or filename[0] == '/' or filename[0] == '.':
        directory = ''

    try:
        if isinstance(filename, str):
            dataset = iris.load_cube(directory + filename, callback=sort_time)
        else:
            dataset = iris.load_cube(
                directory + filename[0], filename[1], callback=sort_time
            )
    except Exception as expt:
        try:
            dataset = iris.load_cube(directory + filename)
        except Exception as expt_:
            logger.debug(
                'read_variable_from_netcdf errors\n%s\n%s', expt, expt_)
            logger.debug(
                "Data cannot be opened." +
                "Check directory %s, filename %s or file format",
                directory, filename
            )
    coord_names = [coord.name() for coord in dataset.coords()]

    if dataset is None:
        return None

    if time_points is not None:
        if 'time' in coord_names:
            dataset = dataset.interpolate(
                [('time', time_points)], iris.analysis.Linear()
            )
        else:
            def addtime(time_point):
                time = iris.coords.DimCoord(
                    np.array([time_point]), standard_name='time',
                    units='days since 1661-01-01 00:00:00'
                )
                dataset_cp = dataset.copy()
                dataset_cp.add_aux_coord(time)
                return dataset_cp
            dataset_time = [addtime(time_point) for time_point in time_points]
            dataset = iris.cube.CubeList(dataset_time).merge_cube()
            # Unused? dataset0 = dataset.copy()

    if extent is not None:
        dataset = dataset.regrid(extent, iris.analysis.Linear())
        # Unused? dataset0 = dataset.copy()

    if units is not None:
        dataset.units = units

    if subset_function is not None:
        if isinstance(subset_function, list):
            for func, args in zip(subset_function, subset_function_args):
                try:
                    dataset = func(dataset, **args)
                except Exception as expt:
                    logger.debug('read_variable_from_netcdf error %s', expt)
                    logger.debug(
                        "Warning! function: %s not applied to file: %s",
                        func.__name__, directory + filename
                    )
        else:
            dataset = subset_function(dataset, **subset_function_args)

    if return_time_points:
        time_points = dataset.coord('time').points

    if return_extent:
        extent = dataset[0]

    if make_flat:
        if time_series is not None:
            years = dataset.coord('year').points
        try:
            dataset = dataset.data.flatten()
        except Exception as expt:
            logger.debug('read_variable_from_netcdf error %s', expt)
        if time_series is not None:
            if not years[0] == time_series[0]:
                dataset = np.append(
                    np.repeat(np.nan, years[0]-time_series[0]), dataset
                )
            if not years[-1] == time_series[1]:
                dataset = np.append(
                    dataset, np.repeat(np.nan, time_series[1]-years[-1])
                )

    if return_time_points:
        dataset = (dataset, time_points)

    if return_extent:
        dataset += (extent,)

    return dataset


def read_all_data_from_netcdf(
        y_filename, x_filename_list, *args, ca_filename=None,
        add_1s_columne=False, y_threshold=None, x_normalise01=False,
        scalers=None, check_mask=True, frac_random_sample=1.0,
        min_data_points_for_sample=None, **kw):
    """Read data from netCDF files

    Arguments:
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
        see read_variable_from_netcdf comments for *arg and **kw.
    Returns:
        y: np.array
            a numpy array of the target variable.
        x: np.array
            an n-D numpy array of the feature variables.
    """
    y_var, time_points, extent = read_variable_from_netcdf(
        y_filename, make_flat=True, *args,
        return_time_points=True, return_extent=True, **kw
    )

    if ca_filename is not None:
        ca_var = read_variable_from_netcdf(
            ca_filename, make_flat=True,
            time_points=time_points, extent=extent,
            *args, **kw
        )

    # Create a new categorical variable based on the threshold
    if y_threshold is not None:
        y_var = np.where(y_var >= y_threshold, 0, 1)

    x_var = np.zeros([len(y_var), len(x_filename_list)])

    for i, filename in enumerate(x_filename_list):
        x_var[:, i] = read_variable_from_netcdf(
            filename, make_flat=True,
            time_points=time_points, extent=extent,
            *args, **kw
        )

    if add_1s_columne:
        # add a column of ones to x
        x_var = np.column_stack((x_var, np.ones(len(x_var))))

    if check_mask:
        if ca_filename is not None:
            cells_we_want = np.array([
                np.all(rw > -9e9) and np.all(rw < 9e9)
                for rw in np.column_stack((x_var, y_var, ca_var))])
            ca_var = ca_var[cells_we_want]
        else:
            cells_we_want = np.array([
                np.all(rw > -9e9) and np.all(rw < 9e9)
                for rw in np.column_stack((x_var, y_var))])
        y_var = y_var[cells_we_want]
        x_var = x_var[cells_we_want, :]

    if x_normalise01 and scalers is None:
        try:
            scalers = np.array([np.min(x_var, axis=0), np.max(x_var, axis=0)])
        except Exception as expt:
            logger.debug('read_all_data_from_netcdf error %s', expt)
        squidge = (scalers[1, :] - scalers[0, :]) / (x_var.shape[0])
        scalers[0, :] = scalers[0, :] - squidge
        scalers[1, :] = scalers[1, :] + squidge

        test = scalers[1, :] == scalers[0, :]
        scalers[0, test] = 0.0
        scalers[1, test] = 1.0

    if frac_random_sample is None:
        frac_random_sample = 1000
    else:
        if min_data_points_for_sample is not None:
            min_data_frac = min_data_points_for_sample / len(y_var)
            frac_random_sample = max(frac_random_sample, min_data_frac)

    if frac_random_sample < 1:
        m_x = x_var.shape[0]
        selected_rows = np.random.choice(
            m_x, size=int(m_x * frac_random_sample), replace=False
        )
        y_var = y_var[selected_rows]
        x_var = x_var[selected_rows, :]
        if ca_filename is not None:
            ca_var = ca_var[selected_rows]

    if scalers is not None:
        x_var = (x_var - scalers[0, :]) / (scalers[1, :] - scalers[0, :])
        if check_mask:
            if ca_filename is not None:
                return y_var, x_var, ca_var, cells_we_want, scalers
        return y_var, x_var, cells_we_want, scalers

    if check_mask or frac_random_sample:
        if ca_filename is not None:
            return y_var, x_var, ca_var, cells_we_want

    return y_var, x_var, cells_we_want


# /fire_models/ConFire.py
class ConFire():
    """Class for a ConFire model object.
    This class contains all the necessary parameters to initialize a ConFire
    model instance and functions to run the evaluation.
    To initialize an object, you need to provide the parameters of the model.
    """
    def __init__(self, params, inference=False):
        """
        Initalise parameters and calculates the key variables needed to
        calculate burnt area.

        Arguments:
            params: dict or list of dict
            inference: bool
                Flag indicating to run inference or not.
        """
        self.inference = inference
        if self.inference:
            self.numpck = __import__('pytensor').tensor
        else:
            self.numpck = __import__('numpy')

        self.params = params

        def select_param_or_default(*args, **kw):
            return select_key_or_default(
                self.params, numpck=self.numpck, *args, **kw
            )

        self.controlid = self.params['controlID']
        self.control_direction = self.params['control_Direction']
        self.x0_param = select_param_or_default('x0', [0])
        self.log_control = select_param_or_default(
            'log_control', [False] * len(self.control_direction)
        )
        self.betas = select_param_or_default('betas', [[0]], stack=False)
        self.powers = select_param_or_default('powers', None, stack=False)
        self.driver_direction = self.params['driver_Direction']
        self.fmax = select_param_or_default('Fmax', None, stack=False)

    def burnt_area(
            self, data, return_controls=False, return_limitations=False):
        """Compute burnt area.

        Arguments:
            data: numpck instance
                Input drivers.
            return_controls: bool
            return_limitation: bool
        Returns:
            ba: numpck instance
                Burnt area data.
        """
        # finds controls
        def cal_control(cid=0):
            ids = self.controlid[cid]
            betas = self.betas[cid] * self.driver_direction[cid]

            x_i = data[:, ids]
            if self.powers is not None:
                x_i = self.numpck.power(x_i, self.powers[cid])

            out = self.numpck.sum(x_i * betas[None, ...], axis=-1)
            if self.log_control[cid]:
                out = self.numpck.log(out)
            out = out + self.x0_param[cid]
            return out

        controls = [cal_control(i) for i in range(len(self.controlid))]
        if return_controls:
            return controls

        def sigmoid(data, factor):
            """Compute sigmoid.

            Arguments:
                y: numpck instance
                    Input data.
                k: float
                    Exponential factor.
            Returns:
                Applied sigmoid function value.
            """
            if factor == 0:
                return None
            return 1.0 / (1.0 + self.numpck.exp(-data * factor))

        limitations = [
            sigmoid(y, k) for y, k in zip(controls, self.control_direction)
        ]

        if return_limitations:
            return limitations

        limitations = [lim for lim in limitations if lim is not None]

        b_a = self.numpck.prod(limitations, axis=0)
        if self.fmax is not None:
            b_a = sigmoid(self.fmax, 1.0) * b_a

        return b_a

    def emc_weighted(self, emc, precip, wd_pg):
        """Compute weighted Event Mean Concentration (EMC).

        Arguments:
            emc: iris.cube
                Cube containing the EMC.
            precip: iris cube
                Cube containing precipitation data.
            wd_pg: iris cube
                Cube containing wet days.
        Returns:
            emcw: iris cube
                Weighted cube of EMC.
        """
        try:
            wet_days = 1.0 - self.numpck.exp(-wd_pg * precip)
            emcw = (1.0 - wet_days) * emc + wet_days
        except Exception as expt:
            logger.debug('emc_weighted error %s', expt)
            emcw = emc.copy()
            emcw.data = 1.0 - self.numpck.exp(-wd_pg * precip.data)
            emcw.data = emcw.data + (1.0 - emcw.data) * emc.data
        return emcw

    def list_model_params(self, params, varnames=None):
        """Summary model parameters.

        Arguments:
            params: list
            varnames: list
        Returns:
            full_df: pd.DataFrame
                Dataframe containing the parameters values in each experiment.
        """
        controlid = params[0]['controlID']

        def list_one_line_of_parmas(param):
            def select_param_or_default(*args, **kw):
                return select_key_or_default(
                    param, numpck=__import__('numpy'), *args, **kw
                )

            control_direction = param['control_Direction']
            x0s = select_param_or_default('x0', [0])
            betas = select_param_or_default('betas', [[0]], stack=False)
            powers = select_param_or_default('powers', None, stack=False)
            driver_direction = param['driver_Direction']
            fmax = select_param_or_default('fmax', 1.0, stack=False)

            directions = [
                np.array(driver) * control for control, driver in
                zip(control_direction, driver_direction)
            ]
            betas = [
                direction * beta for direction, beta in zip(directions, betas)
            ]

            def mish_mash(beta, power, x_0):
                return np.append(np.column_stack((beta, power)).ravel(), x_0)
            varps = [
                mish_mash(beta, power, x0) for beta, power, x0 in zip(
                    betas, powers, x0s)
            ]

            return np.append(fmax, np.concatenate(varps))

        params_sorted = np.array([
            list_one_line_of_parmas(param) for param in params
        ])

        param = ['fmax']
        controln = ['']
        variable_name = ['']

        for ids, idn in zip(controlid, range(len(controlid))):
            for i_d in ids:
                param.append('beta')
                controln.append(str(idn))
                variable_name.append(varnames[i_d])
                param.append('power')
                controln.append(str(idn))
                variable_name.append(varnames[i_d])
            param.append('beta')
            controln.append(str(idn))
            variable_name.append('beta0')

        header_df = pd.DataFrame([param, controln, variable_name])
        index_labels = ["Parameter", "Control", "Variable"] + \
            list(map(str, range(1, params_sorted.shape[0] + 1)))

        data_df = pd.DataFrame(params_sorted)
        full_df = pd.concat([header_df, data_df], ignore_index=True)

        full_df.index = index_labels

        return full_df


def get_parameters(config):
    """Get parameters necessary for ConFire run.

    Arguments:
        config: dict
            Dictionary from ESMValTool recipe.
    Returns:
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
    work_dir = config['work_dir']
    confire_param = config['confire_param_dir']
    # **Define Paths for Parameter Files  and for outputs**
    output_dir = work_dir + '/ConFire_outputs/'
    # Parameter files (traces, scalers, and other model parameters)
    param_file_trace = glob.glob(confire_param + "trace*.nc")[0]
    param_file_none_trace = glob.glob(
        confire_param + "none_trace-params*.txt")[0]
    scale_file = glob.glob(confire_param + "scalers*.csv")[0]
    # **Load Variable Information and NetCDF Files**
    # Replace these lines with user-specified NetCDF files, ensuring variable
    # order is the same.
    nc_files = config['files_input']
    nc_dir = ''
    # **Load Driving Data and Land Mask**
    logger.info('Loading data for ConFire model...')
    scalers = pd.read_csv(scale_file).values
    _, driving_data, lmask, _ = read_all_data_from_netcdf(
        nc_files[0], nc_files, scalers=scalers, directory=nc_dir)
    # Load a sample cube (used for inserting data)
    eg_cube = read_variable_from_netcdf(nc_files[0], directory=nc_dir)
    # **Extract Model Parameters**
    logger.info('Loading ConFire model parameters...')
    params, params_names = select_post_param(param_file_trace)
    extra_params = read_variables_from_namelist(param_file_none_trace)
    control_direction = extra_params['control_Direction'].copy()
    return output_dir, params, params_names, extra_params, driving_data, \
        lmask, eg_cube, control_direction


def diagnostic_run_confire(config, model_name='model', timerange='none'):
    """Run ConFire as a diagnostic.
    The outputs from the model run are saved and the plots are returned.

    Arguments:
        config: dict
            Dictionary containing the ESMValTool recip configuration.
        model_name: str
            Model name to include in plots.
        timerange: str
            Time range of the input data to include in plots.
    Returns:
        figures: list
            List of matplotlib figures produced for the burnt area results.
    """
    # --------------------------------------------------------
    # This script runs the ConFire model using pre-generated parameter files.
    # It calculates burnt area under different control conditions and saves
    # results.
    # --------------------------------------------------------
    output_dir, params, params_names, extra_params, driving_data, lmask, \
        eg_cube, control_direction = get_parameters(config)

    # Number of samples to run from the trace file
    nsample_for_running = 100
    nexp = len(control_direction)
    # **Sample Iterations from Trace**
    nits = len(params[0])
    idx = range(0, nits, int(np.floor(nits / nsample_for_running)))
    # **Storage for Model Outputs (Including full model, physical &
    # stochastic controls)**
    out_cubes = [[] for _ in range(nexp + 2)]

    def run_model_into_cube(param_in, coord):
        """
        Runs the ConFire model with given parameters.
        Results are inserted into an iris cube.

        Arguments:
            param_in: dict
                Dictionary of input parameters for the run.
            coord: iris coord
                Coordinate to add to the output cube.
        Returns:
            cube: iris cube
                ConFire model output cube.
        """
        out = ConFire(param_in).burnt_area(driving_data)
        cube = insert_data_into_cube(out, eg_cube, lmask)
        cube.add_aux_coord(coord)
        return cube

    # **Run ConFire Model with Different Control Scenarios**
    logger.info('Running ConFire model...')
    for index, i in zip(idx, range(len(idx))):
        coord = iris.coords.DimCoord(i, "realization")
        param_in = construct_param_comb(
            index, params, params_names, extra_params
        )

        # **Run Full Model**
        out_cubes[0].append(run_model_into_cube(param_in, coord))

        # **Run Model with Individual Controls Turned On**
        for exp in range(nexp):
            param_in['control_Direction'][:] = [0] * nexp
            param_in['control_Direction'][exp] = control_direction[exp]
            param_in = construct_param_comb(
                exp, params, params_names, extra_params
            )
            out_cubes[exp + 1].append(run_model_into_cube(param_in, coord))

        # **Run Model with All Controls Off (Stochastic Control)**
        param_in['control_Direction'][:] = [0] * nexp
        out_cubes[exp + 2].append(run_model_into_cube(param_in, coord))

    # **Save Output Cubes**
    logger.info('Saving ConFire output cubes...')
    filenames_out = ['burnt_area_model'] + \
        ['burnt_area_control_' + str(i) for i in range(nexp)] + \
        ['burnt_area_control_stochastic']
    os.makedirs(output_dir, exist_ok=True)
    timerange = timerange.replace('/', '-')

    for i, o_c in enumerate(out_cubes):
        cubes = iris.cube.CubeList(o_c).merge_cube()
        iris.save(cubes, os.path.join(
            output_dir, f'{model_name}_{filenames_out[i]}_{timerange}.nc'
        ))

    # --------------------------------------------------------
    # **Visualization: Plot Resultant Maps**
    # --------------------------------------------------------
    # **Load and Plot Output Data**
    logger.info('Plotting model diagnostic outputs...')
    figures = []
    for i, filename in enumerate(filenames_out):
        filepath = os.path.join(
            output_dir, f'{model_name}_{filename}_{timerange}.nc'
        )
        fig, axes = plt.subplots(
            nrows=1, ncols=2, figsize=(12, 8),
            subplot_kw={'projection': ccrs.PlateCarree()})
        axes = axes.flatten()
        plotn = 0
        try:
            # Load saved NetCDF file
            cube = iris.load_cube(filepath)
            for pct in [5, 95]:
                pc_cube = 100 * cube.collapsed(
                    'realization', iris.analysis.PERCENTILE, percent=pct)[0]
                img = iris.quickplot.pcolormesh(
                    pc_cube, axes=axes[plotn], colorbar=False,
                    vmin=0., vmax=100., cmap='Oranges'
                )
                axes[plotn].set_title(
                    filename.replace("_", " ").capitalize() +
                    f' - {timerange}\n{model_name} - [{str(pct)}% percentile]'
                )
                axes[plotn].coastlines()
                # Add colorbar
                fig.colorbar(
                    img, ax=axes[plotn], orientation='vertical',
                    extend='neither', label='Burnt area [%]',
                    fraction=0.046, pad=0.04, shrink=0.3,
                )
                # Iterate plot number
                plotn = plotn + 1
            figures.append(fig)

        except Exception as expt:
            logger.info("Skipping %s due to error: %s", filename, expt)
            logger.debug("Skipping %s due to error: %s", filename, expt)

    return figures
