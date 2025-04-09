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

import arviz as az
import cartopy.crs as ccrs
import cf_units
import iris
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import glob

logger = logging.getLogger(os.path.basename(__file__))


# /libs/select_key_or_default.py
def select_key_or_default(
        dirc, key, default=None, stack=True,
        numPCK=__import__('numpy')):
    dirc = dict(sorted(dirc.items()))
    out = [dirc[name] for name in dirc if key in name]
    
    if len(out) == 0: 
        out = default 
    elif len(out) == 1: 
        out = out[0]
    else: 
        if stack:
            out = numPCK.stack([i[0] for i in out])
        
    if type(out) is list:         
        try:
            if stack:
                out = numPCK.stack(out)[:, 0]
        except:
            pass
       
    return(out)


# /libs/iris_plus.py
def sort_time(cube, field, filename):
    
    cube.coord("time").bounds = None
    tcoord = cube.coord("time")
    tcoord.units = cf_units.Unit(tcoord.units.origin, calendar="gregorian")
    tcoord.convert_units("days since 1661-01-01 00:00:00")
    tcoord.units = cf_units.Unit(
        tcoord.units.origin, calendar="proleptic_gregorian"
    )
    cube.remove_coord("time")
    cube.add_dim_coord(tcoord, 0) # might need to find this dimension

    try:
        iris.coord_categorisation.add_year(cube, 'time')
    except:
        pass
    
    try:               
        try:
            cube.remove_coord("month")
        except:
            pass
        iris.coord_categorisation.add_month_number(cube, 'time', name='month')
    except:
        pass

    try:
        del cube.attributes["history"]
    except:
        pass

    return(cube)


def insert_data_into_cube(x, eg_cube, mask =None):
    """ insert data into cube.
    Arguments:
        x -- np array of data that we want to insert into the cube. 
            Should have same shape of eg_cube, or same length as eg_cube
            or length equal to Trues in mask.
	    eg_cube -- The cube we want to insert data into.
	    mask -- Boolean array of shape or length x.
            Where True, will inster data.
            Defaulk of None which means True for all points in eg_cube.
    Returns:
        eg_cube with data replaced by x
    """

    Pred = eg_cube.copy()
    pred = Pred.data.copy().flatten()

    if mask is None:
        pred[:] = x
    else:
        pred[mask] = x

    Pred.data = pred.reshape(Pred.data.shape)
    return(Pred)


# /libs/namelist_functions.py
def read_variables_from_namelist(file_name):
    """Read variables from a file and create them with their original names.
    Inputs:
        file_name (str): The name of the file containing the variables.
    Returns:
        dict: A dictionary of variable names and their values.
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
        return eval(fun.split('function ')[1].split(' at ')[0])
    
    with open(file_name, 'r') as file:
        for line in file:
            parts = line.strip().split("::")
            if len(parts) == 2:
                variable_name = parts[0].strip()
                variable_value = parts[1].strip()
                if callable(variable_value):
                    # If the variable is a function, save its name
                    variable_value_set = variable_value
                elif variable_value.startswith('"') \
                    and variable_value.endswith('"'):
                    # If the variable is a string, remove the quotes
                    variable_value_set = variable_value[1:-1]
                elif variable_value.startswith('[') \
                    and variable_value.endswith(']'):
                    # If the variable is a list, parse it
                    try:        
                        variable_value_set = eval(variable_value)
                    except:
                        functions = variable_value.split(', ')
                        variable_value_set = [
                            define_function(fun) for fun in functions
                        ]
                else:
                    try:
                        # Try to parse the variable as a dictionary
                        variable_value_set = eval(variable_value)
                    except (SyntaxError, NameError):
                        # If parsing fails, assume it's a non-string,
                        # or non-list variable
                        variable_value_set = variable_value
                #set_trace()
                if variable_name in variables:
                    if type(variables[variable_name]) is list:
                        variables[variable_name].append(variable_value_set)
                    else:
                        variables[variable_name] = [
                            variables[variable_name], variable_value_set]
                else:
                    variables[variable_name] = variable_value_set
    return variables


def read_variable_from_namelist_with_overwite(file_name, **kwargs):

    def merge_variables(dict1):
        merged = dict1.copy()
        merged.update(**kwargs)
        return merged
    
    read_variables = read_variables_from_namelist(file_name)
    return merge_variables(read_variables)


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
        A = out.shape[0]
        B = out.shape[1]
        new_shape = ((A * B), *out.shape[2:])
        return np.reshape(out, new_shape)
    try:        
        trace = az.from_netcdf(trace, engine='netcdf4')
    except:
        pass
    params = trace.to_dict()['posterior']
    params_names = params.keys()
    params = [select_post_param_name(var) for var in params_names]
    return params, [var for var in params_names]


def construct_param_comb(i, params, params_names, extra_params):
    param_in = [
        param[i] if param.ndim == 1 else param[i,:] for param in params
    ]
    param_in = dict(zip(params_names, param_in))
    param_in.update(extra_params)
    return param_in


# /libs/read_variable_from_netcdf.py
def read_variable_from_netcdf(
        filename, dir='', subset_function=None,
        make_flat=False, units=None,
        subset_function_args=None,
        time_series=None, time_points=None, extent=None,
        return_time_points=False, return_extent=False,
        *args, **kw):
    """Read data from a netCDF file
        Assumes that the variables in the netcdf file all have the name
        "variable". Assumes that values < -9E9, you dont want. This could
        be different in some circumstances.
    Arguments:
        filename -- a string with filename or two element python list
            containing the name of the file and the target variable name.
            If just the sting of "filename" assumes variable name is "variable"
        dir -- directory file is in. Path can be in "filename" and None means
            no additional directory path needed.
        subset_function -- a function or list of functions to be applied to
            each data set.
        subset_function_args -- If subset_function is a function, dict
            arguments for that function. If subset_function is a list,
            a list or dict constaining arguments for those functions in turn.
        make_flat -- should the output variable to flattened or remain cube
        time_series -- list comtaining range of years. If making flat and
            returned a time series, checks if that time series contains year.
    Returns:
        Y - if make_flat, a numpy vector of the target variable, otherwise
        returns iris cube.
    """

    logger.info("Opening:")
    logger.info(filename)
    if filename[0] == '~' or filename[0] == '/' or filename[0] == '.':
        dir = ''
    try:
        if isinstance(filename, str):        
            dataset = iris.load_cube(dir + filename, callback=sort_time)
        else:
            dataset = iris.load_cube(
                dir + filename[0], filename[1], callback=sort_time
            )
    except:
        try:
            dataset = iris.load_cube(dir + filename)
        except:
            print("==============\nERROR!")
            print("can't open data.")
            print("Check directory (''" + dir + "''), filename (''" + \
                filename + "'') or file format")
            print("==============")
    coord_names = [coord.name() for coord in dataset.coords()]
    if dataset is None:
        return None
    if time_points is not None:     
        if 'time' in coord_names:
            dataset = dataset.interpolate(
                [('time', time_points)], iris.analysis.Linear()
            )
        else:   
            def addTime(time_point):
                time = iris.coords.DimCoord(
                    np.array([time_point]), standard_name='time',
                    units = 'days since 1661-01-01 00:00:00'
                )
                dataset_cp = dataset.copy()
                dataset_cp.add_aux_coord(time)
                return dataset_cp

            dataset_time = [addTime(time_point) for time_point in time_points]
            dataset = iris.cube.CubeList(dataset_time).merge_cube()
            dataset0 = dataset.copy()
    if extent is not None:
        dataset = dataset.regrid(extent, iris.analysis.Linear())
        dataset0 = dataset.copy()
    if units is not None:
        dataset.units = units
    if subset_function is not None:
        if isinstance(subset_function, list):
            for FUN, args in zip(subset_function, subset_function_args):
                try:    
                    dataset = FUN(dataset, **args)
                except:
                    print(
                        "Warning! function: " + FUN.__name__ + \
                        " not applied to file: " + dir + filename
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
        except:
            pass
            
        if time_series is not None:
            if not years[ 0] == time_series[0]:
                dataset = np.append(
                    np.repeat(np.nan, years[ 0]-time_series[0]), dataset
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
        y_filename, x_filename_list, CA_filename=None, add_1s_columne=False,
        y_threshold=None, x_normalise01=False, scalers=None, check_mask=True,
        frac_random_sample=1.0, min_data_points_for_sample=None, *args, **kw
    ):
                              
    """Read data from netCDF files 
        
    Arguments:
        y_filename -- a two element python list containing the name of the
            file and the target variable name.
        x_filename_list -- a python list of filename containing the feature
            variables.
        CA_filename -- a python list of filename containing the area of the
            cover type.
        y_threshold -- if converting y into boolean, the threshold we use to
            split into 0's and 1's.
        add_1s_columne -- useful for if using for regressions. Adds a variable
            of just 1's t rperesent y = SUM(a_i * x_i) + c
        x_normalise01 -- Boolean. If True, then X's are normalised between 0
            and 1.
        scalers -- None or numpy array of shape 2 by n. columns of X.
            Defines what scalers (min and max) to apply to each X column. 
            If None, doesn't apply anything.
        check_mask -- Boolean. If True, simple checks if there are any large
            negtaive numbers and makes them out. Assunes that values < -9E9,
            you dont want. This could be different in some circumstances.
        frac_random_sample -- fraction of data to be returned
        see read_variable_from_netcdf comments for *arg and **kw.
    Returns:
        Y - a numpy array of the target variable
        X - an n-D numpy array of the feature variables 
    """
    
    Y, time_points, extent = read_variable_from_netcdf(
        y_filename, make_flat=True, *args, 
        return_time_points=True, return_extent=True, **kw
    )
    
    if CA_filename is not None:
        CA = read_variable_from_netcdf(
            CA_filename, make_flat=True,
            time_points=time_points, extent=extent,
            *args, **kw
        )
   
    # Create a new categorical variable based on the threshold
    if y_threshold is not None:
        Y = np.where(Y >= y_threshold, 0, 1)
        # count number of 0 and 1 
        counts = np.bincount(Y)
    
    n = len(Y)
    m = len(x_filename_list)
    
    X = np.zeros([n,m])
    
    for i, filename in enumerate(x_filename_list):
        X[:, i] = read_variable_from_netcdf(
            filename, make_flat=True,
            time_points=time_points, extent=extent,
            *args, **kw
        )
    
    if add_1s_columne: 
        X = np.column_stack((X, np.ones(len(X)))) # add a column of ones to X 
    
    if check_mask:
        if CA_filename is not None:
            cells_we_want = np.array([
                np.all(rw > -9e9) and np.all(rw < 9e9) \
                    for rw in np.column_stack((X, Y, CA))
            ])
            CA = CA[cells_we_want]
        else:
            cells_we_want = np.array([
                np.all(rw > -9e9) and np.all(rw < 9e9) \
                    for rw in np.column_stack((X,Y))
            ])
        Y = Y[cells_we_want]
        X = X[cells_we_want, :]
        
    if x_normalise01 and scalers is None: 
        try:
            scalers = np.array([np.min(X, axis=0), np.max(X, axis=0)])
        except:
            pass
        squidge = (scalers[1,:]-scalers[0,:])/(X.shape[0])
        scalers[0,:] = scalers[0,:] - squidge
        scalers[1,:] = scalers[1,:] + squidge
        
        test = scalers[1,:] == scalers[0,:]
        scalers[0,test] = 0.0
        scalers[1,test] = 1.0

    if frac_random_sample is None: 
        frac_random_sample = 1000
    else:
        if min_data_points_for_sample is not None:
            min_data_frac = min_data_points_for_sample / len(Y)
            if min_data_frac > frac_random_sample:
                frac_random_sample = min_data_frac
    
    if frac_random_sample < 1:
        M = X.shape[0]
        selected_rows = np.random.choice(
            M, size=int(M * frac_random_sample), replace=False
        )
        Y = Y[selected_rows]
        X = X[selected_rows, :]
        if CA_filename is not None:
            CA = CA[selected_rows]
    
    if scalers is not None:
        X = (X - scalers[0, :]) / (scalers[1, :] - scalers[0, :])
        if check_mask: 
            if CA_filename is not None:
                return Y, X, CA, cells_we_want, scalers
        return Y, X, cells_we_want, scalers
        
    if check_mask or frac_random_sample: 
        if CA_filename is not None:
            return Y, X, CA, cells_we_want
    
    return Y, X, cells_we_want


# /fire_models/ConFire.py
class ConFire(object):
    def __init__(self, params, inference=False):
        """
        Initalise parameters and calculates the key variables needed to
        calculate burnt area.
        """
        self.inference = inference
        if self.inference:
            self.numPCK =  __import__('pytensor').tensor
        else:
            self.numPCK =  __import__('numpy')
        
        self.params = params

        def select_param_or_default(*args, **kw):
            return select_key_or_default(
                self.params, numPCK=self.numPCK, *args, **kw
            ) 

        self.controlID = self.params['controlID']
        self.control_Direction = self.params['control_Direction']
        self.x0 = select_param_or_default('x0', [0])
        self.log_control = select_param_or_default(
            'log_control', [False] * len(self.control_Direction)
        )
        
        self.betas = select_param_or_default('betas', [[0]], stack=False)
        self.powers = select_param_or_default('powers', None, stack=False)
        self.driver_Direction = self.params['driver_Direction']
        self.Fmax = select_param_or_default('Fmax', None, stack=False)

    def burnt_area(self, X, return_controls=False, return_limitations=False):
        ## finds controls        
        def cal_control(cid=0):
            ids = self.controlID[cid]
            betas =  self.betas[cid] * self.driver_Direction[cid]
            
            X_i = X[:,ids]
            if self.powers is not None:
                X_i = self.numPCK.power(X_i, self.powers[cid])
            
            out = self.numPCK.sum(X_i * betas[None, ...], axis=-1)
            if self.log_control[cid]:
                out = self.numPCK.log(out)
            out = out + self.x0[cid]
            return(out)
            
        
        controls = [cal_control(i) for i in range(len(self.controlID))]
        if return_controls:
            return controls

        def sigmoid(y, k):
            if k == 0:
                return None
            return 1.0 / (1.0 + self.numPCK.exp(-y * k))
        
        
        limitations = [
            sigmoid(y, k) for y, k in zip(controls, self.control_Direction)
        ]
        
        if return_limitations:
            return limitations

        limitations = [lim for lim in limitations if lim is not None]
        
        
        BA =  self.numPCK.prod(limitations, axis=0)
        if self.Fmax is not None:
            sigmoid(self.Fmax, 1.0) * BA

        return BA
    
    
    def emc_weighted(self, emc, precip, wd_pg):
        
        try:
            wet_days = 1.0 - self.numPCK.exp(-wd_pg * precip)
            emcw = (1.0 - wet_days) * emc + wet_days
        except:
            emcw = emc.copy()
            emcw.data  = 1.0 - self.numPCK.exp(-wd_pg * precip.data)
            emcw.data = emcw.data + (1.0 - emcw.data) * emc.data
        return(emcw)

    def list_model_params(self, params, varnames=None):
        controlID = params[0]['controlID']
        def list_one_line_of_parmas(param):
            def select_param_or_default(*args, **kw):
                return select_key_or_default(
                    param, numPCK=__import__('numpy'), *args, **kw
                )
            
            control_Direction = param['control_Direction']
            x0s = select_param_or_default('x0', [0])
            betas = select_param_or_default('betas', [[0]], stack=False)
            powers = select_param_or_default('powers', None, stack=False)
            driver_Direction = param['driver_Direction']
            Fmax = select_param_or_default('Fmax', 1.0, stack=False)

            directions = [
                np.array(driver)*control for control, driver in \
                zip(control_Direction, driver_Direction)
            ]
            betas = [
                direction * beta for direction, beta in zip(directions, betas)
            ]

            def mish_mash(beta, power, x0):
                return np.append(np.column_stack((beta, power)).ravel(), x0)
            varps = [
                mish_mash(beta, power, x0) for beta, power, x0 in zip(
                    betas, powers, x0s)
            ]
            
            return np.append(Fmax, np.concatenate(varps))
            
        params_sorted = np.array([
            list_one_line_of_parmas(param) for param in params
        ])
    
        param = ['Fmax']
        controlN = ['']
        variable_name = ['']
        
        for IDs, IDn in zip(controlID, range(len(controlID))):
            for ID in IDs:
                param.append('beta')
                controlN.append(str(IDn))
                variable_name.append(varnames[ID])
                param.append('power')
                controlN.append(str(IDn))
                variable_name.append(varnames[ID])
            param.append('beta')
            controlN.append(str(IDn))
            variable_name.append('beta0')

        header_df = pd.DataFrame([param, controlN, variable_name])   
        index_labels = ["Parameter", "Control", "Variable"] + \
            list(map(str, range(1, params_sorted.shape[0] + 1)))

        data_df = pd.DataFrame(params_sorted)
        full_df = pd.concat([header_df, data_df], ignore_index=True) 

        full_df.index = index_labels

        return full_df


def diagnostic_run_ConFire(config, model_name, timerange):
    # --------------------------------------------------------
    # This script runs the ConFire model using pre-generated parameter files.
    # It calculates burnt area under different control conditions and saves
    # results.
    # --------------------------------------------------------

    work_dir = config['work_dir']
    diag_dir = config['diag_dir']
    confire_param = config['confire_param_dir']

    # **Define Paths for Parameter Files  and for outputs**
    output_dir = work_dir + '/ConFire_outputs/'

    # Parameter files (traces, scalers, and other model parameters)
    param_file_trace = glob.glob(confire_param + "trace*.nc")[0]
    param_file_none_trace = glob.glob(
        confire_param + "none_trace-params*.txt")[0]
    scale_file = glob.glob(confire_param + "scalers*.csv")[0]

    # Number of samples to run from the trace file
    nsample_for_running = 100

    # **Load Variable Information and NetCDF Files**
    # Replace these lines with user-specified NetCDF files, ensuring variable
    # order is the same.
    nc_files = config['files_input']
    nc_dir = ''

    # **Load Driving Data and Land Mask**
    logger.info('Loading data for ConFire model...')
    scalers = pd.read_csv(scale_file).values
    obs_data, driving_data, lmask, scalers = read_all_data_from_netcdf(
        nc_files[0], nc_files, scalers=scalers, dir=nc_dir)

    # Load a sample cube (used for inserting data)
    eg_cube = read_variable_from_netcdf(nc_files[0], dir=nc_dir)

    # **Extract Model Parameters**
    logger.info('Loading ConFire model parameters...')
    params, params_names = select_post_param(param_file_trace)
    extra_params = read_variables_from_namelist(param_file_none_trace)
    control_direction = extra_params['control_Direction'].copy()
    Nexp = len(control_direction)

    # **Sample Iterations from Trace**
    nits = len(params[0])
    idx = range(0, nits, int(np.floor(nits/nsample_for_running)))

    # **Storage for Model Outputs (Including full model, physical &
    # stochastic controls)**
    out_cubes = [[] for _ in range(Nexp + 2)]

    def run_model_into_cube(param_in, coord):
        """
        Runs the ConFire model with given parameters and inserts results into
        an Iris cube.
        """
        out = ConFire(param_in).burnt_area(driving_data)
        cube = insert_data_into_cube(out, eg_cube, lmask)
        cube.add_aux_coord(coord)
        return cube

    # **Run ConFire Model with Different Control Scenarios**
    logger.info('Running ConFire model...')
    for id, i in zip(idx, range(len(idx))):
        coord = iris.coords.DimCoord(i, "realization")
        param_in = construct_param_comb(
            id, params, params_names, extra_params
        )
        
        # **Run Full Model**
        out_cubes[0].append(run_model_into_cube(param_in, coord))

        # **Run Model with Individual Controls Turned On**
        for cn in range(Nexp):
            param_in['control_Direction'][:] = [0] * Nexp
            param_in['control_Direction'][cn] = control_direction[cn]
            param_in = construct_param_comb(
                cn, params, params_names, extra_params
            )
            out_cubes[cn+1].append(run_model_into_cube(param_in, coord))

        # **Run Model with All Controls Off (Stochastic Control)**
        param_in['control_Direction'][:] = [0] * Nexp
        out_cubes[cn+2].append(run_model_into_cube(param_in, coord))

    # **Save Output Cubes**
    logger.info('Saving ConFire output cubes...')
    filenames_out = ['burnt_area_model'] + \
        ['burnt_area_control_' + str(i) for i in range(Nexp)] + \
        ['burnt_area_control_stochastic']
    os.makedirs(output_dir, exist_ok=True)

    for i in range(len(out_cubes)):
        cubes = iris.cube.CubeList(out_cubes[i]).merge_cube()
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
        filepath = os.path.join(output_dir, filename + ".nc")

        fig, axes = plt.subplots(
            nrows=1, ncols=2, figsize=(12, 8),
            subplot_kw={'projection': ccrs.PlateCarree()})
        axes = axes.flatten()
        plotN = 0
        
        try:
            cube = iris.load_cube(filepath)   # Load saved NetCDF file
            annual_avg = cube.collapsed('time', iris.analysis.MEAN)
            for pc in [5, 95]:
                pc_cube = 100 * cube.collapsed(
                    'realization', iris.analysis.PERCENTILE, percent=pc)[0]
                ax = axes[plotN]
                img = iris.quickplot.pcolormesh(
                    pc_cube, axes=ax, colorbar=False,
                    vmin=0., vmax=100., cmap='Oranges')
                ax.set_title(
                    filename.replace("_", " ").capitalize() + \
                    f' - {timerange}\n{model_name} - [{str(pc)}% percentile]'
                )
                ax.coastlines()
                # Add colorbar
                fig.colorbar(
                    img, ax=ax, orientation='vertical', extend='neither',
                    label='Burnt area [%]',
                    fraction=0.046, pad=0.04, shrink=0.3,
                )
                # Iterate plot number
                plotN = plotN + 1
            figures.append(fig)
        
        except Exception as e:
            logger.info(f"Skipping {filename} due to error: {e}")
            logger.debug(f"Skipping {filename} due to error: {e}")

    return figures
