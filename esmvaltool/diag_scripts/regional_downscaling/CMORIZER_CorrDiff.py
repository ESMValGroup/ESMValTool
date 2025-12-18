import netCDF4 as nc
import xarray as xr
import numpy as np
from datetime import datetime
import os
import iris
from ncdata.iris_xarray import cubes_from_xarray, cubes_to_xarray
from esmvaltool.cmorizers.data import utilities
import warnings
warnings.filterwarnings('ignore')
# import cftime  # optional, only if you need non-standard calendars

# /!\ IMPORTANT
# Define time frequency and lat x lon bounds of the domain
time_freq_in_hours = 3
lon_values = np.linspace(1, 19, 320)
lat_values = np.linspace(42, 60, 320)
# Define output folder where CMORized data is saved. Add this path to your ESMValTool config file for OBS6
cmor_output_path = '/work/...data/ml_model_cmor_ouput/'
# Absolute path to CorrDiff/CorrDiff++ output"
input_file = "output_CorrPlusPlus_multi_logpr_UNet_025lossweight_xcond_ConFIG_4M_test.nc"
# Name you want to give to the ML-based dataset
method_name = "CorrDiff++"
# List of methods to keep from CorrDiff original output: can only contain "truth" and "prediction"
# If "truth" is in the list it will create a HIGHRES-REF dataset
method_to_keep = ["prediction"]

def fix_time_values(time_vals, freq_hours=3):
    """
    Reconstructs a regularly spaced time coordinate (e.g. 3-hourly)
    from the first valid time entry in a dataset.
    
    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with a possibly corrupted 'time' variable.
    freq_hours : int, optional
        Desired time step in hours (default: 3).
    
    Returns
    -------
    np.ndarray
        Fixed array of time values (same length as ds.time).
    """
    n_times = len(time_vals)

    # Extract the first time value (the correct one)
    first_time = time_vals[0]

    # Create fixed time sequence in hours since the first timestamp
    fixed_time = first_time + np.arange(n_times) * freq_hours

    return fixed_time

def create_cmor_compliant_files(input_file, output_dir='cmor_output', method_name="CORRDIFF", groups_to_keep=["truth","prediction"]):
    """
    Convert NetCDF file to CMOR-compliant format.
    Creates separate files for each variable, group, and ensemble member.
    TO DO: Integrate more existing CMORIZER methods from ESMValCore
    """
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Open the original file
    ds_orig = nc.Dataset(input_file, 'r')
    
    # Read common coordinates
    lat = lat_values
    lon = lon_values
    
    time_vals = ds_orig.variables['time'][:]
    time_vals = fix_time_values(time_vals, freq_hours=time_freq_in_hours)
    
    # Variable name mapping to CMOR standard names
    # Add variables here if necessary or modify the code to use ESMValCore utilities
    var_mapping = {
        'uas': {'standard_name': 'eastward_wind', 'long_name': 'Eastward Near-Surface Wind', 'units': 'm s-1'},
        'vas': {'standard_name': 'northward_wind', 'long_name': 'Northward Near-Surface Wind', 'units': 'm s-1'},
        'tas': {'standard_name': 'air_temperature', 'long_name': 'Near-Surface Air Temperature', 'units': 'K'},
        'logpr': {'standard_name': 'log_of_precipitation_flux', 'long_name': 'Log of Precipitation', 'units': 'kg m-2 s-1'},
        '2rootpr': {'standard_name': 'square_root_of_precipitation_flux', 'long_name': 'Square Root of Precipitation', 'units': 'kg m-2 s-1'},
        'pr': {'standard_name': 'precipitation_flux', 'long_name': 'Precipitation', 'units': 'kg m-2 s-1'},
        'hus_sfc': {'standard_name': 'specific_humidity', 'long_name': 'Near-Surface Specific Humidity', 'units': '1'},
        'pres_sfc': {'standard_name': 'surface_air_pressure', 'long_name': 'Surface Air Pressure', 'units': 'Pa'},
    }

    rename_mapping = {
        '2rootpr': "pr",
        'logpr': "pr",
        'hus_sfc': "huss",
        'pres_sfc': "ps",
    }

    # Process each group
    for group_name in groups_to_keep:
        print(f"\nProcessing group: {group_name}")
        group = ds_orig.groups[group_name]
        
        # Get variables in this group
        for var_name in group.variables.keys():
            if var_name in var_mapping:
                print(f"  Processing variable: {var_name}")
                
                # Handle the case of square root of precipitation or log pr: transform it back
                if var_name == "2rootpr":
                    original_var_name = var_name
                    var_name = "pr"
                    print(f"  Renamed variable to: {var_name}")
                    var_data = np.power(np.clip(group.variables["2rootpr"][:],0, None), 2)
                    var_info = var_mapping[var_name]
                elif var_name == "logpr":
                    original_var_name = var_name
                    var_name = "pr"
                    print(f"  Renamed variable to: {var_name}")
                    var_data = np.clip(np.exp(group.variables["logpr"][:])- 10e-6, 0, None)
                    var_info = var_mapping[var_name]
                else:
                    var_data = group.variables[var_name][:]
                    var_info = var_mapping[var_name]
                    if var_name in rename_mapping: 
                        var_name = rename_mapping[var_name]

                # For prediction group with ensembles
                if group_name == 'prediction' and len(var_data.shape) == 4:  # (ensemble, time, y, x)
                    num_ensembles = var_data.shape[0]
                    
                    # Create separate file for each ensemble member
                    for ens_idx in range(num_ensembles):
                        dataset = method_name
                        filename = f"OBS6_{dataset}_reanaly_{ens_idx+1}_3hr_{var_name}.nc"
                        
                        # Create proper directory structure
                        tier_dir = os.path.join(output_dir, "Tier1", dataset)
                        os.makedirs(tier_dir, exist_ok=True)
                        filepath = os.path.join(tier_dir, filename)
                        
                        create_cmor_file(
                            filepath=filepath,
                            var_name=var_name,
                            var_data=var_data[ens_idx, :, :, :],  # (time, y, x)
                            lat=lat,
                            lon=lon,
                            time_vals=time_vals,
                            var_info=var_info,
                            group_name=group_name,
                            ensemble_member=ens_idx + 1
                        )
                        print(f"    Created: {filepath}")
                
                # For truth group (no ensemble dimension)
                else:
                    dataset = "HIGHRES-REF"
                    filename = f"OBS6_{dataset}_reanaly_1_3hr_{var_name}.nc"
                    
                    # Create proper directory structure
                    tier_dir = os.path.join(output_dir, "Tier1", dataset)
                    os.makedirs(tier_dir, exist_ok=True)
                    filepath = os.path.join(tier_dir, filename)
                    
                    create_cmor_file(
                        filepath=filepath,
                        var_name=var_name,
                        var_data=var_data,  # (time, y, x)
                        lat=lat,
                        lon=lon,
                        time_vals=time_vals,
                        var_info=var_info,
                        group_name=group_name,
                        ensemble_member=None
                    )
                    print(f"    Created: {filepath}")
    
    ds_orig.close()
    print(f"\nCMOR-compliant files created in: {output_dir}")

def create_cmor_file(filepath, var_name, var_data, lat, lon, time_vals, 
                     var_info, group_name, ensemble_member=None):
    """
    Create a single CMOR-compliant NetCDF file with proper height coordinates.
    """
    
    # Create xarray dataset for easier metadata handling
    ds = xr.Dataset()
    #Fix for CorrDiff to rotate lat lon (and be align with genfocal)
    var_data = np.rot90(var_data, k=-1, axes=(1, 2))

    # Determine if this variable needs a height coordinate
    needs_height2m = var_name in ['tas', 'huss']  # 2m height variables
    needs_height10m = var_name in ['uas', 'vas']   # 10m height variables
    
    dims = ['time', 'lon', 'lat']
    var_data_with_height = var_data
    
    # Add variable data
    ds[var_name] = xr.DataArray(
        var_data_with_height,
        dims=dims,
        attrs={
            'standard_name': var_info['standard_name'],
            'long_name': var_info['long_name'],
            'units': var_info['units'],
        }
    )
    
    ds.attrs['coordinates'] = 'lon lat'
    ds[var_name].attrs['coordinates'] = 'lon lat'
    
    # Time coordinate
    ds['time'] = xr.DataArray(
        time_vals,
        dims=['time'],
        attrs={
            'standard_name': 'time',
            'long_name': 'time',
            'axis': 'T',
            'calendar': 'standard',
            'units': 'hours since 1990-01-01 00:00:00'
        }
    )

    # Longitude (X)
    ds['lon'] = xr.DataArray(
        lon,
        dims=['lon'],
        attrs={
            'standard_name': 'longitude',
            'long_name': 'longitude',
            'units': 'degrees_east',
            'axis': 'X',
        }
    )

    # Latitude (Y)
    ds['lat'] = xr.DataArray(
        lat,
        dims=['lat'],
        attrs={
            'standard_name': 'latitude',
            'long_name': 'latitude',
            'units': 'degrees_north',
            'axis': 'Y',
        }
    )

    # Add global attributes (CMOR-compliant)
    ds.attrs.update({
        'Conventions': 'CF-1.8',
        'title': f'{var_info["long_name"]} - {group_name}',
        'institution': 'DLR',
        'source': 'CorrDiff generative model',
        'history': f'Created on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}',
        'references': 'https://doi.org/your-reference',
        'comment': f'Data from {group_name} group',
        'frequency': '3hrPt',
        'realm': 'atmos',
        'product': 'model-output',
        'variable_id': var_name,
        'grid': 'native grid',
        'grid_label': 'gn',
        'nominal_resolution': '~6 km',
    })
    
    # Add ensemble-specific attributes
    if ensemble_member is not None:
        ds.attrs['realization_index'] = ensemble_member
        ds.attrs['variant_label'] = f'r{ensemble_member}i1p1f1'
        ds.attrs['ensemble_member'] = ensemble_member

    cube = cubes_from_xarray(ds)[0]
    if needs_height2m:
        utilities.add_height2m(cube)
    elif needs_height10m:
        utilities.add_height10m(cube)
    fixed_cube = utilities.fix_coords(cube)
    # fixed_cube= utilities.fix_var_metadata(fixed_cube, var_info)
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    iris.save(fixed_cube, filepath, fill_value=1e20, unlimited_dimensions=["time"])

# CMORize CorrDiff++ output
print("Creating CMOR-compliant files...")
create_cmor_compliant_files(input_file, output_dir=cmor_output_path, method_name=method_name, groups_to_keep=["prediction"])
print("\nDone!")