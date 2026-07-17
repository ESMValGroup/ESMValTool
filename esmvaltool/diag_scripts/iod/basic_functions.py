import logging
from pathlib import Path
import iris # type: ignore
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import scipy # type: ignore
import os

from esmvaltool.diag_scripts.shared import select_metadata, save_data  # type: ignore


logger = logging.getLogger(Path(__file__).stem)
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)

def get_provenance_record(caption, ancestor_files):
    """create a provenance record describing the diagnostic data and plot."""
    # associated recipe uses contains a caption string with placeholders
    # like {long_name} that are now populated from attributes dictionary.
    # note that for simple recipes, caption can be set here as a simple string
    
    logger.debug("Creating provenance record.")
    record = {
        'caption': caption,
        'ancestors': ancestor_files,
        'authors': ['righi_mattia']
    }
    return record

def load_data(group_md, variable_group, get_filenames=None):
    """
    Load data for a given variable group and return a dictionary with dataset names as keys.
    """
    logger.info(f"Loading data for {variable_group}...")

    # I had been using this but I can't seem to get a single dict entry out of the list
    # selected_meta = select_metadata(group_md, variable_group=variable_group)
    var_dict = next((item for item in group_md if item["variable_group"] == variable_group), None)
    if var_dict == None:
        return None
    cube_dict = {}

    dataset = var_dict['dataset']
    filename = var_dict['filename']
    cube = iris.load_cube(filename)
    new_cube = iris.util.squeeze(cube)

    if get_filenames:
        cube_dict[dataset] = {'cube': new_cube, 'filename': filename}
    else:
        cube_dict[dataset] = new_cube

    return cube_dict

def to_set(value):
    return set(value) if isinstance(value, list) else {value}

def get_prefix(file_or_files):
    files = file_or_files if isinstance(file_or_files, list) else [file_or_files]
    return os.path.basename(files[0]).split("_")[0]

def compute_cube_diff(cfg, dict1, dict2, output_basename, skew=None):
    """
    Compute the difference between two cubes or skewness of the result and return a dictionary of results.
    1 minus 2.
    """
    diff_results = {}

    (dataset1, info1), = dict1.items()
    (dataset2, info2), = dict2.items()

    if dataset1 != dataset2:
        logger.warning(f"Datasets do not match: {dataset1} vs {dataset2}. Skipping cube_diff.")
        return None
    
    # Compute DMI for each dataset

    cube1 = info1['cube']
    cube2 = info2['cube']
    files = list((info1['filename'], info2['filename']))
    diff_cube = cube1 - cube2
    logger.info(f"Computed difference for {dataset1}")

    if skew:
        skew_value = scipy.stats.skew(diff_cube.data, bias=False)
        logger.info(f"Computed skewness for {dataset1}")

        # Create a scalar cube to hold the skewness value
        skew_cube = iris.cube.Cube(
            np.array(skew_value),
            long_name="Skewness of difference",
            var_name="skewness",
            units="1",
            attributes=diff_cube.attributes)

        result_cube = skew_cube
    else:
        result_cube = diff_cube

    diff_results[dataset1] = {
        'cube': result_cube,
        'filename': files
    }

    # Save output
    save_prefix = get_prefix(files)
    save_string = f"{save_prefix}_{output_basename}_{dataset1}"
    provenance_record = get_provenance_record(save_string, files)
    save_data(save_string, provenance_record, cfg, result_cube)

    return diff_results

def iso_depth_3d(CT_cube, t):
    """
    Compute the depth of a given isotherm for a 3D Iris cube.

    Parameters
    ----------
    CT_cube : iris.cube.Cube
        Conservative Temperature cube with dimensions (depth, latitude, longitude)
    p : array_like
        1D array of depth levels (must match the depth coordinate of CT_cube)
    t : float
        Temperature of the isotherm

    Returns
    -------
    iso_depth_array : np.ndarray
        2D array (latitude, longitude) of isotherm depths
    """

    # Extract data and mask
    CT = CT_cube.data  # shape: (depth, lat, lon)
    mask = ma.getmask(CT)

    p_coord = CT_cube.coord('depth')
    p = p_coord.points

    # Apply mask to pressure array
    p = np.broadcast_to(p[:, np.newaxis, np.newaxis], CT.shape)  # Expand depth to match cube shape
    p = ma.masked_array(p, mask=mask)

    # Create output array for isotherm depth
    iso_depth_array = np.full(CT.shape[1:], np.nan)  # Shape (lat, lon)

    # Loop through each latitude and longitude point
    for j in range(CT.shape[1]):  # Latitude
        for i in range(CT.shape[2]):  # Longitude
            # Extract 1D vertical profile at (lat, lon)
            CT_profile = CT[:, j, i]
            p_profile = p[:, j, i]

            # Check for missing data
            if ma.is_masked(CT_profile) or ma.is_masked(p_profile) or np.all(CT_profile.mask):
                continue  # Skip if all values are masked

            # Find the index where CT crosses the isotherm t
            idx_mld = CT_profile > t

            if not np.any(idx_mld):
                continue  # Skip if no crossing found

            idx_max = p_profile[idx_mld].argmax()

            # Ensure valid interpolation index
            if idx_max >= len(p_profile) - 1:
                continue  # Skip if no valid interpolation

            # Linear interpolation between depth levels
            len_inter = ((CT_profile[idx_max] - t) / 
                         (CT_profile[idx_max] - CT_profile[idx_max + 1])) * (p_profile[idx_max + 1] - p_profile[idx_max])

            # Compute isotherm depth
            iso_depth_array[j, i] = p_profile[idx_max] + len_inter

    iso_depth_array = np.nan_to_num(iso_depth_array, nan=1e+20)

    # Create new cube with same lat/lon coordinates
    if np.size(CT_cube.coord('latitude').shape) == 2:
        iso_cube = iris.cube.Cube(
            iso_depth_array,
            long_name=f'Isotherm Depth ({t} $^\\circ$C)',
            units=p_coord.units,
            aux_coords_and_dims=[(CT_cube.coord('latitude'), (0,1)), (CT_cube.coord('longitude'), (0,1))]
        )
    else:
        iso_cube = iris.cube.Cube(
            iso_depth_array,
            long_name=f'Isotherm Depth ({t} $^\\circ$C)',
            units=p_coord.units,
            aux_coords_and_dims=[(CT_cube.coord('latitude'), 0), (CT_cube.coord('longitude'), 1)]
        )

    return iso_cube

def iso_depth_4d(CT_cube, t, time_measure):
    """
    Compute the depth of a given isotherm for a 4D Iris cube with time.

    Parameters
    ----------
    CT_cube : iris.cube.Cube
        Conservative Temperature cube with dimensions (time, depth, latitude, longitude).
    t : float
        Temperature of the isotherm.
    time_measure: str
        'season_number' or 'month_number'

    Returns
    -------
    iso_cube : iris.cube.Cube
        3D cube (time, latitude, longitude) of isotherm depths.
    """

    # Extract data and mask
    CT = CT_cube.data  # shape: (time, depth, lat, lon)
    mask = ma.getmask(CT)

    p_coord = CT_cube.coord('depth')
    p = p_coord.points  # shape: (depth,)

    # Expand depth array to match CT shape: (time, depth, lat, lon)
    p = np.broadcast_to(p[np.newaxis, :, np.newaxis, np.newaxis], CT.shape)
    p = ma.masked_array(p, mask=mask)

    # Create output array for isotherm depth (time, lat, lon)
    iso_depth_array = np.full(CT[:,0,:,:].shape, np.nan)  # (time, lat, lon)

    # Loop through each time, latitude, and longitude point
    for t_idx in range(CT.shape[0]):  # Time loop
        for j in range(CT.shape[2]):  # Latitude loop
            for i in range(CT.shape[3]):  # Longitude loop
                # Extract vertical profile at (time, lat, lon)
                CT_profile = CT[t_idx, :, j, i]
                p_profile = p[t_idx, :, j, i]

                # Check for missing data
                if ma.is_masked(CT_profile) or ma.is_masked(p_profile) or np.all(CT_profile.mask):
                    continue  # Skip if all values are masked

                # Find the index where CT crosses the isotherm t
                idx_mld = CT_profile > t
                if not np.any(idx_mld):
                    continue  # Skip if no crossing found

                idx_max = p_profile[idx_mld].argmax()

                # Ensure valid interpolation index
                if idx_max >= len(p_profile) - 1:
                    continue  # Skip if no valid interpolation

                # Linear interpolation between depth levels
                len_inter = ((CT_profile[idx_max] - t) /
                             (CT_profile[idx_max] - CT_profile[idx_max + 1])) * (p_profile[idx_max + 1] - p_profile[idx_max])

                # Compute isotherm depth
                iso_depth_array[t_idx, j, i] = p_profile[idx_max] + len_inter

    # Create new cube with the same lat/lon and time coordinates
    iso_cube = iris.cube.Cube(
        iso_depth_array,
        long_name=f'Isotherm Depth ({t} $^\circ$C)',
        units=p_coord.units,
        dim_coords_and_dims=[(CT_cube.coord(time_measure), 0)],
        aux_coords_and_dims=[(CT_cube.coord('latitude'), (1,2)), (CT_cube.coord('longitude'), (1,2))]
    )

    return iso_cube