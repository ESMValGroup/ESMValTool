"""
CMORizer for EN4 dataset.

This script processes EN4 ocean temperature and salinity data to CMOR-compliant format for use in ESMValTool.

Tier
    Tier 2: other freely-available dataset.

Source
    https://www.metoffice.gov.uk/hadobs/en4/download-en4-2-2.html
    
Last access
    2025-06-13

Info
    EN4: quality controlled subsurface ocean temperature and salinity objective analyses.
    Script tested using analyses with Gouretski and Reseghetti (2010) XBT corrections
    and Gouretski and Cheng (2020) MBT corrections applied.
    
Download instructions      
    - Edit the text file for your chosen years, https://www.metoffice.gov.uk/hadobs/en4/EN.4.2.2.analyses.g10.download-list.txt
    - Save .txt file in directory for data to be downloaded to.
    - Run 'wget -i EN.4.2.2.profiles.g10.download-list.txt' in same directory.
    - Unzip files prior to running the cmorizer script.

"""
import logging
from pathlib import Path

import iris

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)

def load_and_prepare_cube(fullpath, var, var_info, glob_attrs, cmor_table):
    """
    Load and prepare a data cube for CMORization. Fix attributes, coordinates, and units.

    Parameters
    ----------
    fullpath : str
        Path to the input file.
    var : str
        Name of the variable to save in output cube.
    var_info : dict
        Variable information dictionary.
    glob_attrs : dict
        Global attributes to set on the cube.
    cmor_table : object
        CMOR table.

    Returns
    -------
    iris.cube.Cube
        The prepared data cube.
    """

    glob_attrs["mip"] = var_info["mip"]
    raw_var = var_info["raw_var"]
    cmor_info = cmor_table.get_variable(var_info["mip"], var)

    cubes = iris.load(fullpath, raw_var)
    iris.util.equalise_attributes(cubes)
    cube = cubes.concatenate_cube()

    if cube.units == "K":
        cube.convert_units("degC")

    cube.coord("depth").units = "m"   
    cube = utils.fix_coords(cube)
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, glob_attrs)

    return cube

def extract_surface_var(cube, cmor_info):
    """
    Extract the surface level variable from a data cube.

    Parameters
    ----------
    cube : iris.cube.Cube
        Input data cube.
    cmor_info : object
        CMOR table object for the surface variable.

    Returns
    -------
    iris.cube.Cube
        The extracted surface level cube.
    """
    logger.info("Extracting surface level")

    depth0 = iris.Constraint(depth=cube.coord("depth").points[0])
    surface_cube = cube.extract(depth0)

    utils.fix_var_metadata(surface_cube, cmor_info)

    return surface_cube

def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """
    CMORization main function call.
    """
    cmor_table = cfg["cmor_table"]
    glob_attrs = cfg["attributes"]
    fullpath = str(Path(in_dir) / cfg["filename"])

    for var, var_info in cfg["variables"].items():
        logger.info("Loading %s", fullpath)
        
        srf_var = var_info["srf_var"]
        cmor_info_srf = cmor_table.get_variable(var_info["mip"], srf_var)

        cube = load_and_prepare_cube(fullpath, var, var_info, glob_attrs, cmor_table)
        surface_cube = extract_surface_var(cube, cmor_info_srf)
        logger.info("Saving for %s", var)
        utils.save_variable(cube, var, out_dir, glob_attrs)

        logger.info("Saving for %s", srf_var)
        utils.save_variable(surface_cube, srf_var, out_dir, glob_attrs)

    
