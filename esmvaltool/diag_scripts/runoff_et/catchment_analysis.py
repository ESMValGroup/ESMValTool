"""
Look at this module for guidance how to write your own.

Module for personal diagnostics (example).
Internal imports from exmvaltool work e.g.:

from esmvaltool.diag_scripts.shared._supermeans import get_supermean

"""
import os
import logging
import pdb

import matplotlib
# use this everytime you import matplotlib
# modules; some machines dont have graphical interface (X)
matplotlib.use('Agg')  # noqa

import iris
import matplotlib.pyplot as plt

import esmvaltool.diag_scripts.shared as diag
from esmvaltool.preprocessor._regrid  import regrid
from esmvaltool.preprocessor._area_pp import area_average

logger = logging.getLogger(os.path.basename(__file__))

class defaults(object):
    """Class containing default dictionaries for predefined catchments

    The properties are used in the routine analysecatchments.
    Properties are
        catchments
        runoffrefdata
        preciprefdata
        ETrefdata
    """
    catchments = {
        # Catchments with name as used in REFFILE as key and the
        # catchment number as used in pcatchment as value
        "Amazon": 94,
        "Parana": 98,
        "Mackenzie": 76,
        "Mississippi": 86,
        "Danube": 14,
        "Congo": 68,
        "Niger": 65,
        "Nile": 60,
        "Lena": 40,
        "Yangtze-Kiang": 52,
        "Ganges-Brahmaputra": 54,
        "Murray": 100
        }

    runoffrefdata = {
        'Amazon': {'data': 1195.4477, 'unit': 'mm a-1'},
        'Congo': {'data': 365.6980, 'unit': 'mm a-1'},
        'Danube': {'data': 250.9211, 'unit': 'mm a-1'},
        'Ganges-Brahmaputra': {'data': 672.5738, 'unit': 'mm a-1'},
        'Lena': {'data': 197.3081, 'unit': 'mm a-1'},
        'Mackenzie': {'data': 173.9881, 'unit': 'mm a-1'},
        'Mississippi': {'data': 182.2420, 'unit': 'mm a-1'},
        'Murray': {'data': 8.2041, 'unit': 'mm a-1'},
        'Niger': {'data': 31.5160, 'unit': 'mm a-1'},
        'Nile': {'data': 48.7528, 'unit': 'mm a-1'},
        'Parana': {'data': 203.0060, 'unit': 'mm a-1'},
        'Yangtze-Kiang': {'data': 531.6936, 'unit': 'mm a-1'}
        }

    preciprefdata = {
        'Amazon': {'data': 2253.61, 'unit': 'mm a-1'},
        'Congo': {'data': 1539.98, 'unit': 'mm a-1'},
        'Danube': {'data': 809.11, 'unit': 'mm a-1'},
        'Ganges-Brahmaputra': {'data': 1387.95, 'unit': 'mm a-1'},
        'Lena': {'data': 399.146, 'unit': 'mm a-1'},
        'Mackenzie': {'data': 445.342, 'unit': 'mm a-1'},
        'Mississippi': {'data': 890.034, 'unit': 'mm a-1'},
        'Murray': {'data': 530.441, 'unit': 'mm a-1'},
        'Niger': {'data': 436.907, 'unit': 'mm a-1'},
        'Nile': {'data': 673.565, 'unit': 'mm a-1'},
        'Parana': {'data': 1311.22, 'unit': 'mm a-1'},
        'Yangtze-Kiang': {'data': 1032.84, 'unit': 'mm a-1'}
        }

    ETrefdata = {
        'Amazon': {'data': 1014.4023, 'unit': 'mm a-1'},
        'Congo': {'data': 1203.182, 'unit': 'mm a-1'},
        'Danube': {'data': 554.5999, 'unit': 'mm a-1'},
        'Ganges-Brahmaputra': {'data': 722.5962, 'unit': 'mm a-1'},
        'Lena': {'data': 187.4469, 'unit': 'mm a-1'},
        'Mackenzie': {'data': 269.2429, 'unit': 'mm a-1'},
        'Mississippi': {'data': 712.192, 'unit': 'mm a-1'},
        'Murray': {'data': 465.1909, 'unit': 'mm a-1'},
        'Niger': {'data': 402.23, 'unit': 'mm a-1'},
        'Nile': {'data': 602.1752, 'unit': 'mm a-1'},
        'Parana': {'data': 1085.554, 'unit': 'mm a-1'},
        'Yangtze-Kiang': {'data': 538.0664, 'unit': 'mm a-1'}
        }


def main(cfg):
    """Run the diagnostic.

    Parameters
    ----------
    cfg : dict
    Configuration dictionary of the recipe.

    """

    # Get dataset and variable information
    datalist = diag.Datasets(cfg)
    logging.debug("Found datasets in recipe:\n%s", datalist)
    varlist  = diag.Variables(cfg)
    logging.debug("Found variables in recipe:\n%s", varlist)

    # Check for correct variables
    if not varlist.vars_available('pr','mrro','evspsbl'):
        raise ValueError("This diagnostic requires input of precipitation, surface runoff and evaporation")

    # Read catchmentmask
    # to check: Correct way to read auxillary data using recipes?
    catchment_filepath = cfg.get('catchmentmask')
    catchment_cube     = iris.load_cube(catchment_filepath)

    # Read data and compute long term means
    # to check: Shouldn't this be part of preprocessing?
    # to check: How to regrid onto catchment_cube grid with preproc recipe statements
    #           instead of using regrid here?
    allcubes = []
    for dataset_path in datalist:
        new_cube        = iris.load(dataset_path, varlist.standard_names())[0]
        # Check for expected unit
        if new_cube.units != 'kg m-2 s-1':
            raise ValueError('Unit [kg m-2 s-1] is expected for ',new_cube.long_name.lower(),' flux')
        # Convert to unit mm/d
        new_cube.data *=  86400.0
        # Aggregate over year
        year_cube = new_cube.aggregated_by('year', iris.analysis.SUM)
        # Compute long term mean
        mean_cube       = year_cube.collapsed([diag.names.TIME], iris.analysis.MEAN)
        mean_cube.units = "mm a-1"
        # Regrid to catchment data grid --> maybe use area_weighted instead?
        mean_cube_regrid = regrid(mean_cube, catchment_cube, 'linear')
        # Doing something as yet unclear
        datalist.set_data(mean_cube_regrid.data, dataset_path)
        allcubes.append(mean_cube_regrid)

    # Write regridded data files
    # to do: update attributes
    filepath = os.path.join(cfg[diag.names.WORK_DIR], cfg.get('output_name', 'pp_runoff_et') + '.nc')
    if cfg[diag.names.WRITE_PLOTS]:
        iris.save(allcubes, filepath)
        logger.info("Writing %s", filepath)


    # Plot catchment data


if __name__ == '__main__':

    with diag.run_diagnostic() as config:
        main(config)
