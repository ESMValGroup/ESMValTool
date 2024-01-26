"""Python example diagnostic."""
import logging
from pathlib import Path
from pprint import pformat

import iris
import numpy as np

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
    select_metadata,
    sorted_metadata,
)
from esmvaltool.diag_scripts.shared.plot import quickplot
import matplotlib.pyplot as plt

logger = logging.getLogger(Path(__file__).stem)


def get_provenance_record(caption, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['reg'],
        'plot_types': ['times'],
        'authors': ['bock_lisa'],
        'references': ['cmug'],
        'ancestors': ancestor_files,
    }
    return record


def read_data(filename):
    """Compute an example diagnostic."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    logger.debug("Running example computation")
    cube = iris.util.squeeze(cube)
    return cube


def plot_diagnostic(cube, basename, provenance_record, cfg):
    """Create diagnostic data and plot it."""

    # Save the data used for the plot
    #save_data(basename, provenance_record, cfg, cube)

    #if cfg.get('quickplot'):
    # Create the plot
    quickplot(cube, **cfg['plot_timeseries'])

    #fig, axs = plt.subplots(25, sharex=True)
    #fig.suptitle(cube.var_name)

    #ipoint = 0
    #time = cube.coord('time').points

    #for ilat, lat in enumerate(cube.coord('latitude')):
    #    for ilon, lon in enumerate(cube.coord('longitude')):
    #        #if cube.data[0,ilat,ilon] > 0.:
    #        if cube.data[0,ilat,ilon] == cube.data[0,ilat,ilon]:
    #            print(ipoint, lat.points, lon.points, cube.data[:,ilat,ilon])
    #            axs[ipoint].plot(time, cube.data[:,ilat,ilon], linewidth=0.5)
    #            axs[ipoint].set_ylim(bottom=np.min(cube.data[:,:,:]), top=np.max(cube.data[:,:,:]))
    #            axs[ipoint].tick_params(axis="y", labelsize=5)
    #            ipoint = ipoint + 1

    # And save the plot
    save_figure(basename, provenance_record, cfg)


def potential_temperature(temperature, plev):
    """Compute potential temperature.

    Parameters
    ----------
    temperature: iris.cube.Cube
        Cube of air temperature ta.
    plev: iris.coords.Coord
        Pressure level coordinates

    Returns
    -------
    theta: iris.cube.Cube
        Cube of potential temperature theta.
    """

    """Reference Pressure [hPa]."""
    reference_pressure = 1000.
    pressure = (reference_pressure / plev)**(2 / 7)
    theta = temperature * pressure
    theta.long_name = 'potential_air_temperature'

    return theta
    

def calculate_theta(var, input_data, cfg):
    # This function calculates the potential temperature for one level.

    print(input_data)
    t_data = select_metadata(input_data, short_name='ta'+var[5:])
    print(t_data)
    t_file = t_data[0]['filename']
    temperature = read_data(t_file)

    pressure = float(var[5:])

    theta_x = potential_temperature(temperature, pressure)  

    # Write output
    logger.info("Saving data.")
    basename = var
    if "caption" not in t_data[0]:
        t_data[0]['caption'] = t_file
    provenance_record = get_provenance_record(
        t_data[0]['caption'], ancestor_files=[t_file])
    save_data(basename, provenance_record, cfg, theta_x)

    return theta_x, t_file


def calculate_lts(var, input_data, cfg):
    # This function calculates the potential temperature for one level.

    theta_700, file_1  = calculate_theta('theta700', input_data, cfg)

    theta_1000, file_2  = calculate_theta('theta1000', input_data, cfg)

    lts = theta_700 - theta_1000

    lts.long_name = 'lower tropospheric stability'

    # Write output
    logger.info("Saving data.")
    basename = var
    caption = 'lower_tropospheric_stability' 
    provenance_record = get_provenance_record(
        caption, ancestor_files=[file_1, file_2])
    save_data(basename, provenance_record, cfg, lts)
    if cfg.get('plot_timeseries'):
        plot_diagnostic(lts, basename, provenance_record, cfg)

    return lts


def main(cfg):
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    ## Demonstrate use of metadata access convenience functions.
    #selection = select_metadata(input_data, short_name='tas', project='CMIP5')
    #logger.info("Example of how to select only CMIP5 temperature data:\n%s",
    #            pformat(selection))

    #selection = sorted_metadata(selection, sort='dataset')
    #logger.info("Example of how to sort this selection by dataset:\n%s",
    #            pformat(selection))

    #grouped_input_data = group_metadata(input_data,
    #                                    'variable_group',
    #                                    sort='dataset')
    #logger.info(
    #    "Example of how to group and sort input data by variable groups from "
    #    "the recipe:\n%s", pformat(grouped_input_data))

    ## Example of how to loop over variables/datasets in alphabetical order
    groups = group_metadata(input_data, 'variable_group', sort='dataset')
    for group_name in groups:
        logger.info("Processing variable %s", group_name)
        for attributes in groups[group_name]:
            logger.info("Processing dataset %s", attributes['dataset'])
            input_file = attributes['filename']
            cube = read_data(input_file)

            output_basename = Path(input_file).stem
            if group_name != attributes['short_name']:
                output_basename = group_name + '_' + output_basename
            if "caption" not in attributes:
                attributes['caption'] = input_file
            provenance_record = get_provenance_record(
                attributes, ancestor_files=[input_file])
            if cfg.get('plot_timeseries'):
                plot_diagnostic(cube, output_basename, provenance_record, cfg)

    if cfg.get('compute'):
        for var in cfg['compute']:
            if var not in ['LTS', 'theta700', 'theta850']:
                logger.error('Computation of %s is not available...', var)
            if 'theta' in var:
                theta, file_theta = calculate_theta(var, input_data, cfg)
            if var == 'LTS':
                for ivar in ['ta700', 'ta1000']:
                    if ivar not in groups:
                        logger.error('Variable %s is needed to calculate LTS but is not available', ivar)
                lts = calculate_lts(var, input_data, cfg)
                #lts = calculate_lts(var, groups, cfg)





if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
