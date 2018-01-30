"""Python example diagnostic."""
import yaml
import sys
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

import iris
import iris.plot as iplt
import iris.quickplot as qplt
import os
import logging

logger = logging.getLogger(__name__)


def get_cfg():
    """Read diagnostic script configuration from settings.yml."""
    settings_file = sys.argv[1]
    with open(settings_file) as file:
        cfg = yaml.safe_load(file)
    return cfg


def get_input_files(cfg, index=0):
    """Get a dictionary with input files from metadata.yml files."""
    metadata_file = cfg['input_files'][index]
    with open(metadata_file) as file:
        metadata = yaml.safe_load(file)
    return metadata


def plot2d(cube, filename):
    logger.info("Creating %s", filename)
    fig = plt.figure()
    qplt.pcolormesh(cube)
    plt.gca().coastlines()
    fig.savefig(filename)


def main():

    cfg = get_cfg()
    logger.setLevel(cfg['log_level'].upper())

    input_files = get_input_files(cfg)
    os.makedirs(cfg['plot_dir'])

    for variable_name, filenames in input_files.items():
        logger.info("Processing variable %s", variable_name)
        for filename, attributes in filenames.items():
            plot_filename = os.path.join(
                cfg['plot_dir'],
                os.path.splitext(os.path.basename(filename))[0] + '.png',
            )
            cube = iris.load_cube(filename)
            cube = cube.collapsed('time', iris.analysis.MEAN)
            plot2d(cube, plot_filename)


if __name__ == '__main__':
    iris.FUTURE.netcdf_promote = True
    logging.basicConfig(
        format=
        "%(asctime)s [%(process)d] %(levelname)-8s %(name)s,%(lineno)s\t%(message)s"
    )
    main()
