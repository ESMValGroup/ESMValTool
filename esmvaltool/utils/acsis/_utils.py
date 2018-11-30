"""
Perform basic functionalities common to acsis routines.

This utility module contains such functionalities
as setting logging, saving cubes to file etc.
"""
import os
import numpy as np
import iris


def _set_logger(logging, out_dir, log_file, log_level):
    # set logging for screen and file output
    root_logger = logging.getLogger()
    out_fmt = "%(asctime)s %(levelname)-8s %(name)s,%(lineno)s\t%(message)s"
    logging.basicConfig(
        filename=os.path.join(out_dir, log_file),
        filemode='a',
        format=out_fmt,
        datefmt='%H:%M:%S',
        level=logging.DEBUG)
    root_logger.setLevel(log_level.upper())
    logfmt = logging.Formatter(out_fmt)
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logfmt)
    root_logger.addHandler(console_handler)


def _save_cubes(time_axis, data_struct, variable, output_dir, cube_name):
    """Save generic cube."""
    if time_axis == 'years':
        times = iris.coords.DimCoord(
            np.array([f for f in data_struct.keys()]),
            standard_name='time',
            units='years')
        cube_data = np.array([f for f in data_struct.values()])
        cspec = [(times, 0)]
    elif time_axis == 'seasons':
        cube_data = data_struct
        cspec = []
    cube = iris.cube.Cube(
        cube_data,
        dim_coords_and_dims=cspec,
        long_name=variable)
    iris.save(cube, os.path.join(output_dir, cube_name))
