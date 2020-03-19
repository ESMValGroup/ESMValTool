#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Calculates the ETCCDI Climate Change Indices."""
import logging
import os
from pprint import pformat
import numpy as np
import iris
import iris.coord_categorisation

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger, get_diagnostic_filename, get_plot_filename)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(attributes, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = ("Average {long_name} between {start_year} and {end_year} "
               "according to {dataset}.".format(**attributes))

    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_type': 'zonal',
        'authors': [
            'andela_bouwe',
            'righi_mattia',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestor_files,
    }
    return record

def _frost_days(cubes):
    logger.info("Computing the annual number of frost days.")
    required_variables = ['tasmin']
    missing = [item for item in required_variables if item not in cubes.keys()]
    if len(missing):
        raise Exception(f"Missing required varaible {' and '.join(missing)}.")
    cube = cubes['tasmin']
    iris.coord_categorisation.add_year(cube, 'time', name='year')
    annual_count = cube.aggregated_by(['year'], iris.analysis.COUNT, function=lambda values: values < 273.15)
    annual_count.rename('number_of_days_with_air_temperature_below_freezing_point')
    return annual_count


def plot_diagnostic(cube, basename, provenance_record, cfg):
    """Create diagnostic data and plot it."""
    diagnostic_file = get_diagnostic_filename(basename, cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)
    iris.save(cube, target=diagnostic_file)

    if cfg['write_plots'] and cfg.get('quickplot'):
        plot_file = get_plot_filename(basename, cfg)
        logger.info("Plotting analysis results to %s", plot_file)
        provenance_record['plot_file'] = plot_file
        quickplot(cube, filename=plot_file, **cfg['quickplot'])

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)

def write_netcdf(cubes, cfg, filename='test.nc', netcdf_format='NETCDF4'):
    outdir = cfg['work_dir']
    with iris.fileformats.netcdf.Saver(os.path.join(outdir, filename), netcdf_format) as sman:
        for cube in cubes:
            sman.write(cube)

def main(cfg):
    """Compute Indices."""
    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(
        input_data, 'alias', sort='alias')

    cubes  = {}
    for alias in grouped_input_data:
        logger.info("Processing alias %s", alias)
        cubes[alias] = {}
        for attributes in grouped_input_data[alias]:
            logger.info("Processing dataset %s", attributes['dataset'])
            input_file = attributes['filename']
            output_basename = os.path.splitext(
                    os.path.basename(input_file))[0]
            cubes[alias][attributes['short_name']] = iris.load_cube(input_file)
            provenance_record = get_provenance_record(
                    attributes, ancestor_files=[input_file]) # TODO: this is wrong here
        for index in cfg['indices']:
            if index == 'annual_number_of_frost_days':
                write_netcdf([_frost_days(cubes[alias])], cfg, filename=f'{alias}_frost-days.nc')
            elif index == 'annual_number_of_summer_days':
                pass
            else:
                logger.info("Index %s not implemented!", index)

if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
