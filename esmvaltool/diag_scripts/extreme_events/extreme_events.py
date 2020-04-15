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
from esmvaltool.diag_scripts.shared._base import (ProvenanceLogger,
                                                  get_diagnostic_filename,
                                                  get_plot_filename)
from esmvaltool.diag_scripts.shared.plot import quickplot
import operator
import yaml

logger = logging.getLogger(os.path.basename(__file__))

index_definition = yaml.load(
"""
annual_number_of_frost_days:
    name: frost days
    required:
        - tasmin
    threshold:
        value: 273.15
        unit: K
        logic: lt
    cf_name: number_of_days_with_air_temperature_below_freezing_point
annual_number_of_summer_days:
    name: summer days
    required:
        - tasmax
    threshold:
        value: 298.15
        unit: K
        logic: gt
    cf_name: number_of_days_with_air_temperature_above_25_degree_Celsius
""")
print("INDEX_DEFINITION:")
print(yaml.dump(index_definition))

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


def _check_required_variables(required, available):
    missing = [item for item in required if item not in available]
    if len(missing):
        raise Exception(f"Missing required varaible {' and '.join(missing)}.")


def _count_days_by_threshold_annually(cubes, index):
    logger.info(f"Computing the annual number of {index['name']}.")
    _check_required_variables(index['required'], [c.var_name for c in cubes])
    # Here we assume that only one variable is required.
    cube = [item for item in cubes if item.var_name in index['required']].pop()
    iris.coord_categorisation.add_year(cube, 'time', name='year')
    annual_count = cube.aggregated_by(
        ['year'],
        iris.analysis.COUNT,
        function=lambda values: getattr(operator, index['threshold']['logic'])
        (values, index['threshold']['value']))
    annual_count.rename(index['cf_name'])
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
    with iris.fileformats.netcdf.Saver(os.path.join(outdir, filename),
                                       netcdf_format) as sman:
        for cube in cubes:
            sman.write(cube)


def main(cfg):
    """Compute Indices."""
    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(input_data, 'alias', sort='alias')
    for alias in grouped_input_data:
        cubes = []
        for attributes in grouped_input_data[alias]:
            cubes.append(iris.load_cube(attributes['filename']))
            #provenance_record = get_provenance_record(
            #    attributes,
            #    ancestor_files=['/to/do.nc'])  # TODO: this is wrong here
        for index_name in cfg['indices']:
            if index_name not in index_definition.keys():
                logger.info("Index %s not implemented!", index_name)
                continue
            if index_name in ['annual_number_of_frost_days', 'annual_number_of_summer_days']:
                write_netcdf([
                    _count_days_by_threshold_annually(cubes,
                                                  index_definition[index_name])
                ],
                         cfg,
                         filename=f'{alias}_{index_name}.nc')


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
