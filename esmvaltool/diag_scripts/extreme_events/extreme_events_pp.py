#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculates the ETCCDI Climate Change Indices."""
import logging
import os
from pprint import pformat
import numpy as np
import iris

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger, get_diagnostic_filename, get_plot_filename)
from esmvaltool.diag_scripts.shared.plot import quickplot

import extreme_events_indices

from extreme_events_utils import convert_ETCCDI_units

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


def write_netcdf(cube, cfg, filename='test.nc', netcdf_format='NETCDF4'):
    outdir = cfg['work_dir']
    with iris.fileformats.netcdf.Saver(os.path.join(outdir, filename),
                                       netcdf_format) as sman:
        sman.write(cube)


def plot_diagnostic(cube, basename, provenance_record, cfg):
    """Create diagnostic data and plot it."""
    diagnostic_file = get_diagnostic_filename(basename, cfg)

    logger.info('Saving analysis results to %s', diagnostic_file)
    iris.save(cube, target=diagnostic_file)

    if cfg['write_plots'] and cfg.get('quickplot'):
        plot_file = get_plot_filename(basename, cfg)
        logger.info('Plotting analysis results to %s', plot_file)
        provenance_record['plot_file'] = plot_file
        quickplot(cube, filename=plot_file, **cfg['quickplot'])

    logger.info('Recording provenance of %s:\n%s', diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)


def main(cfg):
    """Compute the time average for each input dataset."""
    # Get the data that we will use as input.
    input_data = cfg.pop('input_data').values()
    logger.info("cfg:\n%s", pformat(cfg))

    # Grouped and sorted the data by alias.
    grouped_input_data = group_metadata(
        input_data, 'alias', sort='alias')

    # Loop over variables/datasets in alphabetical order of alias
    for alias, alias_data in grouped_input_data.items():
        logger.info('Processing alias %s', alias)

        alias_cubes = {}
        # Save all alias related cubes in dictionary
        for attributes in alias_data:
            logger.info('Processing dataset {}, {}.'.format(attributes['dataset'], attributes['short_name']))
            input_file = attributes['filename']
            alias_cubes[attributes['short_name']] = convert_ETCCDI_units(iris.load_cube(input_file))
            local_atts = attributes.copy()
            base_fname = "_{}_{}_{}_{}".format(local_atts.pop("dataset", "NULL"),
                                          local_atts.pop("exp", "NULL"),
                                          local_atts.pop("ensemble", "NULL"),
                                          "{}-{}".format(
                                                  local_atts.pop("start_year", "NULL"),
                                                  local_atts.pop("end_year", "NULL"),
                      ))

        # get index function
        for index_name in cfg['indices']:
            logger.info('Computing index %s', index_name)
            etccdi_index = getattr(extreme_events_indices,
                                   extreme_events_indices.index_method[
                                           index_name],
                                   None)

            # if there is no index available: report and break
            if etccdi_index is None:
                logger.error('There is no index available called: {}'.format(index_name))
                logger.error('*** calculation failed ***')
                return

#            logger.info(alias_cubes)
#            logger.info(etccdi_index(alias_cubes, cfg))
            # calculate and save cube
            etccdi_cube = etccdi_index(
                {k:v.copy() for k,v in alias_cubes.items()}, cfg=cfg)
            fname = etccdi_cube.attributes["FI_index"] + base_fname
            iris.save(etccdi_cube,
                      cfg['work_dir'] + os.sep + fname + '.nc')

            logger.info('Finalized computation for %s', ', '.join([alias, attributes['dataset'], index_name]))

    #output_basename = os.path.splitext(
    #        os.path.basename(input_file))[0] + '_mean'
    #provenance_record = get_provenance_record(
    #        attributes, ancestor_files=[input_file])
    #call plot idices
    #plot_diagnostic(cube, output_basename, provenance_record, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
