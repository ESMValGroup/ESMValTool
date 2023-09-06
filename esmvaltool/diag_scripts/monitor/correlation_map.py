"""Python diagnostic computing temporal correlation of variable(s) to
reference variable on a global map."""

import logging
import os

import iris
import matplotlib.pyplot as plt
from iris.analysis.stats import pearsonr

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
)
from esmvaltool.diag_scripts.shared.plot._plot import global_pcolormesh

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(attributes, ancestor_files, plot_type, ref_variable):
    """Create a provenance record describing the diagnostic data and plot."""
    names = {'long_name': attributes['long_name'],
             'ref_variable': ref_variable,
             'dataset': attributes['dataset']}
    caption = ("Correlation of {long_name} and "
               "{ref_variable} for {dataset}.".format_map(names))

    record = {
        'caption': caption,
        'statistics': ['corr'],
        'domains': ['global'],
        'plot_type': plot_type,
        'authors': [
            'malles_jan-hendrik',
            'weigel_katja',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestor_files,
    }
    return record


def plot_correlation_map(cube, basename, title, provenance_record, cfg):
    """Create diagnostic data and plot it."""

    # Save the data used for the plot
    save_data(basename, provenance_record, cfg, cube)

    if cfg.get('plot_config') or cfg.get('colorbar_config') or \
       cfg.get('plot_title'):
        # Create the plot
        if cfg.get('colorbar_config'):
            cbar_config = cfg.get('colorbar_config')
            global_pcolormesh(cube, cbar_label=cbar_config.get('label',
                                                               'Correlation'),
                              cbar_ticks=cbar_config.get('ticks'),
                              cbar_center=cbar_config.get('center'),
                              **cfg['plot_config'])
        elif cfg.get('plot_config'):
            global_pcolormesh(cube, **cfg['plot_config'])
        else:
            global_pcolormesh(cube)
        plt.gca().set_title(title)
        # And save the plot
        save_figure(basename, provenance_record, cfg)


def main(cfg):
    """Compute the time average for each input dataset."""
    input_data = group_metadata(
        cfg['input_data'].values(),
        'dataset',
        sort='short_name',
    )

    for dataset in input_data:
        logger.info("Processing dataset %s", dataset)
        # Load reference dataset
        for attributes in input_data[dataset]:
            if cfg['reference_variable'] == attributes['short_name']:
                ref_variable = attributes['short_name']
                reference_name = attributes['long_name']
                logger.info("Using %s as a reference dataset",
                            ref_variable)
                reference_filename = attributes['filename']
                reference = iris.load_cube(reference_filename)
                logger.info("Reference cube:\n%s\n%s", reference_filename,
                            reference)
                break
        else:
            raise ValueError("No reference_variable defined in recipe.")
        # Compute and plot correlation
        for attributes in input_data[dataset]:
            if attributes['short_name'] == ref_variable:
                continue
            logger.info("Processing dataset %s", attributes['dataset'])

            filename = attributes['filename']
            dataset = iris.load_cube(filename)
            kwargs = cfg.get('pearsonr', {})
            logger.info(
                "Computing correlation with settings %s between "
                "reference and cube:\n%s\n%s", kwargs, filename, dataset)

            cube = pearsonr(dataset, reference, **kwargs)
            names = {'long_name': attributes['long_name'],
                     'ref_variable': reference_name,
                     'dataset': attributes['dataset']}
            if cfg.get('plot_title'):
                title = cfg.get('plot_title')
            else:
                title = ("{long_name} vs. {ref_variable} "
                         "for {dataset}".format_map(names))
            name = '{}_correlation_with_{}'.format(
                os.path.splitext(os.path.basename(filename))[0],
                ref_variable)
            provenance_record = get_provenance_record(
                attributes,
                ancestor_files=[reference_filename, filename],
                plot_type=cfg.get('plot_type', 'pcolormesh'),
                ref_variable=cfg['reference_variable'])
            plot_correlation_map(cube, name, title, provenance_record, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
