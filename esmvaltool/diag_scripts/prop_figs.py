"""Python example diagnostic."""
import logging
import os

import iris
from iris.analysis import MEAN
from iris.analysis.stats import pearsonr

from diagnostic import plot_diagnostic
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(attributes, ancestor_files, plot_type):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = ("Correlation of {long_name} between {dataset} and "
               "{reference_dataset}.".format(**attributes))

    record = {
        'caption': caption,
        'statistics': ['corr'],
        'domains': ['global'],
        'plot_type': plot_type,
        'authors': [
            'andela_bouwe',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestor_files,
    }
    return record


def main(cfg):
    """Compute the time average for each input dataset."""
    input_data = group_metadata(
        cfg['input_data'].values(),
        'standard_name',
        sort='dataset',
    )

    for standard_name in input_data:
        logger.info("Processing variable %s", standard_name)
        # Load reference dataset
        for attributes in input_data[standard_name]:
            if attributes['reference_dataset'] == attributes['dataset']:
                reference_name = attributes['dataset']
                logger.info("Using %s as a reference dataset", reference_name)
                reference_filename = attributes['filename']
                reference = iris.load_cube(reference_filename)
                reference = reference.collapsed('time', MEAN)
                logger.info("Reference cube:\n%s\n%s", reference_filename,
                            reference)
                break
        else:
            raise ValueError("No reference_dataset defined in recipe.")

        # Compute and plot correlation
        for attributes in input_data[standard_name]:
            if attributes['dataset'] == reference_name:
                continue
            logger.info("Processing dataset %s", attributes['dataset'])

            filename = attributes['filename']
            dataset = iris.load_cube(filename)
            kwargs = cfg.get('pearsonr', {})
            logger.info(
                "Computing correlation with settings %s between "
                "reference and cube:\n%s\n%s", kwargs, filename, dataset)
            dataset = dataset.collapsed('time', MEAN)
            # Fix issue with losing vertical bounds in extract_level
            # preprocessor
            if reference.coords(axis='Z'):
                ref_coord = reference.coord(axis='Z')
                coord = dataset.coord(ref_coord)
                if not coord.has_bounds():
                    coord.bounds = ref_coord.bounds
            cube = pearsonr(dataset, reference, **kwargs)

            name = '{}_correlation_with_{}'.format(
                os.path.splitext(os.path.basename(filename))[0],
                reference_name)
            provenance_record = get_provenance_record(
                attributes,
                ancestor_files=[reference_filename, filename],
                plot_type=cfg['plot_type'])
            plot_diagnostic(cube, name, provenance_record, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
