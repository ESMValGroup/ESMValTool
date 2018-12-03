"""Python example diagnostic."""
import logging
import os

import iris
from iris.analysis import MEAN
from iris.analysis.stats import pearsonr

from diagnostic import plot_results
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Compute the time average for each input dataset."""
    input_data = group_metadata(
        cfg['input_data'].values(), 'standard_name', sort='dataset')

    for standard_name in input_data:
        logger.info("Processing variable %s", standard_name)
        # Load reference dataset
        for attributes in input_data[standard_name]:
            if attributes['reference_dataset'] == attributes['dataset']:
                reference_name = attributes['dataset']
                logger.info("Using %s as a reference dataset", reference_name)
                filename = attributes['filename']
                reference = iris.load_cube(filename)
                reference = reference.collapsed('time', MEAN)
                logger.info("Reference cube:\n%s\n%s", filename, reference)
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
            cube = pearsonr(dataset, reference, **kwargs)

            name = '{}_correlation_with_{}'.format(
                os.path.splitext(os.path.basename(filename))[0],
                reference_name)
            plot_results(cube, name, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
