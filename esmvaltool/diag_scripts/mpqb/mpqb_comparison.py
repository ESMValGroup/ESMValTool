#!/usr/bin/env python
"""Python example diagnostic."""
import logging
import os
import warnings
from pprint import pformat

import iris
import numpy as np
import yaml

from esmvaltool.diag_scripts.shared.trend_mpqb_common.diag1d import (absdiffaxismean1d, pearsonr1d, reldiffaxismean1d, rmsd1d,
                    ubrmsd1d)
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import (ProvenanceLogger,
                                                  get_diagnostic_filename,
                                                  get_plot_filename)
from mpqb_plots import get_ecv_plot_config, mpqb_mapplot
from esmvaltool.diag_scripts.shared.trend_mpqb_common.sharedutils import parallel_apply_along_axis

logger = logging.getLogger(os.path.basename(__file__))


def _array2cube(array_in, cube_template):
    newcube = cube_template.copy()
    newcube.data = np.ma.fix_invalid(array_in)
    return newcube


class MPQBpair:
    """Class for calculating metrics on one pair of datasets."""

    def __init__(self, ds_cfg, ds1name, ds2name):
        """ds_cfg should be input data grouped on dataset."""
        self.ds_cfg = ds_cfg
        self.ds1 = ds1name
        self.ds2 = ds2name
        self.metrics = {}
        self.template = None
        self.ds1dat = None
        self.ds2dat = None

    def load(self):
        """load."""
        self.template = iris.load_cube(
            self.ds_cfg[self.ds1][0]['filename']).\
            collapsed('time', iris.analysis.MEAN)
        self.ds1dat = iris.load_cube(
            self.ds_cfg[self.ds1][0]['filename']).data.filled(np.nan)
        self.ds2dat = iris.load_cube(
            self.ds_cfg[self.ds2][0]['filename']).data.filled(np.nan)
        # Assert that shapes are the same of the two datasets
        if not self.ds1dat.shape == self.ds2dat.shape:
            logger.error("Shapes of datasets should be the same")
            logger.error("%s has shape:\n", self.ds1)
            logger.error("%s", self.ds1dat.shape)
            logger.error("%s has shape:\n", self.ds2)
            logger.error("%s", self.ds2dat.shape)
            raise ValueError

    def pearsonr(self):
        """pearsonr."""
        self.metrics['pearsonr'], _ = parallel_apply_along_axis(
            pearsonr1d, 0, (self.ds1dat, self.ds2dat))

    def rmsd(self):
        """rmsd."""
        self.metrics['rmsd'] = parallel_apply_along_axis(
            rmsd1d, 0, (self.ds1dat, self.ds2dat))

    def ubrmsd(self):
        """ubrmsd."""
        self.metrics['ubrmsd'] = parallel_apply_along_axis(
            ubrmsd1d, 0, (self.ds1dat, self.ds2dat))

    def absdiff(self):
        """absdiff."""
        with warnings.catch_warnings():  # silence the mean of empty
                # slice warnings this is expected behaviour
            warnings.simplefilter("ignore", category=RuntimeWarning)
            self.metrics['absdiff'] = parallel_apply_along_axis(
                absdiffaxismean1d, 0, (self.ds1dat, self.ds2dat))

    def reldiff(self):
        """reldiff."""
        self.metrics['reldiff'] = parallel_apply_along_axis(
            reldiffaxismean1d, 0, (self.ds1dat, self.ds2dat))

    def results2cube(self):
        """results2cube."""
        for metric in self.metrics:
            self.metrics[metric] = _array2cube(
                self.metrics[metric], self.template)

    def plot(self):
        """plot."""
        for metricname, cube in self.metrics.items():
            basename = '{0}_{1}_{2}'.format(metricname, self.ds1, self.ds2)
            diagnostic_file = get_diagnostic_filename(basename, cfg)
            logger.info("Saving analysis results to %s", diagnostic_file)
            iris.save(cube, target=diagnostic_file)

            plot_file = get_plot_filename(basename, cfg)
            logger.info("Plotting analysis results to %s", plot_file)

            # Create provenance record
            provenance_record = {
                'caption': "{0} between {1} and {2}".format(
                    metricname, self.ds1, self.ds2),
                'plot_file': plot_file
            }

            metrics_plot_dictionary = get_ecv_plot_config(
                self.ds_cfg[self.ds1][0]['short_name'])
            plot_kwargs = metrics_plot_dictionary[metricname]
            # Overwrite plot title to be dataset name
            plot_kwargs['title'] = self.ds1
            mpqb_mapplot(cube, plot_file, **
                         plot_kwargs)

            logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                        pformat(provenance_record))
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(diagnostic_file, provenance_record)


def compute_diagnostic(filename):
    """Compute an example diagnostic."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    logger.debug("Running example computation")
    return cube.collapsed('time', iris.analysis.MEAN)


def main():
    """MPQB diagnostics."""
    # read referenece data set
    reference_dataset = cfg['reference_dataset']

    # The metrics to be calculated.
    metrics_to_calculate = ['pearsonr', 'rmsd', 'absdiff', 'reldiff', 'ubrmsd']

    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(
        input_data, 'dataset', sort='dataset')
    logger.info(
        "Starting MPQB comparison script."
    )

    # Create a pair of two datasets for inter-comparison
    for dataset in grouped_input_data.keys():
        if dataset != reference_dataset:
            logger.info("Opening dataset: %s", dataset)
            # Opening the pair
            pair = MPQBpair(grouped_input_data,
                            dataset, reference_dataset)
            pair.load()
            # Execute the requested metrics
            for metricname in metrics_to_calculate:
                try:
                    getattr(pair, metricname)()
                except AttributeError:
                    logger.error("Metric %s is not defined. ", metricname)
            # Put the results back into cubes for plotting
            pair.results2cube()
            # Plot the results (if configured to plot)
            if cfg['write_plots']:
                pair.plot()
            logger.info(
                "Finished comparison to ref for dataset: %s", dataset)
    logger.info("Finished!")


if __name__ == '__main__':
    with run_diagnostic() as cfg:
        main()
