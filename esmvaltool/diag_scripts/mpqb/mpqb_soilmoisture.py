"""Python example diagnostic."""
import logging
import os
from pprint import pformat
from sharedutils import parallel_apply_along_axis
from diag1d import *
import numpy as np

from mpqb_plots import mpqb_mapplot, get_ecv_plot_config
import iris
import itertools as it
import warnings

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger, get_diagnostic_filename, get_plot_filename)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))



def array2cube(array_in,cube_template):
    newcube = cube_template.copy()
    newcube.data = np.ma.fix_invalid(array_in)
    return newcube

class mpqb_pair:
    def __init__(self,cfg,ds_cfg,ds1name,ds2name):
        '''  ds_cfg should be input data grouped on dataset '''
        self.cfg = cfg
        self.ds_cfg = ds_cfg
        self.ds1 = ds1name
        self.ds2 = ds2name
        self.metrics = {}
    def load(self):
        self.ds1cube = iris.load_cube(self.ds_cfg[self.ds1][0]['filename'])
        self.template = self.ds1cube.collapsed('time', iris.analysis.MEAN)
        self.ds1dat = self.ds1cube.data.filled(np.nan)
        self.ds2dat = iris.load_cube(self.ds_cfg[self.ds2][0]['filename']).data.filled(np.nan)
        # Assert that shapes are the same of the two datasets
        if not len(set([self.ds1dat.shape,self.ds2dat.shape])) == 1:
            logger.error("Shapes of datasets should be the same")
            logger.error("{0} has shape: {1}".format(self.ds1,self.ds1dat.shape))
            logger.error("{0} has shape: {1}".format(self.ds2,self.ds2dat.shape))
            raise ValueError
    def pearsonr(self):
        self.metrics['pearsonr'],_ = parallel_apply_along_axis(pearsonr1d, 0, (self.ds1dat,self.ds2dat))
    def rmsd(self):
        self.metrics['rmsd'] = parallel_apply_along_axis(rmsd1d, 0, (self.ds1dat,self.ds2dat))
    def ubrmsd(self):
        self.metrics['ubrmsd'] = parallel_apply_along_axis(ubrmsd1d, 0, (self.ds1dat,self.ds2dat))
    def absdiff(self):
        with warnings.catch_warnings(): # silence the mean of empty slice warnings
                                        # this is expected behaviour
            warnings.simplefilter("ignore", category=RuntimeWarning) 
            self.metrics['absdiff'] = parallel_apply_along_axis(absdiffaxismean1d, 0, (self.ds1dat,self.ds2dat))
    def reldiff(self):
        self.metrics['reldiff'] = parallel_apply_along_axis(reldiffaxismean1d, 0, (self.ds1dat,self.ds2dat))
    def results2cube(self):
        for metric in self.metrics:
            self.metrics[metric] = array2cube(self.metrics[metric],self.template)
    def plot(self):
        for metricname,cube in self.metrics.items():
            basename = '{0}_{1}_{2}'.format(metricname,self.ds1,self.ds2)
            diagnostic_file = get_diagnostic_filename(basename, cfg)
            logger.info("Saving analysis results to %s", diagnostic_file)
            iris.save(cube, target=diagnostic_file)
            
            plot_file = get_plot_filename(basename, cfg)
            logger.info("Plotting analysis results to %s", plot_file)

            # Create provenance record
            provenance_record = {
                #TODO complete provenance record
                'caption' : "{0} between {1} and {2}".format(metricname,self.ds1,self.ds2),
                'plot_file' : plot_file
            }

            metrics_plot_dictionary = get_ecv_plot_config(self.ds_cfg[self.ds1][0]['short_name'])
            mpqb_mapplot(cube, plot_file, **metrics_plot_dictionary[metricname])

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

def main(cfg):
    """Compute the time average for each input dataset."""

    #TODO move these parameters to config file
    metrics_to_calculate = ['pearsonr', 'rmsd', 'absdiff', 'reldiff', 'ubrmsd']
    reference_dataset = 'cds-era5-land-monthly'

    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(
        input_data, 'dataset', sort='dataset')
    logger.info(
        "Example of how to group and sort input data by standard_name:"
        "\n%s", pformat(grouped_input_data))


    # Create a pair of two datasets for inter-comparison
    for dataset in grouped_input_data.keys():
        if dataset != reference_dataset:
            logger.info("Opening dataset: {0}".format(dataset))
            # Opening the pair 
            pair = mpqb_pair(cfg, grouped_input_data, dataset, reference_dataset)
            pair.load()
            # Execute the requested metrics
            for metricname in metrics_to_calculate:
                try:
                    getattr(pair,metricname)()
                except AttributeError:
                    logger.error("Metric %s is not defined. ", metricname)
            # Put the results back into cubes for plotting
            pair.results2cube()
            # Plot the results (if configured to plot)
            if cfg['write_plots']:
                pair.plot()
            logger.info("Finished comparison to ref for dataset: {0}".format(dataset))
    logger.info("Finished!")


if __name__ == '__main__':
    with run_diagnostic() as cfg:
        main(cfg)
