"""Python example diagnostic."""
import logging
import os
from pprint import pformat
from sharedutils import parallel_apply_along_axis
from diag1d import *
import numpy as np

import iris

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger, get_diagnostic_filename, get_plot_filename)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))


#def array2cube(array_in,cube_template):
#    newcube = cube_template.collapsed('time', iris.analysis.MEAN)
#    newcube.data = np.ma.fix_invalid(array_in)
#    return newcube

#TODO: put results of functions into CUBE
#TODO: loop over pairs of mpqb datasets


class mpqb_pair:
    def __init__(self,ds_cfg,ds1name,ds2name):
        '''  ds_cfg should be input data grouped on dataset '''
        self.ds_cfg = ds_cfg
        self.ds1 = ds1name
        self.ds2 = ds2name
        self.metrics = {}
    def load(self):
        self.ds1cube = iris.load_cube(self.ds_cfg[self.ds1][0]['filename'])
        self.ds1dat = self.ds1cube.data.filled(np.nan)
        self.ds2dat = iris.load_cube(self.ds_cfg[self.ds2][0]['filename']).data.filled(np.nan)
    def pearsonr(self):
        self.metrics['pearsonr'] = parallel_apply_along_axis(pearsonr1d, 0, (self.ds1dat,self.ds2dat))
    def rmsd(self):
        self.metrics['rmsd'] = parallel_apply_along_axis(rmsd1d, 0, (self.ds1dat,self.ds2dat))
    def absdiff(self):
        self.metrics['absdiff'] = parallel_apply_along_axis(absdiffaxismean1d, 0, (self.ds1dat,self.ds2dat))
    def reldiff(self):
        self.metrics['reldiff'] = parallel_apply_along_axis(reldiffaxismean1d, 0, (self.ds1dat,self.ds2dat))
    def template(self):
        self.template = self.ds1cube.collapsed('time', iris.analysis.MEAN)
#    def results2cube(self):
#        for metric in self.metrics:
#            self.metrics[metric] = array2cube(self.metrics[metric],self.ds1cube)




def compute_diagnostic(filename):
    """Compute an example diagnostic."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    logger.debug("Running example computation")
    return cube.collapsed('time', iris.analysis.MEAN)


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


def main(cfg):
    """Compute the time average for each input dataset."""

    #TODO move these parameters to config file
    metrics_to_calculate = ['pearsonr', 'rmsd', 'absdiff', 'reldiff']

    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(
        input_data, 'dataset', sort='dataset')
    logger.info(
        "Example of how to group and sort input data by standard_name:"
        "\n%s", pformat(grouped_input_data))

    # Just create a sample plot for now, to test provenance logging.
    for key in grouped_input_data:
        inputfile = grouped_input_data[key][0]['filename']
        output_basename = os.path.splitext(
            os.path.basename(inputfile))[0] + '_' + key


        provenance_record = {'caption' : 'this is a test record'}
        resultcube = compute_diagnostic(inputfile)
        plot_diagnostic(resultcube, output_basename, provenance_record, cfg)

    # Example of how to loop over variables/datasets in alphabetical order
    pair = mpqb_pair(grouped_input_data, 'cds-era5-monthly', 'cds-era5-monthly')
    pair.load()
   
    # Just as an example for now.
    # Touch the coords, otherwise core dumped 
    pair.ds1cube.coords()
    meancube = pair.ds1cube.collapsed('time', iris.analysis.MEAN)

    # Now go through metrics
    for metricname in metrics_to_calculate:
        try:
            getattr(pair,metricname)()
        except AttributeError:
            logger.error("Metric %s is not defined. ", metricname)

    logger.info("Finished!")


if __name__ == '__main__':
    with run_diagnostic() as cfg:
        main(cfg)
