"""Diagnostic for plotting a taylor diagram."""
import logging
import os
from collections import OrderedDict

import iris
import matplotlib.pyplot as plt
import numpy as np
import skill_metrics as sm
from scipy.stats import pearsonr

from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import get_plot_filename
from mpqb_utils import get_mpqb_cfg

logger = logging.getLogger(os.path.basename(__file__))


def _organize_datasets_for_mpqb(cfg):

    reference_dataset = cfg['reference_dataset']

    # Read all datasets that are provided.
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(input_data, 'alias', sort='alias')

    all_datasets = list(grouped_input_data.keys())

    # Put reference to the front of list
    all_datasets.insert(
        0, all_datasets.pop(all_datasets.index(reference_dataset)))
    ordered = OrderedDict()
    for k in all_datasets:
        ordered[k] = grouped_input_data[k]
    return ordered


def main(cfg):
    """Plot taylor diagram."""
    # read reference data set
    logger.info("Starting MPQB taylordiagram script.")

    grouped_input_data = _organize_datasets_for_mpqb(cfg)

    # Read from mpqb cfg
    reference_dataset = cfg['reference_dataset']

    # Create a pair of two datasets for inter-comparison
    refdat = iris.load_cube(
        grouped_input_data[reference_dataset][0]['filename'])
    rvalue_list = []
    rmsd_list = []
    std_list = []
    labels = []
    for alias in grouped_input_data.keys():
        label = get_mpqb_cfg('datasetname',alias)
        dataset = grouped_input_data[alias][0]['dataset']
        dat = iris.load_cube(grouped_input_data[alias][0]['filename'])
        # Create common masking
        refdat.data.mask |= dat.data.mask
        dat.data.mask |= refdat.data.mask
        adat = refdat.data.compressed()
        bdat = dat.data.compressed()
        rvalue_list.append(pearsonr(adat, bdat)[0])
        rmsd_list.append(sm.centered_rms_dev(adat, bdat))
        std_list.append(np.std(bdat))
        labels.append(label)

    obslabel = get_mpqb_cfg('datasetname', cfg['reference_dataset'])
    sm.taylor_diagram(np.array(std_list),
                      np.array(rmsd_list),
                      np.array(rvalue_list),
                      markerLabel=labels,
                      markerLegend='on',
                      markerColor='r',
                      markerSize=7,
                      rmsLabelFormat='0:.2f',
                      colObs='k',
                      markerObs='x',
                      titleOBS=obslabel,
                      checkstats='on')
    plot_filename = get_plot_filename('taylordiagram', cfg)
    logger.info("Writing plot to: %s", plot_filename)
    plt.savefig(plot_filename)
    logger.info("Finished!")


if __name__ == '__main__':
    with run_diagnostic() as config:
        if config['write_plots']:
            main(config)
        else:
            logger.warning("This diagnostic wants to plot,\
                            but isn't allowed to")
