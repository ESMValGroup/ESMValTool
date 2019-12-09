#!/usr/bin/env python
import logging
import os
import warnings
from pprint import pformat

import yaml
import skill_metrics as sm
import matplotlib.pyplot as plt
import iris
import numpy as np
from scipy.stats import pearsonr
from collections import OrderedDict

from esmvaltool.diag_scripts.shared.trend_mpqb_common.diag1d import rmsd1d
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import (ProvenanceLogger,
                                                  get_diagnostic_filename,
                                                  get_plot_filename)
from mpqb_plots import get_ecv_plot_config, mpqb_mapplot
from esmvaltool.diag_scripts.shared.trend_mpqb_common.sharedutils import parallel_apply_along_axis

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    # read referenece data set
    reference_dataset = cfg['reference_dataset']
    
    # Read all datasets that are provided.
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(
        input_data, 'dataset', sort='dataset')
    logger.info(
        "Starting MPQB taylordiagram script."
    )

    all_datasets = list(grouped_input_data.keys())

    # Put reference to the front of list
    all_datasets.insert(0, all_datasets.pop(all_datasets.index(reference_dataset)))
    ordered = OrderedDict()
    for k in all_datasets:
        ordered[k] = grouped_input_data[k]

    # Create a pair of two datasets for inter-comparison
    refdat = iris.load_cube(grouped_input_data[reference_dataset][0]['filename'])
    rvalue_list = []
    rmsd_list = []
    std_list = []
    for dataset in ordered.keys():
        dat = iris.load_cube(grouped_input_data[dataset][0]['filename'])
        # Create common masking
        refdat.data.mask |= dat.data.mask
        dat.data.mask |= refdat.data.mask
        a = refdat.data.compressed()
        b = dat.data.compressed()
        rvalue_list.append(pearsonr(a,b)[0])
        rmsd_list.append(sm.centered_rms_dev(a,b)) # replace rmsd1d with skillmetrics.centered_rmsd
                                      # https://github.com/PeterRochford/SkillMetrics/
                                      # blob/master/Examples/target3.py
                                      # and switch on checkstat in taylordiagram plot
        std_list.append(np.std(b))


    sm.taylor_diagram(np.array(std_list),
                                 np.array(rmsd_list),
                                 np.array(rvalue_list),
                                 markerLabel=all_datasets,
                                 markerLegend = 'on',
                                 markerColor = 'r',
                                 checkstats='on')
    if cfg['write_plots']:
        plot_filename = get_plot_filename('taylordiagram',cfg)
        logger.info(
                "Writing plot to: %s", plot_filename)
        plt.savefig(plot_filename)

    logger.info("Finished!")





if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
