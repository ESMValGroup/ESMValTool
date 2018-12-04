"""
;;#############################################################################
;; Ensemble Clustering Diagnostics
;; Author: Irene Mavilia (ISAC-CNR, Italy)
;; Copernicus C3S 34a lot 2 (MAGIC)
;;#############################################################################
;; Description
;;    Cluster analysis tool based on the k-means algorithm
;;    for ensembles of climate model simulations
;; Modification history
;;    20181002-A_arno_ea: updating to version2_develpment (recipe/dataset)
;;    20170710-A_mavi_ir: Routines written.
;;
;;#############################################################################
"""

# Basic Python packages

import yaml
import sys
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import os
import logging
from esmvaltool.diag_scripts.shared import\
     (plot, run_diagnostic, save_iris_cube, variables_available,
      extract_variables)
import warnings

plt.switch_backend('Agg')

logger = logging.getLogger(__name__)

if not sys.warnoptions:
    warnings.simplefilter("ignore")


def get_cfg():
    """Read diagnostic script configuration from settings.yml."""
    settings_file = sys.argv[1]
    with open(settings_file) as file:
        cfg = yaml.safe_load(file)
    return cfg


def get_input_files(cfg, index=0):
    """Get a dictionary with input files from metadata.yml files."""
    metadata_file = cfg['input_files'][index]
    with open(metadata_file) as file:
        metadata = yaml.safe_load(file)
    return metadata


def main(cfg):
    logger.setLevel(cfg['log_level'].upper())

    input_files = get_input_files(cfg)
    os.makedirs(cfg['plot_dir'], exist_ok=True)
    os.makedirs(cfg['work_dir'], exist_ok=True)

    out_dir = cfg['work_dir']
    write_plots = cfg['write_plots']

    # Import full diagnostic routines
    from ens_anom import ens_anom
    from ens_eof_kmeans import ens_eof_kmeans
    from ens_plots import ens_plots

    filenames_cat = []
    numens = len(input_files.keys())
    for element in input_files.values():
        logger.info("Processing file %s", element['filename'])
        filenames_cat.append(element['filename'])
    name_outputs = element['short_name'] + '_' + str(numens) + \
        'ens_' + cfg['season'] + '_' + cfg['area'] + \
        '_' + element['project'] + '_' + element['exp']
    variable_name = element['short_name']

    # Building the name of output files
    print('The name of the output files will be <variable>_{0}.txt'
          .format(name_outputs))

    # ###################### PRECOMPUTATION #######################
    # ____________run ens_anom as a module
    ens_anom(filenames_cat, out_dir, name_outputs, variable_name,
             numens, cfg['season'], cfg['area'], cfg['extreme'])

    # ###################### EOF AND K-MEANS ANALYSES #######################
    # ____________run ens_eof_kmeans as a module
    ens_eof_kmeans(out_dir, name_outputs, numens, cfg['numpcs'],
                   cfg['perc'], cfg['numclus'])

    # ###################### PLOT AND SAVE FIGURES ##########################
    # ____________run ens_plots as a module
    if write_plots:
        ens_plots(out_dir, cfg['plot_dir'], name_outputs, cfg['numclus'],
                  'anomalies')  # cfg['file_to_plot']

    print('\n>>>>>>>>>>>> ENDED SUCCESSFULLY!! <<<<<<<<<<<<\n')
    print('')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
