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
#import imp

import yaml
import sys
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

import iris
import iris.plot as iplt
import iris.quickplot as qplt
import os
import logging

logger = logging.getLogger(__name__)

import warnings
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


def main():
    cfg = get_cfg()
    logger.setLevel(cfg['log_level'].upper())

    input_files = get_input_files(cfg)
    os.makedirs(cfg['plot_dir'])
    os.makedirs(cfg['work_dir'])
    
    out_dir=cfg['work_dir']

    sys.path.append(cfg['path_diag_aux'])
    # Import full diagnostic routines
    from ens_anom import ens_anom
    from ens_eof_kmeans import ens_eof_kmeans
    from ens_plots import ens_plots

    for variable_name, filenames in input_files.items():
        logger.info("Processing variable %s", variable_name)

    filenames_cat=[]
    print('_____________________________\n{0} INPUT FILES:'.format(len(filenames)))
    for i in filenames:
        print(i)
        filenames_cat.append(i)
    print('_____________________________\n')

    #____________Building the name of output files
    name_outputs=variable_name+'_'+str(cfg['numens'])+'ens_'+cfg['season']+'_'+cfg['area']+'_'+cfg['kind']
    print('The name of the output files will be <variable>_{0}.txt'.format(name_outputs))
   
    ####################### PRECOMPUTATION #######################
    #____________run ens_anom as a module
    ens_anom(filenames_cat,out_dir,name_outputs,variable_name,cfg['numens'],cfg['season'],cfg['area'],cfg['extreme'])
    
    ####################### EOF AND K-MEANS ANALYSES #######################
    #____________run ens_eof_kmeans as a module
    ens_eof_kmeans(out_dir,name_outputs,cfg['numens'],cfg['numpcs'],cfg['perc'],cfg['numclus'])
    
    ####################### PLOT AND SAVE FIGURES ################################
    #____________run ens_plots as a module
    ens_plots(out_dir,cfg['plot_dir'],name_outputs,cfg['numclus'],cfg['field_to_plot'])
    
    print('\n>>>>>>>>>>>> ENDED SUCCESSFULLY!! <<<<<<<<<<<<\n')
    print('')


if __name__ == '__main__':
    iris.FUTURE.netcdf_promote = True
    logging.basicConfig(
        format="%(asctime)s [%(process)d] %(levelname)-8s "
               "%(name)s,%(lineno)s\t%(message)s"
    )
    main()
