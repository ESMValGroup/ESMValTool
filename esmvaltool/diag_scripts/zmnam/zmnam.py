"""
;;#############################################################################
;; Zonal mean Northern Annular Mode Diagnostics
;; Author: Federico Serva (ISAC-CNR, Italy)
;; Copernicus C3S 34a lot 2 (MAGIC)
;;#############################################################################
;; Description
;;    Evaluation of stratosphere-troposphere coupling
;;    based on EOF/PC analysis of the geopotential height field
;;    
;; Modification history
;;    20180510-A_serv_fe: Routines written.
;;
;;#############################################################################
"""

import yaml
import sys
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

import iris
import iris.plot as iplt
import iris.quickplot as qplt
import os
import logging

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(__name__)

import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")


# Note: this can be removed
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

    plot_dir = cfg['plot_dir']
    out_dir = cfg['work_dir']
    write_plots = cfg['write_plots']

    # Import full diagnostic routines
    #sys.path.append(cfg['path_diag_aux'])
    from zmnam_calc import zmnam_calc
    from zmnam_plot import zmnam_plot
    from zmnam_preproc import zmnam_preproc

    # List of files to be compared as last step
    sim_list = []
    ref_list = []

    filenames_cat = []
    fileprops_cat = []

    # Loop over input cfg
    for key, value in input_files.items():

        # Collect file names
        filenames_cat.append(key)
 
        # Collect relevant information for outputs naming
        fileprops_cat.append([value['project'],
                             value['dataset'],
                             value['exp'],
                             value['ensemble'],
                             str(value['start_year'])+'-'
                             +str(value['end_year'])])

    os.chdir(out_dir)

    # Process list of input files
    for indfile in range(len(filenames_cat)):

        ifile = filenames_cat[indfile]
        ifile_props = fileprops_cat[indfile]

        # Call diagnostics functions
        zmnam_preproc(ifile, [20, 90])  # area selection should be removed
        zmnam_calc(out_dir + '/', out_dir + '/', ifile_props)
        zmnam_plot(out_dir + '/', plot_dir + '/', ifile_props)

    """

    filenames_cat = []
    for filenames in list(input_files.keys()):
        filenames_cat.append(filenames)

        os.chdir(out_dir)

        # List of model simulations

        # List of reference datasets (if any)

        for ifile in filenames_cat:

            # Get 6 model properties: stream, name, exp, ensemble member, period
            ifile_props = ifile.rsplit('/', 1)[1].rsplit('_', 7)
            cmor_table = ifile_props[0]  # project!!!
            dataset = ifile_props[1]
            exp = ifile_props[3]
            ensemble = ifile_props[4]
            period = ifile_props[7].replace('.nc', '')
            ifile_props = [cmor_table, dataset, exp, ensemble, period]

            # Diagnostics calculation. Input parameters for
            # files and plots naming
            zmnam_preproc(ifile, [20, 90])  # area selection should be removed
            zmnam_calc(out_dir + '/', out_dir + '/', ifile_props)
            zmnam_plot(out_dir + '/', plot_dir + '/', ifile_props)

            # len(OBS)>0 and len(mod)>0 do the difference,
            # all-mod (interpolated) minus OBS
            #zmnam_diff()

    """    

# Run the diagnostics
if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
