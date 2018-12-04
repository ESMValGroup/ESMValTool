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

    #print(cfg.keys())
    #print (cfg)

    #cfg = get_cfg()

    logger.setLevel(cfg['log_level'].upper())

    input_files = get_input_files(cfg)
    #os.makedirs(cfg['plot_dir'])
    #os.makedirs(cfg['work_dir'])

    #input_cmor = get_input_cmor_table(cfg)
    #print(input_cmor)
    #stop 

    plot_dir=cfg['plot_dir']
    out_dir=cfg['work_dir']

    # Import full diagnostic routines
    sys.path.append(cfg['path_diag_aux'])
    from zmnam_calc import zmnam_calc
    from zmnam_plot import zmnam_plot
    from zmnam_preproc import zmnam_preproc

    # List of files to be compared as last step
    sim_list = []
    ref_list = []

    filenames_cat = []
    for filenames in list(input_files.keys()): 
        #logger.info("Processing variable %s", variable_name)
        #filenames=list(filenames) 
        #sys.exit() 
        #filenames_cat=[]
        filenames_cat.append(filenames)
        """  
        print('_____________________________\n{0} INPUT FILES:'.format(len(filenames)))
        for i in filenames:
            print(i)
            filenames_cat.append(i)
        #filenames_cat.append(filenames)
        print('_____________________________\n')

        #____________Building the name of output files
        """
        os.chdir(out_dir)

        # List of model simulations

        # List of reference datasets (if any)
       

    #print(filenames_cat)
    #if 1 == 0:

        for ifile in filenames_cat:

            # Get 6 model properties: stream, name, exp, ensemble member, period
            ifile_props = ifile.rsplit('/',1)[1].rsplit('_',7) 
            cmor_table = ifile_props[0]
            dataset = ifile_props[1]
            exp = ifile_props[3]
            ensemble = ifile_props[4]
            period = ifile_props[7].replace('.nc','')
            """
            ifile_props = [ifile_props[0],ifile_props[1],\
                          ifile_props[3],ifile_props[4],\
                          ifile_props[7].replace('.nc','')]
            """ 
            ifile_props = [cmor_table,dataset,exp,ensemble,period]

            # Diagnostics calculation. Input parameters for 
            # files and plots naming
            zmnam_preproc(ifile,[20,90]) # area selection should be removed
            zmnam_calc(out_dir+'/',out_dir+'/',ifile_props) 
            zmnam_plot(out_dir+'/',plot_dir+'/',ifile_props)   
            
            # len(OBS)>0 and len(mod)>0 do the difference, 
            # all-mod (interpolated) minus OBS
            #zmnam_diff()

"""
if __name__ == '__main__':
    iris.FUTURE.netcdf_promote = True
    logging.basicConfig(
        format="%(asctime)s [%(process)d] %(levelname)-8s "
               "%(name)s,%(lineno)s\t%(message)s"
    )
    main()

"""
# added 20180601h1717
if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)

"""
if __name__ == '__main__':
    iris.FUTURE.netcdf_promote = True
    logging.basicConfig(
        format="%(asctime)s [%(process)d] %(levelname)-8s "
               "%(name)s,%(lineno)s\t%(message)s"
    )
    main()

"""

