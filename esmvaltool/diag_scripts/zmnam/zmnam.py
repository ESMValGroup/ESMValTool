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
#import cartopy.crs as ccrs
import matplotlib.pyplot as plt

import iris
#import iris.plot as iplt
#import iris.quickplot as qplt
import os
import logging

from esmvaltool.diag_scripts.shared import run_diagnostic # 20180601h1712

logger = logging.getLogger(__name__)

import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")

def get_input_files(cfg, index=0):
    """Get a dictionary with input files from metadata.yml files."""
    metadata_file = cfg['input_files'][index]
    with open(metadata_file) as file:
        metadata = yaml.safe_load(file)
    return metadata

def get_cfg():
    """Read diagnostic script configuration from settings.yml."""
    settings_file = sys.argv[1]
    print('settings file: ', settings_file)
    with open(settings_file) as file:
        cfg = yaml.safe_load(file)
    print('Configuration: ', cfg)
    return cfg


def main(config):
    #cfg = get_cfg()
    cfg = config
    """
    print('(1)@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')  
    print(cfg)
    print('(1)@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')  
    """
    logger.setLevel(cfg['log_level'].upper())

    input_files = get_input_files(cfg)
    #os.makedirs(cfg['plot_dir'])
    #os.makedirs(cfg['work_dir'])
    """
    print('(2)@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')  
    print(input_files) 
    print(list(input_files.keys()))
    print('(2)@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')  
    print('(3)@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')  
    print(input_files.items())
    print(list(input_files.items()))
    print('(3)@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')  
    """

    #sys.exit() 

    plot_dir=cfg['plot_dir']
    out_dir=cfg['work_dir']

    sys.path.append(cfg['path_diag_aux'])
    # Import full diagnostic routines
    from zmnam_calc import zmnam_calc
    from zmnam_plot import zmnam_plot
    from zmnam_preproc import zmnam_preproc

    #for variable_name, filenames in input_files.items():
    #for filenames,variable_name in input_files.items():
    filenames_cat = []
    for filenames in list(input_files.keys()): 
        #logger.info("Processing variable %s", variable_name)
        print(filenames) #
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

    #print(filenames_cat)
    #if 1 == 0:

        for ifile in filenames_cat:
             
            # Get model properties (stream, name, exp, ensemble member, period)
            ifile_props = ifile.rsplit('/',1)[1].rsplit('_',7) 
            ifile_props = [ifile_props[0],ifile_props[1],\
                          ifile_props[3],ifile_props[4],\
                          ifile_props[7].replace('.nc','')]

            # Diagnostics calculation. Input parameters for 
            # files and plots naming
            zmnam_preproc(ifile,[20,90])
            zmnam_calc(out_dir+'/',out_dir+'/',ifile_props) 
            zmnam_plot(out_dir+'/',plot_dir+'/',ifile_props)   

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

