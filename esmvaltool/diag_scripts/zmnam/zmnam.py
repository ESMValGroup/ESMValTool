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

import logging
import os
from pprint import pformat

import iris

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))

def main(cfg):
    """Compute the time average for each input dataset."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    # Demonstrate use of metadata access convenience functions.
    selection = select_metadata(input_data, short_name='zg', project='CMIP5')
    logger.info("Example of how to select only CMIP5 gh data:\n%s",
                pformat(selection))

    selection = sorted_metadata(selection, sort='dataset')
    logger.info("Example of how to sort this selection by dataset:\n%s",
                pformat(selection))

    grouped_input_data = group_metadata(
        input_data, 'standard_name', sort='dataset')
    logger.info(
        "Example of how to group and sort input data by standard_name:"
        "\n%s", pformat(grouped_input_data))

    # Example of how to loop over variables/datasets in alphabetical order
    for standard_name in grouped_input_data:
        logger.info("Processing variable %s", standard_name)
        for attributes in grouped_input_data[standard_name]:
            logger.info("Processing dataset %s", attributes['dataset'])

            filename = attributes['filename']
            logger.debug("Loading %s", filename)
            cube = iris.load_cube(filename)

            logger.debug("Running example computation")
            cube = cube.collapsed('time', iris.analysis.MEAN)

            name = os.path.splitext(os.path.basename(filename))[0] + '_mean'
            plot_results(cube, name, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)


if 1 == 0:
        
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

        print('**********************')
        print('before meain')
        print('**********************')
        
        cfg = get_cfg()

        logger.setLevel(cfg['log_level'].upper())

        input_files = get_input_files(cfg)
        #os.makedirs(cfg['plot_dir'])
        #os.makedirs(cfg['work_dir'])

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

