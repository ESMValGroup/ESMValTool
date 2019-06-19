"""
Zonal-mean Northern Annular Mode main routine.

Author: Federico Serva (ISAC-CNR & ISMAR-CNR, Italy)
Copernicus C3S 34a lot 2 (MAGIC)

Description:
Evaluation of stratosphere-troposphere coupling
based on EOF/PC analysis of the geopotential height field.

Modification history
20180512-serva_federico: Added output netCDFs, more use of preprocessor.
20180510-serva_federico: Routines written.

"""

import os
import logging

from esmvaltool.diag_scripts.shared import run_diagnostic, ProvenanceLogger

# Import zmnam diagnostic routines
from zmnam_calc import zmnam_calc
from zmnam_plot import zmnam_plot
from zmnam_preproc import zmnam_preproc

logger = logging.getLogger(__name__)


def get_provenance_record(vatt, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = ("Compute Zonal-mean Northern Annular Modes between "
               "{start_year} and {end_year} ".format(**vatt))
    record = {
        'caption': caption,
        'authors': [
            'serva_federico',
            'vonhardenberg_jost',
            'arnone_enrico',
            'cagnazzo_chiara',
        ],
        'projects': ['c3s-magic'],
        'references': ['baldwin09qjrms'],
        'plot_types': ['polar', 'zonal'],
        'realms': ['atmos'],
        'domains': ['polar'],
        'ancestors': ancestor_files,
    }
    return record


def main(cfg):
    """
    Run the zonal-mean NAM diagnostic.

    Calling in order:
    - preprocessing
    - index calculation
    - regression and plot
    """
    logger.setLevel(cfg['log_level'].upper())

    input_files = cfg['input_data']

    plot_dir = cfg['plot_dir']
    out_dir = cfg['work_dir']
    write_plots = cfg['write_plots']
    fig_fmt = cfg['output_file_type']

    filenames_cat = []
    fileprops_cat = []

    # Loop over input cfg
    for key, value in input_files.items():

        # Collect file names
        filenames_cat.append(key)

        # Collect relevant information for outputs naming
        fileprops_cat.append([
            value['project'], value['dataset'], value['exp'],
            value['ensemble'],
            str(value['start_year']) + '-' + str(value['end_year'])
        ])

    # Go to work_dir for running
    os.chdir(out_dir)

    # Process list of input files
    for indfile, ifile in enumerate(filenames_cat):

        ifile_props = fileprops_cat[indfile]

        # Call diagnostics functions
        print("prepro")
        (file_da_an_zm, file_mo_an) = zmnam_preproc(ifile)
        print("calc")
        outfiles = zmnam_calc(file_da_an_zm, out_dir + '/', ifile_props)
        provenance_record = get_provenance_record(
            list(input_files.values())[0], ancestor_files=ifile)
        if write_plots:
            print("plot_files")
            plot_files = zmnam_plot(file_mo_an, out_dir + '/', plot_dir +
                                    '/', ifile_props, fig_fmt, write_plots)
        else:
            plot_files = []
        for file in outfiles + plot_files:
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(file, provenance_record)


# Run the diagnostics
if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
