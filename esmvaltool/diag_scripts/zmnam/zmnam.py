"""

Zonal-mean Northern Annular Mode main routine.

Author: Federico Serva (ISAC-CNR & ISMAR-CNR, Italy)
Copernicus C3S 34a lot 2 (MAGIC)

Description:
Evaluation of stratosphere-troposphere coupling
based on EOF/PC analysis of the geopotential height field.

Modification history
20180512-A_serv_fe: Added output netCDFs, more use of preprocessor.
20180510-A_serv_fe: Routines written.

"""

import os
import logging

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger

# Import zmnam diagnostic routines
from zmnam_calc import zmnam_calc
from zmnam_plot import zmnam_plot
from zmnam_preproc import zmnam_preproc
from zmnam_clean import zmnam_clean

logger = logging.getLogger(__name__)


def get_provenance_record(vatt, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = ("Compute Zonal-mean Northern Annular Modes between "
               "{start_year} and {end_year} ".format(**vatt))
    record = {
        'caption': caption,
        'authors': ['serv_fe', 'hard_jo', 'arno_en', 'cagn_ch'],
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

    Run the zonal-mean NAM diagnostic, calling in order:

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
        zmnam_preproc(ifile)
        outfiles = zmnam_calc(out_dir + '/', out_dir + '/', ifile_props)
        plot_files = zmnam_plot(out_dir + '/', plot_dir +
                                '/', ifile_props, fig_fmt,
                                write_plots)
        provenance_record = get_provenance_record(
            list(input_files.values())[0], ancestor_files=ifile)
        if write_plots:
        # plot_file cannot be an array, so only the first plot is provided
            provenance_record['plot_file'] = plot_files[0]
        zmnam_clean()
        for file in outfiles:
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(file, provenance_record)


# Run the diagnostics
if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
