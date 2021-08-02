"""
Zonal mean annular mode main routine.

Author: Federico Serva (ISAC-CNR & ISMAR-CNR, Italy)
Copernicus C3S 34a lot 2 (MAGIC)

Description:
Evaluation of stratosphere-troposphere coupling
based on EOF/PC analysis of the geopotential height field.

Modification history
202107-serva_federico: Update to include hemisphere selection.
20180512-serva_federico: Added output netCDFs, more use of preprocessor.
20180510-serva_federico: Routines written.

"""

import logging
import os

from esmvaltool.diag_scripts.shared import ProvenanceLogger, run_diagnostic
# Import zmnam diagnostic routines
from zmnam_calc import zmnam_calc
from zmnam_plot import zmnam_plot
from zmnam_preproc import zmnam_preproc

logger = logging.getLogger(__name__)


def get_provenance_record(caption, vatt, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = (caption + " between "
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
    Run the zonal mean annular mode diagnostic.

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

    hemispheres = cfg['hemisphere']

    # Go to work_dir for running
    os.chdir(out_dir)

    # Loop over input cfg
    for ifile, props in input_files.items():

        # Collect relevant information for outputs naming
        ifile_props = [
            props['project'], props['dataset'], props['exp'],
            props['ensemble'],
            str(props['start_year']) + '-' + str(props['end_year'], )
        ]

        for hemisphere in hemispheres:

            # Call diagnostics functions
            print("prepro")
            (file_da_an_zm, file_mo_an) = zmnam_preproc(ifile, hemisphere)
            print("calc")
            outfiles = zmnam_calc(file_da_an_zm, out_dir + '/', ifile_props)

            if write_plots:
                print("plot_files")
                plot_files = zmnam_plot(file_mo_an, out_dir + '/',
                                        plot_dir + '/', ifile_props,
                                        fig_fmt, write_plots, hemisphere)
            else:
                plot_files = []
            for file in outfiles + plot_files:

                if 'pc_da' in file:
                    caption = 'Daily principal component timeseries'
                elif 'pc_mo' in file or 'mo_ts' in file:
                    caption = 'Monthly principal component timeseries'
                elif 'mo_reg' in file:
                    caption = 'Monthly regression maps'
                elif 'da_pdf' in file:
                    caption = 'Daily probability distribution function'
                elif 'eofs' in file:
                    caption = 'Empirical orthogonal functions'
                else:  # please add additional diagnostic description
                    caption = 'Not specified diagnostic'

                provenance_record = get_provenance_record(caption, props,
                                                          ancestor_files=
                                                          [ifile])
                with ProvenanceLogger(cfg) as provenance_logger:
                    provenance_logger.log(file, provenance_record)


# Run the diagnostics
if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
