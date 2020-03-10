"""
Ensemble Clustering Diagnostics.

Author: Irene Mavilia (ISAC-CNR, Italy)
Copernicus C3S 34a lot 2 (MAGIC)

Description
    Cluster analysis tool based on the k-means algorithm
    for ensembles of climate model simulations
Modification history
    20181202-vonhardenberg_jost: cleanup, style, provenance and finalising
    20181002-arnone_enrico: updating to version2_develpment (recipe/dataset)
    20170710-mavilia_irene: routines written.
"""

import os
import logging
import numpy as np
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared import ProvenanceLogger, sorted_metadata

# Import user diagnostic routines
from ens_anom import ens_anom
from ens_eof_kmeans import ens_eof_kmeans
from ens_plots import ens_plots

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(gatt, vatt, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = ("Ensemble Clustering Diagnostics of extreme {extreme} of "
               .format(**gatt) + "variable {long_name} between "
               "{start_year} and {end_year} ".format(**vatt))
    print(gatt)
    record = {
        'caption': caption,
        'authors': ['vonhardenberg_jost', 'arnone_enrico', 'mavilia_irene'],
        'projects': ['c3s-magic'],
        'references': ['straus07jcli'],
        'plot_types': ['other'],
        'realms': ['atmos'],
        'domains': ['reg'],
        'ancestors': ancestor_files,
    }
    return record


def main(cfg):
    """Ensemble Clustering Diagnostics."""
    out_dir = cfg['work_dir']
    write_plots = cfg['write_plots']
    input_data = cfg['input_data'].values()
    input_data = sorted_metadata(input_data, sort='recipe_dataset_index')
    files_dict = group_metadata(input_data, 'filename',
                                sort=False)
    numens = len(files_dict)
    logger.info('numens=%d', numens)

    # Building the name of output files
    element = list(files_dict.values())[0][0]
    name_outputs = (element['short_name'] + '_' + str(numens) +
                    'ens_' + cfg['season'] + '_' + cfg['area'] +
                    '_' + element['project'] + '_' + element['exp'])
    logger.info('The name of the output files will be <variable>_%s.txt',
                name_outputs)
    variable_name = element['short_name']
    max_plot_panels = cfg.get('max_plot_panels', 72)
    numpcs = cfg.get('numpcs', 0)
    perc = cfg.get('numpcs', 80)

    filenames_cat = []
    legend_cat = []
    for value in files_dict.values():
        logger.info("Processing file %s", value[0]['filename'])
        filenames_cat.append(value[0]['filename'])
        leg = (value[0]['project'] + " " +
               value[0]['dataset'] + " " +
               value[0]['exp'] + " " +
               value[0]['mip'] + " " +
               value[0]['short_name'] + " " +
               value[0]['ensemble'] + " " +
               str(value[0]['start_year']) + "-" +
               str(value[0]['end_year']))
        legend_cat.append(leg)
        logger.info('Processing: %s', leg)
    namef = os.path.join(out_dir, 'legend_{0}.txt'.format(name_outputs))
    np.savetxt(namef, legend_cat, fmt='%s')

    # ###################### PRECOMPUTATION #######################
    outfiles = ens_anom(filenames_cat, out_dir, name_outputs, variable_name,
                        numens, cfg['season'], cfg['area'], cfg['extreme'])

    # ###################### EOF AND K-MEANS ANALYSES #######################
    outfiles2 = ens_eof_kmeans(out_dir, name_outputs, numens, numpcs,
                               perc, cfg['numclus'])

    outfiles = outfiles + outfiles2
    provenance_record = get_provenance_record(
        cfg, list(files_dict.values())[0][0], ancestor_files=filenames_cat)

    # ###################### PLOT AND SAVE FIGURES ##########################
    if write_plots:
        plotfiles = ens_plots(out_dir, cfg['plot_dir'], name_outputs,
                              cfg['numclus'], 'anomalies',
                              cfg['output_file_type'], cfg['season'],
                              cfg['area'], cfg['extreme'], max_plot_panels)
    else:
        plotfiles = []

    for file in outfiles + plotfiles:
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(file, provenance_record)

    logger.info('\n>>>>>>>>>>>> ENDED SUCCESSFULLY!! <<<<<<<<<<<<\n')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
