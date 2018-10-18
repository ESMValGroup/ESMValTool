"""Run the CMORization module."""
import logging
import fnmatch
import os
from subprocess import call

from ._task import write_ncl_settings


logger = logging.getLogger(__name__)


def _write_ncl_infofile(project_info, dataset, output_dir):
    """Write the information needed by the ncl reformat script."""
    info = {
        'input_file_path': project_info[dataset]['infile'],
        'output_file_path': project_info[dataset]['outfile'],
    }

    filename = os.path.join(output_dir, dataset + '_info.ncl')
    if not os.path.isfile(filename):
        write_ncl_settings(info, filename)


def cmor_reformat(config):
    """Run the oldskool v1 cmorization scripts."""
    logger.info("Running the CMORization scripts.")

    # set the reformat scripts dir
    reformat_scripts = os.path.join(os.path.dirname(__file__),
                                    'cmor/cmorize_obs')

    # assemble i/o information
    project_info = {}
    raw_obs = config["rootpath"]["RAWOBS"][0]
    for _, datasets, _ in os.walk(raw_obs, followlinks=True):
        for dataset in datasets:
            project_info[dataset] = {}
            reformat_script = os.path.join(reformat_scripts,
                                           'reformat_obs_' + dataset + '.ncl')
            logger.info("Attempting to CMORize using script %s, if it exists",
                        reformat_script)
            if os.path.isfile(reformat_script):
                logger.info("Script exists, will proceed to CMORize...")
                data_dir = os.path.join(raw_obs, dataset)
                out_data_dir = os.path.join(config['output_dir'], dataset)
                if not os.path.isdir(out_data_dir):
                    os.makedirs(out_data_dir)
                # copy over the reformat script and get in the dir
                call(['cp', reformat_script, out_data_dir])
                os.chdir(out_data_dir)
                for path, _, files in os.walk(data_dir, followlinks=True):
                    files = fnmatch.filter(files, '*' + dataset + '*')
                    in_fullpaths = [os.path.join(path, fil) for fil in files]
                    out_fullpaths = [
                        os.path.join(out_data_dir,
                                     dataset, fil) for fil in files
                    ]
                    for in_fil, out_fil in zip(in_fullpaths, out_fullpaths):
                        project_info[dataset]['infile'] = in_fil
                        project_info[dataset]['outfile'] = out_fil
                        _write_ncl_infofile(project_info,
                                            dataset, out_data_dir)
                        ncl_call = ['ncl', os.path.basename(reformat_script)]
                        logger.info("Executing cmd: %s", ' '.join(ncl_call))
                        call(ncl_call)
            else:
                logger.info("No need to CMORize, no CMOR reformat script.")
