"""Run the CMORization module."""
import logging
import fnmatch
import os
import sys
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
    return filename


def _run_ncl_script(in_fullpaths,
                    out_fullpaths,
                    dataset,
                    reformat_script,
                    project_info,
                    out_data_dir):
    """Run the NCL cmorization mechanism."""
    for in_fil, out_fil in zip(in_fullpaths, out_fullpaths):
        project_info[dataset]['infile'] = in_fil
        project_info[dataset]['outfile'] = out_fil
        info_file = _write_ncl_infofile(project_info,
                                        dataset, out_data_dir)
        ncl_call = ['ncl', os.path.basename(reformat_script)]
        logger.info("Executing cmd: %s", ' '.join(ncl_call))
        call(ncl_call)
        to_rename_file = os.path.join(
            os.path.dirname(info_file),
            os.path.basename(in_fil).split('_')[0] + '_' +
            os.path.basename(info_file)
        )
        call(['mv', info_file, to_rename_file])


def _run_pyt_script(in_fullpaths, out_fullpaths):
    """Run the Python cmorization mechanism."""
    import py_cmor
    for in_fil, out_fil in zip(in_fullpaths, out_fullpaths):
        py_cmor.cmorization(in_fil, out_fil)


def _get_all_files(files, path, out_data_dir, dataset):
    """Find the I/O files paths."""
    files = fnmatch.filter(files, '*' + dataset + '*')
    in_fullpaths = [os.path.join(path, fil) for fil in files]
    out_fullpaths = [
        os.path.join(out_data_dir,
                     dataset, fil) for fil in files
    ]
    return in_fullpaths, out_fullpaths


def cmor_reformat(config):
    """Run the oldskool v1 cmorization scripts."""
    logger.info("Running the CMORization scripts.")

    # set the reformat scripts dir
    reformat_scripts = os.path.join(os.path.dirname(__file__),
                                    'cmor/cmorize_obs')

    # assemble i/o information
    project_info = {}
    raw_obs = config["rootpath"]["RAWOBS"][0]

    # loop through main obs directory
    for _, datasets, _ in os.walk(raw_obs, followlinks=True):
        for dataset in datasets:
            project_info[dataset] = {}
            reformat_script_root = os.path.join(reformat_scripts,
                                                'cmorize_obs_' + dataset)

            # build directory tree
            data_dir = os.path.join(raw_obs, dataset)
            out_data_dir = os.path.join(config['output_dir'], dataset)
            if not os.path.isdir(out_data_dir):
                os.makedirs(out_data_dir)
            # all operations are done in the working dir now
            os.chdir(out_data_dir)

            # figure out what language the script is in
            if os.path.isfile(reformat_script_root + '.ncl'):
                reformat_script = reformat_script_root + '.ncl'
                logger.info("Attempting to CMORize using NCL script %s",
                            reformat_script)
                # copy over the reformat script
                call(['cp', reformat_script, out_data_dir])
                for path, _, files in os.walk(data_dir, followlinks=True):
                    in_fullpaths, out_fullpaths = _get_all_files(files, path,
                                                                 out_data_dir,
                                                                 dataset)
                    _run_ncl_script(in_fullpaths,
                                    out_fullpaths,
                                    dataset,
                                    reformat_script,
                                    project_info,
                                    out_data_dir)
            elif os.path.isfile(reformat_script_root + '.py'):
                py_reformat_script = reformat_script_root + '.py'
                logger.info("Attempting to CMORize using Python script %s",
                            py_reformat_script)
                # copy over the reformat script
                call(['cp', py_reformat_script,
                      os.path.join(out_data_dir, 'py_cmor.py')])
                sys.path.append(out_data_dir)
                for path, _, files in os.walk(data_dir, followlinks=True):
                    in_fullpaths, out_fullpaths = _get_all_files(files, path,
                                                                 out_data_dir,
                                                                 dataset)
                    _run_pyt_script(in_fullpaths, out_fullpaths)
            else:
                logger.info("No need to CMORize, could not find CMOR script.")
