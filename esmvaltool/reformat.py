"""Port reformat.py from v1 to v2 (life hack)."""
import fnmatch
import os
from subprocess import call


def cmor_reformat(config, logvar):
    """Run the oldskool v1 reformat scripts."""
    logvar.info("Running the oldschool reformatting scripts.")

    # set the reformat scripts dir
    reformat_scripts = os.path.join(os.path.dirname(__file__),
                                    'reformat_scripts')

    # assemble i/o information
    project_info = {}
    raw_obs = config["rootpath"]["RAWOBS"][0]
    for _, datasets, _ in os.walk(raw_obs, followlinks=True):
        for dataset in datasets:
            project_info[dataset] = {}
            reformat_script = os.path.join(reformat_scripts,
                                           'reformat_obs_' + dataset + '.ncl')
            logvar.info("Attempting to reformat using script %s, if it exists",
                        reformat_script)
            if os.path.isfile(reformat_script):
                logvar.info("Script exists, will proceed to reformatting...")
                data_dir = os.path.join(raw_obs, dataset)
                out_data_dir = os.path.join(config['output_dir'], dataset)
                if not os.path.isdir(out_data_dir):
                    os.makedirs(out_data_dir)
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
                        ncl_call = ['ncl', reformat_script]
                        logvar.info("Executing cmd: %s", ' '.join(ncl_call))
                        call(ncl_call)
            else:
                logvar.info("No need to reformat, no reformat script.")
