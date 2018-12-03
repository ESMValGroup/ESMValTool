"""wrapper diagnostic for the NCAR CVDP (p)ackage."""
import logging
import os
import re
import subprocess

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata)
from esmvaltool._task import DiagnosticError

logger = logging.getLogger(os.path.basename(__file__))


def setup_driver(cfg):
    """Write the driver.ncl file of the cvdp package."""
    cvdp_root = os.path.join(os.path.dirname(__file__), 'cvdp')
    if not os.path.isdir(cvdp_root):
        raise DiagnosticError("CVDP is not available.")

    settings = {
        'outdir': "{0}/".format(cfg['work_dir']),
        'obs': 'False',
        'zp': os.path.join(cvdp_root, "ncl_scripts/"),
        'run_style': 'serial',
        'webpage_title': 'CVDP run via ESMValTool'
    }
    settings['output_data'] = "True" if _nco_available() else "False"

    def _update_settings(line):

        for key, value in settings.items():
            pattern = r'\s*{0}\s*=.*\n'.format(key)
            search_results = re.findall(pattern, line)
            if search_results == []:
                continue
            return re.sub(
                r'".+?"', '"{0}"'.format(value), search_results[0], count=1)

        return line

    content = []
    driver = os.path.join(cvdp_root, "driver.ncl")

    with open(driver, 'r') as driver_file:
        for line in driver_file:
            content.append(_update_settings(line))

    new_driver = os.path.join(cfg['work_dir'], "driver.ncl")

    with open(new_driver, 'w') as new_driver_file:
        new_driver_file.write("".join(content))


def create_link(cfg, inpath):
    """Create link for the input file.

    The link matches the naming convention of the cvdp package.
    Returns the path to the link.

    cfg: configuration dict
    inpath: path to infile
    """
    def _create_link_name(inpath):
        tail = os.path.split(inpath)[1]
        search_result = re.search(r'[0-9]{4}-[0-9]{4}', tail).group(0)
        return tail.replace(search_result,
                            "{0}01-{1}12".format(*search_result.split('-')))

    if not os.path.isdir(inpath):
        raise DiagnosticError("Path {0} does not exist".format(inpath))

    lnk_dir = os.path.join(cfg['work_dir'], "links")

    if not os.path.isdir(lnk_dir):
        os.mkdir(lnk_dir)

    link = os.path.join(lnk_dir, _create_link_name(inpath))
    os.symlink(inpath, link)

    return link


def setup_namelist(cfg):
    """Set the namelist file of the cvdp package."""
    input_data = cfg['input_data'].values()
    selection = select_metadata(input_data, project='CMIP5')
    grouped_selection = group_metadata(selection, 'dataset')

    content = []
    for key, value in grouped_selection.items():
        links = [create_link(cfg, item["filename"]) for item in value]
        head, tail = os.path.split(links[0])
        head, tail = os.path.split(head)
        tail = "_".join(tail.split('_')[:-1])
        ppath = "{}*/".format(os.path.join(head, tail))
        content.append("{0} | {1} | {2} | {3}\n".format(
            key, ppath, value[0]["start_year"], value[0]["end_year"]))

    namelist = os.path.join(cfg['work_dir'], "namelist")

    with open(namelist, 'w') as namelist_file:
        namelist_file.write("\n".join(content))


def log_functions(func):
    """Decorater to check functions."""
    def inner():
        """Inner function."""
        ret = func()
        logger.debug("Function %s returns %s", func.__name__, str(ret))
        return ret

    return inner


@log_functions
def _nco_available():
    """Check if nco is available."""
    try:
        retcode = subprocess.call("which ncks", shell=True)
        if retcode < 0:
            ret = False
        else:
            ret = True
    except OSError:
        ret = False
    return ret


def main(cfg):
    """Set and execute the cvdp package."""
    setup_driver(cfg)
    setup_namelist(cfg)
    subprocess.run(["ncl", "driver.ncl"], cwd=os.path.join(cfg['work_dir']))


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
