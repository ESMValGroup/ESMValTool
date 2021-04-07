"""wrapper diagnostic for the NCAR CVDP (p)ackage."""
import logging
import os
import re
import shutil
import subprocess

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic)
from esmvaltool.diag_scripts.shared import ProvenanceLogger

logger = logging.getLogger(os.path.basename(__file__))


class DiagnosticError(Exception):
    """Error in diagnostic."""


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
            return re.sub(r'".+?"',
                          '"{0}"'.format(value),
                          search_results[0],
                          count=1)

        return line

    content = []
    driver = os.path.join(cvdp_root, "driver.ncl")

    with open(driver, 'r') as driver_file:
        for line in driver_file:
            content.append(_update_settings(line))

    new_driver = os.path.join(cfg['run_dir'], "driver.ncl")

    with open(new_driver, 'w') as new_driver_file:
        new_driver_file.write("".join(content))


def create_link(cfg, inpath, _name):
    """Create link for the input file.

    The link matches the naming convention of the cvdp package.
    Returns the path to the link.

    cfg: configuration dict
    inpath: path to infile
    """
    def _create_link_name():
        tail = os.path.split(inpath)[1]
        search_result = re.search(r'[0-9]{4}-[0-9]{4}', tail).group(0)
        return _name + "_" + tail.replace(
            search_result, "{0}01-{1}12".format(*search_result.split('-')))

    if not os.path.isfile(inpath):
        raise DiagnosticError("Path {0} does not exist".format(inpath))

    lnk_dir = cfg['lnk_dir']

    if not os.path.isdir(lnk_dir):
        os.mkdir(lnk_dir)

    link = os.path.join(lnk_dir, _create_link_name())
    if not os.path.exists(link):
        os.symlink(inpath, link)

    return link


def setup_namelist(cfg):
    """Set the namelist file of the cvdp package."""
    input_data = cfg['input_data'].values()
    grouped_selection = group_metadata(input_data, 'alias')

    content = []
    for _, attributes in grouped_selection.items():
        for item in attributes:
            create_link(cfg, item["filename"], item['alias'])
        ppath = "{0}/".format(cfg['lnk_dir'])
        content.append("{0} | {1}{0} | {2} | {3}\n".format(
            attributes[0]["alias"], ppath, attributes[0]["start_year"],
            attributes[0]["end_year"]))

    namelist = os.path.join(cfg['run_dir'], "namelist")

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
        if shutil.which("ncks") is None:
            ret = False
        else:
            ret = True
    except OSError:
        ret = False
    return ret


def _is_png(path):
    exclude = ['cas-cvdp.png']
    filename = os.path.basename(path)
    return filename.endswith('.png') and filename not in exclude


def _get_caption(filename):
    caption = []
    stat = _get_stat(filename)
    if stat is not None and stat != "other":
        caption.append(stat)
    season = _get_season(filename)
    if season is not None:
        caption.append(season)
    long_name = _get_long_name(filename)
    if long_name is not None:
        caption.append(long_name)
    mode = _get_mode(filename)
    if mode is not None:
        caption.append(mode)
    return " ".join(caption)


def _get_plot_type(filename):
    plot_type = {
        'timeseries': "times",
        'mean': "other",
        'stddev': "other",
        'trends': "other",
        'eight_yr_runtrend': "other",
        'sixteen_yr_runtrend': "other",
        'fourteen_yr_runtrend': "other",
        'twelve_yr_runtrend': "other",
        'ten_yr_runtrend': "other",
        'powspec': "other",
        'reg': "other",
        'hov': "other",
        'monstddev': "other",
        'runstddev': "other",
        'za': "zonal",
    }
    ans = _get_info(filename, plot_type)
    return ans if ans is not None else 'other'


def _get_stat(filename):
    stat = {
        'timeseries': "other",
        'mean': "mean",
        'stddev': "stddev",
        'trends': "trend",
        'eight_yr_runtrend': "trend",
        'sixteen_yr_runtrend': "trend",
        'fourteen_yr_runtrend': "trend",
        'twelve_yr_runtrend': "trend",
        'ten_yr_runtrend': "trend",
        'powspec': "spectrum",
        'reg': "other",
        'hov': "other",
        'monstddev': "stddev",
        'runstddev': "stddev",
        'za': "mean",
    }
    ans = _get_info(filename, stat)
    return ans if ans is not None else 'other'


def _get_season(filename):
    season = {
        'ann': "Annual",
        'djf': "DJF",
        'mam': "MAM",
        'jja': "JJA",
        'son': "SON",
    }
    return _get_info(filename, season)


def _get_long_name(filename):
    variable = {
        'pr': "Precipitation",
        'tas': "Surface temperature",
        'psl': "Sea level pressure",
        'sst': "Sea surface temperature",
    }
    return _get_info(filename, variable)


def _get_mode(filename):
    mode = {
        'iod': "iod",
        'ipo': "ipo",
        'nam': "nam",
        'nao': "nao",
        'lanina': "La nina",
        'nino12': "El nino 12",
        'nino3': "El nino 3",
        'nino34': "El nino 34",
        'nino4': "El nino 4",
        'npi': "npi",
        'npo': "npo",
        'pdo': "pdo",
        'pna': "pna",
        'psa1': "psa1",
        'psa2': "psa2",
        'sam': "sam",
        'socn': "socn",
        'tio': "tio",
        'tna': "tna",
        'tsa': "tsa",
    }
    return _get_info(filename, mode)


def _get_info(filename, dictionary):
    intersection = list(
        set(os.path.basename(filename).split('.')).intersection(
            dictionary.keys()))
    if len(intersection) != 1:
        return None
    return dictionary[intersection[0]]


def _get_global_ancestors(cfg):
    input_data = cfg['input_data'].values()
    grouped_selection = group_metadata(input_data, 'dataset')
    ancestor = []
    for _, attributes in grouped_selection.items():
        ancestor += [item['filename'] for item in attributes]
    return ancestor


def set_provenance(cfg):
    """Add provenance to all image files that the cvdp package creates."""
    def _get_provenance_record(filename, ancestors):
        return {
            'caption': _get_caption(filename),
            'statistics': [_get_stat(filename)],
            'domain': 'global',
            'plot_type': _get_plot_type(filename),
            'plot_file': filename,
            'authors': [
                'phillips_adam',
            ],
            'references': [
                'acknow_project',
                'phillips14eos',
            ],
            'ancestors': ancestors,
        }

    ancestors = _get_global_ancestors(cfg)
    logger.info("Path to work_dir: %s", cfg['work_dir'])
    with ProvenanceLogger(cfg) as provenance_logger:
        for root, _, files in os.walk(cfg['work_dir']):
            for datei in files:
                path = os.path.join(root, datei)
                if _is_png(path):
                    logger.info("Name of file: %s", path)
                    provenance_record = _get_provenance_record(path, ancestors)
                    logger.info("Recording provenance of %s:\n%s", path,
                                provenance_record)
                    provenance_logger.log(path, provenance_record)


def _execute_cvdp(cfg):
    subprocess.check_call(["ncl", "driver.ncl"],
                          cwd=os.path.join(cfg['run_dir']))


def main(cfg):
    """Set and execute the cvdp package."""
    cfg['lnk_dir'] = os.path.join(cfg['run_dir'], "links")
    setup_driver(cfg)
    setup_namelist(cfg)
    _execute_cvdp(cfg)
    set_provenance(cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
