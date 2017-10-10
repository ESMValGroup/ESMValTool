# Data finder module for the ESMValTool
# Authors:
# Valeriu Predoi (URead, UK - valeriu.predoi@ncas.ac.uk)
# Mattia Righi (DLR, Germany - mattia.righi@dlr.de)

import logging
import os
import subprocess
import yaml
import re
from datetime import datetime
import iris
import iris.coord_categorisation
import six

logger = logging.getLogger(__name__)


def cmip5_model2inst(model):
    """ Return the institute given the model name in CMIP5 """

    instdict = {
        'HadGEM2-CC': 'MOHC',
        'HadGEM2-A': 'MOHC',
        'HadCM3': 'MOHC',
        'HadGEM2-ES': 'MOHC',
        'FIO-ESM': 'FIO',
        'fio-esm': 'FIO',
        'CCSM4': 'NCAR',
        'GEOS-5': 'NASA-GMAO',
        'inmcm4': 'INM',
        'CanESM2': 'CCCma',
        'CanCM4': 'CCCma',
        'CanAM4': 'CCCma',
        'GISS-E2-R': 'NASA-GISS',
        'GISS-E2-R-CC': 'NASA-GISS',
        'GISS-E2-H-CC': 'NASA-GISS',
        'GISS-E2-H': ' NASA-GISS',
        'CNRM-CM5': 'CNRM-CERFACS',
        'CNRM-CM5-2': 'CNRM-CERFACS',
        'NICAM-09': 'NICAM',
        'IPSL-CM5A-LR': 'IPSL',
        'IPSL-CM5A-MR': 'IPSL',
        'IPSL-CM5B-LR': 'IPSL',
        'CSIRO-Mk3-6-0': 'CSIRO-QCCCE',
        'CESM1-CAM5': 'NSF-DOE-NCAR',
        'CESM1-CAM5-1-FV2': 'NSF-DOE-NCAR',
        'CESM1-BGC': 'NSF-DOE-NCAR',
        'CESM1-WACCM': 'NSF-DOE-NCAR',
        'CESM1-FASTCHEM': 'NSF-DOE-NCAR',
        'NorESM1-M': 'NCC',
        'NorESM1-ME': 'NCC',
        'CFSv2-2011': 'NOAA-NCEP',
        'ACCESS1-3': 'CSIRO-BOM',
        'ACCESS1-0': 'CSIRO-BOM',
        'CMCC-CM': 'CMCC',
        'CMCC-CESM': 'CMCC',
        'CMCC-CMS': 'CMCC',
        'FGOALS-g2': 'LASG-CESS',
        'FGOALS-s2': 'LASG-IAP',
        'FGOALS-gl': 'LASG-IAP',
        'GFDL-HIRAM-C180': 'NOAA-GFDL',
        'GFDL-ESM2G': 'NOAA-GFDL',
        'GFDL-CM2p1': 'NOAA-GFDL',
        'GFDL-CM3': 'NOAA-GFDL',
        'GFDL-ESM2M': 'NOAA-GFDL',
        'GFDL-HIRAM-C360': 'NOAA-GFDL',
        'EC-EARTH': 'ICHEC',
        'BNU-ESM': 'BNU',
        'CFSv2-2011': 'COLA-CFS',
        'HadGEM2-AO': 'NIMR-KMA',
        'MIROC4h': 'MIROC',
        'MIROC5': 'MIROC',
        'MIROC-ESM': 'MIROC',
        'MIROC-ESM-CHEM': 'MIROC',
        'bcc-csm1-1': 'BCC',
        'bcc-csm1-1-m': 'BCC',
        'HadGEM2-ES': 'INPE',
        'MPI-ESM-LR': 'MPI-M',
        'MPI-ESM-MR': 'MPI-M',
        'MPI-ESM-P': 'MPI-M',
        'MRI-AGCM3-2H': 'MRI',
        'MRI-CGCM3': 'MRI',
        'MRI-ESM1': 'MRI',
        'MRI-AGCM3-2S': 'MRI',
    }

    if model in instdict:
        return instdict[model]

    raise KeyError("CMIP5: cannot map model {} to institute".format(model))


def cmip5_mip2realm_freq(mip):
    """ Returns realm and frequency given the mip in CMIP5 """

    mipdict = {
        'Amon': ['atmos', 'mon'],
        'Omon': ['ocean', 'mon'],
        'Lmon': ['land', 'mon'],
        'LImon': ['landIce', 'mon'],
        'OImon': ['seaIce', 'mon'],
        'aero': ['aerosol', 'mon'],
        #'3hr': ???
        'cfDay': ['atmos', 'day'],
        'cfMon': ['atmos', 'mon'],
        'day': ['atmos', 'day'],
        'fx': ['*', 'fx']
    }

    if mip in mipdict:
        return mipdict[mip]

    raise KeyError("CMIP5: cannot map mip {} to realm".format(mip))


def replace_tags(path, model, var):
    """ Replaces tags in the config-developer's file with actual values """

    path = path.strip('/')

    tlist = re.findall(r'\[([^]]*)\]', path)

    for tag in tlist:

        if tag == 'var':
            replacewith = var['name']
        elif tag == 'field':
            replacewith = var['field']
        elif tag in ('institute', 'freq', 'realm'):
            if tag in model:
                replacewith = str(model[tag])
            else:
                if tag == 'institute':
                    replacewith = cmip5_model2inst(model['name'])
                elif tag == 'freq':
                    replacewith = cmip5_mip2realm_freq(model['mip'])[1]
                elif tag == 'realm':
                    replacewith = cmip5_mip2realm_freq(model['mip'])[0]
        elif tag == 'latestversion':  # handled separately later
            pass
        elif tag == 'tier':
            replacewith = ''.join(('Tier', str(model['tier'])))
        elif tag == 'model':
            replacewith = model['name']
        else:  # all other cases use the corrsponding model dictionary key
            if tag in model:
                replacewith = str(model[tag])
            else:
                raise KeyError('Model key %s must be specified for project %s, check your namelist entry' % (tag, model['project']))

        path = path.replace('[' + tag + ']', replacewith)

    return path


def read_config_file(project, cfg_file=None):
    """ Parses the developer's configuration file and returns the dictionary
        for the given project
    """

    dict = {}
    if (cfg_file is None):
        cfg_file = os.path.join(os.path.dirname(__file__), '../config-developer.yml')
        dict = yaml.load(file(cfg_file, 'r'))

    if project in dict:
        return dict[project]

    raise KeyError('Specifications for {} not found in config-developer file'.format(project))


def get_input_filelist(project_info, model, var):
    """ Returns the full path to input files
    """

    project = model['project']

    project_config = read_config_file(project)

    # Apply variable-dependent model keys
    if 'mip' in var:
        model['mip'] = var['mip']
    if 'ensemble' in var:
        model['ensemble'] = var['ensemble']
    if 'exp' in var:
        model['exp'] = var['exp']

    # Set the rootpath
    dir1 = _get_option_with_default(project_info['GLOBAL'], 'rootpath', project, 'user config')
    if not os.path.isdir(dir1):
        raise OSError('directory not found', dir1)

    # Set the drs
    if project in project_info['GLOBAL']['drs']:
        drs = project_info['GLOBAL']['drs'][project]
    else:
        drs = 'default'

    input_folder = _get_option_with_default(project_config, 'input_dir', drs, 'developer config for %s' % project)
    dir2 = replace_tags(input_folder, model, var)
    dirname = os.path.join(dir1, dir2)

    # Find latest version if required
    if '[latestversion]' in dirname:
        part1 = dirname.split('[latestversion]')[0]
        part2 = dirname.split('[latestversion]')[1]
        list_versions = os.listdir(part1)
        list_versions.sort()
        latest = os.path.basename(list_versions[-1])
        dirname = os.path.join(part1, latest, part2)

    if not os.path.isdir(dirname):
        raise OSError('directory not found', dirname)

    # Set the filename
    input_file = _get_option_with_default(project_config, 'input_file', drs, 'developer config for %s' % project)
    filename = replace_tags(input_file, model, var)

    # Full path to files
    files = veto_files(model, dirname, filename)

    return files

def _get_option_with_default(config, key, option, config_name):
    config_value = config[key]

    if isinstance(config_value, six.string_types):
        value = config_value
    elif option in config_value:
        value = config_value[option]
    elif 'default' in config_value:
        value = config_value['default']
    else:
        raise KeyError('Option %s not specified and no default provided for %s in %s' % (option, key, config_name))
    return value

def get_output_file(project_info, model, var):
    """ Returns the full path to the output (preprocessed) file
    """

    dict = read_config_file(model['project'])

    outfile = os.path.join(project_info['GLOBAL']['preproc_dir'],
                           model['project'],
                           replace_tags(dict['output_file'], model, var))
    outfile = ''.join((outfile, '.nc'))

    return outfile


def find_files(dirname, filename):
    """
    Function that performs local search for files using `find'
    The depth is as high as possible so that find is fast
    """
    flist = []

    # work only with existing dirs or allowed permission dirs
    strfindic = 'find ' + dirname +' -follow -type f -iname ' + '*' + filename + '*'
    proc = subprocess.Popen(strfindic, stdout=subprocess.PIPE, shell=True)
    out, err = proc.communicate()
    if err:
        logger.warning("'%s' says:\n%s", strfindic, err)
    for t in out.split('\n')[0:-1]:
        flist.append(t)
    return flist


def veto_files(model, dirname, filename):
    """
    Function that does direct parsing of available datasource files and establishes
    if files are the needed ones or not

    """

    arname = find_files(dirname, filename)
    fs = []

    if len(arname) > 0:
        yr1 = int(model['start_year'])
        yr2 = int(model['end_year'])
        for s in arname:
            tc = time_check(s, yr1, yr2)
            if tc is True:
                fs.append(s)

    return fs


def time_handling(year1, year1_model, year2, year2_model):
    """
    This function is responsible for finding the correct
    files for the needed timespan:

    year1 - the start year in files
    year1_model - the needed start year of data
    year2 - the last year in files
    year2_model - the needed last year of data
    WARNINGS:
    we reduce our analysis only to years

    """
    # model interval < data interval / file
    # model requirements completely within data stretch
    if year1 <= int(year1_model) and year2 >= int(year2_model):
        return True
    # model interval > data interval / file
    # data stretch completely within model requirements
    elif year1 >= int(year1_model) and year2 <= int(year2_model):
        return True
    # left/right overlaps and complete misses
    elif year1 <= int(year1_model) and year2 <= int(year2_model):
        # data is entirely before model
        if year2 < int(year1_model):
            return False
        # edge on
        elif year2 == int(year1_model):
            return True
        # data overlaps to the left
        elif year2 > int(year1_model):
            return True
    elif year1 >= int(year1_model) and year2 >= int(year2_model):
        # data is entirely after model
        if year1 > int(year2_model):
            return False
        # edge on
        elif year1 == int(year2_model):
            return True
        # data overlaps to the right
        elif year1 < int(year2_model):
            return True

# ---- function to handle various date formats
def date_handling(time1,time2):
    """
    This function deals with different input date formats e.g.
    time1 = 198204 or
    time1 = 19820422 or
    time1 = 198204220511 etc
    More formats can be coded in at this stage.
    Returns year 1 and year 2
    """
    # yyyymm
    if len(list(time1)) == 6 and len(list(time2)) == 6:
        y1 = datetime.strptime(time1, '%Y%m')
        year1 = y1.year
        y2 = datetime.strptime(time2, '%Y%m')
        year2 = y2.year
    else:
        # yyyymmdd
        if len(list(time1)) == 8 and len(list(time2)) == 8:
            y1 = datetime.strptime(time1, '%Y%m%d')
            year1 = y1.year
            y2 = datetime.strptime(time2, '%Y%m%d')
            year2 = y2.year
        # yyyymmddHHMM
        if len(list(time1)) == 12 and len(list(time2)) == 12:
            y1 = datetime.strptime(time1, '%Y%m%d%H%M')
            year1 = y1.year
            y2 = datetime.strptime(time2, '%Y%m%d%H%M')
            year2 = y2.year
    return year1,year2

# ---- function that does time checking on a file
def time_check(fpath, yr1, yr2):
    """
    fpath: full path to file
    yr1, yr2: model['start_year'], model['end_year']
    """
    try:
        ssp = fpath.split('/')
        av = ssp[-1]
        time_range = av.split('_')[-1].strip('.nc')
        time1 = time_range.split('-')[0]
        time2 = time_range.split('-')[1]
        year1 = date_handling(time1,time2)[0]
        year2 = date_handling(time1,time2)[1]
        if time_handling(year1, yr1, year2, yr2) is True:
            return True
        else:
            return False
    except Exception:
        cubes = iris.load(fpath)
        for cube in cubes:
            if cube.coords('time'):
                iris.coord_categorisation.add_year(cube, 'time')
                extracted = cube.extract(iris.Constraint(year=lambda year: yr1 <= year <= yr2))
                if extracted is None:
                    return False
                else:
                    return True
        return False

