"""
TO BE WRITTEN...
Authors: Valeriu Predoi, University of Reading, valeriu.predoi@ncas.ac.uk
         Mattia Righi, DLR, mattia.righi@dlr.de
"""

# ---- Import standard modules to the python path.
import logging
import os
import subprocess
from datetime import datetime

logger = logging.getLogger(__name__)


def cmip5_model2inst(model):
    """
    Return the institute given the model name in CMIP5
    """

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
    """
    Returns realm and frequency given the mip in CMIP5
    """

    mipdict: {
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
    

def get_input_filelist(project_info, model, var): ## FIX-ME


    # project_info is supposed to contain all the information from the config file
    # Here in particular rootpath and drs are needed (also as dictionaries)
    # Order of setting variables:
    # - first get the drs, so we condition on it if at all
    # - second, set the root depending on drs

    ## Directory structure is defined in config and also project-depentend (at present only used for CMIP5)
    key_drs = 'drs_' + model['project']
    drs = project_info['GLOBAL'].get(key_drs)

    ## The rootpath is defined in config but now is project-dependent!
    ## Define as dictionary (is this possible in config.ini? Otherwise switch to yaml...)
    # no dictionary at the moment, just individual subscripts

    # if drs at all
    if drs == 'BADC' or drs == 'DKRZ':
        root = project_info['GLOBAL']['host_root']
    # FIXME check cases ETHZ, SMHI...
    else:
        precise_rootpaths = ['rootpath_CMIP5', 'rootpath_OBS', 'rootpath_obs4mips']
        string = 'rootpath_' + model['project']
        if string in precise_rootpaths:
            root = project_info['GLOBAL'][string]
        elif 'rootpath_default' in project_info['GLOBAL']:
            root = project_info['GLOBAL']['rootpath_default']
        else:
            raise ValueError("rootpath for project {} not defined in config".format(model['project']))
    
    logger.info("Root path set to: %s", root)

    ## Trying to implement the variable-dependent model keys (see section 4.1 yaml document)
    ## Basic idea: same model used by different variable with different mip/ensemble/exp. In
    ##             this case the mip/exp/ensemble is defined in the variable dictionary, while
    ##             the model dictionary contains a wildcard (mip: *)
    if 'mip' in var:
        model['mip'] = var['mip']
    if 'ensemble' in var:
        model['ensemble'] = var['ensemble']
    if 'exp' in var:
        model['exp'] = var['exp']

    module = globals()[model['project']]()
    files = module.infile_path(root, model, var, drs)

    return files


def get_output_file(project_info, model, var):
    
    module = globals()[model['project']]()
    files = module.outfile_path(project_info['GLOBAL']['preproc_dir'], model, var)

    return files

class DataFinder:
    """
    Class for input and output file name based on the project class
    rootpath is the machine-dependent root path specified in config.ini
    model is the model dictionary in main namelist
    var is the variable dictionary in the main namelist
    drs is a flag in config.ini (only used in for CMIP5)
    """

    # The input filename is constructed as
    # rootpath (project dependent, explicitly defined in config) /
    # drs (project dependent, machine dependent, set with label in config) /
    # file name (project dependent, fully determined by the model and variable dicts)
    def infile_path(self, rootpath, model, var, drs):
        return None

    # The output filename is constructed as
    # rootpath (preproc_dir, defined in config) /
    # project name /
    # file name (project dependent, fully determined by the model and variable dicts)
    def outfile_path(self, rootpath, model, var):  # this shall replace get_cf_outfile in preprocess.py
        return None

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

class CMIP5(DataFinder):

    def infile_path(self, rootpath, model, var, drs):

        if (drs == 'BADC'):

            logger.info('Getting CMIP5 data from BADC')

            # Latest version is always called 'latest' at BADC
            latest = 'latest'

            dirname = '/'.join([rootpath,              # CHECK-ME
                                cmip5_model2inst(model['name']),
                                model['name'],
                                model['exp'],
                                cmip5_mip2realm_freq(model['mip'])[1],
                                cmip5_mip2realm_freq(model['mip'])[0],
                                model['mip'],
                                model['ensemble'],
                                latest,
                                var['name']])

        elif (drs == 'DKRZ'):

            logger.info('Getting CMIP5 data from DKRZ')

            ## host_root = '/mnt/lustre01/work/kd0956/CMIP5/data/cmip5/output1/'  ## must be given in config, not here

            # Automatically find the latest version by sorting
            dirname = '/'.join([rootpath,
                                cmip5_model2inst(model['name']),
                                model['name'],
                                model['exp'],
                                cmip5_mip2realm_freq(model['mip'])[0],
                                cmip5_mip2realm_freq(model['mip'])[1],
                                model['mip'],
                                model['ensemble']])

            list_versions = os.listdir(dirname)
            list_versions.sort()
            latest = os.path.basename(list_versions[-1])

            dirname = '/'.join([dirname,
                                latest,
                                var['name']])

        elif (drs == 'ETHZ'):

            logger.info('Getting CMIP5 data with ETHZ root format')

            dirname = '/'.join([rootpath, 
                                model['exp'], 
                                model['mip'], 
                                var, model['name'], 
                                model['ensemble']]) + '/'

        elif (drs == 'SMHI'):
            dirname = '/'.join([rootpath,
                                model["name"],
                                model["ensemble"],
                                model["exp"],
                                cmip5_mip2realm_freq(model['mip'])[1]])

        elif drs is None:

            logger.info('No DRS rootpath specified')

            dirname = rootpath
            if (not dirname.endswith('/')):
                dirname = dirname + '/'

        else:

            raise ValueError("Cannot establish path to CMIP5 files")


        # The CMIP5 filename is always the same, only the drs changes
        filename = '_'.join([var['name'],
                             model['mip'],
                             model['name'],
                             model['exp'],
                             model['ensemble']]) + '*.nc'
        
        # use find here
        in_files_list = veto_files(model, dirname, filename)

        return in_files_list

    def outfile_path(self, rootpath, model, var):
        out_files_path =  rootpath + '/' + model['project'] + '/' + '_'.join([model['project'],
                                     model['name'],
                                     model['mip'],
                                     model['exp'],
                                     model['ensemble'],
                                     var['field'],
                                     var['name'],
                                     str(model['start_year']) + '-' + str(model['end_year'])]) + '.nc'

        return out_files_path

class EMAC(DataFinder):

    def infile_path(self, rootpath, model, var, drs):
        path = '/'.join(rootpath, model['name'])
        return path

    def outfile_path(self, rootpath, model, var):
        path = rootpath + '/' + model['project'] + '/' + '_'.join([model['project'],
                                    model['name'],
                                    model['ensemble'],
                                    var['field'],
                                    var['name'],
                                    str(model['start_year']) + '-' + str(model['end_year'])]) + '.nc'
        return path


class OBS(DataFinder):

    def infile_path(self, rootpath, model, var, drs):
        dirname = '/'.join([rootpath, 
                            'Tier' + str(model['tier']),
                            model['name']])
        filename = '_'.join([model['project'], 
                         model['name'],
                         str(model['type']), 
                         str(model['version']),
                         var['field'],
                         var['name']]) + '*.nc'

        # use find here
        in_files_list = veto_files(model, dirname, filename)

        return in_files_list

    def outfile_path(self, rootpath, model, var):
        out_files_path = rootpath + '/' + model['project'] + '/' + '_'.join([model['project'],
                                     model['name'],
                                     str(model['type']),
                                     str(model['version']),
                                     var['field'],
                                     var['name'],
                                     str(model['start_year']) + '-' + str(model['end_year'])]) + '.nc'
        return out_files_path


## TO BE ADDED ##############
# class obs4mips(DataFinder):
# class ana4mips(DataFinder):
# class MiKlip(DataFinder):
# class CCMval1(DataFinder):
# class CCMval2(DataFinder):
# class ECEARTH(DataFinder):
# class GFDL(DataFinder):
#############################

### FIX-ME: below is the original code, I didn't change anything there.
   


# ---- handling the years for files
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
