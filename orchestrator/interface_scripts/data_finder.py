"""
TO BE WRITTEN...
Authors: Valeriu Predoi, University of Reading, valeriu.predoi@ncas.ac.uk
         Mattia Righi, DLR, mattia.righi@dlr.de
"""

# ---- Import standard modules to the python path.
import sys, os
import subprocess
from datetime import datetime

def cmip5_model2inst(model): ## CHECK-ME: A dictionary is preferred to avoid using find, which causes some issues on some machines in the past (too slow)
    """
    Return the institute given the model name in CMIP5
    """

    instdict = {}
    instdict['HadGEM2-CC'] = 'MOHC'
    instdict['HadGEM2-A'] = 'MOHC'
    instdict['HadCM3'] = 'MOHC'
    instdict['HadGEM2-ES'] = 'MOHC'
    instdict['FIO-ESM'] = 'FIO'
    instdict['fio-esm'] = 'FIO'
    instdict['CCSM4'] = 'NCAR'
    instdict['GEOS-5'] = 'NASA-GMAO'
    instdict['inmcm4'] = 'INM'
    instdict['CanESM2'] = 'CCCma'
    instdict['CanCM4'] = 'CCCma'
    instdict['CanAM4'] = 'CCCma'
    instdict['GISS-E2-R'] = 'NASA-GISS'
    instdict['GISS-E2-R-CC'] = 'NASA-GISS'
    instdict['GISS-E2-H-CC'] = 'NASA-GISS'
    instdict['GISS-E2-H'] = ' NASA-GISS'
    instdict['CNRM-CM5'] = 'CNRM-CERFACS'
    instdict['CNRM-CM5-2'] = 'CNRM-CERFACS'
    instdict['NICAM-09'] = 'NICAM'
    instdict['IPSL-CM5A-LR'] = 'IPSL'
    instdict['IPSL-CM5A-MR'] = 'IPSL'
    instdict['IPSL-CM5B-LR'] = 'IPSL'
    instdict['CSIRO-Mk3-6-0'] = 'CSIRO-QCCCE'
    instdict['CESM1-CAM5'] = 'NSF-DOE-NCAR'
    instdict['CESM1-CAM5-1-FV2'] = 'NSF-DOE-NCAR'
    instdict['CESM1-BGC'] = 'NSF-DOE-NCAR'
    instdict['CESM1-WACCM'] = 'NSF-DOE-NCAR'
    instdict['CESM1-FASTCHEM'] = 'NSF-DOE-NCAR'
    instdict['NorESM1-M'] = 'NCC'
    instdict['NorESM1-ME'] = 'NCC'
    instdict['CFSv2-2011'] = 'NOAA-NCEP'
    instdict['ACCESS1-3'] = 'CSIRO-BOM'
    instdict['ACCESS1-0'] = 'CSIRO-BOM'
    instdict['CMCC-CM'] = 'CMCC'
    instdict['CMCC-CESM'] = 'CMCC'
    instdict['CMCC-CMS'] = 'CMCC'
    instdict['FGOALS-g2'] = 'LASG-CESS'
    instdict['FGOALS-s2'] = 'LASG-IAP'
    instdict['FGOALS-gl'] = 'LASG-IAP'
    instdict['GFDL-HIRAM-C180'] = 'NOAA-GFDL'
    instdict['GFDL-ESM2G'] = 'NOAA-GFDL'
    instdict['GFDL-CM2p1'] = 'NOAA-GFDL'
    instdict['GFDL-CM3'] = 'NOAA-GFDL'
    instdict['GFDL-ESM2M'] = 'NOAA-GFDL'
    instdict['GFDL-HIRAM-C360'] = 'NOAA-GFDL'
    instdict['EC-EARTH'] = 'ICHEC'
    instdict['BNU-ESM'] = 'BNU'
    instdict['CFSv2-2011'] = 'COLA-CFS'
    instdict['HadGEM2-AO'] = 'NIMR-KMA'
    instdict['MIROC4h'] = 'MIROC'
    instdict['MIROC5'] = 'MIROC'
    instdict['MIROC-ESM'] = 'MIROC'
    instdict['MIROC-ESM-CHEM'] = 'MIROC'
    instdict['bcc-csm1-1'] = 'BCC'
    instdict['bcc-csm1-1-m'] = 'BCC'
    instdict['HadGEM2-ES'] = 'INPE'
    instdict['MPI-ESM-LR'] = 'MPI-M'
    instdict['MPI-ESM-MR'] = 'MPI-M'
    instdict['MPI-ESM-P'] = 'MPI-M'
    instdict['MRI-AGCM3-2H'] = 'MRI'
    instdict['MRI-CGCM3'] = 'MRI'
    instdict['MRI-ESM1'] = 'MRI'
    instdict['MRI-AGCM3-2S'] = 'MRI'


def cmip5_mip2realm_freq(mip): ## CHECK-ME: Same as above
    """
    Returns realm and frequency given the mip in CMIP5
    """

    mipdict = {} ## CHECK-ME: I wrote this based on the DKRZ data, some entries might be missing
    mipdict['Amon'] = ['atmos', 'mon']
    mipdict['Omon'] = ['ocean', 'mon']
    mipdict['Lmon'] = ['land', 'mon']
    mipdict['LImon'] = ['landIce', 'mon']
    mipdict['OImon'] = ['seaIce', 'mon']
    mipdict['aero'] = ['aerosol', 'mon']
#    mipdict['3hr'] = ???
    mipdict['cfDay'] = ['atmos', 'day'] 
    mipdict['cfMon'] = ['atmos', 'mon'] 
    mipdict['day'] = ['atmos', 'day'] 
    mipdict['fx'] = ['*', 'fx']

    if mip in mipdict.keys():
        return mipdict[mip]
    else:
        print "ERROR"
    

def get_input_filelist(project_info, model, var): ## FIX-ME


    # project_info is supposed to contain all the information from the config file
    # Here in particular rootpath and drs are needed (also as dictionaries)

    ## The rootpath is defined in config but now is project-dependent!
    ## Define as dictionary (is this possible in config.ini? Otherwise switch to yaml...)
    ## rootpath['CMIP5'] = /.../.../
    ## rootpath['OBS'] = /.../.../
    ## rootpath['default'] = /../.../ ## default is used to set all non-defined rootpaths
    if model['project'] in rootpath.keys():
        root = rootpath[model['project']]
    else:
        if 'default' in rootpath.keys():
            root = rootpath['default']
        else:
            print "ERROR: rootpath for project " + model['project'] + " not defined in config"

    ## Directory structure is defined in config and also project-depentend (at present only used for CMIP5)
    if model['project'] in drs.keys():
        drs = drs[model['project']]
    else:
        if 'default' in drs.keys():
            drs = drs['default']
        else:
            drs = None

    ## Trying to implement the variable-dependent model keys (see section 4.1 yaml document)
    ## Basic idea: same model used by different variable with different mip/ensemble/exp. In
    ##             this case the mip/exp/ensemble is defined in the variable dictionary, while
    ##             the model dictionary contains a wildcard (mip: *)
    if 'mip' in var.keys():
        model['mip'] = var['mip']
    if 'ensemble' in var.keys():
        model['ensemble'] = var['ensemble']
    if 'exp' in var.keys():
        model['exp'] = var['exp']

    ## TO BE ADDED!!
    ## Possibility to concatenate the results from two different experiments which cover a long time range
    ## Example: 1980-2030, with 1980-2005 from exp=historical and 2006-2030 from exp=rcp26 (future)
    ## Need to think how to best realize this (maybe mip as list? mip: [historical, rcp26]). To be discussed.
    
    ## not sure this is the right syntax to call the function depending on the class name
    ## I guess something like getattr is needed
    files = model['project'].infile_path(root, model, var, drs) ## This return just a string with the full path, including wildcards
    
    ## Missing steps:
    ## --> apply ls / find to convert the string into an actual list of files (like in find_local_files / find_file)
    ## --> make sure that links are replaced with their targets
    ## --> apply the time selection functions time_handling / date_handling to keep only the relevant files

    return filelist


def get_output_file(project_info, model, var): ## FIX-ME
    
    ## not sure this is the right syntax to call the function depending on the class name
    ## I guess something like getattr is needed
    file = model['project'].outfile_path(project_info['GLOBAL']['preproc_dir'], model, var)

    return file



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


class CMIP5(DataFinder):

    def infile_path(self, rootpath, model, var, drs):

        if (drs == 'BADC'):

            ## host_root = '/badc/cmip5/data/cmip5/output1/' ## must be given in config, not here

            # Latest version is always called 'latest' at BADC
            latest = 'latest'

            dirname = '/'.join([rootpath,              # CHECK-ME
                                cmip5_model2inst(model['name']),
                                model['name'],
                                model['exp'],
                                cmip5_mip2realm_freq(model['mip'])[0],
                                cmip5_mip2realm_freq(model['mip'])[1],
                                model['mip'],
                                model['ensemble'],
                                latest,
                                var['name']])

        elif (drs == 'DKRZ'):

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

            dirname = '/'.join([rootpath,
                                cmip5_model2inst(model['name']),
                                model['name'],
                                model['exp'],
                                cmip5_mip2realm_freq(model['mip'])[0],
                                cmip5_mip2realm_freq(model['mip'])[1],
                                model['mip'],
                                model['ensemble'],
                                latest,
                                var['name']])

        elif (drs == 'ETHZ'):
            dirname = '/'.join([rootpath, 
                                model['exp'], 
                                model['mip'], 
                                var, model['name'], 
                                model['ensemble']]) + '/'

        elif (drs == 'SMHI'):
            dirname = '/'.join([rootpath,
                                model["name"],
                                model["ensemble"],
                                model["experiment"],
                                cmip5_mip2realm_freq(model['mip'])[1]])

        elif drs is None:
            dirname = rootpath
            if (not dirname.enswith('/')):
                dirname = dirname + '/'

        else:
            print "ERROR"
            ## FIX-ME: error message as appropriate

        # The CMIP5 filename is always the same, only the drs changes
        filename = '_'.join([var['name'],
                             model['mip'],
                             model['name'],
                             model['experiment'],
                             model['ensemble']]) + '*.nc'
        path = dirname + filename
        
        return path

    def outfile_path(self, rootpath, model, var):
        path =  rootpath + '_'.join([model['project'],
                                     model['name'],
                                     model['mip'],
                                     model['exp'],
                                     model['ensemble'],
                                     var['field'],
                                     var['name'],
                                     str(model['start_year']) + '-' + str(model['end_year'])]) + '.nc'


class EMAC(DataFinder):

    def infile_path(self, rootpath, model, var, drs):
        path = '/'.join(rootpath, model['name'])
        return path

    def outfile_path(self, rootpath, model, var):
        path = rootpath + '_'.join([model['project'],
                                    model['name'],
                                    model['ensemble'],
                                    var['field'],
                                    var['name'],
                                    str(model['start_year']) + '-' + str(model['end_year'])]) + '.nc'
        return path


class OBS(DataFinder):

    def infile_path(self, rootpath, model, var, drs):
        dirname = '/'.join([rootpath, 
                            'Tier' + model['tier'],
                            model['name']])
        filename = '_'.([model['project'], 
                         model['name'],
                         model['type'], 
                         model['version'],
                         var['field'],
                         var['name'],
                         '??????-??????']) + '.nc'
        path = dirname + filename
        return path

    def outfile_path(self, rootpath, model, var):
        path = roothpath + '_'.join([model['project'],
                                     model['name'],
                                     model['type'],
                                     model['version'],
                                     var['field'],
                                     var['name'],
                                     str(model['start_year']) + '-' + str(model['end_year'])]) + '.nc'
        return path


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
        if year1 >= int(year2_model):
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

# ---- function that establishes the DRS path for CMIP5 files
def get_cmip5_drs(dir1, sdir, ic, model, var, latest_dir):
    """
    Function that returns DRS.
    dir1: root directory
    sdir: subdir (institution eg MPI-M)
    ic: project name - MPI-ESM-LR
    model -> mip, exp; var: model attrs and variable name
    latest_dir: on badc is /latest/ - this is known in advance
    and is dependant on where the code is run.
    """
    # these are mendatory to establish path


# ---- capture ls in the preferred directory
def lsladir(dirname):
    """
    Calling this function once so we save time; called in root dirname.
    It is needed for generalization and not hardcoding the institutions.
    """
    # capture the ls output
    lsd = 'ls -la ' + dirname
    proc = subprocess.Popen(lsd, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    res = out.split('\n')[3:-1]
    return res

# ---- local file finder
def find_local_files(model, out1, dirname1, latest_dir, var, label):
    """
    Function that performs local search for files using `find'
    The depth is as high as possible so that find is fast.
    latest_dir: latest version directory e.g. /latest/ on badc
    (see above for details)
    """
    flist = []

    for st in out1:
        subdir = st.split()[-1]
        lsd2 = 'ls -la ' + dirname1 + subdir
        proc2 = subprocess.Popen(lsd2, stdout=subprocess.PIPE, shell=True)
        (out2, err2) = proc2.communicate()

        # work only with existing dirs or allowed permission dirs
        if len(out2) > 0:
            for st2 in out2.split('\n')[3:-1]:
                findic = st2.split()[-1]
                if findic == model['name']:
                    if label == 'cmip5':
                        drs = get_cmip5_drs(dirname1, subdir, findic, model, var, latest_dir)

                    # -follow option allows for finding symlinked files
                    strfindic = 'find ' + drs\
                                 +' -follow -type f -iname "*.nc"'
                    proc = subprocess.Popen(strfindic, stdout=subprocess.PIPE, shell=True)
                    (out, err) = proc.communicate()
                    for t in out.split('\n')[0:-1]:
                        flist.append(t)
    return flist

def find_file(model, ldir, rdir, ld, var, label):
    """
    Function that does direct parsing of available datasource files and establishes
    the paths to the needed files; makes use of find_local_files()
    File versioning is controlled by finding the ld = e.g. /latest/ dir 
    in the badc datasource, this may differ on other clusters and should be correctly
    hardcoded in the code when calling the function. For case 'user', this is a
    'jump' directory (it doesn't exist)

    ldir,rdir: args for find_local_files
    rdir: root data directory eg /badc/cmip5/data/cmip5/output1/
          or rootdir in the 'user' case
    ldir: 'ls -la' output of rdir

    label: type of file needed to establish the drs path

    """

    arname = find_local_files(model, ldir, rdir, ld, var, label)
    fs = []

    if len(arname) > 0:
        yr1 = int(model['start_year'])
        yr2 = int(model['end_year'])
        for s in arname:
            ssp = s.split('/')
            av = ssp[-1]
            # FIXME replace this with time_check()
            time_range = av.split('_')[-1].strip('.nc')
            time1 = time_range.split('-')[0]
            time2 = time_range.split('-')[1]
            year1 = date_handling(time1,time2)[0]
            year2 = date_handling(time1,time2)[1]
            if time_handling(year1, yr1, year2, yr2) is True:
                if os.path.exists(s):
                    fs.append(s)

    return fs

def get_single_file(rootdir, yr1, yr2):
    """
    Retrieves a single file given that
    rootdir = model['path'] is a pointer to a single file
    Performs checks on file extension and needed time bracket
    rootdir: model['path'' - an .nc file
    yr1, yr2: model['start_year'], model['end_year']
    """
    # checks
    if os.path.exists(rootdir):
        if rootdir.endswith('.nc') is True:
            print("PY  info:  >>> get_file_from_drs.py >>> data_dir_type set to user_file ")
            print("PY  info:   >>> get_file_from_drs.py >>> Specified path points to file %s" % rootdir)
        else:
            print("PY  info:   >>> get_file_from_drs.py >>> data_dir_type set to user_file ")
            print("PY  info:   >>> get_file_from_drs.py >>> Specified path DOES NOT point to netCDF file %s" % rootdir)
            sys.exit(1)
    else:
        print("PY  info:   >>> get_file_from_drs.py >>> data_dir_type set to user_file ")
        print("PY  info:   >>> get_file_from_drs.py >>> Specified path is non existent %s" % rootdir)
        sys.exit(1)

    # time checks
    tc = time_check(rootdir, yr1, yr2)
    if tc is True:
        files = [rootdir]
    else:
        print("PY  info:   >>> get_file_from_drs.py >>> Specified file failed time checks: %s" % rootdir)

    return files

def get_from_unstructured_dir(rootdir, model, var):
    """
    Looks for matching files in a directory
    containing a bunch of files, no structure given
    """
    files = []

    proj_name = model['project']

    # use standard convention for file naming according to project
    # CMIP5 and any of its derrivatives
    if proj_name.startswith('CMIP5') is True:

        infile_id = '*' + '_'.join([var, model['mip'],
                                    model['name'],
                                    model['exp'],
                                    model['ensemble']]) + '_*.nc*'

    # EMAC and any of its derrivatives (FIXME UNTESTED!!!)
    elif proj_name.startswith('EMAC') is True:
        infile_id = '*' + '_'.join([var,
                                    model['name'],
                                    model['ensemble']]) + '_*.nc*'

    # GFDL and any of its derrivatives (FIXME UNTESTED!!!)
    elif proj_name.startswith('GFDL') is True:
        infile_id = '*' + '_'.join([var,
                                    model['name'], model['realm'],
                                    model['ensemble']]) + '_*.nc*'

    # CCMVal and any of its derrivatives (FIXME UNTESTED!!!)
    elif proj_name.startswith('CCMVal') is True:
        infile_id = '*' + '_'.join([var,
                                    model['name'], model['exp'],
                                    model['ensemble']]) + '_*.nc*'

    # look for files: get paths
    srch = 'ls ' + rootdir + '/' + infile_id
    proc = subprocess.Popen(srch, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    fpaths = out.split('\n')
    

    # loop through filepaths and identify per variable name;
    # last element of ls will always be an empty string
    for fpath in fpaths[0:-1]:
        filepath = fpath.strip()
        if os.path.exists(filepath):
            # time checks
            tc = time_check(filepath, model['start_year'], model['end_year'])
            if tc is True:
                files.append(filepath)
                print("PY  info:   >>> get_file_from_drs.py >>> Using file %s" % filepath)
        else:
            fi = rootdir + '/' + infile_id
            print("PY  info:   >>> get_file_from_drs.py >>> Could not find file type %s" % fi)

    return files

def get_obs(rootdir, obs_var, obs_model):
    """@brief Function that returns the observation file for regridding
       Returns a files list
       namelist field type: {name: ERA-Interim,  project: OBS,  type: reanaly,
       version: 1,  start: 2000,  end: 2002,  path: /obspath/Tier3/ERA-Interim/}
       path type needs to be: FULL e.g. data_dir/OBS/Tier3/ERA-Interim/ where files are
    """

    obsfiles = []

    # check if rootdir is a full path to a file
    if rootdir.endswith('.nc') is True:
        print("PY  info:   >>> get_file_from_drs.py >>> get_file_from_drs.py >>> Specified path points to file " + rootdir)
        # assert true file existence
        srch = 'ls ' + rootdir

    else:

        # use standard convention for OBS file naming
        # FIXME - there may be multiple OBS naming conventions? 
        infile_id = '*' + '_'.join([obs_model['project'], obs_model['name'],
                                    obs_model['type'], str(obs_model['version']),
                                    obs_var['field'], obs_var['name']]) + '_*.nc*'

        # look for files: get paths
        srch = 'ls ' + rootdir + '/' + infile_id

    proc = subprocess.Popen(srch, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    fpaths = out.split('\n')

    # loop through filepaths
    # last element of ls will always be an empty string
    for fpath in fpaths[0:-1]:
        if os.path.exists(fpath.strip()):
            obsfiles.append(fpath.strip())
            print("PY  info:   >>> get_file_from_drs.py >>> get_file_from_drs.py >>> Using OBS file for regridding " + fpath.strip())
        else:
            fi = rootdir + '/' + infile_id
            print("PY  info:   >>> get_file_from_drs.py >>> get_file_from_drs.py >>> Could not find OBS file type " + fi)

    return obsfiles

def get_cmip5(drs, rootdir, var, model):
    """
    drs: input from config (badc, dkrz, ethz, default, None)
    rootdir: input from roject_info[MODEL] via preprocess.py
    var and model: input from project_info
    """

    if drs == 'badc':
        """
        This case together with 'dkrz' allows for an automatic search on 
        local databases at BADC and DKRZ. The search is performed using pre-built
        DRS paths and uses the key directory name latestDir specific to each database
        for file version control (hardcoded here).
        """
        print("PY  info:  >>> get_file_from_drs.py >>> Looking for data on BADC ")
        host_root = '/badc/cmip5/data/cmip5/output1/' #rdir
        ls_host_root = lsladir(host_root) #ldir

        # this is a standard for BADC
        # hence the hardcoding
        latestDir = '/latest/'

        files = find_file(model, ls_host_root, host_root, latestDir, var, label = 'cmip5')
        print("PY  info:  >>> get_file_from_drs.py >>> Found matching files on BADC: ")
        for fi in files:
            print("PY  info:  >>>    file: %s" % fi)

    elif drs == 'dkrz':

        # FIXME this bit needs correct paths
        print("PY  info:  >>> get_file_from_drs.py >>> Looking for data on DKRZ ")
        host_root = 'xxx/xxx' #rdir
        ls_host_root = lsladir(host_root) #ldir

        # this is a standard for DKRZ
        # hence the hardcoding
        latestDir = '/xxx/'

        files = find_file(model, ls_host_root, host_root, latestDir, var, label = 'cmip5')
        print("PY  info:  >>> get_file_from_drs.py >>> Found matching files on DKRZ: ")
        for fi in files:
            print("PY  info:  >>>    file: %s" % fi)

    elif drs == 'ethz':

        print("PY  info:  >>> get_file_from_drs.py >>> User root directory path with ETHZ structure: %s" % rootdir)
        host_root = rootdir + '/ETHZ_CMIP5/'
        ls_host_root = lsladir(host_root) #ldir
       
        # Not needed for ETHZ
        latestDir = '/'

        files = find_file(model, ls_host_root, host_root, latestDir, var, label=None)

        print(files)
        sys.exit(1)

    elif drs == 'default':
        """
        Default structure: DRS enabled
        This case assumes that the model['path'] = rootdir supplied in the namelist
        has a DRS structure e.g. for CMIP5: rootdir/institution/proj_name/experiment/frequency/realm/mip/ensemble/variable_name/
                                            ./test_data_drs/MPI-M/MPI-ESM-LR/historical/mon/atmos/Amon/r1i1p1/ta/
        """
        print("PY  info:  >>> get_file_from_drs.py >>> User root directory path with DRS structure: %s" % rootdir)
        ls_host_root = lsladir(rootdir) #ldir

        # Not needed for default
        latestDir = '/'

        # check path integrity
        if rootdir.endswith('/'):
            host_root = rootdir
        else:
            host_root = rootdir + '/'

        files = find_file(model, ls_host_root, host_root, latestDir, var, label = 'cmip5')
        if len(files) > 0:
            print("PY  info:  >>> get_file_from_drs.py >>> Found matching files in %s: " % rootdir)
            for fi in files:
                print("PY  info:  >>>    file: %s" % fi)
        else:
            print("PY  info:  >>> get_file_from_drs.py >>> No matching files found in %s " % rootdir)
            sys.exit(1)

    elif drs == 'None':
        """
        This case assumes files live in a single directory
        with no structure whatsoever (all jumbled in)
        """
        files = get_from_unstructured_dir(rootdir, model, var)
        if len(files) > 0:
            print("PY  info:  >>> get_file_from_drs.py >>> Found matching files in %s: " % rootdir)
            for fi in files:
                print("PY  info:  >>>    file: %s" % fi)
        else:
            print("PY  info:  >>> get_file_from_drs.py >>> No matching files found in %s " % rootdir)
            sys.exit(1)

    else:
        print(" >>> get_file_from_drs.py >>> Could not establish root directory type: ", drs)    

    return files


def get_emac(drs, rootdir, var, model):
    # FIXME needs implementation
    pass
def get_gfdl(drs, rootdir, var, model):
    # FIXME needs implementation
    pass
def get_ccmval(drs, rootdir, var, model):
    # FIXME needs implementation
    pass
