"""
Tool to establish a structured DRS path to files
Retrieval of data files from either structured databases like BADC or DKRZ
or from DRS-structured ROOTDIR specified by user.

The usage of this tool is goverened by data_dir_type option in config.ini with
the following permitted values and their respective actions:

            Option | For all options just model['path'] from namelist is needed
            --------------------------------------------------------------------
            O1: user_drs: assumes user's 'model['path']' path has a DRS structure of
                          form is model['path']/institution/proj_name/experiment/frequency/realm/mip/ensemble/variable_name/FILE
                          once model['path'] is passed, the code will retrieve the file from the DRS structure
            O2: user_file: user's model['path'] points to a single file (good for small diagnostics)
            O3: user_unstructured: user's model['path'] contains an unstructured collection of files (good for testing or for small diagnostics)
            O4: badc: file search on BADC archive (automatic)
            O5: dkrz: file search on DKRZ archive (automatic) - still to be correctly handled

Author: Valeriu Predoi, University of Reading, valeriu.predoi@ncas.ac.uk
First version: August 2017
"""

# ---- Import standard modules to the python path.
import sys, os
import subprocess
from datetime import datetime
import logging

logger = logging.getLogger(__name__)

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
        if year2 <= int(year1_model):
            return False
        # data overlaps to the left
        elif year2 >= int(year1_model):
            return True
    elif year1 >= int(year1_model) and year2 >= int(year2_model):
        # data is entirely after model
        if year1 >= int(year2_model):
            return False
        # data overlaps to the right
        elif year1 <= int(year2_model):
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
    try:
        mip = model['mip']
        exp = model['exp']
        ens = model['ensemble']
    except KeyError as e:
        logger.error("%s is missing in model %s", e, model)
        sys.exit(1)

    if mip == '3h':
        gdrs = dir1 + sdir + '/' + ic + '/' + exp + '/3h/*/*/' + ens \
               + latest_dir + var + '/'
    elif mip == '6h':
        gdrs = dir1 + sdir + '/' + ic + '/' + exp + '/6h/*/*/' + ens \
               + latest_dir + var + '/'
    elif mip == 'day':
        gdrs = dir1 + sdir + '/' + ic + '/' + exp + '/day/*/*/' + ens \
              + latest_dir + var + '/'
    elif mip == 'cfDay':
        gdrs = dir1 + sdir + '/' + ic + '/' + exp + '/cfDay/*/*/' + ens \
              + latest_dir + var + '/'
    elif mip == 'Amon':
        gdrs = dir1 + sdir\
              + '/' + ic + '/' + exp + '/mon/atmos/Amon/' + ens \
              + latest_dir + var + '/'
    elif mip == 'Omon':
        gdrs = dir1 + sdir\
              + '/' + ic + '/' + exp + '/mon/ocean/Omon/' + ens \
              + latest_dir + var + '/'
    elif mip == 'Lmon':
        gdrs = dir1 + sdir\
              + '/' + ic + '/' + exp + '/mon/land/Lmon/' + ens \
              + latest_dir + var + '/'
    elif mip == 'LImon':
        gdrs = dir1 + sdir\
              + '/' + ic + '/' + exp + '/mon/landIce/LImon/' + ens \
              + latest_dir + var + '/'
    elif mip == 'OImon':
        gdrs = dir1 + sdir\
              + '/' + ic + '/' + exp + '/mon/seaIce/OImon/' + ens \
              + latest_dir + var + '/'
    elif mip == 'aero':
        gdrs = dir1 + sdir\
              + '/' + ic + '/' + exp + '/mon/aerosol/aero/' + ens \
              + latest_dir + var + '/'
    else:
        logger.info('Could not establish custom DRS; trying a general path')
        gdrs = dir1 + sdir + '/' + ic + '/' + exp + '/' + mip + '/*/*/' + ens \
               + latest_dir + var + '/'
        logger.info('Using generalized path: %s', gdrs)
    return gdrs

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
            logger.info("data_dir_type set to user_file ")
            logger.info("Specified path points to file %s", rootdir)
        else:
            logger.info("data_dir_type set to user_file ")
            logger.info("Specified path DOES NOT point to netCDF file %s", rootdir)
            sys.exit(1)
    else:
        logger.info("data_dir_type set to user_file ")
        logger.info("Specified path is non existent %s", rootdir)
        sys.exit(1)

    # time checks
    tc = time_check(rootdir, yr1, yr2)
    if tc is True:
        files = [rootdir]
    else:
        logger.info("Specified file failed time checks: %s", rootdir)

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
                logger.info("Using file %s", filepath)
        else:
            fi = rootdir + '/' + infile_id
            logger.info("Could not find file type %s", fi)

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
        logger.info("Specified path points to file %s", rootdir)
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
            logger.info("Using OBS file for regridding %s", fpath.strip())
        else:
            fi = rootdir + '/' + infile_id
            logger.info("Could not find OBS file type %s", fi)

    return obsfiles

def get_cmip5(drs, rootdir, var, model):
    """
    drs: input from config (badc, dkrz, user_drs)
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
        logger.info("Looking for data on BADC ")
        host_root = '/badc/cmip5/data/cmip5/output1/' #rdir
        ls_host_root = lsladir(host_root) #ldir

        # this is a standard for BADC
        # hence the hardcoding
        latestDir = '/latest/'

        files = find_file(model, ls_host_root, host_root, latestDir, var, label = 'cmip5')
        logger.info("Found matching files on BADC: ")
        for fi in files:
            logger.info("file: %s", fi)

    elif drs == 'dkrz':

        # FIXME this bit needs correct paths
        logger.info("Looking for data on DKRZ ")
        host_root = 'xxx/xxx' #rdir
        ls_host_root = lsladir(host_root) #ldir

        # this is a standard for DKRZ
        # hence the hardcoding
        latestDir = '/xxx/'

        files = find_file(model, ls_host_root, host_root, latestDir, var, label = 'cmip5')
        logger.info("Found matching files on DKRZ: ")
        for fi in files:
            logger.info("file: %s", fi)

    elif drs == 'default':
        """
        Default structure: DRS enabled
        This case assumes that the model['path'] = rootdir supplied in the namelist
        has a DRS structure e.g. for CMIP5: rootdir/institution/proj_name/experiment/frequency/realm/mip/ensemble/variable_name/
                                            ./test_data_drs/MPI-M/MPI-ESM-LR/historical/mon/atmos/Amon/r1i1p1/ta/
        """
        logger.info("User root directory path with DRS structure: %s", rootdir)
        ls_host_root = lsladir(rootdir) #ldir

        # this is a dummy
        latestDir = '/'

        # check path integrity
        if rootdir.endswith('/'):
            host_root = rootdir
        else:
            host_root = rootdir + '/'

        files = find_file(model, ls_host_root, host_root, latestDir, var, label = 'cmip5')
        if len(files) > 0:
            logger.info("Found matching files in %s: ", rootdir)
            for fi in files:
                logger.info("file: %s", fi)
        else:
            logger.info("No matching files found in %s ", rootdir)
            sys.exit(1)

    elif drs == 'None':
        """
        This case assumes files live in a single directory
        with no structure whatsoever (all jumbled in)
        """
        files = get_from_unstructured_dir(rootdir, model, var)
        if len(files) > 0:
            logger.info("Found matching files in %s: ", rootdir)
            for fi in files:
                logger.info("file: %s", fi)
        else:
            logger.info("No matching files found in %s ", rootdir)
            sys.exit(1)

    else:
        logger.info("Could not establish root directory type: ", drs)    

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
