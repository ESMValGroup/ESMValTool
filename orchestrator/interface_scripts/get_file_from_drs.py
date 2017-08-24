"""
Tool to establish a structured DRS path to files
Author: Valeriu Predoi, University of Reading, valeriu.predoi@ncas.ac.uk
First version: August 2017
"""

# ---- Import standard modules to the python path.
import sys, os
import subprocess
from datetime import datetime

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

# ---- function that establishes the DRS path
def get_drs(dir1, sdir, ic, model, var, latest_dir):
    """
    Function that returns DRS.
    dir1: root directory
    sdir: subdir (institution eg MPI-M)
    ic: project name - MPI-ESM-LR
    model -> mip, exp; var: model attrs and variable name
    latest_dir: on badc is /latest/ - this is known in advance
    and is dependant on where the code is run.
    """
    mip = model['mip']
    exp = model['exp']
    ens = model['ensemble']

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
        print('Could not establish custom DRS...')
        gdrs = dir1 + sdir + '/' + ic + '/' + exp + '/' + mip + '/*/*/' + ens \
               + latest_dir + var + '/'
        print('Using generalized path: %s' % gdrs)
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
def find_local_files(model, out1, dirname1, latest_dir, var):
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
                    drs = get_drs(dirname1, subdir, findic, model, var, latest_dir)

                    # -follow option allows for finding symlinked files
                    strfindic = 'find ' + drs\
                                 +' -follow -type f -iname "*.nc"'
                    proc = subprocess.Popen(strfindic, stdout=subprocess.PIPE, shell=True)
                    (out, err) = proc.communicate()
                    for t in out.split('\n')[0:-1]:
                        flist.append(t)
    return flist

def find_file(model, ldir, rdir, ld, var):
    """
    Function that does direct parsing of available datasource files and establishes
    the paths to the needed files; makes use of find_local_files()
    File versioning is controlled by finding the ld = e.g. /latest/ dir 
    in the badc datasource, this may differ on other clusters and should be correctly
    hardcoded in the code!

    ldir,rdir: args for find_local_files
    rdir: root data directory eg /badc/cmip5/data/cmip5/output1/
    ldir: 'ls -la' output of rdir

    """

    arname = find_local_files(model, ldir, rdir, ld, var)
    fs = []

    if len(arname) > 0:
        yr1 = int(model['start_year'])
        yr2 = int(model['end_year'])
        for s in arname:
            ssp = s.split('/')
            av = ssp[-1]
            time_range = av.split('_')[-1].strip('.nc')
            time1 = time_range.split('-')[0]
            time2 = time_range.split('-')[1]
            year1 = date_handling(time1,time2)[0]
            year2 = date_handling(time1,time2)[1]
            if time_handling(year1, yr1, year2, yr2) is True:
                if os.path.exists(s):
                    print("PY  info:  >>> get_file_from_drs.py >>> Found matching file: ", s)
                    fs.append(s)

    return fs


def get_cmip5(drs, rootdir, var, model):
    """
    drs: input from config
    rootdir: input from roject_info[MODEL] via preprocess.py
    var and model: input from project_info
    """

    if drs == 'badc':
        print("PY  info:  >>> get_file_from_drs.py >>> Looking for data on BADC ")
        host_root = '/badc/cmip5/data/cmip5/output1/' #rdir
        ls_host_root = lsladir(host_root) #ldir

        # this is a standard for badc
        latestDir = '/latest/'

        files = find_file(model, ls_host_root, host_root, latestDir, var)
        print("PY  info:  >>> get_file_from_drs.py >>> Found matching files on BADC: ", files)

    elif drs == 'dkrz':
        print(" >>> get_file_from_drs.py >>> Looking for data on DKRZ ")

    elif drs == 'user':
        print(" >>> get_file_from_drs.py >>> User root directory path with DRS structure: %s" % rootdir)
        ls_host_root = lsladir(rootdir) #ldir

        # this is a standard for badc
        latestDir = '/latest/'

        # check path integrity
        if rootdir.endswith('/'):
            host_root = rootdir
        else:
            host_root = rootdir + '/'

        files = find_file(model, ls_host_root, host_root, latestDir, var)
        print("PY  info:  >>> get_file_from_drs.py >>> Found matching files in user data root: ", files)
    else:
        print(" >>> get_file_from_drs.py >>> Could not establish root directory path: ", drs)    

    return files


def get_emac(drs, rootdir, var, model):
    pass
def get_gfdl(drs, rootdir, var, model):
    pass
def get_ccmval(drs, rootdir, var, model):
    pass
