#!/usr/bin/env python
"""
###############################################################################
## REFORMAT SCRIPT FOR THE ESACCI-LANDCOVER SATELLITE DATA
###############################################################################
##
## Tier
##    Tier 2:  other freely-available dataset.
##
## Source
##    ESA Climate Change Initiative, Land Cover
##    http://maps.elie.ucl.ac.be/CCI/viewer/download.php
##
## Download and processing instructions
##    1) Goto http://maps.elie.ucl.ac.be/CCI/viewer/index.php and click on
##       "Download data".
##    2) Select "Climate Research Data Package (CRDP)".
##    3) Register with your name and email.
##    4) Seach for "Land Cover Maps - v1.6.1" and download netCDF:
##        a) 2008-2012 epoch - v1.6.1 (netcdf) - 2.8Go
##        b) 2003-2007 epoch - v1.6.1 (netcdf) - 2.8Go
##        c) 1998-2002 epoch - v1.6.1 (netcdf) - 1.1Go
##    5) The CCI-LC User Tools are available on the same page, search for
##       "User tool".
##    6) Unpack the CCI-LC Tools into the directory specified below
##       (variable path2lctool).
##    7) Please adjust in- and out-path below!
##       Then run: python reformat_obs_ESACCI-LANDCOVER.py
##
## Caveats
##    Requires the CCI-LC User Tools
##    Path to folder must be set in variable path2lctool.
##    The CCI-LC User Tools require an installed Java SE 64Bit JRE version 7
##    or higher.
##    The CCI-LC User Tools are available at:
##    http://maps.elie.ucl.ac.be/CCI/viewer/download.php
##    Produces only single year data expanded to their time range of validity!
##
## Modification history
##    20160714-A_muel_bn: written.
##
###############################################################################
"""
from __future__ import print_function

import datetime
import os
import subprocess
import sys
import tempfile

import numpy as np
from cdo import Cdo

sys.path.append(os.path.join(os.path.dirname(__file__), 'lib', 'python'))  # noqa
from preprocessing_basics import _get_files_in_directory

inpath = "/Work/Reference/OBS_ESACCI_LC/"
outpath = "/Work/Reference/OBS_ESACCI_LC/output/"
path2lctool = "./lib/lc-user-tools-3.10/"

# these LC-classes will be processed
var = ["baresoilFrac", "grassNcropFrac", "shrubNtreeFrac"]

# this is how the LC-classes will be aggregated from LC-User-Tool output
translatorlist = {
    'baresoilFrac': ['Bare_Soil'],
    'grassNcropFrac': ['Managed_Grass', 'Natural_Grass'],
    'shrubNtreeFrac': [
        'Tree_Broadleaf_Evergreen', 'Tree_Needleleaf_Evergreen',
        'Tree_Broadleaf_Deciduous', 'Tree_Needleleaf_Deciduous',
        'Shrub_Broadleaf_Evergreen', 'Shrub_Needleleaf_Evergreen',
        'Shrub_Broadleaf_Deciduous', 'Shrub_Needleleaf_Deciduous'
    ]
}

# the available years. You mayextend this list for further data!
year = [2000, 2005, 2010]

field = "T2Ms"

timestep = "monthly"

pathname = os.path.dirname(sys.argv[0])
os.chdir(os.path.abspath(pathname))


def main():

    timestamp = datetime.datetime.now()

    file_list, list_length = _get_files_in_directory(
        inpath, '*P5Y-*-' + str(year) + '*.nc', False)

    if list_length == 0:
        file_list, list_length = _get_files_in_directory(
            inpath, '*P5Y-' + str(year) + '*.nc', False)
        for locfile in file_list:
            lctool_command = [
                "bash",
                os.path.join(path2lctool, "bin", "aggregate-map.sh"),
                "-PgridName=GEOGRAPHIC_LAT_LON",
                "-PnumMajorityClasses=1",
                "-PoutputAccuracy=false",
                "-PoutputPFTClasses=true",
                "-PoutputLCCSClasses=false",
                locfile,
            ]
            process = subprocess.Popen(lctool_command, stdout=subprocess.PIPE)
            for line in iter(process.stdout.readline, b''):
                print(line, end=' ')
                process.stdout.close()
            process.wait()

    for y in year:

        outfilename = dict([(v, "OBS_" + "ESACCI-LANDCOVER_" +
                             "sat_L4-LCCS-Map-300m-P5Y-aggregated-0.083333Deg_"
                             + field + "_" + v) for v in var])

        #        subprocess.call(["mkdir",outpath + "temp/"])

        fullfilenames = dict((de, outpath + outfilename[de])
                             for de in outfilename.keys())

        for ffn in fullfilenames:

            _preprocess_observations(fullfilenames[ffn], outpath,
                                     translatorlist, timestep, ffn, y, inpath,
                                     False)
        # subprocess.call(["rm -rf",outpath + "temp"])

    print("Time for running preprocessing was: " +
          str(datetime.datetime.now() - timestamp))


###############################################################################


def _preprocess_observations(mainfile,
                             outpath,
                             translist,
                             timestep,
                             var,
                             year,
                             check_folder=None,
                             force=False):
    """
    Preprocess observations to adapt to temporal and spatial resolution needed.

    Parameters:
    -----------
    mainfile : string with path to file
    outpath : string with path to outpath
    translist : translatorlist from land cover classes to new class
    timestep : aggregation timestep
    var : name of the variable within the file
    year : year for dataset
    check_folder : alternativeley check this folder and write mainfile.built.nc
    force : if preprocessing should be forced when data already available
    """

    start_year = year - 2
    stop_year = year + 2

    if check_folder is not None:

        file_list, list_length = _get_files_in_directory(
            check_folder, '*P5Y-*-' + str(year) + '*.nc', False)

        assert list_length == 1, ("There is less or more than one file! "
                                  "This should not happen!")

    ofile = mainfile + "_" + str(start_year) + "01" + "-" + str(
        stop_year) + "12" + '.nc'

    mf_bool = os.path.isfile(mainfile)
    of_bool = os.path.isfile(ofile)

    if (mf_bool or of_bool) and not force:
        print("files already exist!")
        return

    elif ((mf_bool or of_bool) and force) or check_folder is not None:

        tmpfile1 = os.path.join(
            outpath,
            "temp",
            os.path.basename(tempfile.NamedTemporaryFile().name), )
        tmpfile2 = os.path.join(
            outpath,
            "temp",
            os.path.basename(tempfile.NamedTemporaryFile().name), )

        cdo = Cdo()

        # select data depending on translatorlist
        files = _select_given_names(
            outpath, file_list[0], translist[var], remove=False)

        # sum data
        cdo.enssum(input=files, output=tmpfile1, options='-f nc4 -b F32')
        [os.remove(f) for f in files]

        # chain: change name, multiply by 100 for "%", setunit to "%",
        # duplicate for length of time range, set time axis, set day
        # to 15, set reference time, calender and time units

        cdo.settunits(
            "seconds",
            input=' '.join([
                '-setcalendar,standard',
                '-setreftime,1970-01-01,00:00:00',
                '-settaxis,{!s}-01-15,12:00:00,1month'.format(start_year),
                '-duplicate,{!s}'.format((stop_year - start_year + 1) * 12),
                '-mulc,100',
                '-setunit,%',
                '-setmissval,1e20',
                '-setctomiss,0',
                '-setname,{!s}'.format(var),
                tmpfile1,
            ]),
            output=tmpfile2,
            options='-f nc4 -b F32')

        subprocess.call(['cp', tmpfile2, ofile])

        os.remove(tmpfile1)
        os.remove(tmpfile2)

    else:
        print(mainfile)
        assert False, "cannot find any files!"


def _select_given_names(work_dir, infile, translist, remove=True):
    """ select names by list """
    cdo = Cdo()

    tmplist = []

    if not (os.path.exists(work_dir + os.sep + "temp")):
        os.makedirs(work_dir + os.sep + "temp")

    for element in translist:
        tmpfile = os.path.join(
            work_dir,
            "temp",
            os.path.basename(tempfile.NamedTemporaryFile().name), )
        # cdo.selname(
        #     element, input=infile, output=tmpfile, options='-f nc4 -b F32')
        cdo.setmisstoc(
            0,
            input="-selvar," + element + " " + infile,
            output=tmpfile,
            options='-f nc4 -b F32')

        tmplist.append(tmpfile)

    if remove:
        os.remove(infile)

    return tmplist


if __name__ == "__main__":
    main()
