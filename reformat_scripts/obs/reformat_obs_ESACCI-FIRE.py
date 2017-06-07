#!/usr/bin/env python

"""
###############################################################################
## REFORMAT SCRIPT FOR THE ESACCI-FIRE SATELLITE DATA
###############################################################################
##
## Tier
##    Tier 2:  other freely-available dataset.
##
## Source
##    ESA Climate Change Initiative, burned area
##    website
##
## Download and processing instructions
##    1) Goto website and click on "Download data".
##
## Caveats
##
## Modification history
##    20170420-A_muel_bn: written.
##
###############################################################################
"""

#inpath = "/Work/Reference/OBS_ESACCI_FIRE_all/"
inpath = "/media/bmueller/Work/ESMVAL_res/work/Reference/OBS_CCI_FIRE_all/"
#outpath = "/Work/Reference/OBS_ESACCI_FIRE_all/output/"
outpath = "/media/bmueller/Work/ESMVAL_res/work/Reference/OBS_CCI_FIRE_all/output/"

field = "T2Ms"

timestep = "monthly"

basicpath = "./lib/python/"

import sys
import os
import datetime
from cdo import Cdo
import numpy as np

pathname = os.path.dirname(sys.argv[0])
os.chdir(os.path.abspath(pathname))
sys.path.append(basicpath)

from preprocessing_basics import _get_files_in_directory


def main():

    timestamp = datetime.datetime.now()

    try:
        os.mkdir(outpath)

    except:
        print(outpath + " exists!")

    file_list, list_length = \
        _get_files_in_directory(inpath, '*ESACCI*FIRE*fv4.1.nc', False)

    timerange = [int(fname.split("/")[-1].split("-")[0][0:6])
                 for fname in file_list]
    timerange = str(np.min(timerange)) + "-" + str(np.max(timerange))

    outfilename = "OBS_" + "ESACCI-FIRE_" + "sat_L4-BA-MERIS-fv4.1_" + \
        field + "_" + "burntArea" + "_" + timerange + ".nc"

    outfilename = outpath + outfilename

    _preprocess_observations(file_list, outfilename, timestep, True)

    print("Time for running preprocessing was: " +
          str(datetime.datetime.now() - timestamp))


###############################################################################

def _preprocess_observations(infiles, outfile, timestep, force=False):
    """
    preprocess observations to adapt to temporal and spatial resolution needed
    Parameters:
    -----------

    infiles : string with path to files
    outfile : string with path to outfile
    timestep : aggregation timestep
    force : if preprocessing should be forced when data already available
    """

    if not os.path.isfile(outfile) or force:

        if all([os.path.isfile(i) for i in infiles]):

            cdo = Cdo()

            tmp = outpath + "tmp/"

            try:
                os.mkdir(tmp)

            except:
                print(tmp + " exists!")

            for f in infiles:
                f_act = f.split("/")[-1]
                date = f_act.split("-")[0]
                print(date[0:4] + "-" + date[4:6] + "-" + date[6:8])
                cdo.selname("burned_area",
                            input="-setcalendar,standard " +
                            "-setreftime,1970-01-01,00:00:00 -settaxis," +
                            date[0:4] + "-" + date[4:6] + "-" + date[6:8] +
                            ",12:00:00 " + f, output=tmp + f_act,
                            options='-f nc4 -b F32')

            file_list, list_length = \
                _get_files_in_directory(tmp, '*ESACCI*FIRE*fv4.1.nc', False)
            cdo.mergetime(input=file_list, output=tmp + "tmp1",
                          options='-f nc4 -b F32')
            cdo.setctomiss(1.e33, input=tmp + "tmp1", output=tmp + "tmp2")
            cdo.setctomiss(1.e20, input=tmp + "tmp2", output=tmp + "tmp3")
            cdo.setctomiss(1.e38, input=tmp + "tmp3", output=tmp + "tmp4")
            cdo.gridarea(input=tmp + "tmp4", output=tmp + "tmp5",
                         options='-f nc4 -b F32')
            cdo.div(input=[tmp + "tmp4", tmp + "tmp5"], output=tmp + "tmp6",
                    options='-f nc4 -b F32')
            cdo.monsum(input="-chname,burned_area,burntArea " + tmp + "tmp6",
                       output=tmp + "tmp7", options='-f nc4 -b F32')
            cdo.chunit("m2,%", input=tmp + "tmp7", output=tmp + "tmp8",
                       options='-f nc4 -b F32')
            cdo.mulc(100, input=tmp + "tmp8", output=outfile,
                     options='-f nc4 -b F32')

        else:
            for f in infiles:
                print(f)
            assert False, "Cannot find all files!"
    else:
        print("Nothing to do, " + outfile + " already exists!")

    for root, dirs, files in os.walk(tmp, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
        os.rmdir(root)

if __name__ == "__main__":
    main()
