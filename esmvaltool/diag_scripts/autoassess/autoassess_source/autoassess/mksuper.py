#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
(C) Crown Copyright 2017, the Met Office

Create multi-year supermeans from monthly/seasonal/annual data
"""

from __future__ import division, print_function
import argparse
import os
import os.path
import pdb
import sys
import tempfile

from autoassess import argparser

# Create script argument parser
POS_ARGS = [
    ('jobid', dict(help='Job ID for the experiment')),
]
OPT_ARGS = [
    ('nyears', dict(default=10, type=int, help='Number of years (default=10)')),
    ('startyear', dict(default=1982, type=int, help='Start year (default=1982)')),
    ('query', dict(dest='query_file', default=None, help='Specify a query file to only retrieve required fields. This query file should contain the line yr = [STARTYR..ENDYR]')),
    ('noget', dict(default=False, action='store_true', help='By default mksuper will look for existing supermeans in mass and retrieve them instead of making new ones. This can be turned off by including -g at the command line.')),
    ('noarchive', dict(default=False, action='store_true', help='Archiving is ON by default and users who do NOT want to archive supermeans after calculating them should use this option (useful when MASS is busy and causes mksuper to hang).')),
    ('nostdev', dict(default=False, action='store_true', help='Calculating standard deviations is ON by default and users who do NOT want to calculate standard deviations should use this option (useful if you get var_to_stdev.pro errors).')),
    ('vn', dict(default=None, help='The UM version number of your job. Normally this does not matter (and can be omitted) with the exception of UM version 7.9 where --vn=7.9 needs to be set to include _00 at the end of the annual mean file names.')),
    ('oneseason', dict(default=0, help='By default supermeans will be created for all four seasons and using annual means. Use this option to specify only one season. Seasons are: 1=DJF, 2=MAM, 3=JJA, 4=SON, 5=annual. ')),
    ('work-dir', dict(default=os.path.join(os.environ['SCRATCH'], 'makesuper'), help='Directory to which files are restored and means are calculated')),
    ('mean-dir', dict(default=None, help='Directory to which supermean files are written, if not specified writes to current working directory.')),
]
ARG_PARSE = argparser(__doc__, POS_ARGS, OPT_ARGS)

class namelistfile(file):
    '''
    Class for building a namelist file:\n
    Open the file with nmf = namelist.namelistfile(filename).
    Add namelists with nmf.namelist().
    Close it like a normal file: nmf.close().
    '''

    def __init__(self, filename=None):
        if not filename:
            # get a unique filename
            (fd, filename) = tempfile.mkstemp(suffix='.nml')
        file.__init__(self, filename, 'w')

    def namelist(self, name, varvals):
        '''
        Add a namelist to the file.\n
        Varvals is a sequence of (var,value) doubles.
        '''
        namlst = ['&'+name]
        namlst.extend([' %s = %s,' % var_value for var_value in varvals])
        namlst.append('/')
        # Blank line in case of future namelist
        namlst.append('')

        namlst = [line+'\n' for line in namlst]
        self.writelines(namlst)


def tidyup(filearray):
    # Remove the pp files
    for afile in filearray:
        print('Deleting file {0}'.format(afile))
        os.remove(afile)


def ppmnvar(infiles, meanfile, varfile=None, mltann='.true.', days=None,
            run_ppmnvar='/opt/ukmo/utils/bin/run_ppmnvar',
            overwrite=False):

    # Open a file to store your inputs to run_ppmnvar
    (fd, file_list) = tempfile.mkstemp(prefix='ppmnvar_', suffix='.list')
    f = open(file_list, 'w')

    # Populate that file
    lastfile = infiles[-1]
    for file in infiles:
        f.write(file)
        if file == lastfile:
            f.write(' {0} {1}'.format(meanfile, varfile))
        f.write('\n')
    f.close()

    # Call run_ppmnvar
    ret = os.system('{0} -F {1} -V'.format(run_ppmnvar, file_list))
    if ret != 0:
        print('Warning: run_ppmnvar returned code {0}\n'.format(ret))

    return(ret)


def get_yearcode(nyears):
    if nyears <= 9:
        stream = chr(nyears+48)
    elif nyears <= 10:
        stream = 'a'
    elif nyears <= 19:
        stream = 'b'
    elif nyears <= 24:
        stream = 'k'
    elif nyears <= 29:
        stream = 's'
    elif nyears <= 39:
        stream = 't'
    elif nyears <= 49:
        stream = 'q'
    elif nyears <= 99:
        stream = 'l'
    elif nyears <= 249:
        stream = 'u'
    elif nyears <= 499:
        stream = 'w'
    else:
        stream = 'd'
    return stream


def get_seasonstr(season, name_system, source=False):
    if source and season == 'ann':
        if name_system == 1:
            seasonstr = 'c10'
        elif name_system == 2:
            seasonstr = '1201'
        else:
            print('Naming system not defined')
    else:
        seasonstr = season
    return seasonstr


def get_yearstr(year, name_system):
    if name_system == 1:
        (decade, yearindec) = divmod(year, 10)
        cyclic_decade = decade % 36
        if cyclic_decade > 9:
            decadestr = chr(cyclic_decade-10+97)
        else:
            decadestr = chr(cyclic_decade+48)
        yearindecstr = chr(yearindec+48)
        yearstr = decadestr + yearindecstr
    elif name_system == 2:
        yearstr = str(year)
    else:
        print('Naming system not defined')
    return yearstr


def get_seasonstream(season):
    if season == 'ann':
        stream = 'y'
    else:
        stream = 's'
    return stream


def get_moodest(jobid, nyears):
    # Get yearcode for number of years
    yearcode = get_yearcode(nyears)
    return 'moose:/crum/{0}/am{1}.pp'.format(jobid, yearcode)


def get_sourceurl(jobid, season):
    stream = get_seasonstream(season)
    return 'moose:/crum/{0}/ap{1}.pp'.format(jobid, stream)


def def_meanfile(jobid, startyear, nyears, season, name_system, filestr):
    """Generate a mean file name"""

    # startyear = the start year
    # nyears = the number of years
    # season = the season
    # name_system = the naming system to use
    # filestr: m = mean, v = variance, s = standard deviation

    # Get string that represents season in mean file
    seasonstr = get_seasonstr(season, name_system)

    # Get yearcode for number of years
    yearcode = get_yearcode(nyears)

    # Get string for year dependent on naming convention
    yearstr = get_yearstr(startyear+nyears-1, name_system)

    # Define a shorter (5 letter) version of the jobid using the last 5 letters
    jobid5 = jobid[-5:]

    filename = '{0}a.{1}{2}{3}{4}.pp'.format(jobid5, filestr,
                                             yearcode, yearstr,
                                             seasonstr)
    return filename


def def_sourcefile(jobid, startyear, dyear, season, name_system, vn):
    """Generate a source file name"""

    # Get stream that source file resides in
    stream = get_seasonstream(season)

    # Get string that represents season in source file
    seasonstr = get_seasonstr(season, name_system, source=True)
    if stream == 'y' and vn == 7.9:
        seasonstr += '_00'

    # Get string for year dependent on naming convention
    yearstr = get_yearstr(startyear+dyear, name_system)

    # Define a shorter (5 letter) version of the jobid using the last 5 letters
    jobid5 = jobid[-5:]

    filename = '{0}a.p{1}{2}{3}.pp'.format(jobid5, stream, yearstr, seasonstr)

    return filename


def get_namesystem(jobid):

    # Check to see what file name convention you are using
    sourceurl = get_sourceurl(jobid, 'ann')
    command = 'moo ls {0} > inarchive.list'.format(sourceurl)
    returncode = os.system(command)
    if returncode != 0:
        print('No annual means found using command: '+command)
        print('Return code = {0}'.format(returncode))
        sys.exit(20)

    # What files are in the apy stream
    inarchive = open("inarchive.list").read()
    # What is the first file in the apy stream
    firstfile = inarchive.split()[0]
    # What are the four digits that represent the year in the new system
    yearfour = firstfile[-11:-7]
    # Delete the file now that it is no longer needed
    os.remove("inarchive.list")

    # Assign a value to name_system:
    # 1 = old pre vn7.9 file naming
    # 2 = new vn7.9 file naming

    # See if we can get an integer out of yearfour
    try:
        yearfourint = int(yearfour)
    except ValueError:
        name_system = 1
    else:
        name_system = 2
    return name_system


def make_query(query_file, startyear, nyears, season):
    formatted_query = tempfile.mktemp(suffix='.query')
    print('Output query file is {0}'.format(formatted_query))

    month = dict(djf=12, mam=3, jja=6, son=9, ann=12)
    offset = dict(djf=-1, mam=0, jja=0, son=0, ann=-1)

    endyear = startyear + nyears - 1
    new_startyr = str(startyear+offset[season])
    new_endyr = str(endyear+offset[season])
    new_month = str(month[season])

    with open(query_file, 'r') as query_file_f:
        with open(formatted_query, 'w') as formatted_query_f:

            for line in query_file_f:
                corr_start = line.replace('STARTYR', new_startyr)
                corr_end = corr_start.replace('ENDYR', new_endyr)
                corr_month = corr_end.replace('MONTH', new_month)
                formatted_query_f.write(corr_month)

    return formatted_query


# TODO: Write this in Iris
def make_stdev(workdir, varfile, stdevfile):
    idlfile = workdir+'/'+stdevfile+'_idl.pro'
    f = open(idlfile, 'w')
    f.write('status=1\n')
    # This should be using the maverick version, not hadco version.
    f.write('.compile /home/h03/hadco/unix/supermeans/var_to_stdev.pro\n')
    f.write('!PP_LARGEFILE_EXTEND=0\n')
    f.write('var_to_stdev, \''+varfile+'\', \''+stdevfile+'\', status=status\n')
    f.write('exit, status=status'+'\n')
    f.close()
    os.environ['IDL_STARTUP'] = ' '
    os.environ['IDL_PP_PATH'] = '.'
    tempdir = os.environ['LOCALTEMP']
    if not os.path.exists(tempdir):
        os.mkdir(tempdir)
    returncode = os.system('tidl -quiet '+idlfile+' 2>'+tempdir+'/mksuper_idl.log')
    if returncode != 0:
        print('Error converting variances into standard deviations:')
        print('From variance file '+varfile)
        print('To standard deviation file '+stdevfile)
        print('Using temporary idl code '+idlfile)
        print('Check '+tempdir+'/mksuper_idl.log')
        print(returncode)
        sys.exit(9)
    os.remove(idlfile)


def archive_supermean(meanfile, moodest):
    print('Archiving supermean {0} to {1}'.format(meanfile, moodest))
    command = 'moo put {0} {1}'.format(meanfile, moodest)
    returncode = os.system(command)
    if returncode != 0:
        print('Error archiving file {0} to MASS-R'.format(meanfile))
        print('Return code = {0}'.format(returncode))
    return (returncode != 0)


def mksuper(jobid, nyears, startyear, workdir, oneseason=0, query_file=None,
            vn=None, noget=False, nostdev=False, noarchive=False):
    """The main part of the program"""

    # Check that QUERYSPEC environment variable set for MOOSE filtering
    if query_file:
        if not os.path.isfile(query_file):
            print('Query file {0} does not exist'.format(query_file))
            sys.exit(1)
        moo_method = 'moo select'
    else:
        moo_method = 'moo get'

    # Define seasons to loop over
    seasonarr = ['djf', 'mam', 'jja', 'son', 'ann']
    if oneseason:
        seasonarr = seasonarr[oneseason-1]

    # Flag for archive success
    failedarchive = False

    # Check to see what file name convention you are using
    name_system = get_namesystem(jobid)

    # Loop over seasons
    for season in seasonarr:

        sourceurl = get_sourceurl(jobid, season)
        moodest = get_moodest(jobid, nyears)
        meanfile = def_meanfile(jobid, startyear, nyears, season, name_system, 'm')
        varfile = def_meanfile(jobid, startyear, nyears, season, name_system, 'v')
        stdevfile = def_meanfile(jobid, startyear, nyears, season, name_system, 's')

        if not nostdev:
            gotstdev = os.path.exists(stdevfile)
        else:
            gotstdev = True

        # If files already exist then try next season
        if os.path.exists(meanfile) and gotstdev:
            tmpmeanfile = os.path.join(os.getcwd(), meanfile)
            print('{0} already exists. Skipping...'.format(tmpmeanfile))
            continue

        # If the supermean already exists in mass then retrieve it
        if not noget:
            moourl = os.path.join(moodest, meanfile)
            command = 'moo get {0} . 2>/dev/null'.format(moourl)
            returncode = os.system(command)
            if returncode == 0 and not nostdev:
                moourl = os.path.join(moodest, stdevfile)
                command = 'moo get {0} . 2>/dev/null'.format(moourl)
                returncode = os.system(command)
            if returncode == 0:
                print('Supermean {0} already exists in MASS-R and has been retrieved.'.format(meanfile))
                continue

        # These arrays will hold the list of files and urls to mean
        filearray = []
        urlarray = []

        # Make sure the temporary directory exists
        if not os.path.exists(workdir):
            os.mkdir(workdir)

        for dyear in range(0, nyears):

            # Generate the filename
            filename = def_sourcefile(jobid, startyear, dyear, season, name_system, vn)

            # Generate full paths
            tempfilename = os.path.join(workdir, filename)
            tempfileurl = os.path.join(sourceurl, filename)

            # Store the filename in an array
            filearray.append(tempfilename)

            # Look to see if the file has already been downloaded.
            # If not append it to the file url array
            if not os.path.exists(tempfilename):
                urlarray.append(tempfileurl)

        # Only download files if needed
        if len(urlarray) > 0:

            # Print what we are doing
            print('Using {0} to download the following files to {1}'.format(moo_method, workdir))
            print('\n'.join(urlarray))

            # Download the files from MASS-R

            # Use a query get if requiring filtered data
            if query_file is not None:

                # Make a query file
                formatted_query = make_query(query_file, startyear, nyears, season)

                # Get the data
                command = 'moo select {0} {1} {2}'.format(formatted_query, sourceurl, workdir)
                returncode = os.system(command)
                if returncode != 0:
                    print('Error downloading files from {0} using query {1}'.format(sourceurl, formatted_query))
                    print('Return code = {0}'.format(returncode))
                    pdb.set_trace()
                    sys.exit(4)
                else:
                    if os.access(formatted_query, os.F_OK):
                        os.remove(formatted_query)

            else:

                # Use the standard moo get to get the data
                file_str = ' '.join(urlarray)
                command = 'moo get -i {0} {1}'.format(file_str, workdir)
                returncode = os.system(command)
                if returncode != 0:
                    print('Error downloading files from MASS-R.')
                    print('Please check that all the required files exist.')
                    print('Return code = {0}'.format(returncode))
                    sys.exit(5)

        else:
            print('Files already downloaded for making supermeans. Skipping moose download.')

        # Make sure there are no ftn files in your output directory
        if os.path.exists('ftn08'):
            print()
            print('Output directory = {0}'.format(os.getcwd()))
            print('There are ftn* files in your output directory.')
            print('The existance of these files means that either someone')
            print('else is already making supermeans in your output directory')
            print('or you had a previous failed attempt to make supermeans')
            print('here. Fot the latter case please go to your output')
            print('directory and delete all files beginning with ftn before')
            print('re-running this supermean script')
            sys.exit(6)

        # Make a supermean using ppmnvar
        print('Making supermean               {0}'.format(meanfile))
        returncode = ppmnvar(filearray, meanfile, varfile=varfile,
                             mltann='.true.', overwrite=True)
        if returncode != 0:
            print('Error meaning files:')
            print('\n'.join(filearray))
            print(returncode)
            sys.exit(7)

        # Archive the supermean to new MASS
        if not noarchive:
            failedarchive |= archive_supermean(meanfile, moodest)

        # Convert variances to standard deviations for 16222, 5216, 3236
        if not nostdev:
            print('Making standard deviation file '+stdevfile)
            make_stdev(workdir, varfile, stdevfile)
            # Archive the standard deviation file to MASS
            if not noarchive:
                failedarchive |= archive_supermean(stdevfile, moodest)

        # Delete unwanted annual or seasonal means
        print('Deleting temporary files.')
        for myfile in filearray:
            os.remove(myfile)
        os.remove(varfile)

    print('Finished making or downloading supermeans')

    if failedarchive:
        print()
        print('WARNING: mksuper has failed to archive one or more of your supermeans. Please archive them manually to '+moodest)


def main():
    '''Creating multi-year supermeans'''

    # Parse arguments
    args = ARG_PARSE.parse_args(sys.argv[1:])

    # If a means directory is specified then change current working directory
    if args.mean_dir:
        if not os.path.exists(args.mean_dir):
            os.makedirs(args.mean_dir)
        os.chdir(args.mean_dir)

    # Run mksuper
    mksuper(args.jobid, args.nyears, args.startyear, args.work_dir,
            oneseason=args.oneseason, query_file=args.query_file,
            vn=args.vn, noget=args.noget, nostdev=args.nostdev,
            noarchive=args.noarchive)


if __name__ == '__main__':
    main()
