#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
'''
Python script to make source_file.dat for use in the validation note software
'''

import argparse
import os
import os.path


def write_file_location(jobid, mean_root, alphayear, endyears, season, f,
                        average_type, suffix='pp'):
    '''
    Function to write out the location of files in source_file.dat
    jobid = the UMUI jobid or Rose suiteid
    mean_root = the root of the directory where all the supermeans are stored
    alphayear = the character that defines the number of years in your supermean
    endyears = string that represents the final year in the supermean
    season = 3 letter string of the season
    f = file object that you are writing to
    average_type = either 'm' (mean) or 's' (standard deviation)
    suffix = file suffix (default pp)
    '''

    # Define a shorter (5 letter) version of the jobid using the 5 letters
    #  after the -
    if jobid.count('-') >= 1:
        jobid_right = jobid.split('-')[1]
        jobid5 = jobid_right[0:5]
    else:
        jobid5 = jobid[-5:]

    # Write out the main supermean file
    filename = '{0}a.{1}{2}{3}{4}.{5}'.format(jobid5, average_type,
                                              alphayear, endyears,
                                              season, suffix)
    filename = os.path.join(mean_root, jobid, 'supermeans', filename)
    if os.path.exists(filename):
        f.write('   '+season+average_type+':'+filename+'\n')


def get_yearcode(nyear):
    if nyear <= 9:
        cyear = str(nyear)
    elif nyear <= 10:
        cyear = 'a'
    elif nyear <= 19:
        cyear = 'b'
    elif nyear <= 24:
        cyear = 'k'
    elif nyear <= 29:
        cyear = 's'
    elif nyear <= 39:
        cyear = 't'
    elif nyear <= 49:
        cyear = 'q'
    elif nyear <= 99:
        cyear = 'l'
    elif nyear <= 249:
        cyear = 'u'
    elif nyear <= 499:
        cyear = 'w'
    else:
        cyear = 'd'
    return cyear


def run(exper, control, nyears, startyear, exper_name, control_name,
        mean_root, dir_to_write):
    '''Subroutine that generates the source_file'''

    # Set up the character that defines the number of years in your supermean
    alphayear = get_yearcode(nyears)

    # Set the end year
    endyears = str(startyear + nyears - 1)

    # What seasons are we using
    season_array = ['djf', 'mam', 'jja', 'son', 'ann']
    average_type_array = ['m', 's']

    # Make sure write directory exists
    if not os.path.isdir(dir_to_write):
        os.makedirs(dir_to_write)

    # Open a new source_file.dat for writing
    with open(os.path.join(dir_to_write, 'source_file.dat'), 'w') as f:

        # Write out a header for information
        f.write('#\n')
        f.write('# source_file.dat: This file is used by the validation ' +
                'note program to point to your control\n')
        f.write('# and experiment supermeans.\n')
        f.write('#\n')

        # Loop over the experiment and control
        for (model_type, jobid, jobid_name) in (('exper',
                                                 exper,
                                                 exper_name),
                                                ('control',
                                                 control,
                                                 control_name)):
            f.write('start:'+model_type+'\n')
            f.write('   jobid:'+jobid+'\n')
            for average_type in average_type_array:
                for season in season_array:
                    write_file_location(jobid, mean_root, alphayear, endyears,
                                        season, f, average_type)
            f.write('   name:'+jobid_name+'\n')
            f.write('   startyear:'+str(startyear)+'\n')
            f.write('   nyears:'+str(nyears)+'\n')
            f.write('end:'+model_type+'\n')
            f.write('\n')

    print 'Finished writing out source_file.dat'


def parse_args():
    '''Routine to parse the command line arguments'''

    # parse command line arguments
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Add optional arguments to parser
    parser.add_argument('-e', '--exper-name', required=True, help='A short name to label the experiment in the plots')
    parser.add_argument('-c', '--control-name', required=True, help='A short name to label the control in the plots')
    parser.add_argument('-m', '--mean-root', required=True, help='The root location for the supermeans')
    parser.add_argument('-n', '--nyears', default=20, type=int, help='Number of years (default=20)')
    parser.add_argument('-s', '--startyear', default=1989, type=int, help='Start year (default=1989)')

    # Add positional arguments to the parser
    parser.add_argument('exper', help='Job ID for the experiment')
    parser.add_argument('control', help='Job ID for the control')

    # Returned parsed arguments
    return parser.parse_args()


# **************************************************
# Code below is only called if this file is run as a
#  script from the command line
# **************************************************

def main():

    args = parse_args()

    # Run the main part of the program
    # Make the source file in your current directory
    dir_to_write = os.getcwd()
    run(args.exper, args.control, args.nyears, args.startyear,
        args.exper_name, args.control_name, args.mean_root, dir_to_write)

if __name__ == "__main__":
    main()
