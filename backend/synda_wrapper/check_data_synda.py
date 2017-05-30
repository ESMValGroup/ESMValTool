#!/home/valeriu/sdt/bin/python

"""

Script that performs a check on the current/most recent downloads via synda
"""

# -------------------------------------------------------------------------
#      Setup.
# -------------------------------------------------------------------------

# ---- Import standard modules to the python path.
import sys, os, shutil, math, random, copy, getopt, re, string, popen2, time
import numpy as np
from numpy import loadtxt as lt
from numpy import savetxt as st
from xml.dom import minidom
import subprocess
from datetime import datetime
import time

__author__ = "Valeriu Predoi <valeriu.predoi@ncas.ac.uk>"

# ---- Function usage.
# ---- opts parsing
def usage():
  msg = """\
This is a tool to perform a real-time check on the data download via synda.
To use it first have a look and use get_data_synda.py that will start your data download.
For queries, email valeriu.predoi@ncas.ac.uk. Have fun!

Usage:
  check_data_synda.py [options]
  -h, --help                  Display this message and exit
"""
  print >> sys.stderr, msg

########################################
# ---- Operational functions here ---- #
########################################

# ---- get the path to synda executable
def which_synda(synda):
    """

    This function returns the path to the synda exec
    or aborts the whole program if path is not found.

    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(synda)
    if fpath:
        if is_exe(synda):
            #print('We are using the following executable: %s' % synda)
            return synda
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, synda)
            if is_exe(exe_file):
                #print('We are using the following executable: %s' % exe_file)
                return exe_file

    #print('No synda executable found on the system. Aborting program.')
    return None

# ---- synda check download
def synda_check_dll():
    print('Checking if/how your files(s) are being downloaded.')
    print('You can check the download progress yourself with synda queue, see output below')
    synda_queue = which_synda('synda') + ' queue'
    proc = subprocess.Popen(synda_queue, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    print(out)
    statusreport = out.split('\n')
    for entry in statusreport:
        if len(entry)>0:
            if entry.split()[0] == 'waiting':
                print('%i files are waiting, totalling %.2f MB disk' % (int(entry.split()[1]),float(entry.split()[2])))
    synda_watch = which_synda('synda') + ' watch'
    proc = subprocess.Popen(synda_watch, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    print(out)
    return out.strip()

# ---- check the logfile
def synda_check_log():
    synda = which_synda('synda')
    logpath = "/".join(synda.split('/')[0:-2]) + '/log/transfer.log'
    if os.path.isfile(logpath):
        logfile = open(logpath, 'r')
        for line in logfile:
            if line.split()[5] == 'done':
                print(line)
        return logpath
    else:
        print('Can not find synda transfers log file...')
        return None

# -------------------------------------------------------------------------
#      Parse the command line options.
# -------------------------------------------------------------------------
# ---- Initialise command line argument variables.
params_file       = None

# ---- Syntax of options, as required by getopt command.
# ---- Short form.
shortop = "hp:g:r:d:i:n:t:f:m:sc:e:"
# ---- Long form.
longop = [
   "help",
   "params-file="
]

# ---- Get command-line arguments.
try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

# ---- Parse command-line arguments.  Arguments are returned as strings, so
#      convert type as necessary.

command_string = 'check_data_synda.py '
for o, a in opts:
    if o in ("-h", "--help"):
        usage()
        sys.exit(0)
    elif o in ("-p", "--params-file"):
        params_file = a
        command_string = command_string + ' -p ' + a
    else:
        print >> sys.stderr, "Unknown option:", o
        usage()
        sys.exit(1)

# -------------------------------------------------------------------------
#      Status message.  Report all supplied arguments.
# -------------------------------------------------------------------------

print >> sys.stdout
print >> sys.stdout, "####################################################"
print >> sys.stdout, "# Checking Data Search and Download using synda    #"
print >> sys.stdout, "####################################################"
print >> sys.stdout

#--------------------------------------------------------------------------
#     Main code body
#--------------------------------------------------------------------------

# ---- Get the synda path or exit here
print('Looking up synda executable...')
if which_synda('synda') is not None:
    print >> sys.stdout, "Synda found...OK"
    print >> sys.stdout, which_synda('synda')
else:
    print >> sys.stderr, "No synda executable found in path. Exiting."
    sys.exit(1)

# ---- Have us some information from the synda configuration file
# ---- one can add more info if needed, currently just data server
print('\n---------------------------------------------')
print('Information about synda configuration:')
print('---------------------------------------------')
synda_conf_file = which_synda('synda').rsplit('/',2)[0] + '/conf/sdt.conf'
print ('Synda conf file %s' % synda_conf_file)
with open(synda_conf_file, 'r') as file:
    for line in file:
        if line.split('=')[0]=='indexes':
            data_server = line.split('=')[1]
            print('Data server: %s' % line.split('=')[1])
            if data_server is None:
                data_server = 'esgf_generic'

# simple check
if synda_check_dll() == 'No current download':
    print >> sys.stderr, "Currently downloading 0 files. Let's perform a check on what we've downloaded"
    print('Log file: %s' % synda_check_log())
else:
    print >> sys.stderr, "Currently downloading files:"
    synda_check_dll()
    print('Log file: %s' % synda_check_log())
    print >> sys.stderr, "Rerun this script to get updates on your progress."
