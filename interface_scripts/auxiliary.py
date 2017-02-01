import pdb
import sys
import commands
import string


class nclExecuteError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class xmlTagError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class nclExecuteWarning(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class writeProjinfoError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def info(string, verbosity, required_verbosity):
    """ @brief Print an info string to standard out
        @param string the info message to print
        @param verbosity the requested verbosity level
        @param required_verbosity level required to print something
    """
    if verbosity >= required_verbosity:
        print "PY  info: " + str(string)


def error(string):
    """ @brief Print an info string to standard error and exit execution
        @param string the info message to print
    """
    sys.stderr.write("error: " + string + '\n')
    sys.exit(1)

def print_header(projdict, opt):
    """ @brief Print the ESMValTool header
        @param project_info dictionary with the necessary information
        @param opt logical for shorter header for the reformat case
    """
    
    vv = 1
    line1 = 54 * "_"
    line2 = 61 * "_"

    info("", vv, 1)
    info(line1, vv, 1)
    info("  _____ ____  __  ____     __    _ _____           _  ", vv, 1)
    info(" | ____/ ___||  \/  \ \   / /_ _| |_   _|__   ___ | | ", vv, 1)
    info(" |  _| \___ \| |\/| |\ \ / / _` | | | |/ _ \ / _ \| | ", vv, 1)
    info(" | |___ ___) | |  | | \ V / (_| | | | | (_) | (_) | | ", vv, 1)
    info(" |_____|____/|_|  |_|  \_/ \__,_|_| |_|\___/ \___/|_| ", vv, 1)
    info(line1, vv, 1)
    info("", vv, 1)
    info(" http://www.esmvaltool.org/", vv, 1)
    info(line2, vv, 1)
    info("", vv, 1)
    info("CORE DEVELOPMENT TEAM AND CONTACTS:", vv, 1)
    info("  Veronika Eyring (PI; DLR, Germany - veronika.eyring@dlr.de)", vv, 1)
    info("  Bjoern Broetz (DLR, Germany - bjoern.broetz@dlr.de)", vv, 1)
    info("  Axel Lauer (DLR, Germany - axel.lauer@dlr.de)", vv, 1)
    info("  Mattia Righi (DLR, Germany - mattia.righi@dlr.de)", vv, 1)
    info(line2, vv, 1)
    info("", vv, 1)
    if not opt:
        info("NAMELIST = " + projdict['RUNTIME']['xml_name'], vv, 1)
        info("WORKDIR  = " + projdict["GLOBAL"]["wrk_dir"], vv, 1)
        info("CLIMODIR = " + projdict["GLOBAL"]["climo_dir"], vv, 1)
        info("PLOTDIR  = " + projdict["GLOBAL"]["plot_dir"], vv, 1)
        info("LOGFILE  = " + projdict['RUNTIME']['out_refs'], vv, 1)
        info(line2, vv, 1)
    else:
        info("REFORMATTING THE OBSERVATIONAL DATA...", vv, 1)
    info("", vv, 1)

def ncl_version_check():
    """ @brief Check the NCL version
    """
    out = commands.getstatusoutput("ncl -V")

    if out[0] != 0:
        error("NCL not found")

    if out[1] == "6.3.0":
        error("NCL version " + out[1] + 
              " not supported due to a bug " + 
              "(see Known Issues in the ESMValTool user guide)")

    if int(out[1].split(".")[0]) < 6:
        error("NCL version " + out[1] + " not supported, need version 6.2.0 or higher")

    if int(out[1].split(".")[0]) == 6 and int(out[1].split(".")[1]) < 2:
        error("NCL version " + out[1] + " not supported, need version 6.2.0 or higher")
