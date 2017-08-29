import pdb
import sys
import commands
import string
import logging

logger = logging.getLogger(__name__)

class nclExecuteError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class ymlTagError(Exception):
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

def print_header():
    """ @brief Print the ESMValTool header
        @param project_info dictionary with the necessary information
    """

    vv = 1
    line1 = 54 * "_"
    line2 = 70 * "_"

    header = [
        "",
        line1,
        "  _____ ____  __  ____     __    _ _____           _  ",
        " | ____/ ___||  \/  \ \   / /_ _| |_   _|__   ___ | | ",
        " |  _| \___ \| |\/| |\ \ / / _` | | | |/ _ \ / _ \| | ",
        " | |___ ___) | |  | | \ V / (_| | | | | (_) | (_) | | ",
        " |_____|____/|_|  |_|  \_/ \__,_|_| |_|\___/ \___/|_| ",
        line1,
        "",
        " http://www.esmvaltool.org/",
        line2,
        "",
        "CORE DEVELOPMENT TEAM AND CONTACTS:",
        "  Veronika Eyring (PI; DLR, Germany - veronika.eyring@dlr.de)",
        "  Bjoern Broetz (DLR, Germany - bjoern.broetz@dlr.de)",
        "  Nikolay Koldunov (AWI, Germany - nikolay.koldunov@awi.de)",
        "  Axel Lauer (DLR, Germany - axel.lauer@dlr.de)",
        "  Benjamin Mueller (LMU, Germany - b.mueller@iggf.geo.uni-muenchen.de)",
        "  Valeriu Predoi (URead, UK - valeriu.predoi@ncas.ac.uk)",
        "  Mattia Righi (DLR, Germany - mattia.righi@dlr.de)",
        "  Javier Vegas-Regidor (BSC, Spain - javier.vegas@bsc.es)",
        line2,
        "",
    ]

    logger.info("\n".join(header))


def ncl_version_check():
    """ @brief Check the NCL version
    """
    out = commands.getstatusoutput("ncl -V")

    if out[0] != 0:
        logger.error("NCL not found")

    if out[1] == "6.3.0":
        logger.error("NCL version " + out[1] +
              " not supported due to a bug " +
              "(see Known Issues in the ESMValTool user guide)")

    if int(out[1].split(".")[0]) < 6:
        logger.error("NCL version " + out[1] + " not supported, need version 6.2.0 or higher")

    if int(out[1].split(".")[0]) == 6 and int(out[1].split(".")[1]) < 2:
        logger.error("NCL version " + out[1] + " not supported, need version 6.2.0 or higher")
