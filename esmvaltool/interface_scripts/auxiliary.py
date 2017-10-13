import pdb
import sys
import commands
import logging
import warnings

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
    warnings.warn("function info() is deprecated and will be removed "
                  "in the future", DeprecationWarning)
    if verbosity <= 1:
        logger.info("%s", string)
    else:
        logger.debug("%s", string)


def error(string):
    """ @brief Print an info string to standard error and exit execution
        @param string the info message to print
    """
    warnings.warn("function error() is deprecated and will be removed "
                  "in the future", DeprecationWarning)
    logger.error("%s", string)
    sys.exit(1)


def ncl_version_check():
    """ @brief Check the NCL version
    """
    out = commands.getstatusoutput("ncl -V")

    if out[0] != 0:
        logger.error("NCL not found")

    if out[1] == "6.3.0":
        logger.error("NCL version " + out[1] + " not supported due to a bug " +
                     "(see Known Issues in the ESMValTool user guide)")

    if int(out[1].split(".")[0]) < 6:
        logger.error("NCL version " + out[1] +
                     " not supported, need version 6.2.0 or higher")

    if int(out[1].split(".")[0]) == 6 and int(out[1].split(".")[1]) < 2:
        logger.error("NCL version " + out[1] +
                     " not supported, need version 6.2.0 or higher")
