import logging
import subprocess
import sys
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

    try:
        version = subprocess.check_output(['ncl', '-V'])
    except subprocess.CalledProcessError:
        logger.error("NCL not found")

    version = version.decode(sys.stdout.encoding)

    if version == "6.3.0":
        logger.error("NCL version " + version + " not supported due to a bug "
                     + "(see Known Issues in the ESMValTool user guide)")

    if int(version.split(".")[0]) < 6:
        logger.error("NCL version " + version +
                     " not supported, need version 6.2.0 or higher")

    if int(version.split(".")[0]) == 6 and int(version.split(".")[1]) < 2:
        logger.error("NCL version " + version +
                     " not supported, need version 6.2.0 or higher")
