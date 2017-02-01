import pdb
import sys


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
