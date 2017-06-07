import yaml


class Namelist(yaml.YAMLObject):
    """
    Class to hold the information from the namelist
    """
    yaml_tag = u'!Namelist'

    def __init__(self, GLOBAL, PREPROCESS, MODLES, DIAGNOSTICS):
        self.__name__ = "Namelist"
        self.GLOBAL = GLOBAL
        self.MODELS = MODLES
        self.PREPROCESS = PREPROCESS
        self.DIAGNOSTICS = DIAGNOSTICS

    def __repr__(self):
        return '{0}(settings={1}, preprocess={2}, \
               models={3}, diagnostics={4})'.format(
            self.__class__.__name__,
            self.GLOBAL,
            self.PREPROCESS,
            self.MODELS,
            self.DIAGNOSTICS)


class Config(object):
    """
    Class to hold the information from the configfile 
    """

    def __init__(self):
        self.__name__ = "Config"


class Diagnostic(object):
    """
    """

    def __init__(self):
        self.__name__ = "Diagnostic"


class Variable(object):
    """
    """

    def __init__(self):
        self.__name__ = "Variable"
