class Namelist(yaml.YAMLObject):
    """
    Class to hold the information from the namelist
    """
    yaml_tag = u'!Namelist'
    def __init__(self):
        self.__name__       = "Namelist" 
        self.write_plots    = None
        self.verbosity      = None
        self.exit_on_warning = None
        self.ouput_file_type = None
        self.models         = None # sould be []
        self.preprocess     = None # sould be {}
        self.diagnostics    = None # sould be {}

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



