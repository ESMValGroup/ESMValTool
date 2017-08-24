import yaml
import os
import sys

class Namelist(yaml.YAMLObject):
    """
    Class to hold the information from the namelist
    """
    yaml_tag = u'!Namelist'

    def __init__(self, PREPROCESS, MODLES, DIAGNOSTICS, CONFIG):
        self.__name__ = "Namelist"
        self.DIAGNOSTICS = []
        self.MODELS = MODLES
        self.PREPROCESS = PREPROCESS
        self.DIAGNOSTICS.append(DIAGNOSTIC)
        self.CONFIG = CONFIG

    def __repr__(self):
        return '{0}(settings={1}, preprocess={2}, \
               models={3}, diagnostics={4})'.format(
            self.__class__.__name__,
            self.PREPROCESS,
            self.MODELS,
            self.DIAGNOSTIC,
            self.CONFIG)

class Diagnostic(yaml.YAMLObject):
    yaml_tag = u'!Diagnostic'
    def __init__(self, id, description, variables, scripts, additional_models):
        self.id = id
        self.description = description
        self.variables = variables
        self.scripts = scripts
        self.additional_models = additional_models
    def __repr__(self):
        return "%s(id=%r, description=%r, variables=%r, scripts=%r, additional_models=%r)" % (
            self.__class__.__name__, self.id, self.description, self.variables, self.scripts, self.additional_models)

############################################################################################################################
# this class is not currently used but is coded here
# in case we decide to implement Preprocess as a yaml object
class PreprocessItem(yaml.YAMLObject):
    yaml_tag = u'!PreprocessItem'
    def __init__(self, select_level, regrid, target_grid, regrid_scheme, mask_fillvalues, mask_landocean, multimodel_mean):
        self.select_level = select_level
        self.regrid = regrid
        self.target_grid = target_grid
        self.regrid_scheme = regrid_scheme
        self.mask_fillvalues = mask_fillvalues
        self.mask_landocean = mask_landocean
        self.multimodel_mean = multimodel_mean
    def __repr__(self):
        return "%s(select_level=%r, regrid=%r, target_grid=%r, regrid_scheme=%r, mask_fillvalues=%r, mask_landocean=%r, multimodel_mean=%r)" % (
            self.__class__.__name__, self.select_level, self.regrid, self.target_grid, self.regrid_scheme, self.mask_fillvalues, self.mask_landocean, self.multimodel_mean)
##############################################################################################################################

class Parser():
    def load_namelist(self, param_file):
        if not os.path.isfile(param_file):
            print >> sys.stderr,"Error: non existent parameter file: ", \
                param_file
            sys.exit(1)
        s = file(param_file, 'r')
        n = yaml.load(s)
        print('PY  info: >>> yaml_parser.py >>> We have loaded %s parameter file...' % param_file)
        # some conditioning, add more in the future?
        if not isinstance(n, Namelist): raise Exception( "Namelist malformed" )
        if not isinstance(n.MODELS, list): raise Exception( "MODELS is not a list" )
        if not isinstance(n.DIAGNOSTICS, dict): raise Exception( "DIAGNOSTICS is not a dictionary" )
        print('PY  info: >>> yaml_parser.py >>> We have %i diagnostics to do! Starting the main code now...' % len(n.DIAGNOSTICS))
        return n
