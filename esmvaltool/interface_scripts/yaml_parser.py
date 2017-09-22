import yaml

class Namelist(yaml.YAMLObject):
    """
    Class to hold the information from the namelist
    """
    yaml_tag = u'!Namelist'

    def __repr__(self):
        return '{0}(settings={1}, preprocess={2}, \
               models={3}, diagnostics={4})'.format(
            self.__class__.__name__,
            self.PREPROCESS,
            self.MODELS,
            self.DIAGNOSTICS,
            self.CONFIG)

class Diagnostic(yaml.YAMLObject):
    
    yaml_tag = u'!Diagnostic'

    def __repr__(self):
        return "%s(id=%r, description=%r, variables=%r, scripts=%r, additional_models=%r)" % (
            self.__class__.__name__, self.id, self.description, self.variables, self.scripts, self.additional_models)

############################################################################################################################
# this class is not currently used but is coded here
# in case we decide to implement Preprocess as a yaml object
class PreprocessItem(yaml.YAMLObject):
    yaml_tag = u'!PreprocessItem'

    def __repr__(self):
        return "%s(select_level=%r, regrid=%r, target_grid=%r, regrid_scheme=%r, mask_fillvalues=%r, mask_landocean=%r, multimodel_mean=%r)" % (
            self.__class__.__name__, self.select_level, self.regrid, self.target_grid, self.regrid_scheme, self.mask_fillvalues, self.mask_landocean, self.multimodel_mean)
##############################################################################################################################

class Parser():
    def load_namelist(self, param_file):

        with open(param_file, 'r') as file:
            n = yaml.load(file)

        # some conditioning, add more in the future?
        if not isinstance(n, Namelist): raise Exception( "Namelist malformed" )
        if not isinstance(n.MODELS, list): raise Exception( "MODELS is not a list" )
        if not isinstance(n.DIAGNOSTICS, dict): raise Exception( "DIAGNOSTICS is not a dictionary" )
        return n
