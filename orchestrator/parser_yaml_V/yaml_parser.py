import yaml
import os

class Namelist(yaml.YAMLObject):
    """
    Class to hold the information from the namelist
    """
    yaml_tag = u'!Namelist'

    def __init__(self, GLOBAL, PREPROCESS, MODLES, DIAGNOSTICS, CONFIG):
        self.__name__ = "Namelist"
        self.DIAGNOSTICS = []
        self.GLOBAL = GLOBAL
        self.MODELS = MODLES
        self.PREPROCESS = PREPROCESS
        self.DIAGNOSTICS.append(DIAGNOSTIC)
        self.CONFIG = CONFIG

    def __repr__(self):
        return '{0}(settings={1}, preprocess={2}, \
               models={3}, diagnostics={4})'.format(
            self.__class__.__name__,
            self.GLOBAL,
            self.PREPROCESS,
            self.MODELS,
            self.DIAGNOSTIC,
            self.CONFIG)

class Diagnostic(yaml.YAMLObject):
    yaml_tag = u'!Diagnostic'
    def __init__(self, id, description, variable, preprocess, scripts, additional_models):
        self.id = id
        self.description = description
        self.variable = variable
        self.preprocess = preprocess
        self.scripts = scripts
        self.additional_models = additional_models
    def __repr__(self):
        return "%s(id=%r, description=%r, variable=%r, preprocess=%r, scripts=%r, additional_models=%r)" % (
            self.__class__.__name__, self.id, self.description, self.variable, self.preprocess, self.scripts, self.additional_models)

class Parser():
    def load_namelist(self, param_file):
        if not os.path.isfile(param_file):
            print >> sys.stderr,"Error: non existent parameter file: ", \
                param_file
            sys.exit(1)
        s = file(param_file, 'r')
        n = yaml.load(s)
        print('yaml_parser PY: We have loaded %s parameter file...' % param_file)
        assert isinstance(n, Namelist)
        assert n.GLOBAL["write_plots"] == True
        assert n.GLOBAL["write_netcdf"] == True
        assert n.GLOBAL["verbosity"] == 1
        assert n.GLOBAL["exit_on_warning"] == False
        assert n.GLOBAL["output_file_type"] == "ps"
        assert 'select_level' in n.PREPROCESS.keys()
        assert isinstance(n.MODELS, list)
        assert isinstance(n.DIAGNOSTICS, dict)
        print('yaml_parser PY: We have %i diagnostics to do! Starting the main code now...' % len(n.DIAGNOSTICS))
        return n
