import yaml
import os
from Namelist_0 import Namelist

class Parser():
    def load_namelist(self, param_file):
        if not os.path.isfile(param_file):
            print >> sys.stderr,"Error: non existent parameter file: ", \
                param_file
            sys.exit(1)
        s = file(param_file, 'r')
        n = yaml.load(s)
        assert isinstance(n, Namelist)
        assert n.GLOBAL["write_plots"] == True
        assert n.GLOBAL["write_netcdf"] == True
        assert n.GLOBAL["verbosity"] == 1
        assert n.GLOBAL["exit_on_warning"] == False
        assert n.GLOBAL["output_file_type"] == "ps"
        assert 'select_level' in n.PREPROCESS.keys()
        assert isinstance(n.MODELS, list)
        assert isinstance(n.DIAGNOSTICS, dict)
        for k, v in n.DIAGNOSTICS.items():
            assert 'description' in v
            assert 'scripts' in v
        return n
