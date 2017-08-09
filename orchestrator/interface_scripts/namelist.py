import yaml

class NameList(Object):
    """"
    Class representing a namelist containing a number of diagnostics. Each diagnostic will result in a call of some
     diagnostic script. A namelist yaml definition with multiple script lines per diagnostic will therefor lead
     to multiple diagnostics (one per script called) in this namelist object.
    """

    @staticmethod
    def _parse_namelist_file(file):
        with open(file, 'r') as stream:
            namelist_data = yaml.load(stream)


        global_section = namelist_data['GLOBAL']

        logging.debug("loaded global section", global_section)

        global_model_section = namelist_data['MODELS']

        logging.debug("loaded global model section", global_section)

        global_settings_section = namelist_data['OPTIONS']

        logging.debug("loaded global settings section", global_settings_section)


    def __init__(self, file):
        self._diagnostics = []

class Diagnostic(Object):

    def __init__(self, namelist_data, ):
        pass


class Variable(Object):
    pass


class Model(Object):
    pass


