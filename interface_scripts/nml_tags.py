from auxiliary import xmlTagError, error, info
from model import AllModels, Model
from fx_file import AllFXfiles, FX_file, FX_file_exception
from esgf_config import ESGFConfig, ESGFConfigException
from diagdef import AllDiagnostics, Diagnostic
import diagdef
import os
import pdb


class Nml_base:
    def __init__(self):
        pass

    def closing_tag(self, str, attributes):
        return self.project_info

    def add_nml_entry(self, name, str, attributes):
        pass


class GLOBAL(Nml_base):
    def __init__(self):
        Nml_base.__init__(self)
        self.project_info = {}

    def add_nml_entry(self, name, str, attributes):
        project_info = str.strip()
        if attributes:
            if attributes.values()[0] == "boolean":
                if project_info == "True":
                    project_info = True
                elif project_info == "False":
                    project_info = False
                else:
                    raise TypeError("Invalid value for boolean")

            elif attributes.values()[0] == "integer":
                project_info = int(project_info)

            elif attributes.values()[0] == "path":
		if project_info[0] not in "@":
	                project_info = os.path.abspath(project_info)
            else:
                raise TypeError("Invalid value for attributes")

        self.project_info[name] = project_info


class MODELS(Nml_base):
    def __init__(self):
        Nml_base.__init__(self)
        self.project_info = AllModels()

    def add_nml_entry(self, name, str, attributes):
        self.project_info.append(Model(str.strip(),
                                       attributes,
                                       diag_specific_model=False))


class AUXILIARIES(Nml_base):
    def __init__(self):
        Nml_base.__init__(self)
        self.project_info = {"FX_files": AllFXfiles()}
        ## Other types of entry could be added to AUXILIARIES

    def add_nml_entry(self, name, str, attribute):
        if name == "fx_file":
            ## Here assume FX_file ID is stored in an
            ## attribute called 'id' (note parser converts all
            ## attribute names to lower case),
            ## with FX_file fullpath stored in 'str'
            if 'id' in attribute:
                fx = FX_file(attribute['id'], str.strip())
                self.project_info["FX_files"].append(fx)
            else:
                msg = "fx_file entry '" + str.strip() +\
                      "' has no 'id' attribute."
                raise FX_file_exception(msg)


class ESGF(Nml_base):
    """
    This section is for ESGF config information.
    All information is currently held entirely in the ESGF config file,
    so this section may be unneccesary. I'm keeping it, for now, in case
    namelist specific config information is required.
    """
    def __init__(self):
       Nml_base.__init__(self)

       self.project_info = {}
       # self.project_info becomes project_info['ESGF'] in the main program.

    def add_nml_entry(self, name, string, attribute):
        """
        Process namelist entry
        :param name: tag name of XML element
        :param string: content of XML element ('str' in the other classes)
        :param attribute: any attribute(s) given inside start tag
        """
        # Get path to ESGF config file
        if name == 'config_file':
            self.project_info['config_file'] = string.strip()


class DIAGNOSTICS(Nml_base):
    def __init__(self):
        Nml_base.__init__(self)
        self.project_info = AllDiagnostics()
        self.diag_tags = []
        self.reset_temp_diag_storage()

    def reset_temp_diag_storage(self):
        self.variable_def_dir = None
        self.launcher_arguments = None
        self.diag_specific_models = []
        self.diag_script_cfg = []
        self.diag_script = []
        self.field_type = []
        self.variable = []
        self.var_attributes = []
        self.cfg = []

        ## Default value
        self.diag_script_cfg_dir = []

    def add_nml_entry(self, name, string, attributes):
        if name == "diag_script":
            self.diag_script.append(string.strip())
            if attributes:
                self.diag_script_cfg.append(attributes)
            else:
                self.diag_script_cfg.append("default.cfg")

        elif name in ['variable_def_dir']:
            vars(self)[name] = string.strip()

        elif name in ['launcher_arguments']:
            vars(self)[name] = string.strip()

        elif name == 'model':
            self.diag_specific_models.append(Model(string.strip(),
                                                   attributes,
                                                   diag_specific_model=True))
        elif name == 'diag':
            ## These two should arrays of the same length as the
            ## number of variables
            self.extend_array_to_match_length("field_type", len(self.variable))
            self.extend_array_to_match_length("diag_script_cfg_dir",
                                              len(self.diag_script))

            self.diag_tags.append(diagdef.Diag_tag(",".join(self.variable),
                                                   self.variable_def_dir,
                                                   ",".join(self.field_type),
                                                   self.var_attributes,
                                                   self.diag_script,
                                                   self.diag_script_cfg_dir,
                                                   self.diag_script_cfg,
                                                   self.diag_specific_models,
                                                   self.launcher_arguments))
            self.reset_temp_diag_storage()

        elif name == 'description':
            pass

        elif name == 'variable':
            self.variable.append(string.strip())
            self.var_attributes.append(attributes)
        else:
            vars(self)[name].append(string.strip())

    def extend_array_to_match_length(self, extend_var, other_array_length):
        ext_var_len = len(vars(self)[extend_var])

        if ext_var_len == 1:
            vars(self)[extend_var] = vars(self)[extend_var] * other_array_length

        elif ext_var_len == 0:
            vars(self)[extend_var] = [""] * other_array_length

        ext_var_len = len(vars(self)[extend_var])

        if ext_var_len != other_array_length:
            raise xmlTagError("Number of variables vs fields/cfg-files do not match: "
                              + str(ext_var_len)
                              + " != "
                              + str(other_array_length))

    def closing_tag(self, str, attributes):
        all_diags = []
        for curr_diag_tag in self.diag_tags:
            var = curr_diag_tag.get_tag_variable()
            var_attr = curr_diag_tag.get_tag_attr()
            field = curr_diag_tag.get_tag_field_type()
            var_def_dir = curr_diag_tag.get_var_def_dir()
            launch_args = curr_diag_tag.get_launcher_args()

            diags = [Diagnostic(var, var_def_dir, field, var_attr, diag_script, cfg, model, launch_args)
                     for diag_script, cfg, model in curr_diag_tag]

            all_diags.extend(diags)

        self.project_info.extend(all_diags)
        return self.project_info


class namelist_summary(Nml_base):
    def __init__(self):
        self.project_info = ""

    def closing_tag(self, str, attributes):
        return str


class namelist:
    def __init__(self):
        pass

class REFORMAT(Nml_base):
    def __init__(self):
        Nml_base.__init__(self)
        self.project_info = {}

    def add_nml_entry(self, name, string, attributes):
	if attributes['id'] in self.project_info.keys():
		error('Duplicate usage of reformat_script id: {0}'.format(attributes['id']))
        self.project_info[attributes['id']] = string.strip()
