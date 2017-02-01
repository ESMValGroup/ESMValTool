# Diagnostic attribute  Reader

# 2008/06/27 - C. Fischer  taken from namelist.py and changed
# Dependencies:

# Known issues:

# What it does:
#   Reads in a diagnostic attribe files
#
# Usage

#
# Define the classes
#

from auxiliary import info
import re
import copy
import os
import reformat
import operator
import glob
import exceptions
import pdb


class Var:
    """ @brief Class variable the basic (non-derived) variables
    """
    def __init__(self, merged_dict, var0, fld0):
        # Write the variable attributes in 'merged_dict'
        # to class object attributes
        for name, value in merged_dict.iteritems():
            setattr(self, name, value)

        # Special cases, actually not sure what they do
        if (var0 == "none"):
            self.var0 = merged_dict['var']
        else:
            self.var0 = var0
        if (fld0 == "none"):
            self.fld0 = merged_dict['fld']
        else:
            self.fld0 = fld0


class VAR_REQ:
    """ @brief Class holding dependent variables/fields

        For derived variables, this class will holds all
        their dependencies (all required variables/fields)
    """
    def __init__(self, vars, fields):
            self.var = vars
            self.field = fields

    def __iter__(self):
        for i in range(len(self.var)):
            os.environ['__ESMValTool_base_var'] = self.var[i]
            var_mip_exp = [item[i] for item in self.var_attrs.values()]
            key_mip_exp = [re.sub("var_attr_", "", item)
                           for item in self.var_attrs.keys()]
            mip_exp = dict(zip(key_mip_exp, var_mip_exp))
            yield self.var[i], self.field[i], mip_exp


class Diag_tag:
    """ @brief Class to hold an instance for a single Diagnostic xml-file
    """
    def __init__(self, variable, var_def_dir, field, var_attr, diag_scripts,
                 cfg_dir, cfg, diag_specific_models, launch_args):
        self.variable = variable
        self.var_def_dir = var_def_dir
        self.var_attr = var_attr
        self.field = field
        self.diag_scripts = diag_scripts
        self.diag_script_cfg_dir = cfg_dir
        self.diag_script_cfg = cfg
        self.diag_specific_models = diag_specific_models
        self.launcher_arguments = launch_args

    def get_tag_variable(self):
        return self.variable

    def get_tag_attr(self):
        return self.var_attr

    def get_launcher_args(self):
        return self.launcher_arguments

    def get_tag_models(self):
        return self.diag_specific_models

    def get_tag_field_type(self):
        return self.field

    def get_var_def_dir(self):
        return self.var_def_dir

    def get_diag_script_cfg_dir(self):
        return self.diag_script_cfg_dir

    def __iter__(self):
        for idx in range(len(self.diag_scripts)):
            yield self.diag_scripts[idx],\
                os.path.join(self.diag_script_cfg_dir[idx],
                             self.diag_script_cfg[idx]['cfg']),\
                self.diag_specific_models


class Diagnostic:
    """ @brief Class to hold a single Diagnostic instance
    """
    def __init__(self, variable, variable_def_dir, field_type,
                 var_attrs, diag_script, diag_script_cfg, models, launch_args):

        self.diag_script = diag_script
        self.launcher_arguments = launch_args
        self.diag_script_cfg = diag_script_cfg
        variable_list = variable.split(",")
        field_types = field_type.split(",")
        self.variables = []

        # List valid variable attributes
        possible_var_attr_keys = ["mip", "exp", "ref_model", "id", "exclude", "only"]
        attrs = {}
        for key in possible_var_attr_keys:
            tmp_attr = []
            for var_attr in var_attrs:
                if key in var_attr.keys():
                    tmp_attr.append(var_attr[key])
                else:
                    tmp_attr.append("None")
            attrs[key] = tmp_attr

        # Put variable, field + any variable attributes into a dictionary
        for idx in range(len(variable_list)):
            retrieve = operator.itemgetter(idx)
            varfld = dict(zip(['var', 'fld'], [retrieve(mdv) for mdv in [variable_list, field_types]]))
            all_attrs = dict(zip(attrs.keys(), [retrieve(mdv) for mdv in attrs.values()]))

            var_dict = dict(varfld.items() + all_attrs.items())
            self.variables.append(Var(var_dict, "none", "none"))

        self.diag_models = models
        self.variable_def_dir = variable_def_dir

    def __iter__(self):
        for i in range(len(self.variables)):
            os.environ['__ESMValTool_base_var'] = self.variables[i].var
            yield self.variables[i].var,\
                self.variables[i].fld,\
                self.variables[i].mip,\
                self.variables[i].exp, \
                self.variables[i].ref_model, \
                self.variables[i].id, \
                self.variables[i].exclude, \
                self.variables[i].only

    def get_diag_script(self):
        return self.diag_script

    def get_diag_script_cfg(self):
        return self.diag_script_cfg

    def get_launcher_arguments(self):
        return self.launcher_arguments

    def get_variables_list(self):
        return self.variables

    def get_variables(self):
        return [item.var for item in self.variables]

    def get_var_attr_mip(self):
        return [item.mip for item in self.variables]

    def get_var_attr_exp(self):
        return [item.exp for item in self.variables]

    def get_var_attr_ref(self):
        return [item.ref_model for item in self.variables]

    def get_var_attr_exclude(self):
        return [item.exclude for item in self.variables]

    def get_var_attr_only(self):
        return [item.only for item in self.variables]

    def get_var_attrs(self):
        return self.var_attr

    def get_field_types(self):
        return [item.fld for item in self.variables]

    def get_diag_models(self):
        return self.diag_models

    def get_variable_def_dir(self):
        return self.variable_def_dir

    def id_is_explicitly_excluded(self, var, model):
        """ Checks if variable attributes expclitly excludes
            use of this model
        """
        exclude = False
        if hasattr(var, "exclude") and "id" in model.attributes:
            if model.attributes["id"] == var.exclude:
                    exclude = True

        if var.only != "None":
            exclude = True
            if "id" in model.attributes:
                if model.attributes["id"] == var.only:
                    exclude = False

        return exclude

    def find_dimension_entry(self, field):
        """ Checks the number of entries prior to the dimension in
            the 'field_type'-strings, e.g., T2Ms should return 1.
        """
        # Globbed fields defaults to zero, always
        if re.search('[0-9]', field) is None:
            ret_val = 0
        else:
            ret_val = re.search('[0-9]', field).start()
        return ret_val

    def add_base_vars_fields(self, variables, model):
        """ @brief Check for derived variables dependencies
            @param variables List of variables to check for dependencies
            @param model Current model
            @return An instance of the VAR_REQ-class holding all required
            variables/fields

            This function reads the first line of the file
            'var_def/variable.ncl' to check whether the current variable
            is derived from other variables. The syntax is either,

            @code
            ; Requires: var1:field1,var2:field2,...
            @endcode

            or

            @code
            ; Requires: none
            @endcode

            The first and third character in the field can be a
            wild card, '*'.
        """
        dep_vars = []
        model_der = False
        if hasattr(model, "attributes"):
            if "skip_derive_var" in model.attributes.keys():
                model_der = model.attributes["skip_derive_var"]

        for variable in variables:
            keys = [item for item in dir(variable) if re.search("^__", item) is None]
            vars = [getattr(variable, item) for item in dir(variable) if re.search("^__", item) is None]

            var_dict = dict(zip(keys, vars))
            del var_dict['var0']
            del var_dict['fld0']

            f = open(os.path.join(self.variable_def_dir, variable.var
                                  + ".ncl"), 'r')
            for line in f:
                tokens = line.split()

                if "Requires:" in tokens:
                    # If 'none', return orig. field
                    if (tokens[2] == "none" or model_der == "True"):
                        dep_vars.append(Var(var_dict,
                                            "none",
                                            "none"))
                    else:
                        sub_tokens = tokens[2].split(",")
                        for sub in sub_tokens:
                            element = sub.split(":")

                            # Assume first digit is dimension in field type
                            offset = self.find_dimension_entry(element[1]) - 1

                            e_var = element[0]
                            e_fld = element[1]
                            e_fld = (variable.fld[0:offset + 1]
                                     + e_fld[1 + offset]
                                     + variable.fld[2 + offset]
                                     + e_fld[3 + offset:])

                            del keys
                            del vars
                            dep_var = copy.deepcopy(variable)
                            dep_var.var = e_var
                            dep_var.fld = e_fld
                            keys = [item for item in dir(dep_var) if re.search("^__", item) is None]
                            vars = [getattr(dep_var, item) for item in dir(dep_var) if re.search("^__", item) is None]
                            dep_var_dict = dict(zip(keys, vars))

                            dep_vars.append(Var(dep_var_dict,
                                                variable.var,
                                                variable.fld))

        return dep_vars

    def ref_model_lacks_this_variable(self, curr_var, variables, model):
        tmp = dict(zip([var.var for var in variables], [var.ref_model for var in variables]))
        # If we're using a reference model..
        other_ref_vars = [k for k in tmp.keys() if k != curr_var.var]
        other_ref_mods = [tmp[k] for k in other_ref_vars]
        if len(other_ref_mods) > 0:
            if other_ref_mods[0] in model.__dict__['model_line']:
                return True
        return False

    def select_base_vars(self, variables, model, currProject, project_info):
        """ @brief Determine base variables to be read from input file(s)
            @param variables - base variables provided by variable_defs script
            @param model - current model
            @param currProject - current project
            @param project_info - project_info-dictionary
            @return base_vars - variables to be read from input file(s)
        """
        verbosity = project_info['GLOBAL']['verbosity']

        for base_var in variables:
            # Check if this variable should be excluded for this model-id
            if self.id_is_explicitly_excluded(base_var, model):
                continue

            # first try: use base variables provided by variable_defs script
            os.environ['__ESMValTool_base_var'] = base_var.var
            infile = reformat.infile(currProject,
                                     project_info,
                                     base_var,
                                     model)
            precomputed = currProject.get_cf_fullpath(project_info,
                                                      model,
                                                      base_var.fld,
                                                      base_var.var,
                                                      base_var.mip,
                                                      base_var.exp)

            if (len(glob.glob(infile) + glob.glob(precomputed)) == 0):
                info(" No input files found for " + base_var.var +
                     " (" + base_var.fld + ") as " + infile, verbosity, 1)

                base_var.var = base_var.var0
                base_var.fld = base_var.fld0

                # try again with input variable = base variable (non derived)
                infile = reformat.infile(currProject,
                                         project_info,
                                         base_var,
                                         model)
                precomputed = currProject.get_cf_fullpath(project_info,
                                                          model,
                                                          base_var.fld,
                                                          base_var.var,
                                                          base_var.mip,
                                                          base_var.exp)

                if (len(glob.glob(infile) + glob.glob(precomputed)) == 0):
                    raise exceptions.IOError(2, "No input files found in ",
                                             infile)
                else:
                    info(" Using " + base_var.var + " (" + base_var.fld +
                         ") as " + infile, verbosity, 1)
                    base_vars = [base_var]
                    break  # discard other base vars

            else:
                base_vars = variables

        return base_vars

    def get_climate_field_type(self, field):
        """ @brief Update and return the field numbers to be used
            @param field A field number (see the doc/*.pdf:s for
                         further details)
            @return cfield_m - monthly field number
            @return cfield_s - seasonal field number
            @return cfield_a - annual field number
        """
        cfield_m = field.replace("T", "C")
        cfield_s = field.replace("T", "C")
        cfield_a = field.replace("T", "C")
        cfield_s = cfield_s.replace("M", "S")
        cfield_a = cfield_a.replace("M", "A")
        return cfield_m, cfield_s, cfield_a


class AllDiagnostics:
    """ @brief Class to hold a all Diagnostic instances
    """
    def __init__(self):
        self.diagnostics = []

    def extend(self, diagnostic):
        self.diagnostics.extend(diagnostic)

    def __iter__(self):
        for idx in range(len(self.diagnostics)):
            yield self.diagnostics[idx]
