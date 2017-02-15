from auxiliary import writeProjinfoError
from operator import itemgetter
import os
import pdb
import projects
import re
from auxiliary import info


class ESMValTool_interface(object):
    def __init__(self):
        # Array that for standardizing part of figure file name
        # across various diag_scripts
        self.figfiles_suffix = []

        # Arrays that, togehter with additional variables such as the
        # current "variable" and "field_type", define the input file names
        self.infiles_prefix = []
        self.infiles_suffix = []

        self.repackage_these = ["diag_script", "diag_script_cfg",
                                "variables", "field_types",
                                "var_attr_mip", "var_attr_exp",
                                "var_attr_ref", "var_attr_exclude",
                                "variable_def_dir"]

        self.from_proj_info = ["output_file_type"]

    def __iter__(self):
        for interface_object in vars(self).keys():
            yield interface_object


class Data_interface(object):
    def __init__(self, project_info):
        """
            @brief Base class for all the *_data_interface classes
            @param project_info Current namelist in dictionary format

            The init routine repackages all configuration data from xml-
            and variable_definition-files to simpler Python structures,
            typically arrays. These arrays are used by the child
            *_data_interface to rewrite the configuration information
            a format appropriate for the target script language.
        """
        self.project_info = project_info
        self.interface = ESMValTool_interface()

        # Possible variable attributes (in some classes)
        self.mip = '${MIP}'
        self.exp = '${EXP}'

        self.do_not_remove_these =\
            [item for item in os.listdir('interface_data/')
             if re.search('.*_interface_templates', item) is not None]

        self.do_not_remove_these.insert(0, "README")
        self.do_not_remove_these.insert(0, ".svn")

        indata_root = self.get_data_root()
        if 'AUXILIARIES' in project_info:
            # 'CMIP5_fx' is hardcoded as it is (so far) the only class supporting
            # the 'AUXILIARIES'-tag
            fx_project = getattr(globals()['projects'], 'CMIP5_fx')()
            fx_files = fx_project.get_fx_files(project_info)
            self.interface.fx_keys = project_info['AUXILIARIES']['FX_files'].fx_files.keys()
            self.interface.fx_values = [os.path.join(indata_root, item.get_fullpath()) for item in project_info['AUXILIARIES']['FX_files'].fx_files.values()]

        # Repackage the above arrays, i.e., input file names and output figure
        # names (the latter one is optional to use in the diag_script).
        self.interface.figfiles_suffix,\
            self.interface.infiles,\
            self.interface.fullpaths,\
            self.interface.infile_paths\
            = self.get_interface_fileinfo(project_info)

        # Repackage model specific data, i.e. the <model>-tags from the
        # xml-namelist files
        model_specifiers, models, model_attr_id, model_attr_skip = self.get_modelinfo(project_info)
        for modelpart in model_specifiers:
            current_column = map(itemgetter(model_specifiers.index(modelpart)),
                                 models)
            vars(self.interface)["models_" + modelpart] = current_column

        vars(self.interface)["model_attr_id"] = model_attr_id
        vars(self.interface)["model_attr_skip"] = model_attr_skip

        # Extract/repackage some project_info entries directly
        for var in self.interface.repackage_these:
            if "currDiag" in project_info['RUNTIME']:
                currDiag = project_info['RUNTIME']['currDiag']
                curr_entry = getattr(currDiag, "get_" + var)()
                if isinstance(curr_entry, list):
                    vars(self.interface)[var] = curr_entry
                else:
                    vars(self.interface)[var] = [curr_entry]

        for proj_info_key in ['GLOBAL', 'RUNTIME', 'TEMPORARY']:
            if proj_info_key in project_info.keys():
                for var in project_info[proj_info_key].keys():
                    curr_entry = project_info[proj_info_key][var]
                    if isinstance(curr_entry, list):
                        vars(self.interface)[var] = curr_entry
                    else:
                        vars(self.interface)[var] = [curr_entry]

    def get_data_root(self):
        """ @brief Checks if indata root is defined as env. variable
            @return A string (indata root folder)

            This function checks if the environment variable
            'ESMValTool_data_root' is set, if so it is used
            as the root folder (for all) input data sets.
            If not set it defaults to "/".
        """
        if "ESMValTool_data_root" in os.environ.keys():
            indata_root = os.environ["ESMValTool_data_root"]
        else:
            indata_root = ""
        return indata_root

    def get_interface_fileinfo(self, project_info):
        """
            @brief Construct part of the filenames used for input-/output files
            @param project_info Current namelist in dictionary format
        """
        infiles = []
        infile_paths = []
        infile_fullpaths = []
        figfiles_suffix = []

        for model in project_info['MODELS']:
            currProject = getattr(globals()['projects'],
                                  model.split_entries()[0])()

            figfiles_suffix.append(currProject.get_figure_file_names(project_info,
                                                                     model,
                                                                     mip=self.mip,
                                                                     exp=self.exp))

            # Get reformatted infiles
            infile_fullpaths.append(currProject.get_cf_fullpath(project_info, model,
                                                                field="${FIELD}",
                                                                variable="${VARIABLE}",
                                                                mip=self.mip,
                                                                exp=self.exp))

            infiles.append(currProject.get_cf_outfile(project_info, model,
                                                      field="${FIELD}",
                                                      variable="${VARIABLE}",
                                                      mip=self.mip,
                                                      exp=self.exp))
            infile_paths.append(currProject.get_cf_outpath(project_info, model))

        return figfiles_suffix,\
            infiles,\
            infile_fullpaths,\
            infile_paths

    def get_modelinfo(self, project_info):
        """ @brief Extracts the <model>-tag info from the xml-namelists
            @param project_info Current namelist in dictionary format
                                environment variables

            Extracts and collects the model-entries from the <model>-tags
            into Python arrays such that they can be written to file more
            easily.
        """
        # Array to hold the the <model>-tag entry keywords
        model_specifiers = []

        # Collect and extend the model_specifiers array.
        for model in project_info['MODELS']:
            currProject = getattr(globals()['projects'],
                                  model.split_entries()[0])()

            for mspec in currProject.model_specifiers:
                if mspec not in model_specifiers:
                    model_specifiers.append(mspec)
            # Append additional specifiers to list
            if currProject.add_specifier:
                for mspec in currProject.add_specifier:
                    if mspec not in model_specifiers:
                        model_specifiers.append(mspec)

        # Array to hold the entries of the <model>-tag lines
        models = []

        # Collect the values for the different model specifiers. If a specific
        # project lacks a certain specifier it can explicitly added through the
        # add_specifier-array or it is given the value 'No_value'.
        for model in project_info['MODELS']:
            currProject = getattr(globals()['projects'],
                                  model.split_entries()[0])()

            models_tmp = []
            for mspecs in model_specifiers:
                if (currProject.add_specifier
                        and mspecs in currProject.add_specifier):
                    specifier = currProject.add_specifier[mspecs]
                    spec_substitute\
                        = currProject.get_model_subsection(model, specifier)
                    models_tmp.append(spec_substitute)

                elif mspecs in currProject.model_specifiers:
                    models_tmp.append(currProject.get_model_subsection(model,
                                                                       mspecs))
                else:
                    models_tmp.append('No_value')
            models.append(models_tmp)

        model_ids = []
        # Collect ids defined as attributes in the model tags
        for model in project_info['MODELS']:
            if "id" in model.attributes.keys():
                model_ids.append(model.attributes["id"])
            else:
                model_ids.append(currProject.get_model_name(model))

        model_skips = []
        # Collect instances when to skip derived variable calculations
        for model in project_info['MODELS']:
            if "skip_derive_var" in model.attributes.keys():
                model_skips.append(model.attributes["skip_derive_var"])
            else:
                model_skips.append("None")

        return model_specifiers, models, model_ids, model_skips

    def reparse_variable_info(self, variable, variable_def_dir):
        """ @brief Repackages the variable_info to a Python list-list structure
            @param variable List of variables with possible 'variable_info'
            @param variable_def_dir Folder where variable definitions reside

            Repackages the variable specific information that resides in
            the NCL-var_def/-files in the 'variable_info'-attribute.
        """
        variable_info_file = os.path.join(variable_def_dir, variable + ".ncl")

        # Read and parse the existing var_def-file
        variable_info_raw = open(variable_info_file, "r").read()
        var_def = variable_info_raw.split('\n')

        # Check whether the "variable_info"-attribute is True/False
        variable_info_true_regex = re.compile("variable_info\s*=\s*True")
        if variable_info_true_regex.search(variable_info_raw) is not None:
            variable_info_true = True
        else:
            variable_info_true = False

        # Remove comments
        remove_comments = re.compile('(.*?);.*')
        variable_info = [remove_comments.sub(r'\1', entry) for entry in var_def]

        # Remove lines in the "variable_def"-file missing the
        # "variable_info"-attribute
        variable_info_entry_regex = re.compile("variable_info@.*=.*")
        variable_info = [entry for entry in variable_info
                         if variable_info_entry_regex.search(entry) is not None]

        variable_info = [re.sub("variable_info@", "", entry)
                         for entry in variable_info]

        # Split variable_info attribute strings to a [key, value]-list
        regexp_equal = re.compile("\s*=\s*")
        variable_info = [regexp_equal.split(entry) for entry in variable_info]

        return variable_info_true, variable_info

    def write_env_projinfo(self, project_info):
        """ @brief Writes XML-file information to env. variables
            @param project_info Current namelist in dictionary format

            Information between Python and NCL is (partially) exchanged
            through environment variables. This function write the
            project_info-dictionary content to environment variables
            prefixed with "ESMValTool_".
        """
        verbosity = project_info['GLOBAL']['verbosity']

        for section_key in ['GLOBAL', 'RUNTIME', 'TEMPORARY']:
            if section_key in project_info.keys():
                for key in project_info[section_key]:
                    info("writing key to env. variable, key=" + key,
                         verbosity, required_verbosity=11)
                    # Check and fail on duplicate entries
                    if "ESMValTool_" + key in os.environ:
                        raise writeProjinfoError("Environment variable "
                                                 + "'ESMValTool_" + key
                                                 + "' already defined")

                    os.environ["ESMValTool_" + key] = \
                        str(project_info[section_key][key])

    def clean_up_interface_folder(self, exceptions):
        """ @brief Remove files from the interface_data/-folder
            @param exceptions Do not remove these entries
        """
        # Clean up interface-folder
        listdir = os.listdir("./interface_data/")

        # Remove exceptions from listdir
        for listitem in exceptions:
            if listitem in listdir:
                listdir.remove(listitem)

        for direntry in listdir:
            os.remove(os.path.join("./interface_data", direntry))


class Ncl_data_interface(Data_interface):
    """
    This class is responsible for rewriting the information in
    the arrays initialized by the Data_interface base class to files
    that can be loaded into NCL-scripts
    """
    def __init__(self, project_info):
        """ @brief Initiates the NCL data interface class
            @param project_info Current namelist in dictionary format
        """
        Data_interface.__init__(self, project_info)

        # A data structure needed only by NCL, hence defined here instead
        # of in the base class. The dict_keys are are used to keep track
        # of the current data sets
        self.interface.dict_keys = []
        for model in project_info['MODELS']:
            currProject = getattr(globals()['projects'],
                                  model.split_entries()[0])()

            self.interface.dict_keys.append(currProject.get_dict_key(model, self.mip, self.exp))

        # Remove the diag_script_cfg entry if it is not NCL code (e.g., R, python, etc..)
        class_prefix = re.search("([a-zA-Z]*)_.*", self.__class__.__name__).group(1).lower()
        class_regex = re.compile(class_prefix + '$')
        if 'diag_script_cfg' in vars(self.interface).keys():
            if class_regex.search(self.interface.diag_script_cfg[0]) is None:
                del(self.interface.diag_script_cfg)

    def write_data_to_interface(self):
        """ @brief Write the configuration data to NCL format

            This routine writes the configuration data from the xml-,
            diagnostic_def/-files to NCL formatted files the
            interface_data/-folder. These files are then read by the
            NCL diag_scripts.
        """
        self.clean_up_interface_folder(self.do_not_remove_these)

        fsource = open("interface_data/ncl_interface_templates/ncl.tmpl", "r")
        ftarget = open("interface_data/ncl.interface", "w")

        # Replace template file placeholders, <<[A-Z_]+>>, with
        # configuration data
        ncl_src = fsource.readlines()
        for line in ncl_src:
            if re.search("<<[A-Z_]+>>", line):
                # left_hand_side, place_holder, right_hand_side
                lhs, place_holder, rhs \
                    = re.search("(.*)<<([A-Z_]+)>>(.*)", line).group(1, 2, 3)
                if place_holder.lower() in vars(self.interface):
                    place_holder = vars(self.interface)[place_holder.lower()]

                    # Extend 'lhs' array with white space elements
                    len_lhs = len(lhs)
                    lhs = [lhs]
                    lhs.extend([" " * len_lhs] * len(place_holder))

                    # Zip lhs- and placeholder- arrays to a new array. The
                    # whitespace entries in lhs-array provides alignment
                    # if/when data is written to file
                    if isinstance(place_holder[0], int):
                        place_holder = [str(ph) for ph in place_holder]

                    else:  # Assume string
                        place_holder = ['"' + ph + '"' for ph in place_holder]

                    place_holder = [item1 + item2
                                    for item1, item2 in zip(lhs, place_holder)]
                    ftarget.write(', \\\n'.join(place_holder) + rhs + '\n')
                else:
                    pass
            else:
                ftarget.write(line)
        fsource.close()
        ftarget.close()

        # Write proj_info to environment variables
        self.write_env_projinfo(self.project_info)

        # Read and parse the var_def/-file
        if 'variables' in vars(self.interface):
            for curr_var in self.interface.variables:
                variable_def_dir = self.interface.variable_def_dir[0]

                variable_info_true, variable_info\
                    = self.reparse_variable_info(curr_var, variable_def_dir)

                # Write parsed content to temp-file in interface_data
                variable_info_file = "interface_data/"\
                                     + curr_var\
                                     + "_info.tmp"
                fvarinfo = open(variable_info_file, "w")

                if variable_info_true:
                    # A-laue_ax+
                    # Attributes of "variable_info" that are arrays might cause
                    # problems if more than one variable is used by a diagnostic
                    # script as this effectively leads to a redefinition of the
                    # already defined variable attributes. The redefinition will
                    # fail if the number of array elements does not match the
                    # "new" number of array elements.
                    # Work-around: delete "variable_info" if already defined.
                    fvarinfo.write('if (isvar("variable_info")) then\n')
                    fvarinfo.write('    delete(variable_info)\n')
                    fvarinfo.write('end if\n')
                    # A-laue_ax-
                    variable_info = ["variable_info@" + key + "=" + value
                                     for key, value in variable_info]

                    fvarinfo.write('variable_info = True\n' + '\n'.join(variable_info))
                else:
                    fvarinfo.write('variable_info = False')
                fvarinfo.close()


class R_data_interface(Data_interface):
    """
    This class is responsible for rewriting the information in
    the arrays initialized by the Data_interface base class to files
    that can be loaded into R-scripts
    """
    def __init__(self, project_info):
        """ @brief Initiates the R data interface class
            @param project_info Current namelist in dictionary format
        """
        Data_interface.__init__(self, project_info)

    def write_data_to_interface(self):
        """ @brief Write the configuration data to Matlab format
        """
        self.clean_up_interface_folder(self.do_not_remove_these)

        fsource = open("interface_data/r_interface_templates/r.tmpl", "r")
        ftarget = open("interface_data/r.interface", "w")
        line_cont = ', \n'

        # Replace template file placeholders, <<[A-Z_]+>>, with
        # configuration data
        r_src = fsource.readlines()
        for line in r_src:
            if re.search("<<[A-Z_]+>>", line):
                # left_hand_side, place_holder, right_hand_side
                lhs, place_holder, rhs \
                    = re.search("(.*)<<([A-Z_]+)>>(.*)", line).group(1, 2, 3)
                if place_holder.lower() in vars(self.interface):
                    place_holder = vars(self.interface)[place_holder.lower()]

                    # Extend 'lhs' array with white space elements
                    len_lhs = len(lhs)
                    lhs = [lhs]
                    lhs.extend([" " * len_lhs] * len(place_holder))

                    # Zip lhs- and placeholder- arrays to a new array. The
                    # whitespace entries in lhs-array provides alignment
                    # if/when data is written to file
                    place_holder = ['"' + str(ph) + '"' for ph in place_holder]
                    place_holder = [item1 + item2
                                    for item1, item2 in zip(lhs, place_holder)]
                    ftarget.write(line_cont.join(place_holder) + rhs + '\n')
                else:
                    pass
            else:
                ftarget.write(line)
        fsource.close()
        ftarget.close()

        # Write proj_info to environment variables
        self.write_env_projinfo(self.project_info)
        # Read and parse the var_def/-file
        if 'variables' in vars(self.interface):
            for curr_var in self.interface.variables:
                variable_def_dir = self.interface.variable_def_dir[0]

                variable_info_true, variable_info\
                    = self.reparse_variable_info(curr_var, variable_def_dir)

                # Write parsed content to temp-file in data_interface
                variable_info_file = "interface_data/"\
                                     + curr_var\
                                     + "_info.tmp"
                fvarinfo = open(variable_info_file, "w")

                if variable_info_true:
                    variable_info = ["variable_info@" + key + " <- " + value
                                     for key, value in variable_info]

                    fvarinfo.write('variable_info <- True\n' + '\n'.join(variable_info))
                else:
                    fvarinfo.write('variable_info <- False')
                fvarinfo.close()


class Py_data_interface(Data_interface):
    """
    This class is responsible for rewriting the information in
    the arrays initialized by the Data_interface base class to files
    that can be loaded into Python-scripts
    """
    def __init__(self, project_info):
        """ @brief Initiates the Python data interface class
            @param project_info Current namelist in dictionary format
        """
        Data_interface.__init__(self, project_info)

    def write_data_to_interface(self):
        """ @brief Write the configuration data to Matlab format
        """
        self.clean_up_interface_folder(self.do_not_remove_these)
