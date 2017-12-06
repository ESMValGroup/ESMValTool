"""
Completely rewritten tool to be able to deal with
the new yaml parser and simplified interface_scripts
toolbox. Author: Valeriu Predoi, University of Reading,
Initial version: August 2017
contact: valeriu.predoi@ncas.ac.uk
"""
import os
import re
from operator import itemgetter

from .auxiliary import writeProjinfoError
from .data_finder import get_output_file


def get_figure_file_names(project_info, model):
    """ @brief Returns names for plots
        @param project_info Current namelist in dictionary format
        @param some model from namelist
    """
    #     return "_".join([
    #         model['project'],
    #         model['name'],
    #         model['mip'],
    #         model['exp'],
    #         model['ensemble'],
    #         str(model['start_year']) + "-" + str(model['end_year']),
    #     ])
    return "_".join([
        model['project'],
        model['name'],
        str(model['start_year']) + "-" + str(model['end_year']),
    ])


def get_cf_fullpath(project_info, model, variable):
    """ @brief Returns the path (only) to the output file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable (defaults to the variable in

            This function specifies the full output path (directory + file) to
            the outupt file to use in the reformat routines and in climate.ncl
    """
    fullpath = get_output_file(project_info, model, variable)
    return fullpath


def get_cf_outpath(project_info, model):
    """ @brief Returns the path (only) to the output file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (output path)

            Standard path: dir_output/preproc_dir/projectname/
                projectname_expname_ens_field_var_yrstart-yrend.nc
    """
    outdir1 = project_info['GLOBAL']['preproc_dir']
    outdir2 = model['project']
    return os.path.join(outdir1, outdir2)


def get_dict_key(model):
    """ @brief Returns a unique key based on the model entries provided.
            @param model One of the <model>-tags in the yaml namelist file
            @return A string

            This function creates and returns a key used in a number of NCL
            scripts to refer to a specific dataset, see e.g., the variable
            'cn' in 'interface_scripts/read_data.ncl'
    """

    # allow for different projects

    # CMIP5
    if model['project'] == 'CMIP5':
        dict_key = "_".join([
            model['project'],
            model['name'],
            model['mip'],
            model['exp'],
            model['ensemble'],
            str(model['start_year']),
            str(model['end_year']),
        ])
    else:
        dict_key = "_".join([
            model['project'],
            model['name'],
            str(model['start_year']),
            str(model['end_year']),
        ])

    return dict_key


class ESMValTool_interface(object):
    def __init__(self):
        # Array that for standardizing part of figure file name
        # across various diag_scripts
        self.figfiles_suffix = []

        # Arrays that, togehter with additional variables such as the
        # current "variable" and "field_type", define the input file names
        self.infiles_prefix = []
        self.infiles_suffix = []

        self.repackage_these = [
            "diag_script", "diag_script_cfg", "variables", "field_types",
            "var_attr_mip", "var_attr_exp", "var_attr_ref", "var_attr_exclude",
            "variable_def_dir"
        ]

        self.from_proj_info = ["output_file_type"]

    def __iter__(self):
        for interface_object in vars(self).keys():
            yield interface_object


def get_diag_value(di, projinfomodels, projinfoconfig, v):
    if v == 'diag_script':
        currentry = [s['script'] for s in di.scripts]
    elif v == 'diag_script_cfg':
        currentry = [s['cfg_file'] for s in di.scripts]
    elif v == 'variables':
        currentry = [s['name'] for s in di.variables]
    elif v == 'field_types':
        currentry = [s['field'] for s in di.variables]
    elif v == 'var_attr_mip':
        currentry = [s['mip'] for s in projinfomodels if 'mip' in s.keys()]
    elif v == 'var_attr_exp':
        currentry = [s['exp'] for s in projinfomodels if 'exp' in s.keys()]
    elif v == 'var_attr_ref':
        currentry_list = [s['ref_model'] for s in di.variables]
        currentry = [item for sublist in currentry_list for item in sublist]
    elif v == 'var_attr_exclude':
        currentry = ['False' for s in projinfomodels]
    elif v == 'variable_def_dir':
        currentry = projinfoconfig['var_def_scripts']
    return currentry


class Data_interface(object):
    def __init__(self, project_info):
        """
            @brief Base class for all the *_data_interface classes
            @param project_info Current namelist in dictionary format

            The init routine repackages all configuration data from yml-
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

        indata_root = self.get_data_root()
        if 'AUXILIARIES' in project_info:
            # 'CMIP5_fx' is hardcoded as it is (so far) the only class
            # supporting the 'AUXILIARIES'-tag
            fx_project = getattr(globals()['projects'], 'CMIP5_fx')()
            fx_files = fx_project.get_fx_files(project_info)
            self.interface.fx_keys = list(
                project_info['AUXILIARIES']['FX_files'].fx_files.keys())
            self.interface.fx_values = [
                os.path.join(indata_root, item.get_fullpath())
                for item in project_info['AUXILIARIES']['FX_files']
                .fx_files.values()
            ]

        # Repackage the above arrays, i.e., input file names and output figure
        # names (the latter one is optional to use in the diag_script).
        self.interface.figfiles_suffix,\
            self.interface.infiles,\
            self.interface.fullpaths,\
            self.interface.infile_paths\
            = self.get_interface_fileinfo(project_info)

        # Repackage model specific data, i.e. the <model>-tags from the
        # yml-namelist files
        model_specifiers, models, model_attr_id, model_attr_skip = \
            self.get_modelinfo(project_info)
        for modelpart in model_specifiers:
            current_column = list(
                map(itemgetter(model_specifiers.index(modelpart)), models))
            vars(self.interface)["models_" + modelpart] = current_column

        vars(self.interface)["model_attr_id"] = model_attr_id
        vars(self.interface)["model_attr_skip"] = model_attr_skip

        # Extract/repackage some project_info entries directly
        for var in self.interface.repackage_these:
            if "currDiag" in project_info['RUNTIME']:
                currDiag = project_info['RUNTIME']['currDiag']
                curr_entry = get_diag_value(currDiag,
                                            project_info['ALLMODELS'],
                                            project_info['CONFIG'], var)
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

        for model in project_info['ALLMODELS']:
            currProject = model

            figfiles_suffix.append(get_figure_file_names(project_info, model))

            # Get reformatted infiles
            # need the input dict because of old variable derivation
            infile_fullpaths.append(
                get_cf_fullpath(
                    project_info,
                    model,
                    variable={
                        'name': "${VARIABLE}",
                        'field': "${FIELD}"
                    }))
            singfile = get_cf_fullpath(
                project_info,
                model,
                variable={
                    'name': "${VARIABLE}",
                    'field': "${FIELD}"
                }).split('/')[-1]
            infiles.append(singfile)
            infile_paths.append(get_cf_outpath(project_info, model))

        return figfiles_suffix,\
            infiles,\
            infile_fullpaths,\
            infile_paths

    def get_modelinfo(self, project_info):
        """ @brief Extracts the <model>-tag info from the namelists
            @param project_info Current namelist in dictionary format
                                environment variables

            Extracts and collects the model-entries from the <model>-tags
            into Python arrays such that they can be written to file more
            easily.
        """
        # condition: keys for all models are the same !!!
        # this is, in fact, bullshit, because we should
        # allow for a variety of model specifiers. Hacking it
        # here (VP)
        # model_specifiers = project_info['ALLMODELS'][0].keys()

        model_specifiers = [
            'project',
            'start_year',
            'name',
            'exp',
            'mip',
            'end_year',
            'ref',
            'ensemble',
        ]

        # OBS: {'project': 'OBS', 'start_year': 2000, 'version': 1,
        #  'name': 'ERA-Interim', 'ref': 'ERA-Interim', 'tier': 3,
        #  'end_year': 2002, 'type': 'reanaly'}
        # CMIP5 (default) : {'project': 'CMIP5', 'start_year': 2000,
        #  'name': 'bcc-csm1-1', 'exp': 'historical', 'mip': 'Amon',
        #  'end_year': 2002, 'ref': 'ERA-Interim', 'ensemble': 'r1i1p1'}

        # this is a bit hacky
        # but makes sure the key - val order stays fixed
        models = []
        for model in project_info['ALLMODELS']:
            # cmip5
            if model['project'] == 'CMIP5':
                mdls = [
                    model['project'],
                    model['start_year'],
                    model['name'],
                    model['exp'],
                    model['mip'],
                    model['end_year'],
                    model['ref'],
                    model['ensemble'],
                ]

            # obs
            if model['project'] == 'OBS':
                mdls = [
                    model['project'],
                    model['start_year'],
                    model['name'],
                    'exp',
                    'mip',
                    model['end_year'],
                    model['ref'],
                    'ensemble',
                ]

            models.append(mdls)

        model_ids = []
        for model in project_info['ALLMODELS']:
            # print('Model is:')
            # print(model)
            if "id" in model.keys():
                model_ids.append(model['id'])
            if "ref" in model.keys():
                model_ids.append(model['ref'])
            else:
                model_ids.append("None")

        model_skips = []
        # Collect instances when to skip derived variable calculations
        for model in project_info['ALLMODELS']:
            if "skip_derive_var" in model.keys():
                model_skips.append(model["skip_derive_var"])
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
        variable_info = [
            remove_comments.sub(r'\1', entry) for entry in var_def
        ]

        # Remove lines in the "variable_def"-file missing the
        # "variable_info"-attribute
        variable_info_entry_regex = re.compile("variable_info@.*=.*")
        variable_info = [
            entry for entry in variable_info
            if variable_info_entry_regex.search(entry) is not None
        ]

        variable_info = [
            re.sub("variable_info@", "", entry) for entry in variable_info
        ]

        # Split variable_info attribute strings to a [key, value]-list
        regexp_equal = re.compile("\s*=\s*")
        variable_info = [regexp_equal.split(entry) for entry in variable_info]

        return variable_info_true, variable_info

    def write_env_projinfo(self, project_info):
        """ @brief Writes YML-file information to env. variables
            @param project_info Current namelist in dictionary format

            Information between Python and NCL is (partially) exchanged
            through environment variables. This function write the
            project_info-dictionary content to environment variables
            prefixed with "ESMValTool_".
        """

        for section_key in ['GLOBAL', 'RUNTIME', 'TEMPORARY']:
            if section_key in project_info.keys():
                for key in project_info[section_key]:
                    # Check and fail on duplicate entries
                    if "ESMValTool_" + key in os.environ:
                        raise writeProjinfoError("Environment variable " +
                                                 "'ESMValTool_" + key +
                                                 "' already defined")
                    os.environ["ESMValTool_" + key] = \
                        str(project_info[section_key][key])


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
        for model in project_info['ALLMODELS']:
            currProject = model

            self.interface.dict_keys.append(get_dict_key(model))
        # Remove the diag_script_cfg entry if it is not NCL code
        # (e.g., R, python, etc..)
        class_prefix = re.search("([a-zA-Z]*)_.*",
                                 self.__class__.__name__).group(1).lower()
        class_regex = re.compile(class_prefix + '$')
        if 'diag_script_cfg' in vars(self.interface).keys():
            if class_regex.search(self.interface.diag_script_cfg[0]) is None:
                del self.interface.diag_script_cfg

    def write_data_to_interface(self):
        """ @brief Write the configuration data to NCL format

            This routine writes the configuration data from the yml-,
            diagnostic_def/-files to NCL formatted files the
            interface_data/-folder. These files are then read by the
            NCL diag_scripts.
        """

        template = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "interface_scripts",
            "templates",
            "ncl_interface_templates",
            "ncl.tmpl",
        )
        interface_data = self.project_info['RUNTIME']['interface_data']
        fsource = open(template, "r")
        ftarget = open(os.path.join(interface_data, "ncl.interface"), "w")

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
                        # but force it to be a string (VP)
                        place_holder = [
                            '"' + str(ph) + '"' for ph in place_holder
                        ]

                    place_holder = [
                        item1 + item2
                        for item1, item2 in zip(lhs, place_holder)
                    ]
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
                variable_info_file = os.path.join(interface_data,
                                                  curr_var + "_info.tmp")
                fvarinfo = open(variable_info_file, "w")

                if variable_info_true:
                    # A-laue_ax+
                    # Attributes of "variable_info" that are arrays might
                    # cause problems if more than one variable is used by
                    # a diagnostic script as this effectively leads to a
                    # redefinition of the already defined variable attributes.
                    # The redefinition will fail if the number of array
                    # elements does not match the "new" number of array
                    # elements.
                    # Work-around: delete "variable_info" if already defined.
                    fvarinfo.write('if (isvar("variable_info")) then\n')
                    fvarinfo.write('    delete(variable_info)\n')
                    fvarinfo.write('end if\n')
                    # A-laue_ax-
                    variable_info = [
                        "variable_info@" + key + "=" + value
                        for key, value in variable_info
                    ]

                    fvarinfo.write('variable_info = True\n' +
                                   '\n'.join(variable_info))
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
        template = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "interface_scripts",
            "templates",
            "r_interface_templates",
            "r.tmpl",
        )
        interface_data = self.project_info['RUNTIME']['interface_data']
        fsource = open(template, "r")
        ftarget = open(os.path.join(interface_data, "r.interface"), "w")
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
                    place_holder = [
                        item1 + item2
                        for item1, item2 in zip(lhs, place_holder)
                    ]
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
                variable_info_file = os.path.join(interface_data,
                                                  curr_var + "_info.tmp")
                fvarinfo = open(variable_info_file, "w")

                if variable_info_true:
                    variable_info = [
                        "variable_info@" + key + " <- " + value
                        for key, value in variable_info
                    ]

                    fvarinfo.write('variable_info <- True\n' +
                                   '\n'.join(variable_info))
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


def write_data_interface(executable, project_info):
    """ @brief Write Python data structures to target script format interface
        @param executable String pointing to the script/binary to execute
        @param project_info Current namelist in dictionary format

        Data structures in Python are rewritten to the interface folder in
        a format appropriate for the target script/binary
    """
    suffix = os.path.splitext(executable)[1][1:]
    curr_interface = globals()[suffix.title()
                               + '_data_interface'](project_info)
    curr_interface.write_data_to_interface()
