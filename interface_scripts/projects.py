from auxiliary import info
import data_interface
import exceptions
import os
import launchers
import pdb
import re
import glob


class Project:
    """ @brief Base class for all ESMValTool projects
    """
    def __init__(self):
        self.add_specifier = {}

    def get_model_subsection(self, model, model_section):
        """ @brief Retrieve a specific model entry from a <model> tag line
            @param model One of the <model>-tags in the XML namelist file
            @param model_section Which of the model entries to retrieve
        """
        section = \
            model.split_entries()[self.model_specifiers.index(model_section)]
        return section

    def get_model_sections(self, model):
        """ @brief Retrieve all model entries from a <model> tag line
            @param model One of the <model>-tags in the XML namelist file
        """
        model_sections = [self.get_model_subsection(model, modelpart)
                          for modelpart in self.model_specifiers]
        model_sect_dict = dict(zip(self.model_specifiers, model_sections))

        # Replace the ${VARIABLE}-placeholder in infile dir with base var
        if re.search("\$\{VARIABLE\}", model_sect_dict['dir']) is not None:
            base_var = os.environ['__ESMValTool_base_var']
            model_sect_dict['dir'] \
                = re.sub("\$\{VARIABLE\}", base_var, model_sect_dict['dir'])
        return model_sect_dict

    def get_cf_outpath(self, project_info, model):
        """ @brief Returns the path (only) to the output file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (output path)

            This function specifies the output path (only) to the outupt file
            to use in the reformat routines and in climate.ncl
        """
        msd = self.get_model_sections(model)
        outdir = os.path.join(project_info['GLOBAL']['climo_dir'],
                              msd['project'])
        return outdir

    def get_cf_fullpath(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path (only) to the output file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable (defaults to the variable in
                            project_info)
            @return A string (output path)

            This function specifies the full output path (directory + file) to
            the outupt file to use in the reformat routines and in climate.ncl
        """
        outdir = self.get_cf_outpath(project_info, model)
        outfile = self.get_cf_outfile(project_info, model, field, variable,
                                      mip, exp)

        return os.path.join(outdir, outfile)

    def get_fx_files(self, project_info):
        """ @brief Returns the path to all fx files
            @param project_info Current namelist in dictionary format
            @return A dict (fx-ID: fx-path)

            This function looks for the fx files defined in the
            in the AUXILIARIES namelist section.
        """
        if not "fx_file_ID" in self.model_specifiers:
            msg = "This class doesn't support specifing fx-files,"\
                   + " check interface_scripts/projects.py"
            raise exceptions.RuntimeError(msg)

        indata_root = self.get_data_root()

        curr_fx = project_info["AUXILIARIES"]["FX_files"].fx_files
        new_fx = {}

        for key in curr_fx:
            new_fx[key] = os.path.join(indata_root, curr_fx[key].get_fullpath())

        return new_fx

    def get_fx_file(self, project_info, model):
        """ @brief Returns the path to specific fx file
            @param project_info Current namelist in dictionary format
            @return A list (fx-file path)

            This function looks for a specific fx files  defined
            in the AUXILIARIES namelist section.
        """
        if not "fx_file_ID" in self.model_specifiers:
            msg = "This class doesn't support specifing fx-files,"\
                   + " check interface_scripts/projects.py"
            raise exceptions.RuntimeError(msg)

        fx_files = self.get_fx_files(project_info)

        # Get gridfile ID specified in model entry
        fx_ID = self.get_model_subsection(model, "fx_file_ID")

        # Find this grid file within project_info
        if not fx_ID in project_info["AUXILIARIES"]["FX_files"]:
            msg = "fx_file_ID '" + fx_ID +\
                  "' not found in any <fx_file> entry within " +\
                  "<AUXILIARIES> section of namelist."
            raise exceptions.RuntimeError(msg)

        return fx_files[fx_ID]

    def get_cf_areafile(self, project_info, model):
        """ @brief Returns the path to the areacello file used for
                   ocean variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the areafile of the ocean grid
        """
        return None

    def get_cf_porofile(self, project_info, model):
        """ @brief Returns the path to the mrsofc file used for
                   land variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the porosity file
        """
        return None

    def get_cf_lmaskfile(self, project_info, model):
        """ @brief Returns the path to the sftlf file used for masking
                   land variables (regular grid)
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (maskfile path)

            This function looks for the maskfile of the ocean grid
        """
        return None

    def get_cf_omaskfile(self, project_info, model):
        """ @brief Returns the path to the sftof file used for masking
                   ocean variables (irregular grid)
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (maskfile path)

            This function looks for the maskfile of the ocean grid
        """
        return None

    def get_project_basename(self):
        """ Define a function but not the actual variable 'basename'
        """
        return self.basename

    def get_cf_sections(self, model):
        """ @brief Return the sections from the <model> tag needed in reformat
            @param model One of the <model>-tags in the XML namelist file
            @return Six strings

            This function returns the sections from the <model> tag that
            are needed in the 'reformat'-routine.
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        return msd['project'],\
            msd['name'],\
            msd['ensemble'],\
            msd['start_year'],\
            msd['end_year'],\
            msd['dir']

    def get_model_name(self, model):
        return self.get_model_subsection(model, "name")

    def get_project_name(self, model):
        return self.get_model_subsection(model, "project")

    def rewrite_mip_exp(self, msd, mip, exp):
        if mip != "None":
            msd["mip"] = mip
        if exp != "None":
            msd["experiment"] = exp
        return msd

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

    def get_project_variable_name(self, model, variable):
        """ @brief Returns the project specific variable name
            @param model One of the <model>-tags in the XML namelist file
            @param variable The variable (ESMValTool standard)
            @return A string (corresponding project variable name)

            This function looks for a project specific dictionary file
            to translate the standard (ESMValTool) variable name into
            the corresponding variable name of the actual project.
        """

        project = self.get_project_name(model)
        names_file = 'reformat_scripts/fixes/names_' + project + '.dat'

        project_var = variable

        if os.path.isfile(names_file):
            with open(names_file) as f:
                for line in f:
                    first_char = line[:1]
                    if first_char != "#":
                        names = line.split('|')
                        std = names[0].strip()
                        if std == variable:
                            project_var = names[1].strip()

        return project_var


class OBS(Project):
    """ @brief Class defining the specifics for observational data
        case can be insitu, sat, ground, reanaly
    """

    def __init__(self):
        Project.__init__(self)

        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "case_name",
                                 "ensemble",
                                 "start_year",
                                 "end_year",
                                 "dir"]

        """ Define the 'basename'-variable explicitly
        """
        self.basename = "OBS"

    def get_dict_key(self, model, mip, exp):
        """ @brief Returns a unique key based on the model entries provided.
            @param model One of the <model>-tags in the XML namelist file
            @return A string

            This function creates and returns a key used in a number of NCL
            scripts to refer to a specific dataset, see e.g., the variable
            'cn' in 'interface_scripts/read_data.ncl'
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        dict_key = '_'.join([msd['project'],
                             msd['name'],
                             msd['case_name'],
                             msd['ensemble'],
                             msd['start_year'],
                             msd['end_year']])
        return dict_key

    def get_figure_file_names(self, project_info, model, mip, exp):
        """ @brief Returns the full path used for intermediate data storage
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param variable_attributes The variable attributes
            @return A string

            This function creates and returns the full path used in a number
            of NCL scripts to write intermediate data sets.
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        return "_".join([msd['project'],
                         msd['name'],
                         msd['case_name'],
                         msd['ensemble'],
                         msd['start_year']]) + "-" + msd['end_year']

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat routines
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable (defaults to the variable in
                            project_info)
            @param variable_attributes The variable attributes
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """
        indata_root = self.get_data_root()

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        indir = os.path.join(indata_root,
                             msd["dir"])

        variable = self.get_project_variable_name(model, variable)

        infile = '_'.join([msd['project'],
                           msd['name'],
                           msd['case_name'],
                           msd['ensemble'],
                           field,
                           variable]) + '_??????-??????.nc'

        PROJECT = msd['project']
        if ('OBS_gridfile' in msd["project"]):
            PROJECT = 'OBS'
            infile = '_'.join([PROJECT,
                               msd['name'],
                               msd['case_name'],
                               msd['ensemble'],
                               field,
                               variable]) + '_??????-??????.nc'

        return indir, infile

    def get_cf_outfile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path and output file used in reformat
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable (defaults to the variable in
                            project_info)
            @return A string (output path)

            This function specifies the output file to use in reformat
            and in climate.ncl
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        outdir = self.get_cf_outpath(project_info, model)

        outfile = '_'.join([msd['project'],
                            msd['case_name'],
                            msd['name'],
                            msd['ensemble'],
                            field,
                            variable,
                            msd['start_year']]) + '-' + msd['end_year'] + '.nc'

        return os.path.join(outdir, outfile)


class OBS_gridfile(OBS):
    """ @brief Follows the DLR directory structure for OBS data,
               + explicit gridfile
    """
    def __init__(self):
        OBS.__init__(self)

        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "case_name",
                                 "ensemble",
                                 "start_year",
                                 "end_year",
                                 "dir",
                                 "gridfile"]

    def get_cf_areafile(self, project_info, model):
        """ @brief Returns the path to the areacello/areacella file
                   used for variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the areafile of the grid
        """
        return self.get_model_subsection(model, "gridfile")


class obs4mips(Project):
    """ @brief Class defining the specifics of the obs4mips project
    """

    def __init__(self):
        Project.__init__(self)

        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "level",
                                 "ensemble",
                                 "start_year",
                                 "end_year",
                                 "dir"]

        """ Define the 'basename'-variable explicitly
        """
        self.basename = "obs4mips"

    def get_dict_key(self, model, mip, exp):
        """ @brief Returns a unique key based on the model entries provided.
            @param model One of the <model>-tags in the XML namelist file
            @return A string

            This function creates and returns a key used in a number of NCL
            scripts to refer to a specific dataset, see e.g., the variable
            'cn' in 'interface_scripts/read_data.ncl'
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        dict_key = "_".join([msd['project'],
                             msd['name'],
                             msd['level'],
                             msd['ensemble'],
                             msd['start_year'],
                             msd['end_year']])
        return dict_key

    def get_figure_file_names(self, project_info, model, mip, exp):
        """ @brief Returns the full path used for intermediate data storage
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string

            This function creates and returns the full path used in a number
            of NCL scripts to write intermediate data sets.
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        return "_".join([msd['project'],
                         msd['name'],
                         msd['level'],
                         msd['ensemble'],
                         msd['start_year']]) + "-" + msd['end_year']

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        indir = os.path.join(msd["dir"],
                             msd["name"])

        infile = '_'.join([variable,
                           msd['name'],
                           msd['level'],
                           msd['ensemble']]) + '.nc'

        if (not os.path.isfile(os.path.join(indir, infile))):
            infile = '_'.join([variable,
                               msd['name'],
                               msd['level'],
                               msd['ensemble']]) + '*.nc'

        return indir, infile

    def get_cf_outfile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the output file used in the reformat routines
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @return A string (output path)

            This function specifies the output file to use in the reformat
            routines and in climate.ncl
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        outfile = '_'.join([msd['project'],
                            msd['name'],
                            msd['level'],
                            msd['ensemble'],
                            field,
                            variable,
                            msd['start_year']]) + '-' + msd['end_year'] + '.nc'

        return outfile


class ana4mips(Project):
    """ @brief Class defining the specific characteristics of the CMIP5 project
    """
    def __init__(self):
        Project.__init__(self)

        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "mip",
                                 "ensemble",
                                 "start_year",
                                 "end_year",
                                 "dir"]

        """ Define the 'basename'-variable explicitly
        """
        self.basename = "ana4mips"

    def get_dict_key(self, model, mip, exp):
        """ @brief Returns a unique key based on the model entries provided.
            @param model One of the <model>-tags in the XML namelist file
            @return A string

            This function creates and returns a key used in a number of NCL
            scripts to refer to a specific dataset, see e.g., the variable
            'cn' in 'interface_scripts/read_data.ncl'
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        dict_key = "_".join([msd['project'],
                             msd['name'],
                             msd['mip'],
                             msd['ensemble'],
                             msd['start_year'],
                             msd['end_year']])
        return dict_key

    def get_figure_file_names(self, project_info, model, mip, exp):
        """ @brief Returns the full path used for intermediate data storage
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param variable_attributes The variable attributes
            @return A string

            This function creates and returns the full path used in a number
            of NCL scripts to write intermediate data sets.
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        return "_".join([msd['project'],
                         msd['name'],
                         msd['mip'],
                         msd['ensemble'],
                         msd['start_year']]) + "-" + msd['end_year']

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @param variable_attributes The variable attributes
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        indir = os.path.join(msd['dir'],
                             msd['name'])

        variable = self.get_project_variable_name(model, variable)

        infile = '_'.join([variable,
                           msd['mip'],
                           msd['ensemble'],
                           msd['name']]) + '.nc'

        if (not os.path.isfile(os.path.join(indir, infile))):
            infile = '_'.join([variable,
                               msd['mip'],
                               msd['ensemble'],
                               msd['name']]) + '*.nc'

        return indir, infile

    def get_cf_outfile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the output file used in the reformat routines
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @return A string (output path)

            This function specifies the output file to use in the reformat
            routines and in climate.ncl
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        # Overide some model lines with settings from the variable attributes
        msd = self.rewrite_mip_exp(msd, mip, exp)

        outfile = '_'.join([msd['project'],
                            msd['mip'],
                            msd['ensemble'],
                            msd['name'],
                            field,
                            variable,
                            msd['start_year']]) + '-' + msd['end_year'] + '.nc'

        return outfile


class CMIP5(Project):
    """ @brief Class defining the specific characteristics of the CMIP5 project
    """
    def __init__(self):
        Project.__init__(self)

        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "mip",
                                 "experiment",
                                 "ensemble",
                                 "start_year",
                                 "end_year",
                                 "dir"]

        self.add_specifier['case_name'] = 'experiment'

        """ Define the 'basename'-variable explicitly
        """
        self.basename = "CMIP5"

    def get_model_mip(self, model):
        return self.get_model_subsection(model, "mip")

    def get_model_exp(self, model):
        return self.get_model_subsection(model, "experiment")

    def get_model_ens(self, model):
        return self.get_model_subsection(model, "ensemble")

    def get_model_start_year(self, model):
        return self.get_model_subsection(model, "start_year")

    def get_model_end_year(self, model):
        return self.get_model_subsection(model, "end_year")

    def get_model_dir(self, model):
        return self.get_model_subsection(model, "dir")

    def get_dict_key(self, model, mip, exp):
        """ @brief Returns a unique key based on the model entries provided.
            @param model One of the <model>-tags in the XML namelist file
            @return A string

            This function creates and returns a key used in a number of NCL
            scripts to refer to a specific dataset, see e.g., the variable
            'cn' in 'interface_scripts/read_data.ncl'
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        dict_key = "_".join([msd['project'],
                             msd['name'],
                             msd['mip'],
                             msd['experiment'],
                             msd['ensemble'],
                             msd['start_year'],
                             msd['end_year']])
        return dict_key

    def get_figure_file_names(self, project_info, model, mip, exp):
        """ @brief Returns the full path used for intermediate data storage
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param variable_attributes The variable attributes
            @return A string

            This function creates and returns the full path used in a number
            of NCL scripts to write intermediate data sets.
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        return "_".join([msd['project'],
                         msd['name'],
                         msd['mip'],
                         msd['experiment'],
                         msd['ensemble'],
                         msd['start_year']]) + "-" + msd['end_year']

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @param variable_attributes The variable attributes
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """
        indata_root = self.get_data_root()

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        indir = os.path.join(indata_root, msd["dir"])

        variable = self.get_project_variable_name(model, variable)

        infile = '_'.join([variable,
                           msd['mip'],
                           msd['name'],
                           msd['experiment'],
                           msd['ensemble']]) + '.nc'

        if (not os.path.isfile(os.path.join(indir, infile))):
            infile = '_'.join([variable,
                               msd['mip'],
                               msd['name'],
                               msd['experiment'],
                               msd['ensemble']]) + '*.nc'

        return indir, infile

    def get_cf_areafile(self, project_info, model):
        """ @brief Returns the path to the areacello file
                   used for ocean variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the areafile of the ocean grid
        """
        msd = self.get_model_sections(model)

        areadir = msd["dir"]

        areafile = 'areacello_fx_' + msd["name"] + "_" + msd["experiment"] + \
                   "_r0i0p0.nc"

        return os.path.join(areadir, areafile)

    def get_cf_porofile(self, project_info, model):
        """ @brief Returns the path to the mrsofc file used for
                   land variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the porosity file
        """
        msd = self.get_model_sections(model)

        maskdir = msd["dir"]

        maskfile = 'mrsofc_fx_' + msd["name"] + "_" + msd["experiment"] + \
                   "_r0i0p0.nc"

        return os.path.join(maskdir, maskfile)

    def get_cf_lmaskfile(self, project_info, model):
        """ @brief Returns the path to the sftlf file used for masking
                   land variables (regular grid)
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (maskfile path)

            This function looks for the maskfile of the ocean grid
        """
        msd = self.get_model_sections(model)

        maskdir = msd["dir"]

        maskfile = 'sftlf_fx_' + msd["name"] + "_" + msd["experiment"] + \
                   "_r0i0p0.nc"

        return os.path.join(maskdir, maskfile)

    def get_cf_omaskfile(self, project_info, model):
        """ @brief Returns the path to the sftof file used for masking
                   ocean variables (irregular grid)
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (maskfile path)

            This function looks for the maskfile of the ocean grid
        """
        msd = self.get_model_sections(model)

        maskdir = msd["dir"]

        maskfile = 'sftof_fx_' + msd["name"] + "_" + msd["experiment"] + \
                   "_r0i0p0.nc"

        return os.path.join(maskdir, maskfile)

    def get_cf_outfile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the output file used in the reformat routines
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @return A string (output path)

            This function specifies the output file to use in the reformat
            routines and in climate.ncl
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        # Overide some model lines with settings from the variable attributes
        msd = self.rewrite_mip_exp(msd, mip, exp)

        outfile = '_'.join([msd['project'],
                            msd['mip'],
                            msd['experiment'],
                            msd['name'],
                            msd['ensemble'],
                            field,
                            variable,
                            msd['start_year']]) + '-' + msd['end_year'] + '.nc'

        return outfile


class CMIP5_gridfile(CMIP5):
    def __init__(self):
        CMIP5.__init__(self)
        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "mip",
                                 "experiment",
                                 "ensemble",
                                 "start_year",
                                 "end_year",
                                 "dir",
                                 "gridfile"]

        self.add_specifier['case_name'] = 'experiment'

    def get_cf_areafile(self, project_info, model):
        """ @brief Returns the path to the areacello file
                   used for ocean variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the areafile of the ocean grid
        """
        return self.get_model_subsection(model, "gridfile")


class CMIP5_fx(CMIP5):
    """ @brief Class handling a fx-file handler on the <model>-tag line
    """
    def __init__(self):
        CMIP5.__init__(self)
        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "mip",
                                 "experiment",
                                 "ensemble",
                                 "start_year",
                                 "end_year",
                                 "dir",
                                 "fx_file_ID"]

        self.add_specifier['case_name'] = 'fx_file_ID'

    def get_cf_areafile(self, project_info, model):
        """ @brief Returns the path to the areacello file
                   used for ocean variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the areafile of the ocean grid
        """
        return self.get_cf_fx_file(project_info, model)

    def get_cf_fx_file(self, project_info, model):
        """ @brief Returns the path to the areacello file
                   used for ocean variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the areafile of the ocean grid
        """
        fx = self.get_fx_file(project_info, model)
        return os.path.abspath(fx)


class CMIP5_SMHI(CMIP5):
    def __init__(self):
        CMIP5.__init__(self)

        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "mip",
                                 "experiment",
                                 "ensemble",
                                 "start_year",
                                 "end_year",
                                 "freq",
                                 "dir"]

        self.add_specifier['case_name'] = 'experiment'

    def get_cf_areafile(self, project_info, model):
        """ @brief Returns the path to the areacello file
                   used for ocean variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the areafile of the ocean grid
        """
        msd = self.get_model_sections(model)

        areadir = os.path.join(msd["dir"], msd["name"], "fx")
        areafile = 'areacello_fx_' + msd["name"] + "_" + "xxx" + "_r0i0p0.nc"

        return os.path.join(areadir, areafile)

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @param variable_attributes The variable attributes
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        indir = os.path.join(msd["dir"],
                             msd["name"],
                             msd["ensemble"],
                             msd["experiment"],
                             msd["freq"])

        infile = "_".join([variable,
                           msd["mip"],
                           msd["name"],
                           msd["experiment"],
                           msd["ensemble"]]) + ".nc"

        if (not os.path.isfile(os.path.join(indir, infile))):
            infile = "_".join([variable,
                               msd["mip"],
                               msd["name"],
                               msd["experiment"],
                               msd["ensemble"]]) + "*.nc"

        return indir, infile


class CMIP5_ETHZ(CMIP5):
    """ @brief Follows the ETH directory structure for input data
    """
    def __init__(self):
        CMIP5.__init__(self)

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @param variable_attributes The variable attributes
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        indir = os.path.join(msd["dir"],
                             msd["experiment"],
                             msd["mip"],
                             variable,
                             msd["name"],
                             msd["ensemble"])

        infile = "_".join([variable,
                           msd["mip"],
                           msd["name"],
                           msd["experiment"],
                           msd["ensemble"]]) + ".nc"

        if (not os.path.isfile(os.path.join(indir, infile))):
            infile = "_".join([variable,
                               msd["mip"],
                               msd["name"],
                               msd["experiment"],
                               msd["ensemble"]]) + "*.nc"

        return indir, infile

    def get_cf_areafile(self, project_info, model):
        """ @brief Returns the path to the areacello file
                   used for ocean variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the areafile of the ocean grid
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        areadir = os.path.join(msd["dir"],
                               msd["experiment"],
                               'fx',
                               'areacello',
                               msd["name"],
                               'r0i0p0')

        areafile = 'areacello_fx_' + msd["name"] + "_" + msd["experiment"] \
            + "_r0i0p0.nc"

        return os.path.join(areadir, areafile)

    def get_cf_porofile(self, project_info, model):
        """ @brief Returns the path to the mrsofc file used for
                   land variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the porosity file
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        maskdir = os.path.join(msd["dir"],
                               msd["experiment"],
                               'fx',
                               'mrsofc',
                               msd["name"],
                               'r0i0p0')

        maskfile = 'mrsofc_fx_' + msd["name"] + "_" + msd["experiment"] \
            + "_r0i0p0.nc"

        return os.path.join(maskdir, maskfile)

    def get_cf_lmaskfile(self, project_info, model):
        """ @brief Returns the path to the sftlf file used for masking
                   land variables (regular grid)
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (maskfile path)

            This function looks for the areafile of the ocean grid
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        maskdir = os.path.join(msd["dir"],
                               msd["experiment"],
                               'fx',
                               'sftlf',
                               msd["name"],
                               'r0i0p0')

        maskfile = 'sftlf_fx_' + msd["name"] + "_" + msd["experiment"] \
            + "_r0i0p0.nc"

        return os.path.join(maskdir, maskfile)

    def get_cf_omaskfile(self, project_info, model):
        """ @brief Returns the path to the sftof file used for masking
                   ocean variables (irregular grid)
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (maskfile path)

            This function looks for the areafile of the ocean grid
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        maskdir = os.path.join(msd["dir"],
                               msd["experiment"],
                               'fx',
                               'sftof',
                               msd["name"],
                               'r0i0p0')

        maskfile = 'sftof_fx_' + msd["name"] + "_" + msd["experiment"] \
            + "_r0i0p0.nc"

        return os.path.join(maskdir, maskfile)

class MiKlip(Project):
    """ @brief Follows the MiKlip directory structure for input data
               for baseline 1
    """
    def __init__(self):
        Project.__init__(self)

        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "mip",
                                 "experiment",
                                 "ensemble",
                                 "realm",
                                 "start_year",
                                 "end_year",
                                 "dir"]

        self.add_specifier['case_name'] = 'experiment'

        """ Define the 'basename'-variable explicitly
        """
        self.basename = "MiKlip"

    def get_dict_key(self, model, mip, exp):
        """ @brief Returns a unique key based on the model entries provided.
            @param model One of the <model>-tags in the XML namelist file
            @return A string

            This function creates and returns a key used in a number of NCL
            scripts to refer to a specific dataset, see e.g., the variable
            'cn' in 'interface_scripts/read_data.ncl'
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        dict_key = "_".join([msd['project'],
                             msd['name'],
                             msd['mip'],
                             msd['experiment'],
                             msd['ensemble'],
                             msd['realm'],
                             msd['start_year'],
                             msd['end_year']])
        return dict_key

    def get_figure_file_names(self, project_info, model, mip, exp):
        """ @brief Returns the full path used for intermediate data storage
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string

            This function creates and returns the full path used in a number
            of NCL scripts to write intermediate data sets.
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        return "_".join([msd['project'],
                         msd['name'],
                         msd['mip'],
                         msd['experiment'],
                         msd['ensemble'],
                         msd['realm'],
                         msd['start_year']]) + "-" + msd['end_year']

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        if ('6hr' in msd["mip"]):
            mip = '6hr'
        if ('day' in msd["mip"]):
            mip = 'day'
        if ('mon' in msd["mip"]):
            mip = 'mon'
        if ('yr' in msd["mip"]):
            mip = 'yr'
        if ('clim' in msd["mip"]):
            mip = 'monClim'

        indir = os.path.join(msd["dir"],
                             msd["name"],
                             msd["experiment"],
                             mip,
                             msd["realm"],
                             variable,
                             msd["ensemble"])

        infile = "_".join([variable,
                           msd["mip"],
                           msd["name"],
                           msd["experiment"],
                           msd["ensemble"]]) + ".nc"

        if (not os.path.isfile(os.path.join(indir, infile))):
            infile = "_".join([variable,
                               msd["mip"],
                               msd["name"],
                               msd["experiment"],
                               msd["ensemble"]]) + "*.nc"

        return indir, infile

    def get_cf_areafile(self, project_info, model):
        """ @brief Returns the path to the areacello file
                   used for ocean variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the areafile of the ocean grid
        """
        msd = self.get_model_sections(model)

        areadir = os.path.join(msd["dir"],
                               msd["name"],
                               msd["experiment"],
                               'fx',
                               'ocean',
                               'areacello',
                               'r0i0p0')

        areafile = 'areacello_fx_' + msd["name"] + "_" + msd["experiment"] + \
                   "_r0i0p0.nc"

        return os.path.join(areadir, areafile)

    def get_cf_outfile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the output file used in the reformat routines
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @return A string (output path)

            This function specifies the output file to use in the reformat
            routines and in climate.ncl
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        # Overide some model lines with settings from the variable attributes
        msd = self.rewrite_mip_exp(msd, mip, exp)

        outfile = '_'.join([msd['project'],
                            msd['mip'],
                            msd['experiment'],
                            msd['name'],
                            msd['ensemble'],
                            field,
                            variable,
                            msd['start_year']]) + '-' + msd['end_year'] + '.nc'

        return outfile


class MiKlip_baseline0(MiKlip):
    """ @brief Follows the MiKlip directory structure for input data
               + fx file for baseline 0
    """
    def __init__(self):
        MiKlip.__init__(self)

        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "mip",
                                 "experiment",
                                 "ensemble",
                                 "realm",
                                 "start_year",
                                 "end_year",
                                 "dir"]

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        if ('6hr' in msd["mip"]):
            mip = '6hr'
        if ('day' in msd["mip"]):
            mip = 'day'
        if ('mon' in msd["mip"]):
            mip = 'mon'
        if ('yr' in msd["mip"]):
            mip = 'yr'
        if ('clim' in msd["mip"]):
            mip = 'monClim'

        version = 'v20120529'
        if ('historical' in msd["experiment"] or 'rcp' in msd["experiment"]):
            version = 'v20111006'
        if ('MPI-ESM-MR' in msd["name"]):
            version = 'v20120608'
            if ('historical' in msd["experiment"] or 'rcp' in msd["experiment"]):
                if (('r2i1p1' in msd["ensemble"] or 'r3i1p1' in msd["ensemble"]) and 'rcp' in msd["experiment"]):
                    version = 'v20120628'
                else:
                    version = 'v20120503'

        indir = os.path.join(msd["dir"],
                             msd["name"],
                             msd["experiment"],
                             mip,
                             msd["realm"],
                             msd["mip"],
                             msd["ensemble"],
                             version,
                             variable)

        infile = "_".join([variable,
                           msd["mip"],
                           msd["name"],
                           msd["experiment"],
                           msd["ensemble"]]) + ".nc"

        if (not os.path.isfile(os.path.join(indir, infile))):
            infile = "_".join([variable,
                               msd["mip"],
                               msd["name"],
                               msd["experiment"],
                               msd["ensemble"]]) + "*.nc"

        return indir, infile

    def get_cf_areafile(self, project_info, model):
        """ @brief Returns the path to the areacello file
                   used for ocean variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the areafile of the ocean grid
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        areafile = 'areacello_fx_'\
                   + msd["name"]\
                   + "_" + msd["experiment"] \
                   + "_r0i0p0.nc"

        version = 'v20120529'
        if ('historical' in msd["experiment"] or 'rcp' in msd["experiment"]):
            version = 'v20111006'
        if('MPI-ESM-MR' in msd["name"]):
            version = 'v20120608'
            if ('historical' in msd["experiment"] or 'rcp' in msd["experiment"]):
                version = 'v20120503'

        areadir = os.path.join(msd["dir"],
                               msd["name"],
                               msd["experiment"],
                               'fx',
                               'ocean',
                               'fx',
                               'r0i0p0',
                               version,
                               'areacello')

        return os.path.join(areadir, areafile)


class CMIP5_DKRZ(MiKlip):
    """@brief Follows the DKRZ directory structure for input data
              that in turn follows ESGF data organisation.
              + fx file
    """
    def __init__(self):
        MiKlip.__init__(self)

        # The names of the space separated entries in the XML-file <model> tag
        # lines
        self.model_specifiers = ["project",
                                 "institute",
                                 "name",
                                 "mip",
                                 "experiment",
                                 "ensemble",
                                 "realm",
                                 "start_year",
                                 "end_year",
                                 "dir"]

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @return Two strings (input directory and input file)
            This function looks for the input file to use in the
            reformat routines
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        if ('6hr' in msd["mip"]):
            mip = '6hr'
        if ('day' in msd["mip"]):
            mip = 'day'
        if ('mon' in msd["mip"]):
            mip = 'mon'
        if ('yr' in msd["mip"]):
            mip = 'yr'
        if ('clim' in msd["mip"]):
            mip = 'monClim'

        version_path = os.path.join(msd["dir"],
                                    msd["institute"],
                                    msd["name"],
                                    msd["experiment"],
                                    mip,
                                    msd["realm"],
                                    msd["mip"],
                                    msd["ensemble"])
        list_versions = glob.glob(version_path+'/*')
        list_versions.sort()
        version = os.path.basename(list_versions[-1])

        indir = os.path.join(msd["dir"],
                             msd["institute"],
                             msd["name"],
                             msd["experiment"],
                             mip,
                             msd["realm"],
                             msd["mip"],
                             msd["ensemble"],
                             version,
                             variable)
        print(indir)

        infile = "_".join([variable,
                           msd["mip"],
                           msd["name"],
                           msd["experiment"],
                           msd["ensemble"]]) + ".nc"

        if (not os.path.isfile(os.path.join(indir, infile))):
            infile = "_".join([variable,
                               msd["mip"],
                               msd["name"],
                               msd["experiment"],
                               msd["ensemble"]]) + "*.nc"

        return indir, infile

    def get_cf_areafile(self, project_info, model):
        """ @brief Returns the path to the areacello file
                   used for ocean variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the areafile of the ocean grid
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        areafile = 'areacello_fx_'\
                   + msd["name"]\
                   + "_" + msd["experiment"] \
                   + "_r0i0p0.nc"

        # Will sellect the newest version
        version_path = os.path.join(msd["dir"],
                                    msd["institute"],
                                    msd["name"],
                                    msd["experiment"],
                                    'fx',
                                    'ocean',
                                    'fx',
                                    'r0i0p0')
        list_versions = glob.glob(version_path+'/*')
        list_versions.sort()
        version = os.path.basename(list_versions[-1])

        areadir = os.path.join(msd["dir"],
                               msd["name"],
                               msd["institute"],
                               msd["experiment"],
                               'fx',
                               'ocean',
                               'fx',
                               'r0i0p0',
                               version,
                               'areacello')

        return os.path.join(areadir, areafile)


class O3_Cionni(Project):
    """ @brief Class defining the specific characteristics of the CMIP5 project
    """
    def __init__(self):
        Project.__init__(self)

        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project", "name", "ensemble", "experiment", "start_year", "end_year", "dir"]

        """ Define the 'basename'-variable explicitly
        """
        self.basename = "O3_Cionni"

    def get_dict_key(self, model, mip, exp):
        """ @brief Returns a unique key based on the model entries provided.
            @parjam model One of the <model>-tags in the XML namelist file
            @return A string

            This function creates and returns a key used in a number of NCL
            scripts to refer to a specific dataset, see e.g., the variable
            'cn' in 'interface_scripts/read_data.ncl'
        """        

        # msd = model_section_dictionary

        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        dict_key = "_".join(['Ozone_CMIP5_ACC_SPARC',
                             msd['experiment'],
                             msd['start_year'],
                             msd['end_year']])
 
        return dict_key

    def get_figure_file_names(self, project_info, model, mip, exp):
        """ @brief Returns the full path used for intermediate data storage
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param variable_attributes The variable attributes
            @return A string

            This function creates and returns the full path used in a number
            of NCL scripts to write intermediate data sets.
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        return "_".join(['Ozone_CMIP5_ACC_SPARC',
                         msd['name'],
                         msd['experiment'],
                         msd['start_year']]) + "-" + msd['end_year']

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @param variable_attributes The variable attributes
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        indir = msd["dir"]+ 'final_' + msd['experiment']

        variable = self.get_project_variable_name(model, variable)

        infile = '_'.join([variable,
                           msd['experiment'],]) + '.nc'

        if (not os.path.isfile(os.path.join(indir, infile))):
            infile = '_'.join(['Ozone_CMIP5_ACC_SPARC_*',
                               msd['experiment'],
                               field,
                               variable]) + '.nc'

        return indir, infile

    def get_cf_outfile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the output file used in the reformat routines
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @return A string (output path)

            This function specifies the output file to use in the reformat
            routines and in climate.ncl
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        # Overide some model lines with settings from the variable attributes
        msd = self.rewrite_mip_exp(msd, mip, exp)

        outfile = '_'.join(['Ozone_CMIP5_ACC_SPARC',
                            msd['name'],
                            msd['experiment'],
                            field,
                            variable,
                            msd['start_year']]) + '-' + msd['end_year'] + '.nc'

        return outfile

class CCMVal(Project):
    """ @brief Defining the specific characteristics of the CCMVal project
    """
    def __init__(self):
        Project.__init__(self)

        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "case_name",
                                 "ensemble",
                                 "start_year",
                                 "end_year",
                                 "dir"]

    def get_dict_key(self, model, mip, exp):
        """ @brief Returns a unique key based on the model entries provided.
            @param model One of the <model>-tags in the XML namelist file
            @return A string

            This function creates and returns a key used in a number of NCL
            scripts to refer to a specific dataset, see e.g., the variable
            'cn' in 'interface_scripts/read_data.ncl'
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        dict_key = '_'.join([msd['project'],
                             msd['name'],
                             msd['case_name'],
                             msd['ensemble'],
                             msd['start_year'],
                             msd['end_year']])
        return dict_key

    def get_figure_file_names(self, project_info, model, mip, exp):
        """ @brief Returns the full path used for intermediate data storage
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param variable_attributes The variable attributes
            @return A string

            This function creates and returns the full path used in a number
            of NCL scripts to write intermediate data sets.
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        return "_".join([msd['project'],
                         msd['name'],
                         msd['case_name'],
                         msd['ensemble'],
                         msd['start_year']]) + "-" + msd['end_year']

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable Variable (defaults to the variable in project_info)
            @param variable_attributes The variable attributes
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        indir = msd['dir']

        infile = '_'.join([msd['project'],
                           msd['case_name'],
                           msd['name'],
                           msd['ensemble'],
                           field,
                           variable]) + '.nc'
        if (os.path.isfile(os.path.join(indir, infile))):
            return indir, infile

        # Try alternative variable names
        for altvar in find_varname(variable):
            infile = '_'.join([msd['project'],
                               msd['case_name'],
                               msd['name'],
                               msd['ensemble'],
                               field,
                               altvar]) + '.nc'
            if (os.path.isfile(os.path.join(indir, infile))):
                info("  No input files found, trying with the alternative "
                     + "variable name " + altvar,
                     project_info["GLOBAL"]["verbosity"], 1)
                return indir, infile

        raise exceptions.IOError(2, "No input files found in", indir)

    def get_cf_outfile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path and output file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable Variable (defaults to the variable in project_info)
            @return A string (output path)

            This function specifies the output file to use in the
            reformat routines and in climate.ncl
        """
        msd = self.get_model_sections(model)

        outdir = self.get_cf_outpath(project_info, model)

        outfile = '_'.join([msd['project'],
                            msd['case_name'],
                            msd['name'],
                            msd['ensemble'],
                            field,
                            variable,
                            msd['start_year']]) + '-' + msd['end_year'] + '.nc'

        return os.path.join(outdir, outfile)


class CCMVal1(CCMVal):
    """ @brief Defining the specific characteristics of the CCMVal1 project
    """
    def __init__(self):
        CCMVal.__init__(self)

        self.basename = "CCMVal1"


class CCMVal2(CCMVal):
    """ @brief Defining the specific characteristics of the CCMVal2 project
    """
    def __init__(self):
        CCMVal.__init__(self)

        self.basename = "CCMVal2"


class EMAC(CCMVal):
    """ @brief Defining the specific characteristics of the EMAC model output
    """
    def __init__(self):
        CCMVal.__init__(self)

        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "ensemble",
                                 "start_year",
                                 "end_year",
                                 "dir"]

        self.basename = "EMAC"

    def get_dict_key(self, model, mip, exp):
        """ @brief Returns a unique key based on the model entries provided.
            @param model One of the <model>-tags in the XML namelist file
            @return A string

            This function creates and returns a key used in a number of NCL
            scripts to refer to a specific dataset, see e.g., the variable
            'cn' in 'interface_scripts/read_data.ncl'
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        dict_key = '_'.join([msd['project'],
                             msd['name'],
                             msd['ensemble'],
                             msd['start_year'],
                             msd['end_year']])

        return dict_key

    def get_figure_file_names(self, project_info, model, mip, exp):
        """ @brief Returns the full path used for intermediate data storage
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param variable_attributes The variable attributes
            @return A string

            This function creates and returns the full path used in a number
            of NCL scripts to write intermediate data sets.
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        return "_".join([msd['project'],
                         msd['name'],
                         msd['ensemble'],
                         msd['start_year']]) + "-" + msd['end_year']

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable Variable (defaults to the variable in project_info)
            @param variable_attributes The variable attributes
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        indir = msd['dir']

        # reformat_EMAC requires only indir, return empty infile
        return indir, ""

    def get_cf_outfile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path and output file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable Variable (defaults to the variable in project_info)
            @return A string (output path)

            This function specifies the output file to use in the
            reformat routines and in climate.ncl
        """
        msd = self.get_model_sections(model)

        outdir = self.get_cf_outpath(project_info, model)

        outfile = '_'.join([msd['project'],
                            msd['name'],
                            msd['ensemble'],
                            field,
                            variable,
                            msd['start_year']]) + '-' + msd['end_year'] + '.nc'

        return os.path.join(outdir, outfile)


class ECEARTH(CMIP5):
    """ @brief Class defining the specific characteristics of the ECEARTH
               project
    """
    def __init__(self):
        Project.__init__(self)

        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "experiment",
                                 "ensemble",
                                 "start_year",
                                 "end_year",
                                 "dir"]

        self.add_specifier['case_name'] = 'experiment'

        """ Define the 'basename'-variable explicitly
        """
        self.basename = "ECEARTH"

    def get_dict_key(self, model, mip, exp):
        """ @brief Returns a unique key based on the model entries provided.
            @param model One of the <model>-tags in the XML namelist file
            @return A string

            This function creates and returns a key used in a number of NCL
            scripts to refer to a specific dataset, see e.g., the variable
            'cn' in 'interface_scripts/read_data.ncl'
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        dict_key = "_".join([msd['project'],
                             msd['name'],
                             msd['experiment'],
                             msd['start_year'],
                             msd['end_year']])
        return dict_key

    def get_figure_file_names(self, project_info, model, mip, exp):
        """ @brief Returns the full path used for intermediate data storage
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param variable_attributes The variable attributes
            @return A string

            This function creates and returns the full path used in a number
            of NCL scripts to write intermediate data sets.
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        return "_".join([msd['project'],
                         msd['name'],
                         msd['experiment'],
                         msd['start_year']]) + "-" + msd['end_year']

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable (defaults to the variable in
                            project_info)
            @param variable_attributes The variable attributes
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """
        indata_root = self.get_data_root()

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        indir = os.path.join(indata_root,
                             msd["dir"],
                             msd['experiment'])

        # variable and input file pattern, STAGGERGRID is replaced in
        # reformat script
        variable = self.get_project_variable_name(model, variable)
        infile = os.path.join('*',
                              'Output',
                              '_'.join([msd['name'],
                                       'MM',
                                       '*',
                                       '*',
                                       'grid',
                                       '*']) + '.nc')

        return indir, infile

    def get_cf_outfile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the output file used in the reformat routines
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable (defaults to the variable in
                            project_info)
            @return A string (output path)

            This function specifies the output file to use in the reformat
            routines and in climate.ncl
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        # Overide some model lines with settings from the variable attributes
        msd = self.rewrite_mip_exp(msd, mip, exp)

        outfile = '_'.join([msd['project'],
                            msd['experiment'],
                            msd['name'],
                            field,
                            variable,
                            msd['start_year']]) + '-' + msd['end_year'] + '.nc'

        return outfile

    def get_cf_hgridfile(self, project_info, model):
        """ @brief Returns the path to the mesh_hgr file used for ocean
                   variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (mesh_hgr path)

            This function looks for the mesh_hgr file of the ocean grid
        """
        indata_root = self.get_data_root()
        hgrid_path = project_info['AUXILIARIES']['FX_files'].fx_files['nemo_hgrid_file'].get_fullpath()

        return os.path.join(indata_root, hgrid_path)

    def get_cf_zgridfile(self, project_info, model):
        """ @brief Returns the path to the mesh_zgr file used for ocean
                   variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (mesh_zgr path)

            This function looks for the mesh_zgr file of the ocean grid
        """
        indata_root = self.get_data_root()
        zgrid_path = project_info['AUXILIARIES']['FX_files'].fx_files['nemo_zgrid_file'].get_fullpath()

        return os.path.join(indata_root, zgrid_path)

    def get_cf_lsmfile(self, project_info, model, field):
        """ @brief Returns the path to the landmask file used for ocean
                   variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @return A string (landmask path)

            This function looks for the landmask file of the ocean grid
        """
        indata_root = self.get_data_root()
        zgrid_path = project_info['AUXILIARIES']['FX_files'].fx_files['nemo_zgrid_file'].get_fullpath()

        if field == "T3M":
            lsmfile_path = project_info['AUXILIARIES']['FX_files'].fx_files['nemo_lsm3d_file'].get_fullpath()
        else:
            lsmfile_path = project_info['AUXILIARIES']['FX_files'].fx_files['nemo_lsm_file'].get_fullpath()

        return os.path.join(indata_root, lsmfile_path)

    def get_cf_sections(self, model):
        """ @brief Return the sections from the <model> tag needed in reformat
            @param model One of the <model>-tags in the XML namelist file
            @return Six strings

            This function returns the sections from the <model> tag that
            are needed in the 'reformat'-routine.
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        # no support for ensemble in EC-Earth, yet
        return msd['project'],\
            msd['name'],\
            msd['ensemble'],\
            msd['start_year'],\
            msd['end_year'],\
            msd['dir']


class GFDL(Project):
    """ @brief Defining the specific characteristics of the GFDL model output
    """
    def __init__(self):
        Project.__init__(self)

        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "realm",
                                 "ensemble",
                                 "start_year",
                                 "end_year",
                                 "shift_year",
                                 "dir"]

        self.basename = "GFDL"

    def get_dict_key(self, model, mip, exp):
        """ @brief Returns a unique key based on the model entries provided.
            @param model One of the <model>-tags in the XML namelist file
            @return A string

            This function creates and returns a key used in a number of NCL
            scripts to refer to a specific dataset, see e.g., the variable
            'cn' in 'interface_scripts/read_data.ncl'
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        dict_key = '_'.join([msd['project'],
                             msd['name'],
                             msd['realm'],
                             msd['ensemble'],
                             msd['start_year'],
                             msd['end_year'],
                             msd['shift_year']])

        return dict_key

    def get_figure_file_names(self, project_info, model, mip, exp):
        """ @brief Returns the full path used for intermediate data storage
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param variable_attributes The variable attributes
            @return A string

            This function creates and returns the full path used in a number
            of NCL scripts to write intermediate data sets.
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        return "_".join([msd['project'],
                         msd['name'],
                         msd['realm'],
                         msd['ensemble'],
                         msd['start_year']]) + "-" + msd['end_year']

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable Variable (defaults to the variable in project_info)
            @param variable_attributes The variable attributes
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        indir = msd['dir']

        # reformat_GFDL requires only indir, return empty infile
        return indir, ""

    def get_cf_outfile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path and output file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable Variable (defaults to the variable in project_info)
            @return A string (output path)

            This function specifies the output file to use in the
            reformat routines and in climate.ncl
        """
        msd = self.get_model_sections(model)

        outdir = self.get_cf_outpath(project_info, model)

        outfile = '_'.join([msd['project'],
                            msd['name'],
                            msd['realm'],
                            msd['ensemble'],
                            field,
                            variable,
                            msd['start_year']]) + '-' + msd['end_year'] + '.nc'

        return os.path.join(outdir, outfile)


class GO(CCMVal):
    """ @brief Class defining the specific characteristics of the GO project
    """
    def __init__(self):
        CCMVal.__init__(self)

        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "frequency",
                                 "case_name",
                                 "resolution",
                                 "start_year",
                                 "end_year",
                                 "dir",
                                 "gridfile",
                                 "which_reformat"]

        self.add_specifier['experiment'] = 'case_name'
        self.basename = "GO"

    def get_dict_key(self, model, mip, exp):
        """ @brief Returns a unique key based on the model entries provided.
            @param model One of the <model>-tags in the XML namelist file
            @return A string

            This function creates and returns a key used in a number of NCL
            scripts to refer to a specific dataset, see e.g., the variable
            'cn' in 'ncl_code/read_data.ncl'
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        dict_key = "_".join([msd['project'],
                             msd['name'],
                             msd['frequency'],
                             msd['case_name'],
                             msd['resolution'],
                             msd['start_year'],
                             msd['end_year']])
        return dict_key

    def get_cf_sections(self, model):
        """ @brief Return the sections from the <model> tag needed in reformat
            @param model One of the <model>-tags in the XML namelist file
            @return Six strings

            This function returns the sections from the <model> tag that
            are needed in the 'reformat'-routine.
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        return msd['project'],\
            msd['name'],\
            msd['resolution'],\
            msd['start_year'],\
            msd['end_year'],\
            msd['dir']

    def get_figure_file_names(self, project_info, model, mip, exp):
        """ @brief Returns the full path used for intermediate data storage
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string

            This function creates and returns the full path used in a number
            of NCL scripts to write intermediate data sets.
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)
        return "_".join([msd['project'],
                         msd['name'],
                         msd['case_name'],
                         msd['resolution'],
                         msd['start_year']]) + "-" + msd['end_year']

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable Variable (defaults to the variable in project_info)
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        indata_root = self.get_data_root()
        indir = os.path.join(indata_root, msd['dir'])

        info("dir = " + indir, 1, 1)
        which_reformat = msd['which_reformat']
        info("which_reformat = " + which_reformat, 1, 1)

        # reformatGO requires only indir, return empty infile
        if (re.search('use_GO_reformat', which_reformat) is not None):
            return indir, ""

        infile = '_'.join([msd['project'],
                           msd['case_name'],
                           msd['name'],
                           msd['resolution'],
                           field,
                           variable]) + '.nc'

        info("file = " + infile, 1, 1)

        if (not os.path.isfile(os.path.join(indir, infile))):
            infile = '_'.join([msd['project'],
                               msd['case_name'],
                               msd['name'],
                               msd['resolution'],
                               field,
                               variable]) + '*.nc'

        if (len(glob.glob(os.path.join(indir, infile))) == 0):
            raise exceptions.IOError(2, "No input files found in", indir)

        return indir, infile

    def get_cf_outfile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path + output file used in the reformat routines
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable Variable (defaults to the variable in project_info)
            @return A string (output path)

            This function specifies the output file to use in the reformat
            routines and in climate.ncl
        """
        msd = self.get_model_sections(model)
     # Overide some model lines with settings from the variable attributes
        msd = self.rewrite_mip_exp(msd, mip, exp)

        outdir = self.get_cf_outpath(project_info,  model)

        outfile = '_'.join([msd['project'],
                            msd['case_name'],
                            msd['name'],
                            msd['resolution'],
                            field,
                            variable,
                            msd['start_year']]) + '-' + msd['end_year'] + '.nc'

        return os.path.join(outdir, outfile)

    def get_cf_areafile(self, project_info, model):
        """ @brief Returns the path to the areacello file
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A strings (areafile path)
            This function looks for the areafile of the ocean grid
        """
        return self.get_model_subsection(model, "gridfile")


class GO_gridfile(GO):
    def __init__(self):
        GO.__init__(self)
        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project",
                                 "name",
                                 "frequency",
                                 "case_name",
                                 "resolution",
                                 "start_year",
                                 "end_year",
                                 "dir",
                                 "fx_file_ID",
                                 "which_reformat"]

 #       self.add_specifier['case_name'] = 'case_name'

    def get_cf_areafile(self, project_info, model):
        """ @brief Returns the path to the areacello file
                   used for ocean variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the areafile of the ocean grid
        """
        return self.get_cf_fx_file(project_info, model)

    def get_cf_fx_file(self, project_info, model):
        """ @brief Returns the path to the areacello file
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A strings (areafile path)

            This function looks for the areafile of the ocean grid
        """
        fx = self.get_fx_file(project_info, model)
        return os.path.abspath(fx)


class JSBACH(Project):
    def __init__(self):
        #~ super(JSBACH, self).__init__()  # does not work, as Project itself
                                           # is an old-style class (http://stackoverflow.com/questions/9698614/super-raises-typeerror-must-be-type-not-classobj-for-new-style-class)
        Project.__init__(self)
        self.model_specifiers = ["project",
                                 "experiment",
                                 "dir",
                                 "which_reformat"]
        #~ self.model_specifiers = ["project", "name", "case_name", "ensemble",
                                 #~ "start_year", "end_year", "dir"]
        self.basename = 'JSBACH'

    def get_model_name(self, model):
        """
        return model name
        overwrites general Master class function
        """
        sep = '_'
        res = sep.join([self.get_model_subsection(model, "project"),
                        self.get_model_subsection(model, "experiment")])
        return res

    def get_cf_infile(self,
                      project_info,
                      model,
                      field,
                      variable,
                      mip,
                      exp,
                      fmt='sz'):
        """ @brief Returns the path to the input file used in reformat
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @param variable_attributes The variable attributes
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        #~ msd = self.rewrite_mip_exp(msd, mip, exp)

        indir = msd["dir"]
        if indir[-1] != os.sep:
            indir += os.sep
        indir += msd["experiment"]
        if indir[-1] != os.sep:
            indir += os.sep

        infile = indir + '*.' + fmt

        return indir, infile

    def get_cf_outfile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the output file used in the reformat routines
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @return A string (output path)

            This function specifies the output file to use in the reformat
            routines and in climate.ncl
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        # Overide some model lines with settings from the variable attributes
        #~ msd = self.rewrite_mip_exp(msd, mip, exp)

        #~ outfile = '_'.join([msd['project'],
                            #~ msd['mip'],
                            #~ msd['experiment'],
                            #~ msd['name'],
                            #~ msd['ensemble'],
                            #~ field,
                            #~ variable,
                            #~ msd['start_year']]) + '-' + msd['end_year'] + '.nc'
        outfile = 'dummy.nc'

        return outfile

    def get_cf_sections(self, model):
        """
        overwrites same function of master class
        """
        #~ project, name, ensemble, start_year, end_year, dir\
            #~ = currProject.get_cf_sections(model)

        msd = self.get_model_sections(model)

        project = msd["project"]
        name = None
        ensemble = None
        start_year = None
        end_year = None
        directory = msd["dir"]

        return project, name, ensemble, start_year, end_year, directory

    def get_figure_file_names(self, project_info, model, mip, exp):
        return 'dummy_figure_filename'

    def get_dict_key(self, model, mip, exp):
        return "dummy_key"

###added by Cionni Irene,ENEA,Italy irene.cionni@enea.it########################

class CCMI(Project):
    """ @brief Class defining the specific characteristics of the CMIP5 project
    """
    def __init__(self):
        Project.__init__(self)

        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project", "name", "mip", "experiment",
                                 "ensemble", "start_year", "end_year", "dir"]

        self.add_specifier['case_name'] = 'experiment'

        """ Define the 'basename'-variable explicitly
        """
        self.basename = "CCMI"

    def get_dict_key(self, model, mip, exp):
        """ @brief Returns a unique key based on the model entries provided.
            @param model One of the <model>-tags in the XML namelist file
            @return A string

            This function creates and returns a key used in a number of NCL
            scripts to refer to a specific dataset, see e.g., the variable
            'cn' in 'interface_scripts/read_data.ncl'
        """
    # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        dict_key = "_".join([msd['project'],
                             msd['name'],
                             msd['mip'],
                             msd['experiment'],
                             msd['ensemble'],
                             msd['start_year'],
                             msd['end_year']])
        return dict_key

    def get_figure_file_names(self, project_info, model, mip, exp):
        """ @brief Returns the full path used for intermediate data storage
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param variable_attributes The variable attributes
            @return A string

            This function creates and returns the full path used in a number
            of NCL scripts to write intermediate data sets.
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        return "_".join([msd['project'],
                         msd['name'],
                         msd['mip'],
                         msd['experiment'],
                         msd['ensemble'],
                         msd['start_year']]) + "-" + msd['end_year']

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @param variable_attributes The variable attributes
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """
        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        indir = msd["dir"]

        variable = self.get_project_variable_name(model, variable)

        infile = '_'.join([variable,
                           msd['mip'],
                           msd['name'],
                           msd['experiment'],
                           msd['ensemble']]) + '.nc'

        if (not os.path.isfile(os.path.join(indir, infile))):
            infile = '_'.join([variable,
                               msd['mip'],
                               msd['name'],
                               msd['experiment'],
                               msd['ensemble']]) + '*.nc'

        return indir, infile

    def get_cf_areafile(self, project_info, model):
        """ @brief Returns the path to the areacello file
                   used for ocean variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the areafile of the ocean grid
        """
        msd = self.get_model_sections(model)

        areadir = msd["dir"]

        areafile = 'areacello_fx_' + msd["name"] + "_" + msd["experiment"] + \
                   "_r0i0p0.nc"

        return os.path.join(areadir, areafile)

    def get_cf_outfile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the output file used in the reformat routines
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @return A string (output path)

            This function specifies the output file to use in the reformat
            routines and in climate.ncl
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        # Overide some model lines with settings from the variable attributes
        msd = self.rewrite_mip_exp(msd, mip, exp)

        outfile = '_'.join([msd['project'],
                            msd['mip'],
                            msd['experiment'],
                            msd['name'],
                            msd['ensemble'],
                            field,
                            variable,
                            msd['start_year']]) + '-' + msd['end_year'] + '.nc'

        return outfile
#############################
class CCMI_gridfile(CCMI):
    def __init__(self):
        CCMI.__init__(self)
        ## The names of the space separated entries in the XML-file <model> tag
        ## lines
        self.model_specifiers = ["project", "name", "mip", "experiment",
                                 "ensemble", "start_year", "end_year", "dir",
                                 "gridfile"]

        self.add_specifier['case_name'] = 'experiment'

    def get_cf_areafile(self, project_info, model):
        """ @brief Returns the path to the areacello file
                   used for ocean variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the areafile of the ocean grid
        """
        return self.get_model_subsection(model, "gridfile")

############
class CCMI_DLR(CCMI):
    """ @brief Follows the ETH directory structure for input data
    """
    def __init__(self):
        CCMI.__init__(self)

    def get_cf_infile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the path to the input file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @param variable_attributes The variable attributes
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        msd = self.rewrite_mip_exp(msd, mip, exp)

        indir = os.path.join(msd["dir"],
                             msd["experiment"],
                             msd["mip"],
                             variable,
                             msd["name"],
                             msd["ensemble"])

        infile = "_".join([variable,
                           msd["mip"],
                           msd["name"],
                           msd["experiment"],
                           msd["ensemble"]]) + ".nc"

        if (not os.path.isfile(os.path.join(indir, infile))):
            infile = "_".join([variable,
                               msd["mip"],
                               msd["name"],
                               msd["experiment"],
                               msd["ensemble"]]) + "*.nc"

        return indir, infile

    def get_cf_areafile(self, project_info, model):
        """ @brief Returns the path to the areacello file
                   used for ocean variables
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (areafile path)

            This function looks for the areafile of the ocean grid
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        areadir = os.path.join(msd["dir"],
                               msd["experiment"],
                               'fx',
                               'areacello',
                               msd["name"],
                               'r0i0p0')

        areafile = 'areacello_fx_' + msd["name"] + "_" + msd["experiment"] \
            + "_r0i0p0.nc"

        return os.path.join(areadir, areafile)
########end add by Cionni Irene,ENEA,Italy irene.cionni@enea.it#################

class OneFile(Project):
    """simple project class intended to pass the filename of one file directly
    to the diagnostics routine"""
    def __init__(self):
        Project.__init__(self)
        self.model_specifiers = ["project", "experiment", "dir",
                                 "infile", "which_reformat"]
        self.basename = 'OneFile'

    def get_model_name(self, model):
        """
        return model name
        overwrites general Master class function
        """
        sep = '_'
        res = self.get_model_subsection(model, "experiment")
        return res

    def get_cf_infile(self, project_info, model, field, variable, mip, exp,
                      fmt='sz'):
        """ @brief Returns the path to the input file used in reformat
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @param variable_attributes The variable attributes
            @return Two strings (input directory and input file)

            This function looks for the input file to use in the
            reformat routines
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)
        #~ msd = self.rewrite_mip_exp(msd, mip, exp)

        indir = msd["dir"]
        if indir[-1] != os.sep:
            indir += os.sep

        infile = msd['infile']

        return indir, infile

    def get_cf_outfile(self, project_info, model, field, variable, mip, exp):
        """ @brief Returns the output file used in the reformat routines
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable
                            (defaults to the variable in project_info)
            @return A string (output path)

            This function specifies the output file to use in the reformat
            routines and in climate.ncl
        """

        # msd = model_section_dictionary
        msd = self.get_model_sections(model)

        # Overide some model lines with settings from the variable attributes
        #~ msd = self.rewrite_mip_exp(msd, mip, exp)

        #~ outfile = '_'.join([msd['project'],
                            #~ msd['mip'],
                            #~ msd['experiment'],
                            #~ msd['name'],
                            #~ msd['ensemble'],
                            #~ field,
                            #~ variable,
                            #~ msd['start_year']]) + '-' + msd['end_year'] +
                            #'.nc'
        #outfile = 'dummy.nc'
        indir = msd["dir"]
        if indir[-1] != os.sep:
            indir += os.sep

        infile = msd['infile']
        return indir + infile

    def get_cf_sections(self, model):
        """
        overwrites same function of master class
        """
        #~ project, name, ensemble, start_year, end_year, dir\
            #~ = currProject.get_cf_sections(model)

        msd = self.get_model_sections(model)

        project = msd["project"]
        name = None
        ensemble = None
        start_year = None
        end_year = None
        directory = os.path.dirname(msd["infile"])

        return project, name, ensemble, start_year, end_year, directory

    def get_figure_file_names(self, project_info, model, mip, exp):
        return 'dummy_figure_filename'

    def get_dict_key(self, model, mip, exp):
        return "dummy_key"


def find_varname(var):
    """
    @brief Read and return alternative names for the given var
    @param var Variable name according to the CMOR standard
    """

    fn = open("./reformat_scripts/recognized_vars.dat", "r")

    lines = []
    for ln in fn:
        lines.append(ln.strip())

    ii = 0
    while ii < len(lines):
        if ('std_name' in lines[ii] and var in lines[ii] and not '#' in lines[ii]):
            altern = lines[ii + 1]
            altern = altern.replace(" ", "")
            altern = altern.split("=")[1].split(",")
            return altern
        ii += 1

    return []


def add_model(project_info, models_to_add):
    """
    @brief Add a number of models to project_info['MODELS']
    @param project_info Current namelist in dictionary format
    @param models_to_add List of 'Model'-instances to add
    """
    for model in models_to_add:
        project_info['MODELS'].append(model)


def remove_diag_specific_models(project_models):
    """
    @brief Remove a number of models from project_info['MODELS']
    @param project_models Current models from the project_info['MODELS']
    """
    project_models = [model for model in project_models
                      if not model.is_diag_specific()]
    return project_models


def write_data_interface(executable, project_info):
    """ @brief Write Python data structures to target script format interface
        @param executable String pointing to the script/binary to execute
        @param project_info Current namelist in dictionary format

        Data structures in Python are rewritten to the interface folder in
        a format appropriate for the target script/binary
    """
    suffix = os.path.splitext(executable)[1][1:]
    currInterface = vars(data_interface)[suffix.title() + '_data_interface'](project_info)
    currInterface.write_data_to_interface()


def run_executable(string_to_execute,
                   project_info,
                   verbosity,
                   exit_on_warning,
                   launcher_arguments=None,write_di=True):
    """ @brief Executes script/binary
        @param executable String pointing to the script/binary to execute
        @param project_info Current namelist in dictionary format
        @param verbosity The requested verbosity level
        @param exit_on_warning Boolean defining whether the wrapper should
                               crash on warnings

        Check the type of script/binary from the executable string suffix and
        execute the script/binary properly.
    """
    
    if write_di:
	write_data_interface(string_to_execute, project_info)

    suffix = os.path.splitext(string_to_execute)[1][1:]
    currLauncher = vars(launchers)[suffix + '_launcher']()
    if launcher_arguments is not None:
        currLauncher.arguments = launcher_arguments
    currLauncher.execute(string_to_execute,
                         project_info,
                         verbosity,
                         exit_on_warning)
