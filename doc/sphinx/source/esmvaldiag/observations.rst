:: _CMORobs:

************************************************************
Contributing a CMORizing script for an observational dataset
************************************************************

ESMValTool v2 is able to work with observational data sets, however, since many features of this version are based on the Iris library 
(https://scitools.org.uk/iris/docs/latest/) that uses the CF conventions (http://cfconventions.org/), the format of the data sets has to follow very strict CMOR rules. If the original observational data set is not in the right CMOR format, it has to be reformatted before the ESMValTool preprocessor can work with it. In the following, steps are explained that are necessary to prepare an observational data set for the use in ESMValTool v2. 

1) Check if your variable is CMOR standard

Some variables are defined as CMOR standard. You can find several files with CMOR tables in the folder “/esmvaltool/cmor/tables/cmip5/Tables/”, differentiated according to the “MIP” they belong to. If you find the variable in one of these files, you do not have to do anything further. Please note: if your variable is somehow more specific than the entries in the tables (e.g. your variable is called “bihemispherical directional albedo” and in the table only “albedo” is available), please write your own custom variable information
as outlined below. If you don’t find your variable here, you have to provide the variable information in an additional file that will have to be placed in the folder “/esmvaltool/cmor/tables/custom/”. 

The file that you have to create for a custom variable needs to follow these guidelines:

- Provide information about the SOURCE (for all observations and reanalysis it was decided by the ESMValTool development team that this is set to “CMIP5”; it basically describes the framework in which the individual variables are associated with realms; it might not be necessary to provide this information here, but before we know for sure that we can leave it out, please just provide it)

- Provide the “variable_entry”

- Provide the “modeling_realm”

- Provide the variable attributes, BUT IMPORTANT: do not provide a standard name! The custom variables will only be read by Iris if this information is not provided! Necessary variable attributes are: “units”, “cell_methods”, “cell_measures”, “long_name”, “comment”. 

- Provide some additional variable attributes. Necessary additional variable attributes are: “dimensions”, “out_name”, “type”. There are also additional variable attributes that can be defined here (see observation reformatting scripts that are available already).

Use the file structure as given in the files that already exist, and save the file with 
the ending “.dat”.

2) Configuration file

It is important to have your configuration file defined correctly, to be sure that your data files are found by the ESMValTool, and that the reformatted data files are written in the correct folders. Here is an example of a configuration file, that has to be filled out with the paths information of your specific machine or the VM if you will run the ESMValTool there.

3) Store your dataset in the right place

For the ESMValTool to find the data set that it needs, it needs a specific folder structure. You have your path that points towards the observations. That path is defined in your configuration file (see Section 2) under “RAWOBS”. From there you need the folders "Tier1", "Tier2" and "Tier3". In the folder of your choice (for C3S_511 datasets that would be "Tier3"), you need a folder with the name of the dataset (for the example of the xch4 dataset that folder would be named "CDS-XCH4"), and in that folder you have the netcdf-file that you downloaded from the CDS.

4) Reformat your data set to meet CMOR standard - script

To make sure that your dataset is in the correct format for Iris to work with it, you most likely will have to process it with a reformat script. In the rare, but possible, case that your dataset is in the correct format already, it is advisable to run it through a reformat script anyway, to be sure the file name is correctly defined for the ESMValTool to find the file.

Find here an example of a reformat script, written for the ESACCI XCH4 dataset that is available on the CDS (/esmvaltool/utils/cmorizers/obs/cmorize_obs_CDS-XCH4.ncl).

The first part of the script collects all the information about the dataset that are necessary to write the filename correctly and to understand which variable is of interest here. Please make sure to pay special attention to the following:

- DIAG_SCRIPT: fill in the name of the current reformatting script

- VAR: here you can list all the different variables you want to store in the reformatted output file. In this example only the variable “xch4” is listed. If more than one variable is supposed to be reformatted, you define VAR as follows: (/"xch4", “xch4stddev”, “xch4_num”/). Note: each variable needs to be stored in a separate file for the ESMValTool to be able to work with the files, e.g. “xch4” and “xch4stddev” would not be stored in the same file, but in two separate files.

- NAME: these are the names of the variables you want to extract out of the original data file. The names do not need to be in CMOR standard therefore there is the distinction between “VAR” and “NAME”

- MIP: this is the “mip” in which the variable is defined (or would be defined if it is not a custom variable or a derived variable) and which describes the realm and the temporal resolution of the dataset; “Amon” from the example stands for an atmospheric variable (“A”) in a monthly resolution (“mon”). 

- Note: The description of the MIP is not necessarily structured the same way as described above. The available choices for MIP are: 3hr, 6hrLev, 6hrPlev, aero, Amon, cf3hr, cfDay, cfMon, cfOff, cfSites, day, fx, grids, LImon, Lmon, Oclim, OImon, Omon, Oyr. See for more details on these different MIPs: “/esmvaltool/cmor/tables/cmip5/Tables/”.

- FREQ: describes the temporal resolution of the dataset

- CMOR_TABLE: provides the link to the CMOR table in which the variable is defined. If the CMOR table is a custom table (like it is here in the example) you need to provide the path and the name of the file in which the definition is stored (here: "/cmor/tables/custom/CMOR_xch4.dat"). The more basic path information is pulled out of the configuration file (see section 2) that you will have to provide to run the reformat script. If your variable is not a custom variable, you would provide here the path to the folder to the table where the variable would be available (see for example “/esmvaltool/utils/cmorizers/obs/cmorize_obs_ERA-Interim.ncl”).

- Note: the fields VAR, NAME, MIP and FREQ all ask for one or more entries. If more than one entry is provided, make sure that the order of the entries is the same for all four fields! (for example, that the first entry in all four fields describe the variable “xch4” that you would like to extract)

- Note: the functions in the script that are written in red are all ncl specific and are available through the loading of the program “interface.ncl”. There are similar functions available for python scripts.

In the second part of the script each variable defined in VAR is separately extracted from the original data file and processed. Most parts of the code are commented, and therefore it should be easy to understand roughly what is happening. This example coded in NCL, since many reformatting scripts that are available so far are written in NCL, and adapting existing code and using existing libraries is easier than writing something totally new. However, in theory this script could also be written in Python but without the help of the Iris package and rather based on netCDF4 or similar packages (this is why we need the cmorizing scripts). There is at least one python-based reformatting script available already (/esmvaltool/utils/cmorizers/obs/cmorize_obs_WOA.py) that can be used as guideline in case you would like to write your reformatting script in python.

For the second part of the program, the following points are important to keep in mind:

- fname: it is the combination of the input path that is defined in the configuration file (see Section 2) that has to be defined to run the reformatting script, and the name of the file with the “raw” data

- “output = f->xch4”: In this line it is hardcoded that the variable with the name “xch4” is processed. If you have defined more than one variable, this statement has to be adjusted, so that the correct variable name is used with each loop of the program.

- “format_coords”: this call is a routine that is available for NCL code already and which takes care of reformatting the coordinates of the current variable if necessary (e.g. longitudes ranging from -180 to 180 degrees instead of 0 to 360 degrees).

- “fout”: the filepath and filename of the output file are set here. The path is taken from the configuration file (see Section 2) that is necessary to run the reformatting script, and the filename is put together from the information given in the first part of the script, following the rules for filenames so that the ESMValTool can read in the files.

The script as it is detailed here would only be able to reformat some minor problems with the coordinates (e.g. latitudes in the wrong order, longitudes in the wrong order, etc.). Everything else will have to be added to the script for it to deal with it. There are many reformat scripts available by now in the folder “/esmvaltool/utils/cmorizers/obs/” where solutions to all kinds of observational problems are provided. Most of these reformat scripts are written in NCL, but there are also a few examples for Python-based reformatting scripts. 

How much reformatting an observational data set needs is strongly dependent on the original netCDF file and how close the original formatting already is to the strict CMOR standard.

5) Run the reformatting script

In order to actually run the reformatting script, you have to use the following statement:

::

    cmorize_obs -c name_of_your_configuration_file -o name_of_your_dataset

This call only works, of course, if you are already in the folder in which also the configuration file “name_of_your_configuration_file” is stored.

Note that the output path given in the configuration file is the path where your reformatted dataset will be stored. The ESMValTool will create a folder with the correct tier information (see Section 2) if that tier folder is not already available, and then a folder named after the data set. In this folder the reformatted data set will be stored as a netCDF file.

Your run was successful if a netCDF file was produced in your output directory, and if at some point the output on your screen shows an info line similar to this and no additional error message is shown (note: this is the example for the example data set “xch4”):

::

    INFO    Processing xch4 (Amon)

6) Naming convention of the observational data files

For the ESMValTool to be able to read the observations from the netCDF file, the file name needs a very specific structure and order of information parts (very similar to the naming convention for observations in version 1). The file name will be automatically correctly created if a reformat script has been used to create the netCDF file. If the file with the observations is already in the correct CMOR format, it is possible to use the data set without using a reformat script. In this case it is important to follow the guidelines below on how the filename should be structured.

The file for the CDS-XCH4 observations in the correct format is named as follows:

::

    OBS_CDS-XCH4_sat_L3_Amon_xch4_200301-201612.nc

The different parts of the name are explained in more detail here:

- OBS: describes what kind of data can be expected in the file, in this case “observations”

- CDS-XCH4: that is the name of the dataset. It has been named this way for illustration purposes (so that everybody understands it is the xch4 dataset downloaded from the CDS), but a better name would indeed be “ESACCI-XCH4” since it is a ESA-CCI dataset

- sat: describes the source of the data, here we are looking at satellite data (therefore “sat”), could also be “reanaly” for reanalyses

- L3: describes the version of the dataset

- Amon: is the information in which “mip” the variable is to be expected, and what kind of temporal resolution it has; here we expect “xch4” to be part of the atmosphere (“A”) and we have the dataset in a monthly resolution (“mon”)

- xch4: Is the name of the variable. Each observational data file is supposed to only include one variable per file.

- 200301-201612: Is the period the dataset spans with “200301” being the start year and month, and “201612” being the end year and month

7) Running a preprocessing test recipe

To verify that the reformatted data file is indeed correctly formatted (following the strict CMOR standard Iris needs), it is good to run a preprocessing test recipe, that does not include any diagnostic, but only reads in the data file and has it processed in the preprocessor.

Our branch “C3S_511_V2” includes such a recipe. It is called “recipe_preprocessor_test.yml” and it can be found in the folder “/esmvaltool/recipes/examples/”. It looks like this:


The following parts are variable specific and will have to be adjusted to use this recipe as a test recipe:

- The information about the data set under the header “datasets” needs to be adjusted. All fields that are given in this line need to be filled out, since these are used to put together the filename with the data that the ESMValTool will look for. The order of the different pieces of information is not important. They will be sorted how they need to be according to their “identifier” before the “:”.

- The name of the variable (in this example it is “xch4”) in the line below “variables” in the diagnostics block needs to be adjusted according to the variable that is supposed to be processed. 

- It is important to add the information about the “mip” under the variables information (here: “Amon”). This field is necessary now for version2 since it is included in the filename (see Section 6). The field “field” might not be necessary anymore, or at least will soon not be necessary anymore. For now it is safer to just provide the information. If it is not needed, it does not hurt to have it there anyway (until we know for sure if it is needed or not).

If the recipe is adjusted as outlined above, run it with the following call:

::

    esmvaltool -c name_of_your_configuration_file recipes/examples/recipe_preprocessor_test.yml

This call only works, of course, if you are already in the folder in which also the configuration file “name_of_your_configuration_file” is stored. If your reformatted data set has the correct format, the ESMValTool will read the data and run it through the preprocessor, basically confirming that Iris can work with the data field. The ESMValTool will then print the statement “Run was successful”. If something is still wrong with the format, the ESMValTool will crash.



 








