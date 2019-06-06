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






