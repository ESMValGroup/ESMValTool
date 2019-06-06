:: _CMORobs:

************************************************************
Contributing a CMORizing script for an observational dataset
************************************************************

ESMValTool v2 is able to work with observational data sets, however, since many features of 
this version are based on the Iris library 
(https://scitools.org.uk/iris/docs/latest/) that uses the CF conventions (http://cfconventions.org/),
the format of the data sets has to follow very strict CMOR rules. If the original observational 
data set is not in the right CMOR format, it has to be reformatted before the ESMValTool 
preprocessor can work with it. In the following, steps are explained that are necessary 
to prepare an observational data set for the use in ESMValTool v2. 

1. Check if your variable is CMOR standard

Some variables are defined as CMOR standard. You
can find several files with CMOR tables in the folder “/esmvaltool/cmor/tables/cmip5/Tables/”,
differentiated according to the “MIP” they belong to. If you find the variable
in one of these files, you do not have to do anything further. Please note: if
your variable is somehow more specific than the entries in the tables (e.g.
your variable is called “bihemispherical directional albedo” and in the table
only “albedo” is available), please write your own custom variable information
as outlined below. If you don’t find your variable here, you have to provide
the variable information in an additional file that will have to be placed in
the folder “/esmvaltool/cmor/tables/custom/”. 


