.. _data_processing:

Guidelines for data processing
******************************

Observations
============

**Basic steps**

#. Start from the original, if possible referenced data product (avoid using unpublished data pre-processed by a third party).
#. Document how to obtain the data (e.g., link to the data source, download instructions, info on contact persons, references, etc.) in as much detail as necessary for another person to successfully obtain the data.
#. Document any pre-processing steps applied to the original data (e.g., conversion of units, interpolation to pressure levels, arithmetic calculations, conversion of file format, etc.). If possible, use only commonly available tools such as, for instance, NCAR Command Language (NCL), Climate Data Operators (cdo), netCDF Operators (NCO), Python, etc.
#. Save the data in netCDF format following the CMOR/CMIP5 or CMOR/CMIP6 conventions (see links below for variable names, units, grid definition, etc.). N.B.: The processed observational data set (step 3) does not follow the CMIP5-DRS but uses a different DRS. When naming the output observational data set file in step 3 above, follow the naming conventions in any of the existing reformat scripts in the folder "reformat_scripts/obs/".
#. Write/provide a script (e.g., shell script, NCL program, Python code, etc.) handling all pre-processing and conversion steps, and including all the necessary information to retrieve the data from the original source. Many examples of such scripts are available as part of the ESMValTool distribution and can be found in the ESMValTool directory reformat_scripts/obs/.

Model data
==========

* If possible, provide model data in netCDF format following the CMOR conventions (see links below, Section :numref:`data_links`).
* In some cases, it might be worth creating reformat routines that read the native model output and convert the data on the fly within the ESMValTool to CMOR compliant files. Examples for such reformatting routines for the models EC-Earth, EMAC, and GFDL are available as part of the ESMValTool distribution (ESMValTool directory reformat_scripts/<model>/.)

.. _data_links:

Links
=====

* CMOR user’s guide: http://www-pcmdi.llnl.gov/software/cmor/cmor_users_guide.pdf
* ESGF CMOR Data Converter: https://cds.nccs.nasa.gov/tools-services/esgf/obs4mips/
* CMIP5 output data requirement (variables and units): http://cmip-pcmdi.llnl.gov/cmip5/docs/standard_output.pdf
* CMIP6 preliminary output data requirement (variables and units): http://clipc-services.ceda.ac.uk/dreq/mipVars.html

Help
====

Support from the ESMValTool Core Development Team (Section :numref:`core_dev_team`) for observational data is available provided that you followed the basic steps outlined above. For model data, help from the core development team is available if you can provide us with:

* netCDF files and a complete list of all details on the model output including e.g., a brief description of the quantity
* native variable name
* units
* vertical grid definition
* other relevant information on the data set

