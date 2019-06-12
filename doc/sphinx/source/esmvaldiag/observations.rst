************************************************************
Contributing a CMORizing script for an observational dataset
************************************************************

ESMValTool is able to process observational data sets, however, since many
features of this version are based on the `Iris library 
<https://scitools.org.uk/iris/docs/latest/>`_ that uses the `CF conventions
<http://cfconventions.org/>`_, the format of the data sets has to follow very
strict CMOR rules. If the original observational data set is not in the right
CMOR format, it has to be CMORized before the ESMValTool preprocessor can
work with it. The following steps are necessary to prepare an observational
data set for the use in ESMValTool.

| `1. Check if your variable is CMOR standard`_
| `2. Edit your configuration file`_
| `3. Store your dataset in the right place`_
| `4. Create a cmorizer for the dataset`_
| `5. Run the cmorizing script`_
| `6. Naming convention of the observational data files`_
| `7. Test the cmorized dataset`_


1. Check if your variable is CMOR standard
==========================================

Most variables are defined in the CMIP data request and can be found in the
CMOR tables in the folder ``/esmvalcore/cmor/tables/cmip5/Tables/``,
differentiated according to the ``MIP`` they belong to. If you find the
variable in one of these tables, you can proceed to the next point.

If your variable is not available in the standard CMOR tables,
you need to write a custom CMOR table for the variable
as outlined below and add it to ``/esmvalcore/cmor/tables/custom/``.

To create a new custom CMOR table you need to follow these
guidelines:

- Provide information about the SOURCE (for all observations and reanalysis it
  was decided by the ESMValTool development team that this is set to ``CMIP5``;
  it basically describes the framework in which the individual variables are
  associated with realms; it might not be necessary to provide this information
  here, but before we know for sure that we can leave it out, please just
  provide it);
- Provide the ``variable_entry``;
- Provide the ``modeling_realm``;
- Provide the variable attributes, but leave `standard_name` blank.
Necessary variable attributes are: ``units``, ``cell_methods``,
  ``cell_measures``, ``long_name``, ``comment``.  
- Provide some additional variable attributes. Necessary additional variable
  attributes are: ``dimensions``, ``out_name``, ``type``. There are also
  additional variable attributes that can be defined here (see the already
  available cmorizers). 

It is recommended to use an existing custom table as a template, to edit the content and save it as
``<short_name>.dat``.

2. Edit your configuration file
=====================

It is important to have your configuration file defined correctly, to be sure
that your data files are found by the ESMValTool, and that the cmorized data
files are written in the correct folders. See Section **link to be added** for
an example of a configuration file, that has to be filled out with the paths
information of your specific machine.

3. Store your dataset in the right place
========================================

For the ESMValTool to find the data set that it needs, a specific folder
structure is defined. You have your path that points towards the
observations. That path is defined in your configuration file (see Section 2)
under ``RAWOBS``. From there you need the subdirectories ``Tier1``, ``Tier2``
and ``Tier3``. In the folder of your choice, you need a folder with the name of
the dataset, and in that folder you put the data files with the raw
observations (as you can download them from the observations data source
webpage).

4. Reformat your data set to meet CMOR standard - script
========================================================

To make sure that your dataset is in the correct format for Iris to work with
it, you most likely will have to process it with a cmorizing script. In the
rare, but possible, case that your dataset is in the correct format already, it
is advisable to run it through a cmorizing script anyway, to be sure the file
name is correctly defined for the ESMValTool to find the file. 

Find here an example of a cmorizing script, written for the ``ESACCI XCH4``
dataset that is available on the Copernicus Climate Data Store: `cmorize_obs_CDS-XCH4.ncl
<https://github.com/ESMValGroup/ESMValTool/blob/version2_development/esmvaltool/utils/cmorizers/obs/cmorize_obs_CDS-XCH4.ncl>`_.

The first part of the script collects all the information about the dataset
that are necessary to write the filename correctly and to understand which
variable is of interest here. Please make sure to pay special attention to the
following: 

- DIAG_SCRIPT: fill in the name of the current cmorizing script;
- VAR: here you can list all the different variables you want to store in the
  cmorized output file. In this example only the variable ``xch4`` is
  listed. If more than one variable is supposed to be cmorized, you define
  VAR as an array, e.g. ``(/"xch4", "xch4stddev", "xch4_num"/)``. **Note:**
  each variable needs to be saved in a separate NetCDF file for the ESMValTool
  to be able to work with the files, e.g. ``xch4`` and ``xch4stddev`` should
  not be stored in the same file, but in two separate files;
- NAME: these are the names of the variables you want to extract out of the
  original data file. The names do not need to be in CMOR standard therefore
  there is the distinction between ``VAR`` and ``NAME``;
- MIP: this is the ``mip`` in which the variable is defined (or would be
  defined if it is not a custom variable or a derived variable) and which
  describes the realm and the temporal resolution of the dataset; ``Amon`` from
  the example stands for an atmospheric variable (``A``) in a monthly
  resolution (``mon``).  **Note:** The description of the MIP is not
  necessarily structured the same  way as described above. The available
  choices for MIP are: 3hr, 6hrLev, 6hrPlev, aero, Amon, cf3hr, cfDay, cfMon,
  cfOff, cfSites, day, fx, grids, LImon, Lmon, Oclim, OImon, Omon, Oyr. See for
  more details on these different MIPs see the 
  `CMOR tables <https://github.com/ESMValGroup/ESMValCore/development/esmvalcore/cmor/tables/cmip5/Tables/>`_;
- FREQ: describes the temporal resolution of the dataset;
- CMOR_TABLE: provides the link to the CMOR table in which the variable is
  defined. If the CMOR table is a custom table (like it is here in the example)
  you need to provide the path and the name of the file in which the definition
  is stored (here: ``/cmor/tables/custom/CMOR_xch4.dat``). The more basic path
  information is pulled out of the configuration file (see section 2) that you
  will have to provide to run the cmorizing script. If your variable is not a
  custom variable, you would provide here the path to the folder to the table
  where the variable is available (see for example `cmorize_obs_ERA-Intermim.ncl
  <https://github.com/ESMValGroup/ESMValTool/blob/version2_development/esmvaltool/utils/cmorizers/obs/cmorize_obs_ERA-Interim.ncl>`_;
- **Note:** the fields ``VAR``, ``NAME``, ``MIP`` and ``FREQ`` all ask for one
  or more entries. If more than one entry is provided, make sure that the order
  of the entries is the same for all four fields! (for example, that the first
  entry in all four fields describe the variable ``xch4`` that you would like
  to extract);
- **Note:** some functions in the script are NCL-specific and are available
  through the loading of the script ``interface.ncl``. There are similar
  functions available for python scripts.

In the second part of the script each variable defined in ``VAR`` is separately
extracted from the original data file and processed. Most parts of the code are
commented, and therefore it should be easy to follow what is happening. This
example coded in NCL, since many cmorizing scripts that are available so far
are written in NCL, and adapting existing code and using existing libraries is
easier than writing something totally new. However, in theory this script could
also be written in Python but without the help of the Iris package and rather
based on NetCDF4 or similar packages (this is why we need the cmorizing
scripts). There is are several python-based cmorizing scripts available
already, that can be used as guideline in case you would like to write your
cmorizing script in python. 

For the second part of the program, the following points are important to keep in mind:

- fname: it is the combination of the input path that is defined in the
  configuration file (see Section 2) that has to be defined to run the
  cmorizing script, and the name of the file with the ``raw`` data; 
- ``output = f->xch4``: In this line it is hardcoded that the variable with the
  name ``xch4`` is processed. If you have defined more than one variable, this
  statement has to be adjusted, so that the correct variable name is used with
  each loop of the program. 
- ``format_coords``: this call is a routine that is available for NCL code
  already and which takes care of cmorizing the coordinates of the current
  variable if necessary (e.g., longitudes ranging from -180 to 180 degrees
  instead of 0 to 360 degrees). 
- ``fout``: the filepath and filename of the output file are set here. The path
  is taken from the configuration file (see Section 2) that is necessary to run
  the cmorizing script, and the filename is put together from the
  information given in the first part of the script, following the rules for
  filenames so that the ESMValTool can read in the files. 

The script as it is detailed here would only be able to correct some minor
problems with the coordinates (e.g. latitudes in the wrong order, longitudes in
the wrong order, etc.). Everything else will have to be added to the script for
it to deal with it. There are many cmorizing scripts available in
`/esmvaltool/utils/cmorizers/obs/
<https://github.com/ESMValGroup/ESMValTool/blob/version2_development/esmvaltool/utils/cmorizers/obs/>`_
where solutions to many kinds of format issues with observational data are
addressed. Most of these scripts are written in NCL, but there more and more
examples for Python-based cmorizing scripts become available.

How much cmorizing an observational data set needs is strongly dependent on
the original NetCDF file and how close the original formatting already is to
the strict CMOR standard. 

5. Run the cmorizing script
===========================

In order to actually run the cmrorizing script, you have to use the following statement:

.. code-block:: console

 cmorize_obs -c <config-user.yml> -o <dataset-name>

This call only works, of course, if you are already in the folder in which also
the ``config-user.yml`` is stored.

.. note::

   The output path given in the configuration file is the path where
   your cmorized dataset will be stored. The ESMValTool will create a folder
   with the correct tier information (see Section 2) if that tier folder is not
   already available, and then a folder named after the data set. In this
   folder the cmorized data set will be stored as a netCDF file. 

Your run was successful if a netCDF file was produced in your output directory,
and if at some point the output on your screen shows an info line similar to
this and no additional error message is shown (note: this is the example for
the example data set ``xch4``): 

.. code-block:: console

  INFO    Processing xch4 (Amon)

6. Naming convention of the observational data files
====================================================

For the ESMValTool to be able to read the observations from the NetCDF file,
the file name needs a very specific structure and order of information parts
(very similar to the naming convention for observations in ESMValTool
v1.0). The file name will be automatically correctly created if a cmorizing
script has been used to create the netCDF file. If the file with the
observations is already in the correct CMOR format, it is possible to use the
data set without using a cmorizing script. In this case it is important to
follow the guidelines below on how the filename should be structured.

The file for the ``CDS-XCH4`` observations example in the correct format is
named as follows:

.. code-block:: console

  OBS_CDS-XCH4_sat_L3_Amon_xch4_200301-201612.nc

The different parts of the name are explained in more detail here:

- OBS: describes what kind of data can be expected in the file, in this case
  ``observations``; 
- CDS-XCH4: that is the name of the dataset. It has been named this way for
  illustration purposes (so that everybody understands it is the xch4 dataset
  downloaded from the CDS), but a better name would indeed be ``ESACCI-XCH4``
  since it is a ESA-CCI dataset; 
- sat: describes the source of the data, here we are looking at satellite data
  (therefore ``sat``), could also be ``reanaly`` for reanalyses;
- L3: describes the version of the dataset:
- Amon: is the information in which ``mip`` the variable is to be expected, and
  what kind of temporal resolution it has; here we expect ``xch4`` to be part
  of the atmosphere (``A``) and we have the dataset in a monthly resolution
  (``mon``);
- xch4: Is the name of the variable. Each observational data file is supposed
  to only include one variable per file; 
- 200301-201612: Is the period the dataset spans with ``200301`` being the
  start year and month, and ``201612`` being the end year and month;

.. note::
   There is a different naming convention for ``obs4mips`` data (see the exact
   specifications for the obs4mips data file naming convention in the
   ``config-developer.yml`` file).

7. Running a preprocessing test recipe
======================================

To verify that the cmorized data file is indeed correctly formatted
(following the strict CMOR standard Iris needs), it is good to run a
preprocessing test recipe, that does not include any diagnostic, but only reads
in the data file and has it processed in the preprocessor. Such a recipe is
called ``recipes/examples/recipe_check_obs.yml``. It can be adjusted to only
define one preprocessor that is then called in the diagnostics part. **Note:**
For the recipe to only run the preprocessor and no diagnostic analyses
whatsoever, make sure to set the "script" key in the diagnostics part to
"null".

If the recipe is adjusted as outlined above, run it with the following call:

.. code-block:: console

  esmvaltool -c *name_of_your_configuration_file* recipes/examples/recipe_preprocessor_test.yml

This call only works, of course, if you are already in the folder in which also
the configuration file ``name_of_your_configuration_file`` is stored. If your
cmorized data set has the correct format, the ESMValTool will read the data and
run it through the preprocessor, basically confirming that Iris can work with
the data field. The ESMValTool will then print the statement ``Run was
successful``. If something is still wrong with the format, the ESMValTool will
crash. 
