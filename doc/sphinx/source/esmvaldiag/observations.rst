************************************************************
Contributing a CMORizing script for an observational dataset
************************************************************

ESMValTool is designed to work with `CF compliant <http://cfconventions.org/>`_
data and follows the CMOR tables from the CMIP data request, therefore
the observational datasets need to be CMORized for usage in ESMValTool.
The following steps are necessary to prepare an observational
data set for the use in ESMValTool.

| `1. Check if your variable is CMOR standard`_
| `2. Edit your configuration file`_
| `3. Store your dataset in the right place`_
| `4. Create a cmorizer for the dataset`_
| `4.1 Cmorizer script written in python`_
| `4.2 Cmorizer script written in NCL`_
| `5. Run the cmorizing script`_
| `6. Naming convention of the observational data files`_
| `7. Test the cmorized dataset`_


1. Check if your variable is CMOR standard
==========================================

Most variables are defined in the CMIP data request and can be found in the
CMOR tables in the folder `/esmvalcore/cmor/tables/cmip6/Tables/
<https://github.com/ESMValGroup/ESMValCore/tree/development/esmvalcore/cmor/tables/cmip6/Tables>`_,
differentiated according to the ``MIP`` they belong to. The tables are a
copy of the `PCMDI <https://github.com/PCMDI>`_ guidelines. If you find the
variable in one of these tables, you can proceed to the next section.

If your variable is not available in the standard CMOR tables,
you need to write a custom CMOR table entry for the variable
as outlined below and add it to `/esmvalcore/cmor/tables/custom/
<https://github.com/ESMValGroup/ESMValCore/tree/development/esmvalcore/cmor/tables/custom>`_.

To create a new custom CMOR table you need to follow these
guidelines:

- Provide the ``variable_entry``;
- Provide the ``modeling_realm``;
- Provide the variable attributes, but leave ``standard_name`` blank. Necessary
  variable attributes are: ``units``, ``cell_methods``, ``cell_measures``,
  ``long_name``, ``comment``.
- Provide some additional variable attributes. Necessary additional variable
  attributes are: ``dimensions``, ``out_name``, ``type``. There are also
  additional variable attributes that can be defined here (see the already
  available cmorizers).

It is recommended to use an existing custom table as a template, to edit the
content and save it as ``CMOR_<short_name>.dat``.

2. Edit your configuration file
===============================

Make sure that beside the paths to the model simulations and observations, also
the path to raw observational data to be cmorized (``RAWOBS``) is present in
your configuration file.

3. Store your dataset in the right place
========================================

The folder ``RAWOBS`` needs the subdirectories ``Tier1``, ``Tier2`` and
``Tier3``. The different tiers describe the different levels of restrictions
for downloading (e.g. providing contact information, licence agreements)
and using the observations. The unformatted (raw) observations
should then be stored then in the appropriate of these three folders.

4. Create a cmorizer for the dataset
====================================

There are many cmorizing scripts available in `/esmvaltool/cmorizers/obs/
<https://github.com/ESMValGroup/ESMValTool/blob/version2_development/esmvaltool/cmorizers/obs/>`_
where solutions to many kinds of format issues with observational data are
addressed. Most of these scripts are written in NCL at the moment, but more
and more examples for Python-based cmorizing scripts become available.

.. note::
  NCL support will terminate soon, so new cmorizer scripts should preferably be
  written in Python.

How much cmorizing an observational data set needs is strongly dependent on
the original NetCDF file and how close the original formatting already is to
the strict CMOR standard.

In the following two subsections two cmorizing scripts, one written in Python
and one written in NCL, are explained in more detail.

4.1 Cmorizer script written in python
-------------------------------------

Find here an example of a cmorizing script, written for the ``MTE`` dataset
that is available at the MPI for Biogeochemistry in Jena: `cmorize_obs_mte.py
<https://github.com/ESMValGroup/ESMValTool/blob/version2_development/esmvaltool/cmorizers/obs/cmorize_obs_mte.py>`_.

All the necessary information about the dataset to write the filename
correctly, and which variable is of interest, is stored in a seperate
configuration file: `MTE.yml
<https://github.com/ESMValGroup/ESMValTool/blob/version2_development/esmvaltool/cmorizers/obs/cmor_config/MTE.yml>`_
in the directory ``ESMValTool/esmvaltool/cmorizers/obs/cmor_config/``. Note
that the name of this configuration file has to be identical to the name of
your data set. It is recommended that you set ``project`` to ``OBS6`` in the
configuration file. That way, the variables defined in the CMIP6 CMOR table,
augmented with the custom variables described above, are available to your script.

The first part of this configuration file defines the filename of the raw
observations file, the second part defines the common global attributes for
the cmorizer output, e.g. information that is needed to piece together the
final observations file name in the correct structure (see Section `6. Naming convention of the observational data files`_). The
third part defines the variables that are supposed to be cmorized.

The actual cmorizing script ``cmorize_obs_mte.py`` consists of a header with
information on where and how to download the data, and noting the last access
of the data webpage.

The main body of the CMORizer script must contain a function called

.. code-block:: python

   def cmorization(in_dir, out_dir, cfg, config_user):

with this exact call signature. Here, ``in_dir`` corresponds to the input
directory of the raw files, ``out_dir`` to the output directory of final
reformatted data set and ``cfg`` to the configuration dictionary given by
the  ``.yml`` configuration file. The return value of this function is ignored. All
the work, i.e. loading of the raw files, processing them and saving the final
output, has to be performed inside its body. To simplify this process, ESMValTool
provides a set of predefined utilities.py_, which can be imported into your CMORizer
by

.. code-block:: python

   from . import utilities as utils

Apart from a function to easily save data, this module contains different kinds
of small fixes to the data attributes, coordinates, and metadata which are
necessary for the data field to be CMOR-compliant.

Note that this specific CMORizer script contains several subroutines in order
to make the code clearer and more readable (we strongly recommend to follow
that code style). For example, the function ``_get_filepath`` converts the raw
filepath to the correct one and the function ``_extract_variable`` extracts and
saves a single variable from the raw data.

.. _utilities.py: https://github.com/ESMValGroup/ESMValTool/blob/version2_development/esmvaltool/cmorizers/obs/utilities.py


4.2 Cmorizer script written in NCL
----------------------------------

Find here an example of a cmorizing script, written for the ``ESACCI XCH4``
dataset that is available on the Copernicus Climate Data Store:
`cmorize_obs_cds_xch4.ncl
<https://github.com/ESMValGroup/ESMValTool/blob/version2_development/esmvaltool/cmorizers/obs/cmorize_obs_cds_xch4.ncl>`_.

The first part of the script collects all the information about the dataset
that are necessary to write the filename correctly and to understand which
variable is of interest here. Please make sure to provide the correct
information for following key words: DIAG_SCRIPT, VAR, NAME, MIP, FREQ,
CMOR_TABLE.

- **Note:** the fields ``VAR``, ``NAME``, ``MIP`` and ``FREQ`` all ask for one
  or more entries. If more than one entry is provided, make sure that the order
  of the entries is the same for all four fields! (for example, that the first
  entry in all four fields describe the variable ``xch4`` that you would like
  to extract);
- **Note:** some functions in the script are NCL-specific and are available
  through the loading of the script interface.ncl_. There are similar
  functions available for python scripts.

.. _interface.ncl: https://github.com/ESMValGroup/ESMValTool/blob/version2_development/esmvaltool/cmorizers/obs/interface.ncl

.. _utilities.ncl: https://github.com/ESMValGroup/ESMValTool/blob/version2_development/esmvaltool/cmorizers/obs/utilities.ncl

In the second part of the script each variable defined in ``VAR`` is separately
extracted from the original data file and processed. Most parts of the code are
commented, and therefore it should be easy to follow. ESMValTool provides a set
of predefined utilities.ncl_, which can be imported into your CMORizer
by

.. code-block:: NCL

   loadscript(getenv("esmvaltool_root") + "/esmvaltool/cmorizers/obs/utilities.ncl")

This module contains different kinds of small fixes to the data attributes,
coordinates, and metadata which are necessary for the data field to be
CMOR-compliant.

5. Run the cmorizing script
===========================

The cmorizing script for the given dataset can be run with:

.. code-block:: console

 cmorize_obs -c <config-user.yml> -o <dataset-name>


.. note::

   The output path given in the configuration file is the path where
   your cmorized dataset will be stored. The ESMValTool will create a folder
   with the correct tier information (see Section `2. Edit your configuration file`_) if that tier folder is not
   already available, and then a folder named after the data set. In this
   folder the cmorized data set will be stored as a netCDF file.

If your run was successful, one or more NetCDF files are produced in your
output directory.


6. Naming convention of the observational data files
====================================================

For the ESMValTool to be able to read the observations from the NetCDF file,
the file name needs a very specific structure and order of information parts
(very similar to the naming convention for observations in ESMValTool
v1.0). The file name will be automatically correctly created if a cmorizing
script has been used to create the netCDF file.

The correct structure of an observational data set is defined in
`config-developer.yml
<https://github.com/ESMValGroup/ESMValCore/blob/development/esmvalcore/config-developer.yml>`_,
and looks like the following:

.. code-block:: console

  OBS_[dataset]_[type]_[version]_[mip]_[short_name]_YYYYMM-YYYYMM.nc

For the example of the ``CDS-XCH4`` data set, the correct structure of the
file name looks then like this:

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

7. Test the cmorized dataset
======================================

To verify that the cmorized data file is indeed correctly formatted, you can
run a dedicated test recipe, that does not include any diagnostic, but only
reads in the data file and has it processed in the preprocessor. Such a recipe
is called ``recipes/examples/recipe_check_obs.yml``. You just need to add a
diagnostic for your dataset following the existing entries.
Only the diagnostic of interest needs to be run, the others should be commented
out for testing.

