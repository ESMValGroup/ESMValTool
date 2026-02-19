.. _new-cmorizer:

Writing a CMORizer script for an additional dataset
***************************************************

ESMValTool is designed to work with `CF compliant <http://cfconventions.org/>`_
data and follows the CMOR tables from the CMIP data request, therefore
the observational datasets need to be CMORized for usage in ESMValTool.
The following steps are necessary to prepare an observational
data set for the use in ESMValTool.

| `1. Check if your variable is CMOR standard`_
| `2. Store your dataset in the right place`_
| `2.1 Downloader script (optional)`_
| `3. Create a cmorizer for the dataset`_
| `3.1 Cmorizer script written in python`_
| `3.2 Cmorizer script written in NCL`_
| `4. Run the cmorizing script`_
| `5. Naming convention of the observational data files`_
| `6. Test the cmorized dataset`_

.. note::
  **CMORization as a fix.** As of early 2020, we've started implementing cmorization as
  *fixes*. As compared to the workflow described below, this has the advantage that
  the user does not need to store a duplicate (CMORized) copy of the data. Instead, the
  CMORization is performed 'on the fly' when running a recipe. **ERA5** is the first dataset
  for which this 'CMORization on the fly' is supported. For more information, see
  :ref:`inputdata_native_datasets`.


1. Check if your variable is CMOR standard
==========================================

Most variables are defined in the CMIP data request and can be found in the
CMOR tables in the folder `/esmvalcore/cmor/tables
<https://github.com/ESMValGroup/ESMValCore/tree/main/esmvalcore/cmor/tables>`_,
differentiated according to the ``MIP`` they belong to. A more extensive
introduction is available in :ref:`esmvalcore:cmor_tables`. If you find the
variable in one of these tables, you can proceed to the next section.

If your variable is not available in the standard CMOR tables,
you need to write a custom CMOR table entry for the variable
as outlined below and add it to `/esmvalcore/cmor/tables/cmip6-custom/
<https://github.com/ESMValGroup/ESMValCore/tree/main/esmvalcore/cmor/tables/cmip6-custom>`_.

To create a new custom CMOR table you need to follow these
guidelines:

- Provide the ``variable_entry`` (for CMIP5-style tables) or the variable name
  (``short_name``) followed by an optional underscore and branding suffix
  as the key for the variable definition (for CMIP6-style tables);
- Provide the ``modeling_realm``;
- Provide the variable attributes, but leave ``standard_name`` blank unless
  the variable has a valid standard name. Necessary
  variable attributes are: ``units``, ``cell_methods``, ``cell_measures``,
  ``long_name``, ``comment``.
- Provide some additional variable attributes. Necessary additional variable
  attributes are: ``dimensions``, ``out_name`` (usually equal to ``short_name``).
  There are also additional variable attributes that can be defined here (see
  the already available custom variables).

It is recommended to extend the file
`esmvalcore/cmor/tables/cmip6-custom/CMIP6_custom.json <https://github.com/ESMValGroup/ESMValCore/blob/main/esmvalcore/cmor/tables/cmip6-custom/CMIP6_custom.json>`__
for CMIP6-style tables. For CMIP5-style projects, use an existing custom table from
`esmvalcore/cmor/tables/cmip5-custom/ <https://github.com/ESMValGroup/ESMValCore/blob/main/esmvalcore/cmor/tables/cmip5-custom>`__
as a template, edit the content and save it in the same directory as ``CMOR_<short_name>.dat``.

2. Store your dataset in the right place
========================================

The folder specified by ``--original-data-dir`` needs the subdirectories ``Tier1``, ``Tier2`` and
``Tier3``. The different tiers describe the different levels of restrictions
for downloading (e.g. providing contact information, licence agreements)
and using the observations. The unformatted (raw) observations
should then be stored in the appropriate of these three folders.

For each additional dataset, an entry needs to be made to the file
`datasets.yml
<https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/data/datasets.yml>`_.
The dataset entry should contain:

- the correct ``tier`` information;
- the ``source`` of the raw data;
- the ``last_access`` date;
- the ``info`` that explain how to download the data.

Note that these fields should be identical to the content of the header
of the cmorizing script (see Section `3. Create a cmorizer for the dataset`_).

2.1 Downloader script (optional)
--------------------------------

A Python script can be written to download raw observations
from source and store the data in the appropriate tier subdirectory of the
folder specified by the ``--original-data-dir`` flag automatically.
There are many downloading scripts available in
`/esmvaltool/cmorizers/data/downloaders/datasets/
<https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/cmorizers/data/downloaders/datasets>`_
where several data download mechanisms are provided:

- A `wget` get based downloader for http(s) downloads, with a specific derivation for NASA datasets.
- A `ftp` downloader with a specific derivation for ESACCI datasets available from CEDA.
- A Climate Data Store downloader based on `cdsapi`.

Note that the name of this downloading script has to be identical to the
name of the dataset.

Depending on the source server, the downloading script needs to contain paths to
raw observations, filename patterns and various necessary fields to retrieve
the data.
Default ``start_date`` and ``end_date`` can be provided in cases where raw data
are stored in daily, monthly, and yearly files.

The downloading script for the given dataset can be run with:

.. code-block:: console

 esmvaltool data download --original-data-dir </path/to/save/data>  <dataset-name>

The options ``--start`` and ``--end`` can be added to the command above to
restrict the download of raw data to a time range. They will be ignored if a specific dataset
does not support it (i.e. because it is provided as a single file). Valid formats are
``YYYY``, ``YYYYMM`` and ``YYYYMMDD``. By default, already downloaded data are not overwritten
unless the option ``--overwrite=True`` is used.

3. Create a cmorizer for the dataset
====================================

There are many cmorizing scripts available in
`/esmvaltool/cmorizers/data/formatters/datasets/
<https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/data/formatters/datasets/>`_
where solutions to many kinds of format issues with observational data are
addressed. These scripts are either written in Python or in NCL.

.. note::
  NCL support will terminate soon, so new cmorizer scripts should preferably be
  written in Python.

How much cmorizing an observational data set needs is strongly dependent on
the original NetCDF file and how close the original formatting already is to
the strict CMOR standard.

In the following two subsections two cmorizing scripts, one written in Python
and one written in NCL, are explained in more detail.

3.1 Cmorizer script written in python
-------------------------------------

Find here an example of a cmorizing script, written for the ``MTE`` dataset
that is available at the MPI for Biogeochemistry in Jena: `mte.py
<https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/data/formatters/datasets/mte.py>`_.

All the necessary information about the dataset to write the filename
correctly, and which variable is of interest, is stored in a separate
configuration file: `MTE.yml
<https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/data/cmor_config/MTE.yml>`_
in the directory ``ESMValTool/esmvaltool/cmorizers/data/cmor_config/``. Note
that both the name of this configuration file and the cmorizing script have to be
identical to the name of your dataset.
It is recommended that you set ``project`` to ``OBS6`` in the
configuration file. That way, the variables defined in the CMIP6 CMOR table,
augmented with the custom variables described above, are available to your script.

The first part of this configuration file defines the filename of the raw
observations file. The second part defines the common global attributes for
the cmorizer output, e.g. information that is needed to piece together the
final observations file name in the correct structure (see Section `5. Naming convention of the observational data files`_).
Another global attribute is ``reference`` which includes a ``doi`` related to the dataset.
Please see the section `adding references
<https://docs.esmvaltool.org/en/latest/community/diagnostic.html#adding-references>`_
on how to add reference tags to the ``reference`` section in the configuration file.
If a single dataset has more than one reference,
it is possible to add tags as a list e.g. ``reference: ['tag1', 'tag2']``.
The third part in the configuration file defines the variables that are supposed to be cmorized.

The actual cmorizing script ``mte.py`` consists of a header with
information on where and how to download the data, and noting the last access
of the data webpage.

The main body of the CMORizer script must contain a function called

.. code-block:: python

   def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):

with this exact call signature. Here, ``in_dir`` corresponds to the input
directory of the raw files, ``out_dir`` to the output directory of final
reformatted data set, ``cfg`` to the dataset-specific configuration file,
``cfg_user`` to the configuration object (which behaves basically like a
dictionary), ``start_date`` to the start
of the period to format, and ``end_date`` to the end of the period to format.
If not needed, the last three arguments can be ignored using underscores.
The return value of this function is ignored.
Note that this function will trigger a Codacy style issue because it has more
than 5 arguments; you can ignore that issue in line with the :ref:`code_quality` guidance.

All the work, i.e. loading of the raw files, processing them and saving the final
output, has to be performed inside its body. To simplify this process, ESMValTool
provides a set of predefined utilities.py_, which can be imported into your CMORizer
by

.. code-block:: python

   from esmvaltool.cmorizers.data import utilities as utils

Apart from a function to easily save data, this module contains different kinds
of small fixes to the data attributes, coordinates, and metadata which are
necessary for the data field to be CMOR-compliant.

Note that this specific CMORizer script contains several subroutines in order
to make the code clearer and more readable (we strongly recommend to follow
that code style). For example, the function ``_get_filepath`` converts the raw
filepath to the correct one and the function ``_extract_variable`` extracts and
saves a single variable from the raw data.

.. _utilities.py: https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/data/utilities.py


3.2 Cmorizer script written in NCL
----------------------------------

.. warning::

  Writing new NCL code is not recommended because the
  `NCL interpreter <https://github.com/NCAR/ncl>`__ is no longer maintained.

Find here an example of a cmorizing script, written for the ``ESACCI XCH4``
dataset that is available on the Copernicus Climate Data Store:
`cds_xch4.ncl
<https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/data/formatters/datasets/cds_xch4.ncl>`_.

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

.. _interface.ncl: https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/data/formatters/interface.ncl

.. _utilities.ncl: https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/data/formatters/utilities.ncl

In the second part of the script each variable defined in ``VAR`` is separately
extracted from the original data file and processed. Most parts of the code are
commented, and therefore it should be easy to follow. ESMValTool provides a set
of predefined utilities.ncl_, which are imported by default into your CMORizer.
This module contains different kinds of small fixes to the data attributes,
coordinates, and metadata which are necessary for the data field to be
CMOR-compliant.

4. Run the cmorizing script
===========================

The cmorizing script for the given dataset can be run with:

.. code-block:: console

 esmvaltool data format --config_dir </path/to/config/dir/> <dataset-name>

The options ``--start`` and ``--end`` can be added to the command above to
restrict the formatting of raw data to a time range. They will be ignored if a specific dataset
does not support it (i.e. because it is provided as a single file). Valid formats are
``YYYY``, ``YYYYMM`` and ``YYYYMMDD``.

.. note::

   The ``output_dir`` path given in the :ref:`configuration file <esmvalcore:config>` is the path where
   your cmorized dataset will be stored. The ESMValTool will create a folder
   with the correct tier information if that tier folder is not
   already available, and then a folder named after the dataset.
   In this folder the cmorized data set will be stored as a NetCDF file.
   The cmorized dataset will be automatically moved to the correct tier
   subfolder of your OBS or OBS6 directory if the option
   ``--install=True`` is used in the command above and no such directory
   was already created.

If your run was successful, one or more NetCDF files are produced in your
output directory.

If a downloading script is available for the dataset, the downloading and
the cmorizing scripts can be run in a single command with:

.. code-block:: console

 esmvaltool data prepare --config_dir </path/to/config/dir/> <dataset-name>

Note that options from the ```esmvaltool data download`` and
``esmvaltool data format`` commands can be passed to the above command.

5. Naming convention of the observational data files
====================================================

For the ESMValTool to be able to read the observations from the NetCDF file,
the file name needs a very specific structure and order of information parts
(very similar to the naming convention for observations in ESMValTool
v1.0). The file name will be automatically correctly created if a cmorizing
script has been used to create the netCDF file.

The correct structure of an observational data set is defined in
`data-local-esmvaltool.yml
<https://github.com/ESMValGroup/ESMValCore/blob/main/esmvalcore/config/configurations/data-local-esmvaltool.yml>`_,
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
- 200301-201812: Is the period the dataset spans with ``200301`` being the
  start year and month, and ``201812`` being the end year and month;

.. note::
   There is a different naming convention for ``obs4MIPs`` data (see the exact
   specifications for the obs4MIPs data file naming convention in the
   `data-local.yml <https://github.com/ESMValGroup/ESMValCore/blob/main/esmvalcore/config/configurations/data-local.yml>`
   file).

6. Test the cmorized dataset
============================

To verify that the cmorized data file is indeed correctly formatted, you can
run a dedicated test recipe, that does not include any diagnostic, but only
reads in the data file and has it processed in the preprocessor. Such a recipe
is called ``recipes/examples/recipe_check_obs.yml``. You just need to add a
diagnostic for your dataset following the existing entries.
Only the diagnostic of interest needs to be run, the others should be commented
out for testing.
