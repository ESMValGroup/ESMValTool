.. _config_developper:

*******************
Configuration files
*******************

There are several configuration files in ESMValTool:

  - config-user.yml
  - config-developer.yml
  - config-references.yml
  - environment.yml

This section describes those files.


User config file
================

This file is described in the user guide section of this document.


Developper config file
======================

This configuration file describes the file system structure for several
key projects (CMIP5, CMIP6) on several key machines (BADC, CP4CDS, DKRZ, ETHZ,
SMHI, BSC).

For instance, the code structure of the CMIP5 project is set up differently
at each site.

.. code-block:: yml

  CMIP5:
    input_dir:
      default: '/'
      BADC: '[institute]/[dataset]/[exp]/[frequency]/[modeling_realm]/[mip]/[ensemble]/latest/[short_name]'
      CP4CDS: '[institute]/[dataset]/[exp]/[frequency]/[modeling_realm]/[mip]/[ensemble]/[short_name]/latest/'
      DKRZ: '[institute]/[dataset]/[exp]/[frequency]/[modeling_realm]/[mip]/[ensemble]/[latestversion]/[short_name]'
      ETHZ: '[exp]/[mip]/[short_name]/[dataset]/[ensemble]/'
      SMHI: '[dataset]/[ensemble]/[exp]/[frequency]'
      BSC: '[project]/[exp]/[dataset.lower]'

So, for instance, a typical file path on BADC would be:

.. code-block:: yml

    MOHC/HadGEM2-CC/rcp85/mon/ocean/Omon/r1i1p1/latest


To be more explicit

.. code-block:: yml

    [institute]/[dataset ]/[exp]/[frequency]/[modeling_realm]/[mip]/[ensemble]/latest/[short_name]
    MOHC       /HadGEM2-CC/rcp85/mon        /ocean           /Omon /r1i1p1    /latest/tos

.. warning::
    Section to be completed

References config file
================================

The ``config-references.yml`` file is the full list of ESMValTool authors,
references and projects. Each author in the header of the recipe needs to
be in the authors list.


Environment config file
================================

This is the conda envronment which informs conda as to which packages
and which versions are needed in order to install ESMValTool.
It is called in the update and install commands:

.. code-block:: bash

    conda env install --name esmvaltool --file ESMValTool/environment.yml
    conda env update --name esmvaltool --file ESMValTool/environment.yml
