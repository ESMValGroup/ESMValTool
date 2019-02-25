.. _config_developper:

*******************
Configuration files
*******************

There are several configuration files in ESMValTool:

  - config-user.yml
  - config-developer.yml
  - config-references.yml
  - environment.yml


User config file
================

The ``config-user.yml`` file is the only user-facing configuration file and is
described in the user guide section of this document.


Developper config file
======================

This configuration file describes the file system structure for several
key projects (CMIP5, CMIP6) on several key machines (BADC, CP4CDS, DKRZ, ETHZ,
SMHI, BSC).

The data directory structure of the CMIP5 project is set up differently
at each site. The following code snipper is an example of several paths
descriptions for the CMIP5 at various sites:

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

As an example, the CMIP5 file path on BADC would be:

.. code-block:: yml

        [institute]/[dataset ]/[exp]/[frequency]/[modeling_realm]/[mip]/[ensemble]/latest/[short_name]

When loading these files, ESMValTool replaces the placeholders with the true
values. The resulting real path would look something like this:

.. code-block:: yml

    MOHC/HadGEM2-CC/rcp85/mon/ocean/Omon/r1i1p1/latest/tos


References config file
================================

The ``config-references.yml`` file is the full list of ESMValTool authors,
references and projects. Each author, project and reference in the documentation
section of a recipe needs to be in this file in the relevant section.

For instance, the recipe ``recipe_ocean_example.yml`` file contains the following
documentation section:

.. code-block:: yml

  documentation
    authors:
      - demo_le

    maintainer:
      - demo_le

    references:
      - demora2018gmd

    projects:
      - ukesm


All four items here are named people, references and projects listed in the 
``config-references.yml`` file.



Environment config file
================================

This is the conda envronment which informs conda as to which packages
and which versions are needed in order to install ESMValTool.
It is called in the update and install commands:

.. code-block:: bash

    conda env install --name esmvaltool --file ESMValTool/environment.yml
    conda env update --name esmvaltool --file ESMValTool/environment.yml
