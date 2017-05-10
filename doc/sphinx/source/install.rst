Software installation
*********************

.. _prerequisites:

Prerequisites
=============

The ESMValTool has the following software requirements (note that specific diagnostics might require additional software packages):

* Unix(-like) operating system
* Python version 2.7.x for running the Python script main.py; most diagnostics written in Python require installation of additional Python packages such as, for instance, Geometry Engine (GEOS), scientificpython, netCDF4, cdo, geoval, cartopy, and iris.

  *The required Python packages can be installed with the following commands:*

  conda install basemap

  conda install --channel https://conda.anaconda.org/Clyde_Fare scientificpython

  conda install netcdf4

  conda install --channel https://conda.anaconda.org/auto cdo

  pip install geoval

  conda install -c scitools cartopy

  conda install -c scitools iris

  .. attention:: It is strongly recommended to use the Python distribution Anaconda (https://www.continuum.io/), as it allows the user to install additional Python libraries and extensions in a simple way and without modifying the installed Python distribution (i.e., without root permissions). The installation instructions for the additional Python packages listed above are given for Anaconda.

* NCAR Command Language (NCL 2014) version 6.2 or higher (note: NCL version 6.3 is not supported, see known issues) to run the quality check and reformat routines processing all input files. See the control flow description on reformat_default in Table S11 for details.
* The statistical computing software R to run diagnostics written in R. A working installation of R and the executable Rscript in the
* default search path are required. In addition, the netCDF for R libraries (ncdf / ncdf4) are needed. Currently, only the diagnostic
* "Standardized Precipitation index (SPI)" (see Annex C) requires R. More diagnostics written in R might be added in the future.
* The sea ice diagnostics (and derived diagnostics such as, for instance, the ESA CCI namelist) require the Climate Data Operators (CDO): https://code.zmaw.de/projects/cdo. The CDO executable has to be in the default search path (callable via the command "cdo").
* Input files in netCDF with required attributes and dimension names. Valid input files are:

  * files in CMIP or similarly standardized format using a CMIP5 table, or with discrepancies that can be handled via the definitions in the files reformat_scripts/recognized_units.dat and reformat_scripts/recognized_units.dat, respectively.
  * any input file with a (user-)supplied reformat routine that converts the input data during run-time, see the control flow description on reformat_EMAC in Table S11 for details

* Common GNU utilities such as "wc", "date", "basename", and "more", which are usually part of the standard Linux distribution.

Obtaining the source code
=========================

The ESMValTool is available on GitHub at https://github.com/ESMValGroup/ESMValTool (see also :numref:`git_repository`). The ESMValTool is released under the Apache License, version 2.0 and citation of the ESMValTool paper ("Software Documentation Paper") is kindly requested upon use alongside with the software doi (doi:10.17874/ac8548f0315) and version number:

  *Eyring et al., ESMValTool (v1.0) - a community diagnostic and performance metrics tool for routine evaluation of Earth System Models in CMIP, Geosci. Model Dev., 9, 1747-1802, 2016.*

Besides the above citation, users are kindly asked to register any journal articles (or other scientific documents) that use the software at the ESMValTool webpage (http://www.esmvaltool.org/). Citing the Software Documentation Paper and registering your paper(s) will serve to document the scientific impact of the Software, which is of vital importance for securing future funding. You should consider this an obligation if you have taken advantage of the ESMValTool, which represents the end product of considerable effort by the development team.

**The ESMValTool will be further developed in a version controlled repository (see** :numref:`git_repository` **for details).** In addition to using the software, we would therefore like to encourage the community to join the Software Development Team and to contribute additional diagnostics and performance metrics or other software improvements. Contributing back the new diagnostics and performance metrics or other software improvements will help to enhance the capability of the Software, which is of vital importance for securing future funding. You should consider this an obligation if you have taken advantage of the Software, which represents a product of considerable effort by the development team.

A wiki page (see :numref:`wiki` for details) that describes ongoing developments is available for the ESMValTool development team. Interested developers are welcome to contact the core development team (see :numref:`core_dev_team`).

Software installation
=====================

The ESMValTool can be downloaded from GitHub (:numref:`git_repository`) to any local directory. While the ESMValTool itself does not need to be installed besides downloading/copying the ESMValTool directories to a local folder, it relies on specific software to be available on your system. Please see :numref:`prerequisites` for details.

