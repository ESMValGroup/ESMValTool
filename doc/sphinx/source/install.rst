Software installation
*********************

.. _prerequisites:

Prerequisites
=============

The ESMValTool has the following software requirements (note that specific diagnostics might require additional software packages):

* Unix(-like) operating system
* NCAR Command Language (NCL 2014, http://www.ncl.ucar.edu/) version 6.4 or higher (note: NCL version 6.3 is not supported, see known issues, Part :numref:`known_issues`).
* Common GNU utilities such as "wc", "date", "basename", and "more", which are usually part of the standard Linux distribution.
* The statistical computing software R (https://www.r-project.org/) to run diagnostics written in R. A working installation of R and the executable Rscript in the default search path are required. In addition, the netCDF for R libraries (ncdf / ncdf4) are needed. Currently, only the diagnostics "Standardized Precipitation index (SPI)" and "Ozone and associated climate impacts (Eyring13, fig. 6)" (see Part :numref:`annex_c`) require R. More diagnostics written in R might be added in the future.
* The sea ice diagnostics (and derived diagnostics such as, for instance, the ESA CCI namelist - see Section :numref:`nml_esacci`) require the Climate Data Operators (cdo): https://code.zmaw.de/projects/cdo. The cdo executable has to be in the default search path (callable via the command "cdo").
* Python version 2.7.x for running the ESMValTool workflow manager main.py; most diagnostics written in Python require installation of additional Python packages such as, for instance, Geometry Engine (GEOS), scientificpython, netCDF4, cdo, geoval, cartopy, and iris.

*The required Python packages can be installed with the following commands:*

1) install *anaconda* (2-5.0.1 or later):

   * download anaconda from https://www.anaconda.com/download/#linux
   * install anaconda, e.g.

   .. code:: bash

      chmod 755 Anaconda2-5.0.1-Linux-x86_64.sh
      Anaconda2-5.0.1-Linux-x86_64.sh

2) install *basemap*

   .. code:: bash

      conda install basemap

3) install netcdf library

   .. code:: bash

      conda install netcdf4

4) install *geoval*

   .. code:: bash

      git clone https://github.com/pygeo/geoval.git
      cd geoval
      python setup.py build
      python setup.py install

   create symbolic link for geoval-lib in anaconda directory, e.g.

   .. code:: bash

      ln -s /home/username/geoval/build/lib.linux-x86_64-2.7/geoval /home/username/anaconda2/lib/

5) install *python cdo*

   .. code:: bash

      conda config --add channels conda-forge
      conda install cdo
      conda install python-cdo

6) it might be needed to replace the cdo executable with more stable version if the version installed in step 5) crashes (e.g. Python diagnostics of namelist "ESA CCI" (Section :numref:`nml_esacci`)).

   * go to https://code.mpimet.mpg.de/projects/cdo/ and download executable or source code and compile your own executable
   * copy the new cdo executable to your anaconda bin directory, e.g. /home/username/anaconda2/bin/

7) install *cartopy*

   .. code:: bash

      conda install cartopy

8) install *gdal*

   .. code:: bash

      conda install gdal
      
9) install *iris*

   .. code :: bash
   
      conda install iris

10) update all conda packages

   .. code:: bash

      conda update --all

.. attention:: It is strongly recommended to use the Python distribution Anaconda (https://www.continuum.io/), as it allows the user to install additional Python libraries and extensions in a simple way and without modifying the installed Python distribution (i.e., without root permissions). The installation instructions for the additional Python packages listed above are given for Anaconda.

Obtaining the source code
=========================

The ESMValTool is available on GitHub at https://github.com/ESMValGroup/ESMValTool (for details see Section :numref:`git_repository`). The ESMValTool is released under the Apache License, version 2.0 and citation of the ESMValTool paper ("Software Documentation Paper") is kindly requested upon use alongside with the software doi (doi:10.17874/ac8548f0315) and version number:

  * Eyring et al., ESMValTool (v1.0) -- a community diagnostic and performance metrics tool for routine evaluation of Earth System Models in CMIP, Geosci. Model Dev., 9, 1747-1802, 2016.*

Besides the above citation, users are kindly asked to register any journal articles (or other scientific documents) that use the software at the ESMValTool webpage (http://www.esmvaltool.org/). Citing the Software Documentation Paper and registering your paper(s) will serve to document the scientific impact of the Software, which is of vital importance for securing future funding. You should consider this an obligation if you have taken advantage of the ESMValTool, which represents the end product of considerable effort by the development team.

**The ESMValTool is developed in a version controlled repository (see Section** :numref:`git_repository` **for details).** In addition to using the software, we would therefore like to encourage the community to join the Software Development Team and to contribute additional diagnostics and performance metrics or other software improvements. Contributing back the new diagnostics and performance metrics or other software improvements will help to enhance the capability of the Software, which is of vital importance for securing future funding. You should consider this an obligation if you have taken advantage of the Software, which represents a product of considerable effort by the development team.

Interested developers are welcome to contact the core development team (see Section :numref:`core_dev_team`).

Software installation
=====================

The ESMValTool can be downloaded from GitHub (for details see Section :numref:`git_repository`) to any local directory. While the ESMValTool itself does not need to be installed besides downloading/copying the ESMValTool directories to a local folder, it relies on specific software to be available on your system. Please see Section :numref:`prerequisites` for details.

Verification of the installation
================================

Once you have ESMValTool installed you can verify your installation following the out-of-the-box steps listed below. These tests will let you execute a few simplified namelists that will verify that the dependencies for the general control flow of ESMValTool are in place and working properly. The tests will not verify more specific dependencies used by some Python and R diagnostics, such dependencies will have to be installed separately. Test procedure:

.. code:: bash

   1. <INSTALL ESMValTool and dependencies>
   2. <cd INTO YOUR INSTALLATION>
   3. wget http://goo.gl/ciHCsO -O test-data.tar
   4. tar xf test-data.tar
   5. wget http://goo.gl/A7pPEz -O test-nml.tar
   6. tar xf test-nml.tar
   7. ./main.py test-nml/namelist_SAMonsoon-and-WAMonsoon.xml
   8. ./main.py test-nml/namelist_SAMonsoon-pr-with-MPI.xml
   9. ./main.py test-nml/namelist_SAMonsoon-pr-with-TRMM.xml
   10. ./main.py test-nml/namelist_SAMonsoon-pr.xml

For each of step 7-10, manually verify that no errors were reported (standard out) and check that diagnostic output figures have been produced in the subfolders "work/plot_*".

