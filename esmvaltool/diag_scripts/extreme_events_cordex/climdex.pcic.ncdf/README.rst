What is climdex.pcic.ncdf?
=====================

* `climdex.pcic.ncdf` is a companion library for `climdex.pcic` which helps in using NetCDF input grids and writing to NetCDF output files when computing the `27 core indices of extreme climate`_. The code allows for parallel computation of indices using either a SOCK or MPI cluster. It was written for the `R statistical programming language`_ by the `Pacific Climate Impacts Consortium`_.

.. _27 core indices of extreme climate: http://etccdi.pacificclimate.org/list_27_indices.shtml
.. _R statistical programming language: http://www.r-project.org/
.. _Pacific Climate Impacts Consortium: http://pacificclimate.org/

Getting Help
============

New to programming or to R?
---------------------------

* Read the the `Software Carpentry`_  `Programming in R`_ lessons
* Read one of the man `R Manuals`_.
* Attend an `R Users Group`_ meeting.

.. _Software Carpentry: http://software-carpentry.org/index.html
.. _Programming in R: http://software-carpentry.org/v5/novice/r/index.html
.. _R Manuals: http://cran.r-project.org/manuals.html
.. _R Users Group: http://r-users-group.meetup.com/

Looking for code?
-----------------

* Get the latest `climdex.pcic.ncdf release from our website`_.
* Explore the `development repository`_.
* Install it with devtools ::

    > library(devtools)
    > install_github('pacificclimate/climdex.pcic.ncdf', ref='release')

.. _climdex.pcic.ncdf release from our website: http://www.pacificclimate.org/sites/default/files/climdex.pcic_.ncdf_0.5-4.tar_.gz
.. _development repository: https://github.com/pacificclimate/climdex.pcic.ncdf/

Need help using the package?
----------------------------

* Read the manual ::

    > library(climdex.pcic.ncdf)
    Loading required package: PCICt
    > ?climdex.pcic.ncdf

* Create a `new issue`_ on the `package issue tracker`_ and label it "help wanted"[1]_.

.. _new issue: https://github.com/pacificclimate/climdex.pcic.ncdf/issues/new

Want to contribute?
-------------------

* To report a bug in pcic.climdex use the `package issue tracker`_ (after you've read the `bug reporting guide`_).
* To help with development read through the `contributor's guide`_

.. _bug reporting guide: https://github.com/pacificclimate/climdex.pcic.ncdf/blob/master/CONTRIBUTING.rst#bug-reports
.. _package issue tracker: https://github.com/pacificclimate/climdex.pcic.ncdf/issues
.. _contributor's guide: https://github.com/pacificclimate/climdex.pcic.ncdf/blob/master/CONTRIBUTING.rst

Still need help?
----------------

* Contact climate@uvic.ca and let us know what we can do.

.. [1] Please know that the pool of people who can provide support for the package is extremely small and time is limited.  We don't necessarily have the capacity for long, open-ended user support. If you keep your questions short, specific and direct, there's a greater probability that someone will take on the ticket.
