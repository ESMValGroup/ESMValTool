.. (C) Crown Copyright 2022, the Met Office.

.. include:: ../common.txt

The workflow
============

An overview of the workflow
---------------------------

The |CAW| performs the following steps:

``install_cold``
  :Description:
     Activates the environment for |ESMValTool|, based on the ``SITE`` provided
  :Runs on:
     Localhost
  :Executes:
     The ``install_cold.sh`` script from the |Rose| app
  :Details:
     Runs once at the start of the workflow

``configure``
  :Description:
     Creates and modifies the |ESMValTool| user configuration file
  :Runs on:
     Localhost
  :Executes:
     The ``configure.py`` script from the |Rose| app
  :Details:
     Runs once at the start of the workflow, immediately after the successful
     completion of the ``install_cold`` job

``process``
  :Description:
     Runs the requested recipes using |ESMValTool|
  :Runs on:
     ``COMPUTE``, which depends on the ``SITE``; at the Met Office, the
     ``process`` jobs will run on SPICE
  :Executes:
     The |ESMValTool| command line script
  :Details:
     Runs for every metric defined in the workflow

Design considerations
---------------------

Portability
~~~~~~~~~~~

The |CAW| is portable; site-specific information can be found in the ``site``
and ``opt`` directories within the workflow. The files required are:

``site/<site>.cylc``
  Contains task definitions specific to the ``SITE``, for example, ``COMPUTE``

``site/<site>-env``
  Contains details on how to set up the environment for ESMValTool at the
  ``SITE``

``opt/rose-suite-<site>.conf``
  Contains configuration items specific to the ``SITE``, including ``SITE``

Metadata
~~~~~~~~

The |CAW| uses Rose metadata. Every item defined in the suite configuration
file (``rose-suite.conf``) will have an entry in the main metadata
configuration file (``meta/rose-meta.conf``).

Resources
~~~~~~~~~

The resources used by the ``process`` jobs are defined in the
``site/<site>.cylc`` file, allowing the jobs to be configured by ``SITE`` as
well as by recipe. This ensures only the required resources are requested when
running each of the ``process`` jobs.
