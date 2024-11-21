.. include:: ../common.txt

The workflow
============

An overview of the workflow
---------------------------

The |RTW| performs the following steps:

``install_env_file``
  :Description:
     Copies the environment file for |ESMValTool|, based on the ``SITE``
     provided
  :Runs on:
     Localhost
  :Executes:
     The ``cp`` command from the |Rose| app
  :Details:
     Runs once at the start of the workflow

``get_esmval``
  :Description:
     Either clones the latest versions of |ESMValTool| and |ESMValCore| from
     GitHub, or gets the latest container image from DockerHub and converts to
     a singularity image, depending on ``SITE``
  :Runs on:
     Localhost (if cloning), or ``COMPUTE`` (if getting container), which
     depends on the ``SITE``; on JASMIN, the ``get_esmval`` jobs will run on
     LOTUS
  :Executes:
     The ``clone_latest_esmval.sh`` script (if cloning), or a
     ``singularity build`` command (if getting container) from the |Rose| app
  :Details:
     Runs at the start of each cycle

``configure``
  :Description:
     Creates the |ESMValTool| user configuration file and validates it
  :Runs on:
     Localhost
  :Executes:
     The ``configure.py`` script from the |Rose| app
  :Details:
     Runs each cycle after ``get_esmval`` has completed

``process``
  :Description:
     Runs the requested recipes using |ESMValTool|
  :Runs on:
     ``COMPUTE``, which depends on the ``SITE``; at the Met Office, the
     ``process`` jobs will run on SPICE
  :Executes:
     The |ESMValTool| command line script from the |Rose| app
  :Details:
     Runs each cycle for every recipe defined in the |RTW| after ``configure``
     has completed

``compare``
  :Description:
     Compares the output from the ``process`` job with |KGOs|
  :Runs on:
     ``COMPUTE``, which depends on the ``SITE``; at the Met Office, the
     ``compare`` jobs will run on SPICE
  :Executes:
     The :ref:`compare.py <compare_recipe_runs>` script from |ESMValTool|
     from the |Rose| app
  :Details:
     Runs each cycle for every recipe defined in the |RTW| after ``process``
     has completed

Design considerations
---------------------

Portability
~~~~~~~~~~~

The |RTW| is portable; site-specific information can be found in the ``site``
and ``opt`` directories within the |RTW|. The files required are:

``site/<site>.cylc``
  Contains task definitions specific to the ``SITE``, for example, ``COMPUTE``

``site/<site>-env``
  Contains details on how to set up the environment for ESMValTool at the
  ``SITE``

``opt/rose-suite-<site>.conf``
  Contains configuration items specific to the ``SITE``, including ``SITE``
  and where to find the required datasets

Metadata
~~~~~~~~

The |RTW| uses Rose metadata. Every item defined in the suite configuration
file (``rose-suite.conf``) will have an entry in the main metadata
configuration file (``meta/rose-meta.conf``).

Resources
~~~~~~~~~

The resources used by the ``process`` jobs are defined in the
``site/<site>.cylc`` file, allowing the jobs to be configured by ``SITE`` as
well as by recipe. This ensures only the required resources are requested when
running each of the ``process`` jobs.
