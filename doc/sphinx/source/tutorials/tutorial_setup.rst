.. _tutorial_setup:

ESMValTool Tutorial: Preparations for participating in the tutorial
===================================================================

This page includes some information on how to prepare for participating in this
tutorial.

.. admonition:: Prerequisites
   :class: important

   *Minimal requirements:*

   * Basic understanding of your preferred command line interface (for example
     a bash terminal)
   * Access to CMIP data

   *Optional, but useful:*

   * Basic understanding of git
   * Access to a suitable computing system (for example CEDA-Jasmin or
     DKRZ-Mistral)
   * GitHub account


Command line & git tutorials
-----------------------------

We typically use the command line to interact with ESMValTool. While most of us
are likely to have experience with the command line, novices may want to work
through the Software Carpentry Unix Shell course.

* Command line: `Software Carpentry Shell Novice
  <https://swcarpentry.github.io/shell-novice/>`_

Git is a distributed version-control system for tracking changes in source code
during software development. It is how we distribute, share, and manage the
ESMValTool code.

* git: `Software Carpentry Git Novice
  <https://swcarpentry.github.io/git-novice/>`_


Access to CMIP and observational data and a suitable compute cluster
---------------------------------------------------------------------

To complete this tutorial and use ESMValTool, you will need access to data in a
reasonable format. Some data will be provided, but there is simply too much
data available for your tutors to make it all available directly.

ESMValTool may be run on multiple platforms, from your local machine to large
computing clusters. The best option is to use a computing cluster with an
`Earth System Grid Federation (ESGF) node <https://esgf.llnl.gov/nodes.html>`_.
The benefit of using a compute cluster with an ESGF node is that the
`Coupled Model Intercomparison Project (CMIP)
<https://en.wikipedia.org/wiki/Coupled_Model_Intercomparison_Project>`_ is
locally stored on disk and accessible directly by the tool. Similarly,
observational data would also be available at these sites.

The ESGF also hosts observations for Model Intercomparison Projects (obs4MIPs)
and reanalyses data (ana4MIPs).

Here are a few options for compute clusters with ESGF nodes:

* :ref:`ceda-jasmin`
* :ref:`dkrz`

For more information see:

* `CMIP5 <https://pcmdi.llnl.gov/mips/cmip5/index.html>`_ and
  `CMIP6 <https://pcmdi.llnl.gov/CMIP6/Guide/dataUsers.html>`_ data obey the
  `CF conventions <http://cfconventions.org/>`_. Available variables can be
  found under the `CMIP5 data request
  <https://pcmdi.llnl.gov/mips/cmip5/docs/standard_output.pdf?id=28>`_ and the
  `CMIP6 Data Request <http://clipc-services.ceda.ac.uk/dreq/index.html>`_.
* List of all `CMIP named variables
  <http://clipc-services.ceda.ac.uk/dreq/index/CMORvar.html>`_.
* List of all `ESGF nodes <https://esgf.llnl.gov/nodes.html>`_.
* A good `tutorial
  <https://esgf.github.io/esgf-user-support/user_guide.html#data-search-and-download>`_
  on how to search and download CMIP data from ESGF nodes.
* `Exploring climate model data
  <https://climate4impact.eu/impactportal/data/esgfsearch.jsp>`_ on
  infrastructure for the European network for Earth system modelling.


.. _ceda-jasmin:

CEDA-Jasmin (UK)
~~~~~~~~~~~~~~~~

Please skip this section if you are not going to use JASMIN and continue to the
:ref:`github-account-advanced` section.

If you do not already have an account on JASMIN, then request an account as
soon as possible. Please follow the steps below (also documented in detail
here: `Essential steps to gain login access to JASMIN
<https://help.jasmin.ac.uk/docs/getting-started/get-started-with-jasmin/#essential-steps>`_).

* Generate an `SSH key
  <https://help.jasmin.ac.uk/docs/getting-started/generate-ssh-key-pair>`_
* Get a `JASMIN portal account
  <https://help.jasmin.ac.uk/docs/getting-started/get-jasmin-portal-account>`_
* Get a `jasmin-login account
  <https://help.jasmin.ac.uk/docs/getting-started/get-login-account>`_

Note that the JASMIN portal is only an account for the web interface. A
jasmin-login account is also required.

Also note that if you are working from home, JASMIN may not be directly
accessible from your home. You may need to use SSH to connect to a machine in
your institute and then on to JASMIN. Please test your connection before the
tutorial.

Here are some further, general resources for `getting started with JASMIN
<https://help.jasmin.ac.uk/article/189-get-started-with-jasmin>`_.

Access to data on JASMIN
^^^^^^^^^^^^^^^^^^^^^^^^

Please request access to the working groups:

* `esmeval working group
  <https://accounts.jasmin.ac.uk/services/group_workspaces/esmeval>`_
* `CMIP5 data
  <https://services.ceda.ac.uk/cedasite/resreg/application?attributeid=cmip5_research>`_

Once you have access to the data archive on CEDA, make sure to link your CEDA
and JASMIN accounts. This can be done by checking the link to CEDA box on `your
JASMIN profile page <https://accounts.jasmin.ac.uk/account/profile/>`_.

The linking may take a few hours to take effect and is necessary for you to
access the BADC archives via JASMIN. Some CMIP5 data sets such as MIROC are not
accessible by default and special permission has to be requested to access them
via `the CEDA catalogue page <https://catalogue.ceda.ac.uk/>`_.

Test your setup
^^^^^^^^^^^^^^^

Log into jasmin-login:

.. code-block:: bash

   ssh -X JASMIN-USERNAME@login1.jasmin.ac.uk

Then log into the sci1 machine:

.. code-block:: bash

   ssh -X jasmin-sci1

Can you see the following locations:

.. code-block:: bash

   ls /gws/ssde/j25a/esmeval/obsdata-v2/
   ls /badc/cmip5/data/cmip5/output1/MOHC/HadGEM2-ES
   ls /badc/cmip6/data/CMIP6/CMIP/*/*/historical/r1i1p1f?/Omon/[ts]os/gn/latest/*.nc

Note that JASMIN is only open to certain locations (mostly universities and
research centres). You may need a VPN if you wish to connect from your home
network.

Congratulations! Please continue to the :ref:`github-account-advanced` section
next.


.. _dkrz:

DKRZ (Germany)
~~~~~~~~~~~~~~

Please skip this section if you are not going to use DKRZ and continue to the
:ref:`github-account-advanced` section.

If you do not already have an account at the DKRZ, then `register
<https://luv.dkrz.de/projects/newuser/>`_ as soon as possible. You could find a
short introduction on how to get started at DKRZ `here
<https://docs.dkrz.de/doc/getting_started/index.html>`_.

There is also a `user manual <https://docs.dkrz.de/doc/levante/index.html>`_
for Levante, which is DKRZ's current supercomputer.

Join a project
^^^^^^^^^^^^^^

To use the resources on DKRZ you have to join a project. One option is to join
an existing project by logging into `https://luv.dkrz.de/
<https://luv.dkrz.de/>`_ with your account and selecting ``Join existing
project``. Once you are accepted by the manager of your chosen project, your
web account will be turned into a full LDAP account which will allow you to log
into and use the DKRZ resources. If you do not have access to an existing
project, another option for you would be to apply for resources at DKRZ. Here
are some instructions on `how to apply for resources
<https://www.dkrz.de/services/bereitstellung-von-rechenleistung?set_language=en&cl=en>`_.

Access to data on DKRZ
^^^^^^^^^^^^^^^^^^^^^^

CMIP5 and CMIP6 data are available in these directories:

* CMIP5: ``/work/kd0956/CMIP5/data/cmip5/output1/``
* CMIP6: ``/work/ik1017/CMIP6/data/CMIP6/CMIP/``

Test your setup
^^^^^^^^^^^^^^^

Log into Levante (DKRZ):

.. code-block:: bash

   ssh -X user-account@levante.dkrz.de

Additional information
^^^^^^^^^^^^^^^^^^^^^^

Login nodes are for compiling and job submission only. For all other tasks, you
can use the `interactive partition
<https://docs.dkrz.de/doc/levante/running-jobs/partitions-and-limits.html#interactive>`_
or start an `interactive session
<https://docs.dkrz.de/doc/levante/data-processing-on-levante.html#>`_.

Data storage:

* Personal data: home directory (30 GiB)
* Project data: ``/work/project_id/username``
* Temporary data: scratch directory on ``/scratch/*/username`` is automatically
  deleted after 14 days (15 TiB). Please use this directory for all your
  testing; do not use the work directory for tests. See also `this page
  <https://docs.dkrz.de/doc/levante/file-systems.html>`_.

Running batch jobs: information and examples on the SLURM job scheduling system
at DKRZ can be found `here
<https://docs.dkrz.de/doc/levante/running-jobs/index.html>`_.

Congratulations! Please continue to the :ref:`github-account-advanced` section
next.


Using your own machine
----------------------

Please skip this section if you are not going to use ESMValTool on your local
machine and continue to the :ref:`github-account-advanced` section.

If you are planning on running ESMValTool on your own machine, please make sure
that you are able to download CMIP data and that you have a few GB of space
available to install conda and ESMValTool, but also enough to make a copy of
some data (about 125 MB) needed for this tutorial.

You can use ESMValTool to automatically download data needed for test recipes.
Please see the Configuration episode or the configuration file documentation
for more information. This is the recommended option as it has the advantage
that data is stored in subdirectories, and features such as wildcards and
recording the version of the data will work automatically.

Alternatively, you can run the following command using `wget
<https://en.wikipedia.org/wiki/Wget>`_:

.. code-block:: shell

   wget --no-clobber --input-file \
     https://github.com/ESMValGroup/ESMValTool_Tutorial/raw/main/data/dataset.urls \
     --directory-prefix $HOME/esmvaltool_tutorial/data/


.. _github-account-advanced:

GitHub account (Advanced)
-------------------------

You do not need a GitHub account to participate in the tutorial. However, if
you want to raise an issue, contribute to the discussions, or share your code,
please `create a GitHub account <https://github.com/>`_.

To learn how to use GitHub, please have a look at this `introduction to GitHub
<https://lab.github.com/githubtraining/introduction-to-github>`_.

You may hear a few of the following phrases during the tutorial. Do not be
alarmed; they will make sense eventually.

GitHub issues
~~~~~~~~~~~~~

Issues are GitHub's ticketing system. They allow users and developers to
discuss problems, identify bugs, or make suggestions. Each issue is assigned a
number and will have its own page on GitHub.

Here is an explanation of the `GitHub issues
<https://guides.github.com/features/issues/>`_.

Raising an issue is the act of creating a new issue. If you are asked to raise
an issue, please follow any instructions that you are given, and also make sure
that you read the default issue text.

GitHub pull requests
~~~~~~~~~~~~~~~~~~~~

A GitHub pull request is the act of requesting that a branch is merged with
another branch.

This is an advanced feature of GitHub, and will generally be performed by the
ESMValTool development team.
