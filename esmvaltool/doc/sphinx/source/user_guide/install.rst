.. _installation_guide:

*********************
Installing ESMValTool
*********************

ESMValTool 2.0 requires a Unix(-like) operating system and Python 2.7+ or 3.6+.

The ESMValTool supports three different installation methods:

* Installation through Conda package manager (see https://www.continuum.io/ )

* Deployment through a Docker container (see https://www.docker.com/ )

* From the source code available at https://github.com/ESMValGroup/ESMValTool

The next sections will detail the procedure to install ESMValTool from each of
this channels


The ESMValTool requires a Unix(-like) operating system


Conda installation
==================

A conda package will be available after the release of ESMValTool 2.0. As
conda has no support for beta releases, it is unlikely that conda packages for
any of the pre-release versions will be provided


Docker installation
===================

.. warning::
    Docker section to be added


Install from source
===================


Obtaining the source code
-------------------------

The ESMValTool source code is available on a public GitHub repository:
https://github.com/ESMValGroup/ESMValTool

The easiest way to obtain it is to clone the repository:

.. code-block:: bash

    git clone https://github.com/ESMValGroup/ESMValTool.git

By default, this command will create a folder called ESMValTool containing the
source code of the tool

.. attention::
    The newly created clone of the git repository will point by default
    to the master branch. To change to another branch or release execute:
    git checkout origin/$BRANCH_OR_RELEASE_NAME, i.e git checkout origin/2.0.0

GitHub also allows to download the source code in as a tar.gz or zip file. If
you choose to use this option, download the compressed file and extract its
contents at the desired location.


Prerequisites
-------------
The ESMValTool requires a Unix(-like) operating system

It is strongly recommended to use conda to manage ESMValTool dependencies.
For a minimal conda installation go to https://conda.io/miniconda.html . To
simplify the process, an environment definition file is provided within the
repository (environment.yml in theroot folder)

ESMValTool's conda environment file can also be used as a requirements list
for those cases in which a conda installation is not possible or advisable.
From now on, we will asume that the installation is going to be done through
conda.

Ideally, you should want to create a conda environment for ESMValTool, so it is
independent from any other Python tools present in the system.

To create a environment using Python 3.x

.. code-block:: bash

    conda env create --file ESMValTool/environment.yml python=3


To create a environment using Python 2.x

.. code-block:: bash

    conda env create --file ESMValTool/environment.yml python=2

The environment is called "esmvaltool" by default, but it is possible to use
the option -n $(ENVIRONMENT_NAME) to use a custom name. You can activate the
environment using the command

.. code-block:: bash

    source activate esmvaltool

It is also possible to update an existing environment from the environment
file. This can be very useful when updating an older installation of ESMValTool

.. code-block:: bash

    conda env update --file environment.yml --name $(ENVIRONMENT_TO_UPDATE)



Software installation
---------------------

Once all prerequesites are fullfilled, ESMValTool 2.0 can be installed using
the following command

.. code-block:: bash

    python ESMValTool/setup.py

.. attention::
    From now on, we assume that the conda environment for ESMValTool is
    activated

The next step is to check that the installation works properly.
To do this, first activate the environment and then run the tool with --version.

.. code-block:: bash

    esmvaltool --version

For a more complete installation verification, run the automated tests and
confirm that no errors are reported.

.. code-block:: bash

    python ESMValTool/setup.py test

