.. _installation_guide:

*********************
Installing ESMValTool
*********************

ESMValTool 2.0 requires a Unix(-like) operating system and Python 2.7+ or 3.6+.
Python 2.7+ will be discontinued in the near future, so we encourage you to use
Python 3.6+ if possible.

The ESMValTool supports three different installation methods:

* Installation through Conda package manager (see https://www.continuum.io/);

* Deployment through a Docker container (see https://www.docker.com/);

* From the source code available at https://github.com/ESMValGroup/ESMValTool.

The next sections will detail the procedure to install ESMValTool for each of
this methods.


Conda installation
==================

A conda package will be available after the release of ESMValTool 2.


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

The easiest way to obtain it is to clone the repository using git
(see https://git-scm.com/):

.. code-block:: bash

    git clone https://github.com/ESMValGroup/ESMValTool.git

By default, this command will create a folder called ESMValTool containing the
source code of the tool.

.. attention::
    The newly created clone of the git repository will point by default
    to the master branch. To change to another branch or release execute:
    git checkout origin/$BRANCH_OR_RELEASE_NAME, i.e git checkout origin/2.0.0.

GitHub also allows to download the source code in as a tar.gz or zip file. If
you choose to use this option, download the compressed file and extract its
contents at the desired location.


Prerequisites
-------------

It is strongly recommended to use conda to manage ESMValTool dependencies.
For a minimal conda installation go to https://conda.io/miniconda.html. To
simplify the process, an environment definition file is provided within the
repository (``environment.yml`` in the root folder).

Note that the standard conda installation has some issues with the ``csh``/``tcsh``
login shell. If you are using such shell, do not prepend the install
location (``<prefix>``) to PATH in your ``.tcshrc`` file (as suggested by
the standard installation procedure which assumes ``bash``). Instead, add
the following line to your ``.cshrc``/``.tcshrc`` file: 

.. code-block:: bash

    source <prefix>/etc/profile.d/conda.csh

The ESMValTool conda environment file can also be used as a requirements list
for those cases in which a conda installation is not possible or advisable.
From now on, we will assume that the installation is going to be done through
conda.

Ideally, you should want to create a conda environment for ESMValTool, so it is
independent from any other Python tools present in the system.

To create a environment using Python 3.x

.. code-block:: bash

    conda create --name esmvaltool python=3
    conda env update --name esmvaltool --file ESMValTool/environment.yml

To create a environment using Python 2.x

.. code-block:: bash

    conda create --name esmvaltool python=2
    conda env update --name esmvaltool --file ESMValTool/environment.yml

The environment is called ``esmvaltool`` by default, but it is possible to use
the option -n $(ENVIRONMENT_NAME) to define a custom name. If you are using the
``bash`` shell, you can activate the environment using the command:

.. code-block:: bash

    source activate esmvaltool

while for the ``csh``/``tcsh`` you need to use:

.. code-block:: bash
    
    conda activate esmvaltool

It is also possible to update an existing environment from the environment
file. This can be very useful when updating an older installation of ESMValTool:

.. code-block:: bash

    conda env update --file environment.yml --name $(ENVIRONMENT_TO_UPDATE)

.. attention::
    From now on, we assume that the conda environment for ESMValTool is
    activated.

Software installation
---------------------

Once all prerequisites are fulfilled, ESMValTool 2.0 can be installed using
the following command:

.. code-block:: bash

    python ESMValTool/setup.py install

The next step is to check that the installation works properly.
To do this, run the tool with --version:

.. code-block:: bash

    esmvaltool --version

If everything was installed properly, ESMValTool should have printed the
version number at the console and exited.

For a more complete installation verification, run the automated tests and
confirm that no errors are reported:

.. code-block:: bash

    python ESMValTool/setup.py test

