.. _installation_guide:

*********************
Installing ESMValTool
*********************

ESMValTool 2.0 requires a Unix(-like) operating system and Python 3.6+.

The ESMValTool supports three different installation methods:

* Installation through Conda package manager (see https://www.continuum.io);

* Deployment through a Docker container (see https://www.docker.com);

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
(see https://git-scm.com/). To clone the public repository:

.. code-block:: bash

    git clone https://github.com/ESMValGroup/ESMValTool.git

It is also possible to work in one of the ESMValTool private repository, e.g.:

.. code-block:: bash

    git clone https://github.com/ESMValGroup/ESMValTool-private.git

By default, this command will create a folder called ESMValTool containing the
source code of the tool.

.. attention::
    The newly created clone of the git repository will point to the default
    branch. To change to another branch or release execute:
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

.. attention::
    Some systems provides a preinstalled version of conda (e.g., via the module environment).
    Several users however reported problems when installing NCL with such versions. It is
    therefore preferable to use a local, fully user-controlled conda installation.
    Using an older version of conda can also be a source of problems, so if you have conda
    installed already, make sure it is up to date by running ``conda update -n base conda``.

To enable the ``conda`` command, please source the appropriate configuration file
from your ``~/.bashrc``  file:

.. code-block:: bash

    source <prefix>/etc/profile.d/conda.sh

or ``~/.cshrc``/``~/.tcshrc`` file:

.. code-block:: bash

    source <prefix>/etc/profile.d/conda.csh

The ESMValTool conda environment file can also be used as a requirements list
for those cases in which a conda installation is not possible or advisable.
From now on, we will assume that the installation is going to be done through
conda.

Ideally, you should create a conda environment for ESMValTool, so it is
independent from any other Python tools present in the system.

Note that it is advisable to update conda to the latest version before
installing ESMValTool. Using the following command:

.. code-block::bash

    conda update -n base conda

To create an environment

.. code-block:: bash

    conda create --name esmvaltool python=3
    conda env update --name esmvaltool --file ESMValTool/environment.yml

The environment is called ``esmvaltool`` by default, but it is possible to use
the option -n $(ENVIRONMENT_NAME) to define a custom name. You can activate
the environment using the command:

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

    cd ESMValTool
    python setup.py install

If you would like to run R diagnostic scripts, you will also need to install the R
dependencies. To do that you will need to install gcc, if gcc is not available on your system or it is too old. On Debian based systems, this can be done by running `apt install build-essential gcc`, on managed systems you can often use the `module avail` command to see what compilers are available (note on gcc version: gcc 7.3.0 works well; gcc 8.2.0 is reported to have issues installing the R packages; example of loading gcc7 on the CEDA Jasmin cluster: `module load contrib/gnu/gcc/7.3.0`, similar module loading functionality should be present at your local cluster as well). Install the R dependency packages:

.. code-block:: bash

    Rscript ESMValTool/esmvaltool/install/R/setup.R

The next step is to check that the installation works properly.
To do this, run the tool with --help:

.. code-block:: bash

    esmvaltool --help

If everything was installed properly, ESMValTool should have printed the
version number at the console and exited.

For a more complete installation verification, run the automated tests and
confirm that no errors are reported:

.. code-block:: bash

    python setup.py test --installation

