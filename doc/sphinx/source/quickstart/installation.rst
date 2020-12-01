.. _install:

************
Installation
************

ESMValTool 2.0 requires a Unix(-like) operating system and Python 3.6+.

The ESMValTool supports five different installation methods:

* Installation through Conda package manager (see https://www.continuum.io);

* Installation with Pip and Conda package manager (see https://pypi.org);

* Deployment through a Docker container (see https://www.docker.com);

* Deployment through a Singularity container (see https://sylabs.io/guides/latest/user-guide/);

* From the source code available at https://github.com/ESMValGroup/ESMValTool.

The next sections will detail the procedure to install ESMValTool for each of
this methods.


Conda installation
==================

In order to install the `Conda <https://docs.conda.io>`_ package, you will need
both Conda and Julia pre-installed, this is because Julia cannot be installed
from Conda.
For a minimal conda installation (recommended) go to https://conda.io/miniconda.html.
It is recommended that you always use the latest version of Conda, as problems
have been reported when trying to use older versions.
Installation instructions for Julia can be found on the
`Julia installation instructions page <https://julialang.org/downloads/platform/>`_.

Once you have installed the above prerequisites, you can install ESMValTool by running:

.. code-block:: bash

    conda install esmvaltool -c esmvalgroup -c conda-forge

Here ``conda`` is the executable calling the Conda package manager to install
``esmvaltool`` and the ``-c`` flag specifies the Conda software channels in which the
``esmvaltool`` package and its dependencies can be found.

It is also possible to create a new
`Conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments>`_
and install ESMValTool into it with a single command:

.. code-block:: bash

    conda create --name esmvaltool -c esmvalgroup -c conda-forge esmvaltool

Don't forget to activate the newly created environment after the installation:

.. code-block:: bash

    conda activate esmvaltool

Of course it is also possible to choose a different name than ``esmvaltool`` for the environment.

.. note::

	  Creating a new Conda environment is often much faster and more reliable than trying to update an existing Conda environment.

Installation of subpackages
---------------------------

The diagnostics bundled in ESMValTool are scripts in four different programming languages: Python, NCL, R, and Julia.

There are four language specific packages available:

* ``esmvaltool-julia``
* ``esmvaltool-ncl``
* ``esmvaltool-python``
* ``esmvaltool-r``

The main ``esmvaltool`` package contains all four subpackages listed above.

If you only need to run a recipe with diagnostics in some of these languages, it is possible to install only the dependencies needed to do just that.

* The diagnostic script(s) used in each recipe, are documented in :ref:`recipes`. The extension of the diagnostic script can be used to see in which language a diagnostic script is written.
* Some of the CMORization scripts are written in Python, while others are written in  NCL. Therefore, both ``esmvaltool-pyhon`` and ``esmvaltool-ncl`` need to be installed in order to be able to run all CMORization scripts.

For example, to only install support for diagnostics written in Python and NCL, run

.. code-block:: bash

    conda install esmvaltool-python esmvaltool-ncl -c esmvalgroup -c conda-forge

Note that it is only necessary to install Julia prior to the conda installation if you are going to install the ``esmvaltool-julia`` package.

Note that the ESMValTool source code is contained in the ``esmvaltool-python`` package, so this package will always be installed as a dependency if you install one or more of the packages for other languages.

There is also a lesson available in the 
`ESMValTool tutorial <https://esmvalgroup.github.io/ESMValTool_Tutorial/>`_
that describes the installation of the ESMValTool in more detail. It can be found
`here <https://esmvalgroup.github.io/ESMValTool_Tutorial/02-installation/index.html>`_.

Pip installation
================

It is also possible to install ESMValTool from `PyPI <https://pypi.org/project/ESMValTool/>`_.
However, this requires first installing dependencies that are not available on PyPI in some other way.
By far the easiest way to install these dependencies is to use conda_.
For a minimal conda installation (recommended) go to https://conda.io/miniconda.html.

After installing Conda, download
`the file with the list of dependencies <https://raw.githubusercontent.com/ESMValGroup/ESMValTool/master/environment.yml>`_:

.. code-block:: bash

    wget https://raw.githubusercontent.com/ESMValGroup/ESMValTool/master/environment.yml

and install these dependencies into a new conda environment with the command

.. code-block:: bash

    conda env create --name esmvaltool --file environment.yml

Finally, activate the newly created environment

.. code-block:: bash

    conda activate esmvaltool

and install ESMValTool as well as any remaining Python dependencies with the command:

.. code-block:: bash

    pip install esmvaltool

If you would like to run Julia diagnostic scripts, you will also need to
`install Julia <https://julialang.org/downloads/platform/>`_ and the Julia dependencies:

.. code-block:: bash

    esmvaltool install Julia

If you would like to run R diagnostic scripts, you will also need to install the R
dependencies:

.. code-block:: bash

    esmvaltool install R

Docker installation
===================

ESMValTool is also provided through `DockerHub <https://hub.docker.com/u/esmvalgroup/>`_
in the form of docker containers.
See https://docs.docker.com for more information about docker containers and how to
run them.

You can get the latest release with

.. code-block:: bash

   docker pull esmvalgroup/esmvaltool:stable

If you want to use the current master branch, use

.. code-block:: bash

   docker pull esmvalgroup/esmvaltool:latest

To run a container using those images, use:

.. code-block:: bash

   docker run esmvalgroup/esmvaltool:stable --help

Note that the container does not see the data or environmental variables
available in the host by default. You can make data available with
``-v /path:/path/in/container`` and environmental variables with ``-e VARNAME``.

For example, the following command would run a recipe

.. code-block:: bash

   docker run -e HOME -v "$HOME":"$HOME" -v /data:/data esmvalgroup/esmvaltool:stable run examples/recipe_python.yml

with the environmental variable ``$HOME`` available inside the container and
the data in the directories ``$HOME`` and ``/data``, so these can be used to
find the configuration file, recipe, and data.

It might be useful to define a `bash alias
<https://opensource.com/article/19/7/bash-aliases>`_
or script to abbreviate the above command, for example

.. code-block:: bash

	 alias esmvaltool="docker run -e HOME -v $HOME:$HOME -v /data:/data esmvalgroup/esmvaltool:stable"

would allow using the ``esmvaltool`` command without even noticing that the
tool is running inside a Docker container.


Singularity installation
========================

Docker is usually forbidden in clusters due to security reasons. However,
there is a more secure alternative to run containers that is usually available
on them: `Singularity <https://sylabs.io/guides/3.0/user-guide/quick_start.html>`_.

Singularity can use docker containers directly from DockerHub with the
following command

.. code-block:: bash

   singularity run docker://esmvalgroup/esmvaltool:stable run examples/recipe_python.yml

Note that the container does not see the data available in the host by default.
You can make host data available with ``-B /path:/path/in/container``.

It might be useful to define a `bash alias
<https://opensource.com/article/19/7/bash-aliases>`_
or script to abbreviate the above command, for example

.. code-block:: bash

	 alias esmvaltool="singularity run -B $HOME:$HOME -B /data:/data docker://esmvalgroup/esmvaltool:stable"

would allow using the ``esmvaltool`` command without even noticing that the
tool is running inside a Singularity container.

Some clusters may not allow to connect to external services, in those cases
you can first create a singularity image locally:

.. code-block:: bash

   singularity build esmvaltool.sif docker://esmvalgroup/esmvaltool:stable

and then upload the image file ``esmvaltool.sif`` to the cluster.
To run the container using the image file ``esmvaltool.sif`` use:

.. code-block:: bash

   singularity run esmvaltool.sif run examples/recipe_python.yml


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

It is also possible to work in one of the ESMValTool private repositories, e.g.:

.. code-block:: bash

    git clone https://github.com/ESMValGroup/ESMValTool-private.git

By default, this command will create a folder called ESMValTool containing the
source code of the tool.

GitHub also allows one to download the source code in as a tar.gz or zip file. If
you choose to use this option, download the compressed file and extract its
contents at the desired location.


Prerequisites
-------------

It is recommended to use conda to manage ESMValTool dependencies.
For a minimal conda installation go to https://conda.io/miniconda.html. To
simplify the installation process, an environment definition file is provided
in the repository (``environment.yml`` in the root folder).

.. attention::
    Some systems provide a preinstalled version of conda (e.g., via the module environment).
    However, several users reported problems when installing NCL with such versions. It is
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

where ``<prefix>`` is the install location of your anaconda or miniconda
(e.g. ``/home/$USER/anaconda3`` or ``/home/$USER/miniconda3``).


.. note::
    Note that during the installation, conda will ask you
    if you want the installation to be automatically sourced from your
    ``.bashrc`` or ``.bash-profile`` files; if you answered yes, then conda
    will write bash directives to those files and every time you get to your
    shell, you will automatically be inside conda's ``(base)`` environment.
    To deactivate this feature, look for the ``# >>> conda initialize >>>``
    code block in your ``.bashrc`` or ``.bash-profile`` and comment the whole block out.


The ESMValTool conda environment file can also be used as a requirements list
for those cases in which a conda installation is not possible or advisable.
From now on, we will assume that the installation is going to be done through
conda.

Ideally, you should create a conda environment for ESMValTool, so it is
independent from any other Python tools present in the system.

Note that it is advisable to update conda to the latest version before
installing ESMValTool, using the command (as mentioned above)

.. code-block:: bash

    conda update --name base conda

To create an environment, go to the directory containing the ESMValTool source
code (called ESMValTool if you did not choose a different name) and run

.. code-block:: bash

    conda env create --name esmvaltool --file environment.yml

This installs the ESMValCore package from conda as a dependency.

The environment is called ``esmvaltool`` by default, but it is possible to use
the option ``--name SOME_ENVIRONMENT_NAME`` to define a custom name. You should then activate
the environment using the command:

.. code-block:: bash

    conda activate esmvaltool

It is also possible to update an existing environment from the environment
file. This may be useful when updating an older installation of ESMValTool:

.. code-block:: bash

    conda env update --name esmvaltool --file environment.yml

but if you run into trouble, please try creating a new environment.

.. attention::
    From now on, we assume that the conda environment for ESMValTool is
    activated.

Software installation
---------------------

Once all prerequisites are fulfilled, ESMValTool can be installed by running
the following commands in the directory containing the ESMValTool source code
(called ESMValTool if you did not choose a different name):

.. code-block:: bash

    pip install -e '.[develop]'

If you would like to run Julia diagnostic scripts, you will also need to
`install Julia <https://julialang.org/downloads/platform/>`_ and the Julia dependencies:

.. code-block:: bash

    esmvaltool install Julia

If you would like to run R diagnostic scripts, you will also need to install the R
dependencies. Install the R dependency packages:

.. code-block:: bash

    esmvaltool install R

The next step is to check that the installation works properly.
To do this, run the tool with:

.. code-block:: bash

    esmvaltool --help

If everything was installed properly, ESMValTool should have printed a
help message to the console.
