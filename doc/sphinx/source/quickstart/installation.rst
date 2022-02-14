.. _install:

************
Installation
************

.. note::
   ESMValTool now uses `mamba` instead of `conda` for the recommended installation.
   For more information about the change, have a look at :ref:`Move to Mamba<move-to-mamba>`.

ESMValTool 2.0 requires a Unix(-like) operating system and Python 3.7+.

The ESMValTool can be installed in multiple ways.

Recommended installation methods:

* On Linux, please install via the :ref:`mamba package manager<install_with_mamba>` (see https://anaconda.com);

* For MacOSX, please follow separate instructions for :ref:`installation on MacOSX<install_on_macosx>`.

Further options for installation are:

* :ref:`Installation with pip and mamba<install_with_pip>` (see https://pypi.org);

* :ref:`Deployment through a Docker container<install_with_docker>` (see https://www.docker.com);

* :ref:`Deployment through a Singularity container<install_with_singularity>` (see https://sylabs.io/guides/latest/user-guide/);

* :ref:`From the source code<install_from_source>` available at https://github.com/ESMValGroup/ESMValTool;

* :ref:`From pre-installed versions on HPC clusters<install_on_hpc>`.

The next sections will detail the procedure to install ESMValTool through each
of these methods.

Note that there is also a
`Tutorial <https://esmvalgroup.github.io/ESMValTool_Tutorial/>`__
available with more explanations and guidance if the installation instructions
here are too concise.
See `common installation issues`_ if you run into trouble.

.. _install_with_mamba:

Mamba/Conda installation
========================

In order to install the `conda <https://docs.conda.io>`_ package, you will need
mamba pre-installed.

For a minimal mamba installation (recommended) go to
https://mamba.readthedocs.io/en/latest/installation.html.
It is recommended that you always use the latest version of mamba, as problems
have been reported when trying to use older versions.

First download the installation file for
`Linux <https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh>`_
or
`MacOSX <https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh>`_.
After downloading the installation file from one of the links above, execute it
by running (Linux example):

.. code-block:: bash

    bash Mambaforge-Linux-x86_64.sh

and follow the instructions on your screen.
Immediately update mamba after installing it:

.. code-block:: bash

    mamba update --name base mamba

You can check that mamba installed correctly by running

.. code-block:: bash

    which mamba

this should show the path to your mamba executable, e.g.
``~/mambaforge/bin/mamba``.

ESMValTool installation
=======================

Once you have installed the above prerequisites, you can install the entire
ESMValTool package by running:

.. code-block:: bash

    mamba create --name esmvaltool esmvaltool python=3.9

Here ``mamba`` is the executable calling the mamba package manager to install
``esmvaltool``. The reason why we are also specifying ``python=3.9`` is that it
will make it easier for mamba to find a working combination of all required
packages, see `Mamba fails to solve the environment`_ in `common installation issues`_
for an in-depth explanation. Python 3.7 and 3.8 are also supported, in case you
prefer to work with an older version of Python.

This will create a new
`conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments>`_
and install ESMValTool into it with a single command.


.. code-block:: bash

    conda activate esmvaltool

Of course it is also possible to choose a different name than ``esmvaltool`` for the environment.

The next step is to check that the installation works properly.
To do this, run the tool with:

.. code-block:: bash

    esmvaltool --help

If everything was installed properly, ESMValTool should have printed a help
message to the console.

.. note::

    Creating a new conda environment is often much faster and more reliable than
    trying to update an existing conda environment.

Julia installation
------------------

If you want to use the ESMValTool Julia functionality, you will also need to
install Julia. If you are just getting started, we suggest that you
come back to this step later when and if you need it.
To perform the Julia installation, make sure that your conda
environment is activated and then execute

.. code-block:: bash

    mamba install julia

.. _conda subpackages:

Installation of subpackages
---------------------------

The diagnostics bundled in ESMValTool are scripts in four different programming
languages: Python, NCL, R, and Julia.

There are three language specific packages available:

* ``esmvaltool-ncl``
* ``esmvaltool-python``
* ``esmvaltool-r``

The main ``esmvaltool`` package contains all four subpackages listed above. If
you only need to run a recipe with diagnostics in some of these languages, it is
possible to install only the dependencies needed to do just that. The diagnostic
script(s) used in each recipe, are documented in :ref:`recipes`. The extension
of the diagnostic script can be used to see in which language a diagnostic
script is written.

To install support for diagnostics written in Python and NCL into an existing
environment, run

.. code-block:: bash

    mamba install esmvaltool-python esmvaltool-ncl

Some of the CMORization scripts are written in Python, while others are written
in NCL. Therefore, both ``esmvaltool-python`` and ``esmvaltool-ncl`` need to be
installed in order to be able to run all CMORization scripts.

Note that the ESMValTool source code is contained in the ``esmvaltool-python``
package, so this package will always be installed as a dependency if you install
one or more of the packages for other languages.

There is also a lesson available in the
`ESMValTool tutorial <https://esmvalgroup.github.io/ESMValTool_Tutorial/>`_
that describes the installation of the ESMValTool in more detail. It can be found
`here <https://esmvalgroup.github.io/ESMValTool_Tutorial/02-installation/index.html>`_.

.. _install_on_macosx:

Installation on MacOSX
======================

The Python diagnostics of the ESMValTool are supported on MacOSX, but Julia, NCL,
and R are not. If any of these are needed, deployment through a :ref:`Docker<install_with_docker>`
container is advised.

The ``esmvaltool-python`` diagnostics can be installed as follows:

First, ensure mamba is pre-installed (see `Mamba installation`_ for more details).

Create a new environment with the ``esmvaltool-python`` package:

.. code-block:: bash

    mamba create --name esmvaltool esmvaltool-python python=3.9

Activate the new environment:

.. code-block:: bash

    conda activate esmvaltool

Confirm that the ESMValTool is working with:

.. code-block:: bash

    esmvaltool --help

Note that some recipes may depend on the OpenMP library, which does not
install via mamba on MacOSX. To install this library, run:

.. code-block:: bash

    brew install libomp

to install the library with Homebrew. In case you do not have Homebrew, follow
installation instructions `here <https://brew.sh/>`__.

.. _install_with_pip:

Pip installation
================

It is also possible to install ESMValTool from `PyPI <https://pypi.org/project/ESMValTool/>`_.
However, this requires first installing dependencies that are not available on PyPI in some other way.
By far the easiest way to install these dependencies is to use `mamba`.
For a minimal mamba installation (recommended) go to https://mamba.readthedocs.io/en/latest/installation.html.

After installing mamba, download
`the file with the list of dependencies <https://raw.githubusercontent.com/ESMValGroup/ESMValTool/main/environment.yml>`_:

.. code-block:: bash

    wget https://raw.githubusercontent.com/ESMValGroup/ESMValTool/main/environment.yml

and install these dependencies into a new conda environment with the command

.. code-block:: bash

    mamba env create --name esmvaltool --file environment.yml


Finally, activate the newly created environment

.. code-block:: bash

    conda activate esmvaltool

and install ESMValTool as well as any remaining Python dependencies with the command:

.. code-block:: bash

    pip install esmvaltool

If you would like to run Julia diagnostic scripts, you will also need to
install the Julia dependencies:

.. code-block:: bash

    esmvaltool install Julia

If you would like to run R diagnostic scripts, you will also need to install the R
dependencies:

.. code-block:: bash

    esmvaltool install R

.. _install_with_docker:

Docker installation
===================

ESMValTool is also provided through `DockerHub <https://hub.docker.com/u/esmvalgroup/>`_
in the form of docker containers.
See https://docs.docker.com for more information about docker containers and how to
run them.

You can get the latest release with

.. code-block:: bash

   docker pull esmvalgroup/esmvaltool:stable

If you want to use the current main branch, use

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

.. _install_with_singularity:

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

.. _install_from_source:

Install from source
===================

Installing the tool from source is recommended if you need the very latest
features or if you would like to contribute to its development.

Obtaining the source code
-------------------------

The ESMValTool source code is available on a public GitHub repository:
https://github.com/ESMValGroup/ESMValTool

The easiest way to obtain it is to clone the repository using git
(see https://git-scm.com/). To clone the public repository:

.. code-block:: bash

    git clone https://github.com/ESMValGroup/ESMValTool

or

.. code-block:: bash

    git clone git@github.com:ESMValGroup/ESMValTool

if you prefer to connect to the repository over SSH.

The command above will create a folder called ``ESMValTool``
containing the source code of the tool in the current working directory.

.. note::
    Using SSH is much more convenient if you push to the repository regularly
    (recommended to back up your work), because then you do not need to type
    your password over and over again.
    See
    `this guide <https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account>`__
    for information on how to set it up if you have not done so yet.
    If you are developing ESMValTool on a shared compute cluster, you can set up
    `SSH agent forwarding <https://docs.github.com/en/free-pro-team@latest/developers/overview/using-ssh-agent-forwarding>`__
    to use your local SSH keys also from the remote machine.

It is also possible to work in one of the ESMValTool private repositories, e.g.:

.. code-block:: bash

    git clone https://github.com/ESMValGroup/ESMValTool-private

GitHub also allows one to download the source code in as a ``tar.gz`` or ``zip``
file.
If you choose to use this option, download the compressed file and extract its
contents at the desired location.


Software installation
---------------------

Once all prerequisites are fulfilled, ESMValTool can be installed by running
the following commands in the directory containing the ESMValTool source code
(called ``ESMValTool`` if you did not choose a different name):

.. code-block:: bash

    pip install --editable '.[develop]'

Using the ``--editable`` flag will cause the installer to create a symbolic link
from the installation location to your source code, so any changes you make to
the source code will immediately be available in the installed version of the
tool.
This command will also install extra development dependencies needed for
building the documentation, running the unit tests, etc.

If you would like to run Julia diagnostic scripts, you will need to
install the ESMValTool Julia dependencies:

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

**MacOSX note:** some recipes may depend on the OpenMP library, which does not
install via mamba on MacOSX. Instead run

.. code-block:: bash

    brew install libomp

to install the library with Homebrew. In case you do not have Homebrew, follow
installation instructions `here <https://brew.sh/>`__.

For a more complete installation verification, run the automated tests and
confirm that no errors are reported:

.. code-block:: bash

    pytest -m "not installation"

or if you want to run the full test suite remove the ``-m "not installation"`` flag;
also if you want to run the tests on multiple threads, making the run faster, use
the `-n N` flag where N is the number of available threads e.g:

.. code-block:: bash

    pytest -n 4


.. _esmvalcore-development-installation:

Using the development version of the ESMValCore package
-------------------------------------------------------

If you need the latest developments of the ESMValCore package, you
can install it from source into the same conda environment.

.. attention::
    The recipes and diagnostics in the ESMValTool repository are compatible
    with the latest released version of the ESMValCore.
    Using the development version of the ESMValCore package is only recommended
    if you are planning to develop new features for the ESMValCore, e.g.
    you want to implement a new preprocessor function.

First follow all steps above.
Next, go to the place where you would like to keep the source code and clone the
ESMValCore github repository:

.. code-block:: bash

    git clone https://github.com/ESMValGroup/ESMValCore

or

.. code-block:: bash

    git clone git@github.com:ESMValGroup/ESMValCore

The command above will create a folder called ``ESMValCore``
containing the source code of the tool in the current working directory.

Go into the folder you just downloaded

.. code-block:: bash

    cd ESMValCore

and then install ESMValCore in development mode

.. code-block:: bash

    pip install --editable '.[develop]'

To check that the installation was successful, run

.. code-block:: bash

    python -c 'import esmvalcore; print(esmvalcore.__path__[0])'

this should show the directory of the source code that you just downloaded.

If the command above shows a directory inside your conda environment instead,
e.g. ``~/mamba/envs/esmvaltool/lib/python3.8/site-packages/esmvalcore``, you
may need to manually remove that directory and run
```pip install -e '.[develop]'``
again.

.. _install_on_hpc:

Pre-installed versions on HPC clusters
======================================

The ESMValTool is also available on the HPC clusters CEDA-JASMIN and DKRZ-MISTRAL and there will be no need
to install it yourself if you are just running diagnostics:

 - CEDA-JASMIN: `esmvaltool` is available on the scientific compute nodes (`sciX.jasmin.ac.uk` where
   `X = 1, 2, 3, 4, 5`) after login and module loading via `module load esmvaltool`; see the helper page at
   `CEDA <https://help.jasmin.ac.uk/article/4955-community-software-esmvaltool>`__ ;
 - DKRZ-Mistral: `esmvaltool` is available on login nodes (`mistral.dkrz.de`) and pre- and post-processing
   nodes (`mistralpp.dkrz.de`) after login and module loading via `module load esmvaltool`; the command
   `module help esmvaltool` provides some information about the module.

.. _common installation issues:

Common installation problems and their solutions
================================================

Mamba fails to solve the environment
------------------------------------
If you see the text ``Solving environment:`` with the characters ``-\|/`` rotating
behind it for more than 10 minutes, mamba may be having problems finding a
working combination of versions of the packages that the ESMValTool depends on.
Because the ESMValTool is a community tool, there is no strict selection of
which tools can be used and installing the ESMValTool requires installing almost
any package that is available for processing climate data.
To help mamba solve the environment, you can try the following.

Always use the latest version of mamba, as problems have been reported by people
using older versions, to update, run:

.. code-block:: bash

    mamba update --name base mamba

Usually mamba is much better at solving new environments than updating older
environments, so it is often a good idea to create a new environment if updating
does not work.

It can help mamba if you let it know what version of certain packages you want,
for example by running

.. code-block:: bash

    mamba create -n esmvaltool esmvaltool python=3.9

you ask for Python 3.9 specifically and that makes it much easier for mamba to
solve the environment, because now it can ignore any packages that were built
for other Python versions. Note that, since the esmvaltool package is built with Python>=3.7,
asking for an older Python version, e.g. `python=3.6`, in this way, it will result in
installation failure.

Problems with proxies
---------------------
If you are installing ESMValTool from source from behind a proxy that does not
trust the usual PyPI URLs you can declare them with the option
``--trusted-host``, e.g.

.. code-block:: bash

    pip install --trusted-host=pypi.python.org --trusted-host=pypi.org --trusted-host=files.pythonhosted.org -e .[develop]

If R packages fail to download, you might be able to solve this by
setting the environment variable ``http_proxy`` to the correct value, e.g.
in bash:

.. code-block:: bash

    export http_proxy=http://user:pass@proxy_server:port

the username and password can be omitted if they are not required. See e.g.
`here <https://support.rstudio.com/hc/en-us/articles/200488488-Configuring-R-to-Use-an-HTTP-or-HTTPS-Proxy>`__
for more information.

Anaconda servers connection issues
----------------------------------
HTTP connection errors (of e.g. type 404) to the Anaconda servers are rather common, and usually a retry
will solve the problem.

Installation of R packages fails
--------------------------------
Problems have been reported if the ``R`` interpreter was made available
through the ``module load`` command in addition to installation from mamba.
If your ESMValTool conda environment is called ``esmvaltool`` and you want to
use the R interpreter installed from mamba, the path to the R interpreter should
end with ``mamba/envs/esmvaltool/bin/R`` or ``conda/envs/esmvaltool/bin/R``.
When the conda environment for ESMValTool is activated, you can check which R
interpreter is used by running

.. code-block:: bash

    which R

The Modules package is often used by system administrators to make software
available to users of scientific compute clusters.
To list any currently loaded modules run ``module list``, run ``module help``
or ``man module`` for more information about the Modules package.

Problems when using ssh
-----------------------
If you log in to a cluster or other device via SSH and your origin
machine sends the ``locale`` environment via the SSH connection,
make sure the environment is set correctly, specifically ``LANG`` and
``LC_ALL`` are set correctly (for GB English UTF-8 encoding these
variables must be set to ``en_GB.UTF-8``; you can set them by adding
``export LANG=en_GB.UTF-8`` and ``export LC_ALL=en_GB.UTF-8``) in your
origin or login machines’ ``.profile``.

Problems when updating the conda environment
--------------------------------------------
Usually mamba is much better at solving new environments than updating older
environments, so it is often a good idea to create a new environment if updating
does not work. See also `Mamba fails to solve the environment`_.

Do not run ``mamba update --update-all`` in the ``esmvaltool``
environment since that will update some packages that are pinned to
specific versions for the correct functionality of the tool.


.. _move-to-mamba:

Move to Mamba
=============

Mamba is a much faster alternative to `conda`, and environment creation and updating
benefits from the use of a much faster (C++ backend) dependency solver; tests have been performed
to verify the integrity of the `esmvaltool` environment built with `mamba`, and we are
now confident that the change will not affect the way ESMValTool is installed and run, whether it be on a Linux or OSX platform.
From the user's perspective, it is a straightforward use change: the CLI (command line
interface) of `mamba` is identical to `conda`: any command that was run with `conda` before
will now be run with `mamba` instead, keeping all the other command line arguments and
flags as they were before. The only place where `conda` should not be replaced with `mamba`
at command line level is at the environment activation point: `conda activate` will still
have to be used.
