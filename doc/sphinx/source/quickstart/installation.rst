.. _install:

************
Installation
************

.. note::
   ESMValTool now uses `mamba` instead of `conda` for the recommended installation.
   For more information about the change, have a look at :ref:`Move to Mamba<move-to-mamba>`.

ESMValTool supports Python 3.11 and later and requires Linux or MacOS.
Successful usage on Windows has been reported by following the Linux
installation instructions with
`WSL <https://learn.microsoft.com/en-us/windows/wsl/install>`__.

ESMValTool can be installed in multiple ways.

Recommended installation method:

Install the :ref:`mamba package manager <install_with_mamba>` and then follow
the instructions for

* :ref:`ESMValTool installation on Linux<install_on_linux>`
* :ref:`ESMValTool installation on MacOS<install_on_macosx>`.

Further options for installation are:


* :ref:`From the source code<install_from_source>` available at https://github.com/ESMValGroup/ESMValTool;
* :ref:`From pre-installed versions on HPC clusters<install_on_hpc>`;
* :ref:`Deployment through a Docker container<install_with_docker>` (see https://www.docker.com);
* :ref:`Deployment through a Singularity container<install_with_singularity>` (see https://sylabs.io/guides/latest/user-guide/);
* :ref:`Installation with pip <install_with_pip>` (see https://pypi.org);
* :ref:`installation_from_the_conda_lock_file`.

The next sections will detail the procedure to install ESMValTool through each
of these methods.

There is also a lesson available in the
`ESMValTool tutorial <https://tutorial.esmvaltool.org/>`_
that describes the installation of the ESMValTool in more detail.
It can be found
`here <https://tutorial.esmvaltool.org/02-installation/index.html>`_.

See `common installation issues`_ if you run into trouble.

.. _install_with_mamba:

Mamba/Conda installation
========================

In order to install ESMValTool and its dependencies from
`conda-forge <https://conda-forge.org/>`__, you will first need to install the
`mamba package manager <https://mamba.readthedocs.io>`__.
We recommend using `mamba <https://mamba.readthedocs.io>`__ as a package manager
for your conda environments instead of
`conda <https://docs.conda.io/projects/conda/en/stable/>`__ because it is
much faster, see `move-to-mamba`_ for more information.

For a minimal mamba installation (recommended) go to
https://mamba.readthedocs.io/en/latest/installation.html.

.. note::
    It is recommended that you always use the latest version of mamba, as
    problems have been reported when trying to use older versions.

.. note::
    Some systems provide a pre-installed version of conda or mamba (e.g. via
    the module environment).
    However, several users reported problems when installing with such versions.
    It is therefore preferable to use a local, fully user-controlled mamba
    installation.

First download the installation file for
`Linux <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh>`_
or
`MacOSX <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh>`_.
After downloading the installation file from one of the links above, execute it
by running (Linux example):

.. code-block:: bash

    bash Miniforge3-Linux-x86_64.sh

and follow the instructions on your screen.

.. note::
    Make sure to choose an installation location where you have at least 10 GB
    of disk space available.

During installation, mamba will ask you if you want ``mamba`` to be
automatically loaded from your ``.bashrc`` or ``.bash-profile`` files.
It is recommended that you answer yes.
If you answered no, you can load the correct paths and environment variables
later by running:

.. code-block:: bash

    source <prefix>/etc/profile.d/conda.sh

where ``<prefix>`` is the installation location of mamba (e.g.
``/home/$USER/miniforge3`` if you chose the default installation path).

If you use another shell than Bash, have a look at the available configurations
in the ``<prefix>/etc/profile.d`` directory.

You can check that mamba installed correctly by running

.. code-block:: bash

    which mamba

this should show the path to your mamba executable, e.g.
``~/miniforge3/bin/mamba``.

It is recommended to update both mamba and conda after installing:

.. code-block:: bash

    mamba update --name base mamba conda

.. _install_on_linux:

ESMValTool installation on Linux
--------------------------------

Once you have installed the mamba package manager, you can install
the entire ESMValTool package by running:

.. code-block:: bash

    mamba create --name esmvaltool esmvaltool

It is also possible to install just a subset of the ESMValTool dependencies
by installing one or more of the :ref:`subpackages <conda subpackages>`
described in the next section.

The command above will create a new
`conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments>`_
called ``esmvaltool``, and install ESMValTool in it.
Of course it is also possible to choose a different name than ``esmvaltool``
for the environment.

.. note::

    Creating a new conda environment is often much faster and more reliable than
    trying to update an existing conda environment.
    Therefore it is recommended that you create a new environment when you
    want to upgrade to the latest version.

The next step is to check that the installation works properly.

First activate the environment with the command:

.. code-block:: bash

    conda activate esmvaltool

and then run the tool with the command:

.. code-block:: bash

    esmvaltool --help

If everything was installed properly, ESMValTool should have printed a help
message to the console.


.. _conda subpackages:

Installation of subpackages
---------------------------

The diagnostics bundled in ESMValTool are scripts in four different programming
languages: Python, NCL, R, and Julia.

There are three language specific packages available:

* ``esmvaltool-ncl``
* ``esmvaltool-python``
* ``esmvaltool-r``

The main ``esmvaltool`` package contains all three subpackages listed above.
For the Julia dependencies, there is no subpackage yet, but there are special
:ref:`installation instructions <install_julia_dependencies>`.
If you only need to run a recipe with diagnostics in some of these languages, it
is possible to install only the dependencies needed to do just that.
The diagnostic script(s) used in each recipe, are documented in :ref:`recipes`.
The extension of the diagnostic script can be used to see in which language a
diagnostic script is written (``.py`` for Python, ``.ncl`` for NCL, ``.R`` for R, and ``.jl`` for Julia diagnostics).

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

.. _install_julia_dependencies:

Installation of Julia dependencies
----------------------------------

If you want to use the ESMValTool Julia functionality, you will also need to
install Julia. If you are just getting started, we suggest that you
come back to this step later when, and if you need it.
To perform the Julia installation, make sure that your conda
environment is activated and then execute

.. code-block:: bash

    curl -fsSL https://install.julialang.org | sh -s -- --yes
    esmvaltool install Julia
.. _install_on_macosx:

ESMValTool installation on MacOS
---------------------------------

The Python diagnostics of the ESMValTool are supported on MacOS, but Julia,
NCL, and R are not.
If any of these are needed, deployment through a
:ref:`Docker<install_with_docker>`
container is advised.

The ``esmvaltool-python`` diagnostics can be installed as follows:

First, ensure mamba is installed (see install_with_mamba_ for more details).

Create a new environment with the ``esmvaltool-python`` package:

.. code-block:: bash

    mamba create --name esmvaltool esmvaltool-python

Activate the new environment:

.. code-block:: bash

    conda activate esmvaltool

Confirm that the ESMValTool is working with:

.. code-block:: bash

    esmvaltool --help

Note that some recipes may depend on the OpenMP library, which does not
install via mamba on MacOS. To install this library, run:

.. code-block:: bash

    brew install libomp

to install the library with Homebrew. In case you do not have Homebrew, follow
installation instructions `here <https://brew.sh/>`__.

.. _install_from_source:

Install from source
===================

Installing the tool from source is recommended if you need the very latest
features or if you would like to contribute to its development.

*Obtaining the source code*

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

*Install dependencies*

It is recommended to use mamba to manage ESMValTool dependencies.
See the :ref:`mamba installation instructions <install_with_mamba>` at the top
of this page for instructions on installing mamba.
To simplify the installation process, an environment definition file is provided
in the repository (``environment.yml`` in the root folder).

The ESMValTool conda environment file can also be used as a requirements list
for those cases in which a mamba installation is not possible or advisable.
From now on, we will assume that the installation is going to be done through
mamba.

Ideally, you should create a separate conda environment for ESMValTool, so it is
independent from any other Python tools present in the system.

To create an environment, go to the directory containing the ESMValTool source
code that you just downloaded. It is called ``ESMValTool`` if you did not
choose a different name.

.. code-block:: bash

    cd ESMValTool

and create a new environment called ``esmvaltool`` with the command (when on
Linux):

.. code-block:: bash

    mamba env create --name esmvaltool --file environment.yml

or (when on MacOS)

.. code-block:: bash

    mamba env create --name esmvaltool --file environment_osx.yml

This will install all of the required development dependencies.
Note that the MacOS environment file contains only Python dependencies,
so you will not be able to run NCL, R, or Julia diagnostics with it.

.. note::
    The environment is called ``esmvaltool`` in the example above, but it is
    possible to use the option ``--name some_environment_name`` to define a
    different name.
    This can be useful when you have an older ESMValTool installation that you
    would like to keep.
    It is recommended that you create a new environment when updating ESMValTool.

.. note::
    There is also a pure-Python environment file ``esmvaltool_python.yml``
    which is a softlink of the ``environment_osx.yml`` file; this one is used
    by any build that needs only Python packages (i.e. no NCL and R), currently
    this is used by our documentation builds, but it could be used by anyone
    needing just the Python dependencies.

Next, activate the environment by using the command:

.. code-block:: bash

    conda activate esmvaltool

.. attention::
    From now on, we assume that the conda environment containing the
    development dependencies for ESMValTool is activated.

*Install ESMValTool*

Once all dependencies have been installed, ESMValTool itself can be installed by
running the following command in the directory containing the ESMValTool source
code (called ``ESMValTool`` if you did not choose a different name):

.. code-block:: bash

    pip install --editable '.[develop]'

Using the ``--editable`` flag will cause the installer to create a symbolic link
from the installation location to your source code, so any changes you make to
the source code will immediately be available in the installed version of the
tool.

If you would like to run Julia diagnostic scripts, you will need to
install Julia and the ESMValTool Julia dependencies:

.. code-block:: bash

    curl -fsSL https://install.julialang.org | sh -s -- --yes
    esmvaltool install Julia

If you are planning to do any coding, install the :ref:`esmvaltool:pre-commit`
hooks by running:

.. code-block:: bash

    pre-commit install

these will make sure that when you commit your changes, they will be formatted
correctly.

The next step is to check that the installation works properly.
To do this, run the tool with:

.. code-block:: bash

    esmvaltool --help

If everything was installed properly, ESMValTool should have printed a
help message to the console.

.. note::
    **MacOS users:** some recipes may depend on the OpenMP library, which does not
    install via mamba on MacOS. Instead run

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

This concludes the installation from source guide. However, if you would like
to do development work on ESMValCore, please read on.

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

First follow the steps in the section above to
:ref:`install ESMValTool from source <install_from_source>`.
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
e.g. ``~/miniforge3/envs/esmvaltool/lib/python3.11/site-packages/esmvalcore``,
you may need to manually remove that directory and run
``pip install --editable '.[develop]'`` again.

Finally, also install the :ref:`esmvaltool:pre-commit` hooks by running:

.. code-block:: bash

    pre-commit install

these will make sure that when you commit your changes, they will be formatted
correctly.

.. _install_on_hpc:

Pre-installed versions on HPC clusters / other servers
======================================================

ESMValTool is available on the HPC clusters CEDA-JASMIN and DKRZ-Levante, and on the Met Office Linux
estate, so there is no need to install ESMValTool if you are just running recipes:

 - CEDA-JASMIN: `esmvaltool` is available on the scientific compute nodes (`sciX.jasmin.ac.uk` where
   `X = 1, 2, 3, 4, 5`) after login and module loading via `module load esmvaltool`; see the helper page at
   `CEDA <https://help.jasmin.ac.uk/article/4955-community-software-esmvaltool>`__ .
 - DKRZ-Levante: `esmvaltool` is available on login nodes (`levante.dkrz.de`) after login and module loading
   via `module load esmvaltool`; the command `module help esmvaltool` provides some information about the module.
   A Jupyter kernel based on the latest module is available from `DKRZ-JupyterHub <https://jupyterhub.dkrz.de/hub/home>`__.
 - Met Office: `esmvaltool` is available on the Linux estate after login and module loading via `module load`;
   see the ESMValTool Community of Practice SharePoint site for more details.
 - NSC-Tetralith and Freja: `esmvaltool` is available after login and module loading via `module load esmvaltool`.

The ESMValTool Tutorial provides a `quickstart guide <https://tutorial.esmvaltool.org/01-quickstart/index.html>`__
that is particularly suited for new users that have an access to pre-installed version of ESMValTool.

Information on how to request an account at CEDA-JASMIN and DKRZ-Levante and to get started with these HPC clusters
can be found on the setup page of the tutorial `here <https://tutorial.esmvaltool.org/setup.html>`__.

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

.. _install_with_pip:

Pip installation
================

It is also possible to install ESMValTool from
`PyPI <https://pypi.org/project/ESMValTool/>`_.
However, this requires first installing dependencies that are not available
on PyPI in some other way.
The list of required dependencies can be found in
:download:`environment.yml <../../../../environment.yml>`.

.. warning::

    It is recommended to use the installation with mamba instead, as it may not
    be easy to install the correct versions of all dependencies.

After installing the dependencies that are not available from PyPI_, install
ESMValTool and any remaining Python dependencies with the command:

.. code-block:: bash

    pip install esmvaltool

If you would like to run Julia diagnostic scripts, you will also need to
install Julia and the ESMValTool Julia dependencies:

.. code-block:: bash

    curl -fsSL https://install.julialang.org | sh -s -- --yes
    esmvaltool install Julia

.. _installation_from_the_conda_lock_file:

Installation from the conda lock file
=====================================

The conda lock file is an alternative to the ``environment.yml`` file used in
the :ref:`installation from source instructions <install_from_source>`.
All other steps in those installation instructions are the same.

The conda lock file can be used to install the dependencies of ESMValTool
whenever the conda environment defined by ``environment.yml`` can not be solved
for some reason.
A conda lock file is a reproducible environment file that contains links to
dependency packages as they are hosted on the Anaconda cloud;
these have frozen version numbers, build hashes, and channel names.
These parameters are established at the time of the conda lock file creation, so
may be outdated after a while.
Therefore, we regenerate these lock files every 10 days through automatic
Pull Requests (or more frequently, since the automatic generator runs on merges
on the ``main`` branch too), to minimize the risk of dependencies becoming
outdated.

Conda environment creation from a lock file is done with the following command:

.. code-block:: bash

   conda create --name esmvaltool --file conda-linux-64.lock

The latest, most up-to-date file can always be downloaded directly from the source code
repository, a direct download link can be found `here <https://raw.githubusercontent.com/ESMValGroup/ESMValTool/main/conda-linux-64.lock>`__.

.. note::
   For instructions on how to manually create the lock file, see
   :ref:`these instructions <esmvalcore:condalock-installation-creation>`.

.. _common installation issues:

Common installation problems and their solutions
================================================

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
does not work.

Do not run ``mamba update --update-all`` in the ``esmvaltool``
environment since that will update some packages that are pinned to
specific versions for the correct functionality of the tool.


.. _move-to-mamba:

Move to Mamba
=============

Mamba is a much faster alternative to `conda`, and environment creation and updating
benefits from the use of a much faster (C++ backend) dependency solver; tests have been performed
to verify the integrity of the `esmvaltool` environment built with `mamba`, and we are
now confident that the change will not affect the way ESMValTool is installed and run, whether it be on a Linux or OS platform.
From the user's perspective, it is a straightforward use change: the CLI (command line
interface) of `mamba` is identical to `conda`: any command that was run with `conda` before
will now be run with `mamba` instead, keeping all the other command line arguments and
flags as they were before. The only place where `conda` should not be replaced with `mamba`
at command line level is at the environment activation point: `conda activate` will still
have to be used.
