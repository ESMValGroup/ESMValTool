.. _install:

************
Installation
************

ESMValTool supports Python 3.12 and later and requires Linux or MacOS.
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
`conda <https://docs.conda.io/projects/conda/en/stable/>`__ because it is faster.

For a minimal mamba installation (recommended) go to
https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html.

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
languages: Python, NCL, and R.

There are three language specific packages available:

* ``esmvaltool-ncl``
* ``esmvaltool-python``
* ``esmvaltool-r``

The main ``esmvaltool`` package contains all three subpackages listed above.
If you only need to run a recipe with diagnostics in some of these languages, it
is possible to install only the dependencies needed to do just that.
The diagnostic script(s) used in each recipe, are documented in :ref:`recipes`.
The extension of the diagnostic script can be used to see in which language a
diagnostic script is written (``.py`` for Python, ``.ncl`` for NCL, ``.R`` for R diagnostics).

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

.. _install_on_macosx:

ESMValTool installation on MacOS
---------------------------------

The Python diagnostics of the ESMValTool are supported on MacOS, but NCL,
and R are not.
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

*Install ESMValTool*

It is recommended to use `pixi <https://pixi.prefix.dev>`__ to manage ESMValTool dependencies.
Follow the `pixi installation instructions <https://pixi.prefix.dev/latest/installation/>`__
to install it.

To create and activate a software environment to develop ESMValTool, go to the directory
containing the ESMValTool source code that you just downloaded. It is called
``ESMValTool`` if you did not choose a different name.

.. code-block:: bash

    cd ESMValTool

and run the command:

.. code-block:: bash

    pixi shell

This will install all of the required dependencies for running and developing
Python diagnostics.

.. tip::

    To exit the pixi environment, run ``exit`` or press ``Ctrl+D``.

If you are on Linux, it is also possible to install the dependencies for NCL
and R diagnostics by running

.. code-block:: bash

    pixi shell -e esmvaltool-dev

or just the dependencies for NCL diagnostics with ``pixi shell -e esmvaltool-ncl-dev``
and for R diagnostics with ``pixi shell -e esmvaltool-r-dev``. There are also
environments with the dependencies for running diagnostics in specific languages
without the development dependencies, i.e. ``pixi shell -e esmvaltool-python``
for running Python diagnostics, ``pixi shell -e esmvaltool-ncl`` for running
NCL diagnostics, ``pixi shell -e esmvaltool-r`` for running R diagnostics, and
``pixi shell -e esmvaltool`` for running any diagnostic. The environment names
correspond to the ESMValTool subpackages described in :ref:`conda subpackages`
and the ``-dev`` suffix indicates that the development dependencies are included.

.. tip::

    If you find that solving the environments (i.e. finding out which
    combination of package versions is compatible and can be installed) is
    slow, you can add the ``--frozen`` flag to the commands above to skip the
    solve step. Add ``export PIXI_FROZEN=true`` to your ``~/.bashrc`` file to
    make this the default behavior.

.. tip::

    If you find that solving the environments uses too much memory, you can
    `set the number of parallel solves to one <https://pixi.prefix.dev/latest/reference/pixi_configuration/#concurrency>`__
    by running ``pixi config set concurrency.solves 1``.

If you are planning to do any coding, install the :ref:`pre-commit`
hooks by running:

.. code-block:: bash

    pre-commit install

these will make sure that when you commit your changes, they will be formatted
correctly.

*Check your ESMValTool installation*

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
can install it from source into the same pixi environment.

.. attention::
    The recipes and diagnostics in the ESMValTool repository are compatible
    with the latest released version of the ESMValCore.
    Using the development version of the ESMValCore package is only recommended
    if you are planning to develop new features for the ESMValCore, e.g.
    you want to implement a new preprocessor function.

First follow the steps in the section above to
:ref:`install ESMValTool from source <install_from_source>`.

If you just want to use the latest features of ESMValCore, you can create a
new pixi environment with the command:

.. code-block:: bash

    pixi shell -e esmvalcore-dev

and verify that you now have a development version of ESMValCore installed by
running

.. code-block:: bash

    esmvaltool version

If you are planning to do development work on ESMValCore, you will need to
clone the ESMValCore repository and tell pixi to use it.

Go to the place where you would like to keep the ESMValCore source code and
clone the ESMValCore github repository:

.. code-block:: bash

    git clone https://github.com/ESMValGroup/ESMValCore

or

.. code-block:: bash

    git clone git@github.com:ESMValGroup/ESMValCore

The command above will create a folder called ``ESMValCore``
containing the source code of the tool in the current working directory.

Edit the ``pyproject.toml`` file in the ESMValTool source code and remove
this line:

.. literalinclude:: ../../../../pyproject.toml
   :start-at: ESMValCore = { git = "https://github.com/ESMValGroup/ESMValCore.git", branch = "main" }
   :end-at: ESMValCore = { git = "https://github.com/ESMValGroup/ESMValCore.git", branch = "main" }

and uncomment this line:

.. literalinclude:: ../../../../pyproject.toml
   :start-at: # ESMValCore = { path = "../ESMValCore", editable = true }
   :end-at: # ESMValCore = { path = "../ESMValCore", editable = true }

Make sure the path in the line you just uncommented points to the ESMValCore
repository you just cloned.

Finally, also install the ESMValCore :ref:`pre-commit` hooks by going into the
directory containing the ESMValCore source code and running:

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
:download:`pyproject.toml <../../../../pyproject.toml>`.

.. warning::

    It is recommended to use the installation with mamba instead, as it may not
    be easy to install the correct versions of all dependencies.

After installing the dependencies that are not available from PyPI_, install
ESMValTool and any remaining Python dependencies with the command:

.. code-block:: bash

    pip install esmvaltool

.. _common installation issues:

Common installation problems and their solutions
================================================

Problems with proxies
---------------------
If you are installing ESMValTool from source from behind a proxy that does not
trust the usual PyPI URLs you can declare them with the option
``--trusted-host``, e.g.

.. code-block:: bash

    pip install --trusted-host=pypi.python.org --trusted-host=pypi.org --trusted-host=files.pythonhosted.org --no-deps -e .[develop]

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
