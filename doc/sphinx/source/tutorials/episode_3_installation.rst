.. _episode_3_installation:

Installation
============

.. admonition:: Overview
   :class: note

   .. grid:: 2
      :gutter: 1
      :margin: 3 3 0 5

      .. grid-item-card:: Timings
         :columns: 12

         * Teaching: 10 min
         * Exercises: 10 min

      .. grid-item-card:: Questions

         * What are the prerequisites for installing ESMValTool?
         * How do I confirm that the installation was successful?

      .. grid-item-card:: Objectives

         * Install ESMValTool.
         * Demonstrate that the installation was successful.

      .. grid-item-card:: Compatibility
         :columns: 12

         ESMValTool v2.16.0


Overview
--------

The instructions help with the installation of ESMValTool on operating systems
like Linux/MacOSX/Windows. We use the `Mamba <https://mamba.readthedocs.io/en/latest/index.html>`_
package manager to install ESMValTool. Other installation methods are also
available; they can be found in the
`ESMValTool documentation <https://docs.esmvaltool.org/en/latest/quickstart/installation.html>`_.
We will first install Mamba, and then ESMValTool. We end this chapter by
testing that the installation was successful.

Before we begin, here are all the possible ways in which you can use ESMValTool
depending on your level of expertise or involvement with ESMValTool and
associated software such as GitHub and Mamba.

1. If you have access to a server where ESMValTool is already installed as a
   module, for example the
   `CEDA JASMIN <https://help.jasmin.ac.uk/article/4955-community-software-esmvaltool>`_
   server, you can simply load the module with the following command:

   .. code-block:: bash

      module load esmvaltool

   After loading ``esmvaltool``, we can start using ESMValTool right away.
   Please see the :doc:`next lesson <episode_4_configuration>`.

2. If you would like to install ESMValTool as a mamba package, then this
   lesson will tell you how!

3. If you would like to start experimenting with existing diagnostics or
   contributing to ESMValTool, please see the instructions for source
   installation in the lesson `Development and contribution <https://tutorial.esmvaltool.org/07-development-setup/index.html>`_
   and in the
   `documentation <https://docs.esmvaltool.org/en/latest/quickstart/installation.html>`_.

.. admonition:: Install ESMValTool on Windows
   :class: note

   ESMValTool does not directly support Windows, but successful usage has been
   reported through the
   `Windows Subsystem for Linux (WSL) <https://docs.microsoft.com/en-us/windows/wsl/>`_,
   available in Windows 10. To install the WSL please follow the instructions
   on the `Windows Documentation page <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_.
   After installing the WSL, installation can be done using the same
   instructions for Linux/MacOSX.

Install ESMValTool on Linux/MacOSX
----------------------------------

Install Mamba
^^^^^^^^^^^^^

ESMValTool is distributed using
`Mamba <https://mamba.readthedocs.io/en/latest/index.html>`_. To install
mamba on ``Linux`` or ``MacOSX``, follow the instructions below:

1. Please download the installation file for the latest Mamba version for
   `Linux <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh>`_
   or for
   `MacOSX <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh>`_.

2. Next, run the installer from the place where you downloaded it.

   On ``Linux`` run:

   .. code-block:: bash

      bash Miniforge3-Linux-x86_64.sh

   On ``MacOSX`` run:

   .. code-block:: bash

      bash Miniforge3-MacOSX-x86_64.sh

3. Follow the instructions in the installer. The defaults should normally
   suffice.

4. You will need to restart your terminal for the changes to have effect.

5. We recommend updating mamba and conda before the esmvaltool installation.
   To do so, run:

   .. code-block:: bash

      mamba update --name base mamba conda

6. Verify you have a working mamba installation by running:

   .. code-block:: bash

      which mamba

   This should show the path to your mamba executable, e.g.
   ``~/miniforge3/bin/mamba``.

For more information about installing mamba, see the
`mamba installation documentation <https://docs.esmvaltool.org/en/latest/quickstart/installation.html#mamba-conda-installation>`_.

Install the ESMValTool package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ESMValTool package contains diagnostic scripts in four languages: R,
Python and NCL. This introduces a lot of dependencies, and therefore the
installation can take quite long. It is, however, possible to install
"subpackages" for each of the languages. The following (sub)packages are
available:

- ``esmvaltool-python``
- ``esmvaltool-ncl``
- ``esmvaltool-r``
- ``esmvaltool`` -> the complete package, i.e. the combination of the above.

For the tutorial, we will install the complete package. Thus, to install the
ESMValTool package, run

.. code-block:: bash

   mamba create --name esmvaltool esmvaltool

On MacOSX ESMValTool functionalities in NCL and R are not supported. To install
a Mamba environment on MacOSX, please refer to specific
`information <https://docs.esmvaltool.org/en/latest/quickstart/installation.html#esmvaltool-installation-on-macos>`_.

This will create a new
`Mamba environment <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_
called ``esmvaltool``, with the ESMValTool package and all of its dependencies
installed in it.

.. admonition:: Common issues
   :class: warning

   You can find a list of common installation problems and their solutions
   `here <https://docs.esmvaltool.org/en/latest/quickstart/installation.html#common-installation-problems-and-their-solutions>`_.

Test that the installation was successful
-----------------------------------------

To test that the installation was successful, run

.. code-block:: bash

   mamba activate esmvaltool

This should activate the environment called ``esmvaltool``. In the shell prompt
the active conda environment should have been changed from ``(base)`` to
``(esmvaltool)``.

Next, run

.. code-block:: bash

   esmvaltool --help

to display the command line help.

.. admonition:: Version of ESMValTool
   :class: admonition-todo

   Can you figure out which version of ESMValTool has been installed?

   .. dropdown:: Solution
      :color: secondary
      :icon: eye

      The ``esmvaltool --help`` command lists ``version`` as a command to get
      the version.

      When you run

      .. code-block:: bash

         esmvaltool version

      The version of ESMValTool installed should be displayed on the screen as:

      .. code-block:: text

         ESMValCore: 2.15.0
         ESMValTool: 2.15.0

      Note that on HPC servers such as JASMIN, sometimes a more recent
      development version may be displayed for ESMValTool, for example
      ``ESMValTool: 2.16.0.dev4+g74d7837f0``.

.. admonition:: Key points
   :class: important

   * All the required packages can be installed using mamba.
   * You can find more information about installation in the
     `documentation <https://docs.esmvaltool.org/en/latest/quickstart/installation.html>`_.
