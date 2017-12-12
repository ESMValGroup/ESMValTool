.. _annex_b:

Annex B - development environment
*********************************

.. _git_repository:

Git repository
==============

The development of the ESMValTool source code is managed by the version control system Git
(https://www.git-scm.com/). The ESMValTool development environment is hosted on GitHub. The following
description gives an overview of the typical workflow and usage for implementing new diagnostics or technical
changes into the ESMValTool. For general information on Git, see e.g. the online documentation at
https://www.git-scm.com/doc.

.. .. _fig_1:
.. .. figure::  /figures/Attention.png
..    :align:   center
.. note::
   .. figure::  /figures/Attention.png

   The OPEN ESMValTool repository is located at https://github.com/ESMValGroup/ESMValTool.
   The PRIVATE ESMValTool repository for the ESMValTool development team is located at https://github.com/ESMValGroup/ESMValTool-private

There are two ESMValTool GitHub repositories available:

#. **OPEN GitHub repository** that is open to the public (https://github.com/ESMValGroup/ESMValTool). The ESMValTool is released as open-source software under the Apache License 2.0, version 2.0. Use of the software constitutes acceptance of this license and terms.
#. **PRIVATE GitHub repository** that is restricted to the ESMValTool Development Team (https://github.com/ESMValGroup/ESMValTool-private). This repository is only accessible to ESMValTool developers that accepted the terms of use (see http://www.esmvaltool.org/license.html) for the ESMValTool development environment. The use of the ESMValTool software and access to the private ESMValTool GitHub repository constitutes acceptance of these terms. **When you fork or copy this repository, you must ensure that you do not copy the PRIVATE repository into the open domain!**

All developments can be made in either of the two repositories. The creation of *FEATURE BRANCHES* (see below),
however, is restricted to registered ESMValTool developers in both repositories. We encourage all developers to
join the ESMValTool development team. Please contact the ESMValTool Core Development Team (:numref:`core_team`) if you
want to join the ESMValTool development team.

The PRIVATE GitHub repository offers a central protected environment for ESMValTool developers who would like to
keep their contributions undisclosed (e.g., unpublished scientific work, work of PhD students in progress) while
at the same time benefiting from the possibilities of collaborating with other ESMValTool developers and having
a backup of their work. *FEATURE BRANCHES* created in the PRIVATE repository are only visible to the ESMValTool
development team but not to the entire public. The concept of a PRIVATE repository has proven to be very useful
to efficiently share code during the development across institutions and projects in a common repository without
having the contributions immediately accessible to the public.

The OPEN and the PRIVATE repository both contain the following kinds of branches (see also :numref:`fig_git`):

* *MASTER BRANCH* (official releases),
* *DEVELOPMENT BRANCH* (includes approved new contributions but version is not yet fully tested),
* *FEATURE BRANCH* (development branches for new features and diagnostics created by ESMValTool developers, the naming convention for *FEATURE BRANCHES* is <Project>_<myfeature>).

**Access rights**

* Write access to the *MASTER* and *DEVELOPMENT BRANCH* in both the OPEN and the PRIVATE GitHub repositories is restricted to the ESMValTool core development team. Both branches have access right READ ALL (noting again that access to the PRIVATE GitHub repository is restricted to the ESMValTool development team).
* *FEATURE BRANCHES* in both the PUBLIC and the PRIVATE repository can be created by all members of the ESMValTool development team (i.e. members team "ESMValTool-DevelopmentTeam" in the GitHub organization "ESMValGroup"). If needed, branches can be individually write-protected within each repository so that other developers cannot accidently push changes to these branches.

Workflow
========

The *MASTER BRANCH* of the PRIVATE repository will be synchronized with the *MASTER BRANCH* of the OPEN repository
in an automated manner. This ensures that they are identical at all times (see schematic in :numref:`fig_git`). The
recommended workflow for members of the ESMValTool development team is to create additional *FEATURE BRANCHES* in
either the PUBLIC or the PRIVATE repository, see further instructions below.

The following description gives an overview of the typical workflow and usage for implementing new diagnostics
or technical changes into the ESMValTool. The description assumes that your local development machine is running
a Unix-like operating system. For a general introduction to Git tutorials such as, for instance,
https://www.git-scm.com/docs/gittutorial are recommended.

.. _fig_git:
.. figure::  /figures/git_diagram.png
   :align:   center

   Schematic diagram of the ESMValTool GitHub repositories. *FEATURE BRANCHES* can only be created by members of the ESMValTool Development Team. The naming convention for *FEATURE BRANCHES* is <Project>_<myfeature>.

**Getting started**

First make sure that you have Git installed on your development machine. You can ask your administrator for
assistance. You can test this with the command:

.. code:: bash

   git --version

In order to properly identify your contributions to the ESMValTool you need to configure your local Git with
some personal data. This can be done with the following commands:

.. code:: bash

   git config --global user.name "YOUR NAME"
   git config --global user.email "YOUR EMAIL"

.. note:: For working on GitHub you need to create an account and login to https://github.com/.

**Option 1: working with the ESMValTool GitHub repository by creating a fork**

In order to start working with the ESMValTool source code you need to get a copy from the OPEN or PRIVATE
repository (e.g., https://github.com/ESMValGroup/ESMValTool). You can fork the ESMValTool repository to your
GitHub account. When you fork or copy the PRIVATE repository, you must ensure that you do not copy it into the
open domain!

* Login to GitHub.com
* On GitHub, go to the website of the ESMValTool repository (e.g., https://github.com/ESMValGroup/ESMValTool) and click on the button "fork"

.. figure::  /figures/git_fork.png

* Choose to create the fork of the ESMValTool repository under your account
* Select the "*DEVELOPMENT BRANCH*" and create a new *FEATURE BRANCH* for the diagnostic/feature you want to implement. Please follow the following naming convention for your new *FEATURE BRANCH*: <Project>_<myfeature>.

.. figure::  /figures/git_branch.png

* On this fork click the button "Clone or Download" and copy the URL shown there
* Open a terminal window and go to the folder where you would like to store your local copy of the ESMValTool source code
* Run git clone with the URL copied:

.. code:: bash

   git clone <URL_OF_YOUR_FORK>

This will clone your fork of the ESMValTool repository at GitHub to a local folder. You can now query the status of your local working copy with:

.. code:: bash

   git status

You will see that you are on a branch called master and your local working copy is up to date with the remote
repository (your fork). With

.. code:: bash

   git branch --all

you can list all available remote and local branches; now switch to your feature branch by:

.. code:: bash

   git checkout <NAME_OF_FEATURE_BRANCH>

You can now start coding. To check your current developments you can use the command

.. code:: bash

   git status

You can add new files and folders that you want to have tracked by Git using:

.. code:: bash

   git add <NEW_FILE|FOLDER>

To simply add all new files use:

.. code:: bash

   git add .

It is recommended to commit your changes to your local working copy often via:

.. code:: bash

   git commit "YOUR COMMIT MESSAGE"

Alternatively, type:

.. code:: bash

   git commit -a

Then an editor window will open, and you can type a long commit message. In order to inspect your changes you
can use the gitk viewer (use man gitk for all options):

.. code:: bash

   gitk

Or if you are in textmode only you can inspect your changes with (use man git-log for all options):

.. code:: bash

   git log

To share your work and to have an online backup, push your local development to your fork at GitHub. **We strongly
recommend doing this on a regular basis:**

.. code:: bash

   git push origin

Once your development is finished, go to the GitHub website of your fork and initiate a pull request to the
ESMValTool Core Development Team by clicking on the button "Pull request". Your changes will then be tested,
discussed and then implemented into the *DEVELPOMENT BRANCH*.

**Option 2: working with the ESMValTool GitHub Repositories without creating a fork**

As a member of the ESMValTool development team you can create *FEATURE BRANCHES* in the OPEN as well as in the
PRIVATE repository. We encourage all ESMValTool developers to use the following workflow for long-lived
developments (>2 weeks).

* Login to GitHub.com
* On GitHub, go to the website of the ESMValTool repository (https://github.com/ESMValGroup/ESMValTool-private or https://github.com/ESMValGroup/ESMValTool)
* Click on the button create *FEATURE BRANCH*
* Select the *"DEVELOPMENT" BRANCH* and create a new feature branch for the diagnostic/feature you want to implement. Please follow the following naming convention for your new *FEATURE BRANCH*: <Project>_<myfeature>.

.. figure::  ./figures/git_branch_2.png

* Click the button “Clone or Download” and copy the URL shown there
* Open a terminal window and go to the folder where you would like to store your local copy of the ESMValTool source
* Type git clone, and paste the URL:

.. code:: bash

   git clone <URL_FROM_CLIPBOARD>

This will clone the ESMValTool repository at GitHub to a local folder.
You can now query the status of your local working copy with:

.. code:: bash

   git status

You will see that you are on a branch called master and your local working copy is up to date with the remote
repository. With

.. code:: bash

   git branch --all

you can list all available remote and local branches; now switch to your feature branch by:

.. code:: bash

   git checkout <NAME_OF_YOUR_FEATURE_BRANCH>

You can now start coding. To check your current developments you can use the command

.. code:: bash

   git status

You can add new files and folders that you want to have tracked by Git using:

.. code:: bash

   git add <NEW_FILE|FOLDER>

To simply add all new files use:

.. code:: bash

   git add .

It is recommended to commit your changes to your local working copy often via:

.. code:: bash

   git commit –am "YOUR COMMIT MESSAGE"

Alternatively, type:

.. code:: bash

   git commit -a

Then an editor window will open, and you can type a long commit message. In order to inspect your changes you
can use the gitk viewer (use man gitk for all options):

.. code:: bash

   gitk

Or if you are in textmode only you can inspect your changes with (use man git-log for all options):

.. code:: bash

   git log

To share your work and to have an online backup, push your local development to your FEATURE BRANCH at GitHub.
**We strongly recommend doing this on a regular basis**:

.. code:: bash

   git push origin <YOUR_FEATURE_BRANCH>

Once your development is finished, go to the GitHub website of the ESMValTool repository and switch to your
*FEATURE BRANCH*. You can then initiate a pull request for the *DEVELPOMENT BRANCH* to the ESMValTool Core
Development Team by clicking on the button "Pull request". Your changes will then be tested, discussed and then
implemented into the *DEVELPOMENT BRANCH*.

General do-s and don't-s
========================

**Do-s**

* Create a *FEATURE BRANCH* (see :numref:`git_repository` for details) for developing the ESMValTool. The naming convention for *FEATURE BRANCHES* is <Project>_<myfeature>.
* Try using self-explanatory names for new branches (avoid things like: "my_branch" or "my_development")
* Comment your code as much as possible.
* Use short but self-explanatory variable names (e.g., model_input and reference_input instead of xm and xr).
* Consider a modular/functional programming style. This often makes code easier to read and deletes intermediate variables from memory immediately. If possible, separate diagnostic calculations from plotting routines.
* Consider reusing or extending existing code (see plotting functions, general calculations). General-purpose code can be found in diag_scripts/lib/ and in plot_scripts/.
* Comment all switches and parameters including a list of all possible settings/options in the header section of your code.
* Use templates for namelists and diagnostics to ensure proper documentation (see :numref:`std_nml`).
* Keep your development branch updated regularly with the master/development branch.

**Don't-s**

* Do not use other programming languages than the ones currently supported (NCL, Python, R). If you want to use a programming language not yet used, please contact the ESMValTool core development team.
* Avoid large (memory, disk space) intermediate results. Delete intermediate files/variables or see modular/functional programming style.
* Do not use hard-coded pathnames or filenames.

.. _wiki:

ESMValTool development team wiki
================================

The latest information on the ESMValTool and diagnostics under development can be found on the wiki of the OPEN
and PRIVATE GitHub repository:

* OPEN GitHub repository: https://github.com/ESMValGroup/ESMValTool/wiki
* PRIVATE GitHub repository: https://github.com/ESMValGroup/ESMValTool-private/wiki

All users and developers are strongly encouraged to frequently check the ESMValTool wiki for new information,
contact data or observational data. Please contact the ESMValTool Core Development Team for access to the wiki
(see :numref:`core_dev_team`).

