.. _preparation-new-release:

Release schedule and procedure for ESMValCore and ESMValTool
============================================================

This document describes the process for the release of ESMValCore
and ESMValTool.
By following a defined process, we streamline the work, reduce
uncertainty about required actions, and clarify the state of the code for the
user.

ESMValTool follows a strategy of timed releases.
That means that we do releases with a regular frequency and all features
that are implemented up to a certain cut-off-point can go
into the upcoming release; those that are not are deferred to the next
release.
This means that generally no release will be delayed due to a pending feature.
Instead, the regular nature of the release guarantees that every feature can be
released in a timely manner even if a specific target release is missed.

Because of limited resources, only the latest released versions of ESMValTool and ESMValCore is maintained.
If your project requires longer maintenance or you have other concerns about
the release strategy, please contact the ESMValTool core development team, see
:ref:`Support-and-Contact`.


Overall Procedure
-----------------

Timeline
~~~~~~~~~

.. figure::  /figures/release-timeline.png
   :align:   center

   Example of a Release Timeline (in this case for 2.1.0)

1. Contributors assign issues (and pull requests) that they intend to finish before the due date, there is a separate milestone for ESMValCore and ESMValTool
2. The ESMValCore feature freeze takes place on the ESMValCore due date
3. Some additional testing of ESMValCore takes place
4. ESMValCore release
5. The ESMValTool feature freeze takes place
6. Some additional testing of ESMValTool takes place
7. ESMValTool release
8. Soon after the release, the core development team meets to coordinate the content of the milestone for the next release

.. _release_schedule:

Release schedule
~~~~~~~~~~~~~~~~

With the following release schedule, we strive to have three releases per year and to avoid releases too close to holidays, as well as avoiding weekends.

Upcoming releases
^^^^^^^^^^^^^^^^^

- 2.13.0 (Release Manager: `Julien Lenhardt`_)

+------------+------------+----------------------------------------+-------------------------------------+
|  Planned   |    Done    |            Event                       |             Changelog               |
+============+============+========================================+=====================================+
| 2025-08-18 |            | ESMValCore `Feature Freeze`_           |                                     |
+------------+------------+----------------------------------------+-------------------------------------+
| 2025-08-29 |            | ESMValCore Release 2.13.0              |                                     |
+------------+------------+----------------------------------------+-------------------------------------+
| 2025-09-01 |            | ESMValTool `Feature Freeze`_           |                                     |
+------------+------------+----------------------------------------+-------------------------------------+
| 2025-09-12 |            | ESMValTool Release 2.13.0              |                                     |
+------------+------------+----------------------------------------+-------------------------------------+

Past releases
^^^^^^^^^^^^^

- 2.12.0 (Release Manager: `Saskia Loosveldt Tomas`_)

+------------+------------+----------------------------------------+-------------------------------------+
|  Planned   |    Done    |            Event                       |             Changelog               |
+============+============+========================================+=====================================+
| 2025-01-13 |            | ESMValCore `Feature Freeze`_           |                                     |
+------------+------------+----------------------------------------+-------------------------------------+
| 2025-01-20 | 2025-02-27 | :esmvalcore-release:`v2.12.0` released | :ref:`esmvalcore:changelog-v2-12-0` |
+------------+------------+----------------------------------------+-------------------------------------+
| 2025-01-27 |            | ESMValTool `Feature Freeze`_           |                                     |
+------------+------------+----------------------------------------+-------------------------------------+
| 2025-02-03 | 2025-03-05 | :release:`v2.12.0` released            | :ref:`changelog-v2-12-0`            |
+------------+------------+----------------------------------------+-------------------------------------+

- 2.11.0 (Release Manager: Met Office: `Emma Hogan`_, `Chris Billows`_, `Ed Gillett`_)

+------------+------------+----------------------------------------+-------------------------------------+
|  Planned   |    Done    |            Event                       |             Changelog               |
+============+============+========================================+=====================================+
| 2024-04-22 |            | ESMValCore `Feature Freeze`_           |                                     |
+------------+------------+----------------------------------------+-------------------------------------+
| 2023-05-03 | 2024-07-03 | :esmvalcore-release:`v2.11.0` released | :ref:`esmvalcore:changelog-v2-11-0` |
+------------+------------+----------------------------------------+-------------------------------------+
| 2023-05-06 |            | ESMValTool `Feature Freeze`_           |                                     |
+------------+------------+----------------------------------------+-------------------------------------+
| 2023-05-17 | 2024-07-04 | :release:`v2.11.0` released            | :ref:`changelog-v2-11-0`            |
+------------+------------+----------------------------------------+-------------------------------------+

- 2.10.0 (Release Manager: `Klaus Zimmermann`_)

+------------+------------+----------------------------------------+-------------------------------------+
|  Planned   |    Done    |            Event                       |             Changelog               |
+============+============+========================================+=====================================+
| 2023-10-02 |            | ESMValCore `Feature Freeze`_           |                                     |
+------------+------------+----------------------------------------+-------------------------------------+
| 2023-10-09 | 2023-12-19 | :esmvalcore-release:`v2.10.0` released | :ref:`esmvalcore:changelog-v2-10-0` |
+------------+------------+----------------------------------------+-------------------------------------+
| 2023-10-16 |            | ESMValTool `Feature Freeze`_           |                                     |
+------------+------------+----------------------------------------+-------------------------------------+
| 2023-10-16 | 2023-12-20 | :release:`v2.10.0` released            | :ref:`changelog-v2-10-0`            |
+------------+------------+----------------------------------------+-------------------------------------+

- 2.9.0 (Release Manager: `Bouwe Andela`_)

+------------+------------+---------------------------------------+-------------------------------------+
|  Planned   |    Done    |            Event                      |             Changelog               |
+============+============+=======================================+=====================================+
| 2023-06-05 |            | ESMValCore `Feature Freeze`_          |                                     |
+------------+------------+---------------------------------------+-------------------------------------+
| 2023-06-12 | 2023-07-04 | :esmvalcore-release:`v2.9.0` released | :ref:`esmvalcore:changelog-v2-9-0`  |
+------------+------------+---------------------------------------+-------------------------------------+
| 2023-06-19 |            | ESMValTool `Feature Freeze`_          |                                     |
+------------+------------+---------------------------------------+-------------------------------------+
| 2023-06-26 | 2023-07-06 | :release:`v2.9.0` released            | :ref:`changelog-v2-9-0`             |
+------------+------------+---------------------------------------+-------------------------------------+

- 2.8.1 (Bugfix, Release Manager: `Valeriu Predoi`_)

+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|    Done    |                                            Event                                            |             Changelog              |
+============+=============================================================================================+====================================+
| 2023-06-02 | `ESMValCore Release 2.8.1 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.8.1>`_ | :ref:`esmvalcore:changelog-v2-8-1` |
+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.8.0 (Release Manager: `Rémi Kazeroni`_)

+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|  Planned   |    Done    |                                            Event                                            |             Changelog              |
+============+============+=============================================================================================+====================================+
| 2023-03-03 |            |                                  ESMValCore Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2023-03-20 | 2023-03-23 | `ESMValCore Release 2.8.0 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.8.0>`_ | :ref:`esmvalcore:changelog-v2-8-0` |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2023-03-17 |            |                                  ESMValTool Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2023-03-27 | 2023-03-28 | `ESMValTool Release 2.8.0 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.8.0>`_ |      :ref:`changelog-v2-8-0`       |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.7.1 (Bugfix, Release Manager: `Valeriu Predoi`_)

+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|    Done    |                                            Event                                            |             Changelog              |
+============+=============================================================================================+====================================+
| 2022-12-12 | `ESMValCore Release 2.7.1 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.7.1>`_ | :ref:`esmvalcore:changelog-v2-7-1` |
+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.7.0 (Release Manager: `Valeriu Predoi`_)

+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|  Planned   |    Done    |                                            Event                                            |             Changelog              |
+============+============+=============================================================================================+====================================+
| 2022-10-03 |            |                                  ESMValCore Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2022-10-10 | 2022-10-13 | `ESMValCore Release 2.7.0 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.7.0>`_ | :ref:`esmvalcore:changelog-v2-7-0` |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2022-10-17 |            |                                  ESMValTool Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2022-10-24 | 2022-10-28 | `ESMValTool Release 2.7.0 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.7.0>`_ |      :ref:`changelog-v2-7-0`       |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.6.0 (Release Manager: `Saskia Loosveldt Tomas`_)

+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|  Planned   |    Done    |                                            Event                                            |             Changelog              |
+============+============+=============================================================================================+====================================+
| 2022-06-06 |            |                                  ESMValCore Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2022-06-13 | 2022-07-15 | `ESMValCore Release 2.6.0 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.6.0>`_ | :ref:`esmvalcore:changelog-v2-6-0` |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2022-06-20 |            |                                  ESMValTool Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2022-06-27 | 2022-07-25 | `ESMValTool Release 2.6.0 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.6.0>`_ |      :ref:`changelog-v2-6-0`       |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.5.0 (Coordinating Release Manager: `Axel Lauer`_, team members: `Manuel Schlund`_, `Rémi Kazeroni`_)

+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|  Planned   |    Done    |                                            Event                                            |             Changelog              |
+============+============+=============================================================================================+====================================+
| 2022-02-07 |            |                                  ESMValCore Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2022-02-14 | 2022-03-14 | `ESMValCore Release 2.5.0 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.5.0>`_ | :ref:`esmvalcore:changelog-v2-5-0` |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2022-02-21 |            |                                  ESMValTool Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2022-02-28 | 2022-03-15 | `ESMValTool Release 2.5.0 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.5.0>`_ |      :ref:`changelog-v2-5-0`       |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.4.0 (Release Manager: `Klaus Zimmermann`_)

+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|  Planned   |    Done    |                                            Event                                            |             Changelog              |
+============+============+=============================================================================================+====================================+
| 2021-10-04 |            |                                  ESMValCore Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-10-11 | 2021-11-08 | `ESMValCore Release 2.4.0 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.4.0>`_ | :ref:`esmvalcore:changelog-v2-4-0` |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-10-18 |            |                                  ESMValTool Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-10-25 | 2021-11-09 | `ESMValTool Release 2.4.0 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.4.0>`_ |      :ref:`changelog-v2-4-0`       |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.3.1 (Bugfix, Release Manager: `Klaus Zimmermann`_)

+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|    Done    |                                            Event                                            |             Changelog              |
+============+=============================================================================================+====================================+
| 2021-07-23 | `ESMValCore Release 2.3.1 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.3.1>`_ | :ref:`esmvalcore:changelog-v2-3-1` |
+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.3.0 (Release Manager: `Klaus Zimmermann`_)

+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|  Planned   |    Done    |                                            Event                                            |             Changelog              |
+============+============+=============================================================================================+====================================+
| 2021-06-07 |            |                                  ESMValCore Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-06-14 | 2021-06-14 | `ESMValCore Release 2.3.0 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.3.0>`_ | :ref:`esmvalcore:changelog-v2-3-0` |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-06-21 |            |                                  ESMValTool Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-06-28 | 2021-07-27 | `ESMValTool Release 2.3.0 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.3.0>`_ |      :ref:`changelog-v2-3-0`       |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.2.0 (Release Manager: `Javier Vegas-Regidor`_)

+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|  Planned   |    Done    |                                            Event                                            |             Changelog              |
+============+============+=============================================================================================+====================================+
| 2021-02-01 |            |                                  ESMValCore Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-02-07 | 2021-02-09 | `ESMValCore Release 2.2.0 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.2.0>`_ | :ref:`esmvalcore:changelog-v2-2-0` |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-02-14 |            |                                  ESMValTool Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-02-21 | 2021-02-25 | `ESMValTool Release 2.2.0 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.2.0>`_ |      :ref:`changelog-v2-2-0`       |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.1.1 (Bugfix, Release Manager: `Valeriu Predoi`_)

+------------+---------------------------------------------------------------------------------------------+-------------------------+
|    Done    |                                            Event                                            |        Changelog        |
+============+=============================================================================================+=========================+
| 2020-12-01 | `ESMValTool Release 2.1.1 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.1.1>`_ | :ref:`changelog-v2-1-1` |
+------------+---------------------------------------------------------------------------------------------+-------------------------+

- 2.1.0 (Release Manager: `Valeriu Predoi`_)

+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|  Planned   |    Done    |                                            Event                                            |             Changelog              |
+============+============+=============================================================================================+====================================+
| 2020-10-05 |            |                                  ESMValCore Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2020-10-12 | 2020-10-12 | `ESMValCore Release 2.1.0 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.1.0>`_ | :ref:`esmvalcore:changelog-v2-1-0` |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2020-10-19 |            |                                  ESMValTool Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2020-10-26 | 2020-10-26 | `ESMValTool Release 2.1.0 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.1.0>`_ |      :ref:`changelog-v2-1-0`       |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.0.0 (Release Manager: `Bouwe Andela`_)

+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|  Planned   |    Done    |                                            Event                                            |             Changelog              |
+============+============+=============================================================================================+====================================+
| 2020-07-01 |            |                                  ESMValCore Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2020-07-20 | 2020-07-20 | `ESMValCore Release 2.0.0 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.0.0>`_ | :ref:`esmvalcore:changelog-v2-0-0` |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2020-07-22 |            |                                  ESMValTool Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2020-08-03 | 2020-08-03 | `ESMValTool Release 2.0.0 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.0.0>`_ |      :ref:`changelog-v2-0-0`       |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+



.. _release_steps:

Detailed timeline steps
~~~~~~~~~~~~~~~~~~~~~~~

These are the detailed steps to take to make a release.

#. Populate the milestone

   - The core development team will make sure it adds issues that it intends to work on as early as possible.
   - Any contributor is welcome to add issues or pull requests that they intend to work on themselves to a milestone.


#. ESMValCore feature freeze, testing, and release candidates

   - A release branch is created and branch protection rules are set up so only the release manager (i.e. the person in charge of the release branch) can push commits to that branch.
   - Make a release candidate with the release branch following the :ref:`ESMValCore release instructions <esmvalcore:how-to-make-a-release>`.
   - Uncomment the release candidate channel item (i.e. ``conda-forge/label/esmvalcore_rc``) in the ``environment.yml`` of ESMValTool to add it to the list of channels used. Adjust the pin on ESMValCore after each release candidate (e.g. ``esmvalcore==2.8.0rc1``). Check that the environment creation of ESMValTool works fine and contains the latest release candidate version.
   - Run all the recipes (optionally with a reduced amount of data) to check that they still work with the release candidate.
   - If a bug is discovered that needs to be fixed before the release, a pull request can be made to the main branch to fix the bug. The person making the pull request can then ask the release manager to cherry-pick that commit into the release branch.
   - Make another release candidate including the bugfix(es) and run the affected recipes again to check for further bugs.
   - Make as many release candidates for ESMValCore as needed in order to fix all the detected bugs.


#. ESMValTool feature freeze

   - A release branch is created and branch protection rules are set up so only the release manager (i.e. the person in charge of the release branch) can push commits to that branch.
   - The creation of the release branch is announced to the ESMValTool development team along with the procedures to use the branch for testing and making last-minute changes (see next step).


#. Some additional testing of ESMValTool

   - :ref:`Run all the recipes to check that they still work and generate the overview HTML pages <detailed_release_procedure>`.
   - Upload the results to the webpage at https://esmvaltool.dkrz.de/shared/esmvaltool/.
   - :ref:`Compare the results to those obtained with the previous release <compare_recipe_runs>`.
   - Create a `GitHub discussion <https://github.com/ESMValGroup/ESMValTool/discussions>`__ to communicate about the results.
   - If there are differences with the previous release, ask recipe maintainers
     or authors to review the plots and NetCDF files of their diagnostics, for
     example by
     `mentioning <https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax#mentioning-people-and-teams>`__
     them in the discussion.
   - If a bug is discovered that needs to be fixed before the release, a pull request can be made to the main branch to fix the bug. The person making the pull request can then ask the release manager to cherry-pick that commit into the release branch.
   - Update the :ref:`list of broken recipes <broken-recipe-list>` with new recipes that could not be run successfully during the testing.
     Open a separate GitHub issue for each failing recipe and assign the next milestone.
     Open an overview issue, see :issue:`3484` for an example, and review past overview issues.
     Take action to ensure that the broken recipe policy is followed.


#. ESMValCore release

   - Make the official ESMValCore release with the last release candidate by following the :ref:`ESMValCore release instructions <esmvalcore:how-to-make-a-release>`.


#. ESMValTool release

   - Pin ESMValCore to the same version as ESMValTool in the ``environment.yml`` and on `conda-forge
     <https://github.com/conda-forge/esmvaltool-suite-feedstock>`__.
     This way, we make sure that ESMValTool uses the ESMValCore version with which it has been tested.
     Make sure to comment again the release candidate channel once ESMValCore has been released.
   - Make the release by following :ref:`How to make a release`.


#. Announce the releases

   - Ask the user engagement team to announce the releases to the user mailing list, the development team mailing list, and on twitter.


#. Core development team meets to coordinate the content of next milestone

   - Create a doodle for the meeting or even better, have the meeting during an ESMValTool workshop
   - Prepare the meeting by filling the milestone
   - At the meeting, discuss

     - If the proposed issues cover everything we would like to accomplish
     - Are there things we need to change about the release process
     - Who will be the release manager(s) for the next release

Bugfix releases
---------------

Next to the feature releases described above, it is also possible to have bugfix releases (2.0.1, 2.0.2, etc). In general bugfix releases will only be done on the latest release, and may include ESMValCore, ESMValTool, or both.


Procedure
~~~~~~~~~

#. One or more issues are resolved that are deemed (by the core development team) to warrant a bugfix release.
#. A release branch is created from the last release tag and the commit that fixes the bug/commits that fix the bugs are cherry-picked into it from the main branch.
#. Some additional testing of the release branch takes place.
#. The release takes place.

Compatibility between ESMValTool and ESMValCore is ensured by the appropriate version pinning of ESMValCore by ESMValTool.

Glossary
--------

Feature freeze
~~~~~~~~~~~~~~
The date on which no new features may be submitted for the upcoming release.
After this date, only critical bug fixes can still be included to the :ref:`release_branch`.
Development work can continue in the main branch.
If you are unsure whether new developments could interfere with the release, check with the :ref:`release_manager`.


Milestone
~~~~~~~~~
A milestone is a list of issues and pull-request on GitHub. It has a due date, this date is the date of the feature freeze. Adding an issue or pull request indicates the intent to finish the work on this issue before the due date of the milestone. If the due date is missed, the issue can be included in the next milestone.

.. _release_manager:

Release manager
~~~~~~~~~~~~~~~
The person in charge of making the release, both technically and organizationally. Appointed for a single release.
Check the :ref:`release_schedule` to see who is the manager of the next release.

.. _release_branch:

Release branch
~~~~~~~~~~~~~~
The release branch can be used to do some additional testing before the release, while normal development work continues in the main branch. It will be branched off from the main branch after the feature freeze and will be used to make the release on the release date. The only way to still get something included in the release after the feature freeze is to ask the release manager to cherry-pick a commit from the main branch into this branch.


.. _How to make a release:

How to make an ESMValTool release
---------------------------------

Before the actual release, a number of tests, and pre-release steps must be performed,
a detailed workflow description can be found here :ref:`detailed_release_procedure`.

The release manager makes the release, assisted by the release manager of the
previous release, or if that person is not available, another previous release
manager.
Perform the steps listed below with two persons, to reduce the risk of
error.

.. note::

   The previous release manager ensures the current release manager has the
   required administrative permissions to make the release.
   Consider the following services:
   `conda-forge <https://github.com/conda-forge/esmvaltool-suite-feedstock>`__,
   `DockerHub <https://hub.docker.com/orgs/esmvalgroup>`__,
   `PyPI <https://pypi.org/project/ESMValTool/>`__, and
   `readthedocs <https://readthedocs.org/dashboard/esmvaltool/users/>`__.

The release of ESMValTool should come after the release of ESMValCore.
To make a new release of the package, follow these steps:

1. Check that all tests and builds work
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Check that the ``nightly``
  `test run on CircleCI <https://circleci.com/gh/ESMValGroup/ESMValTool/tree/main>`__
  was successful.
- Check that the
  `GitHub Actions test runs <https://github.com/ESMValGroup/ESMValTool/actions>`__
  were successful.
- Check that the documentation builds successfully on
  `readthedocs <https://readthedocs.org/projects/esmvaltool/builds/>`__.
- Check that the
  `Docker images <https://hub.docker.com/repository/docker/esmvalgroup/esmvaltool/builds>`__
  are building successfully.

All tests should pass before making a release (branch).

2. Increase the version number
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The version number is automatically generated from the information provided by
git using `setuptools-scm <https://pypi.org/project/setuptools-scm/>`__, but a
static version number is stored in ``CITATION.cff``.
Make sure to update the version number and release date in ``CITATION.cff``.
See https://semver.org for more information on choosing a version number.
Make sure that the ESMValCore version that is being used is set to the latest version.
See the :ref:`dependencies <dependencies>` section in order to find more details on how update the ESMValCore version.
Make a pull request and get it merged into ``main``.

.. _add-release-notes:

3. Add release notes
~~~~~~~~~~~~~~~~~~~~
Use the script :ref:`draft_release_notes.py` to create a draft of the
release notes.
This script uses the titles and labels of merged pull requests since the
previous release.
Open a discussion to allow members of the development team to nominate pull requests
as highlights. Add the most voted pull requests as highlights at the beginning of
changelog.
After the highlights section, list any backward incompatible changes that the
release may include.
The :ref:`backward compatibility policy <guidance-on-releasing-backward-incompatible-changes>`
lists the information that should be provided by the developer of any backward
incompatible change.
Make sure to also list any deprecations that the release may include, as well
as a brief description on how to upgrade a deprecated feature.
Review the results, and if anything needs changing, change it on GitHub and
re-run the script until the changelog looks acceptable.
Copy the result to the file ``doc/sphinx/source/changelog.rst``.
If possible, try to set the script dates to the date of the release
you are managing.
Make a pull request and get it merged into ``main``.

4. Create a release branch
~~~~~~~~~~~~~~~~~~~~~~~~~~
Create a branch off the ``main`` branch and push it to GitHub.
Ask someone with administrative permissions to set up branch protection rules
for it so only you and the person helping you with the release can push to it.
Announce the name of the branch in an issue and ask the members of the
`ESMValTool development team <https://github.com/orgs/ESMValGroup/teams/esmvaltool-developmentteam>`__
to run their favourite recipe using this branch.

5. Make the release on GitHub
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Do a final check that all tests on CircleCI and GitHub Actions completed
successfully.
Then click the
`releases tab <https://github.com/ESMValGroup/ESMValTool/releases>`__
and create the new release from the release branch (i.e. not from ``main``).
The release tag always starts with the letter ``v`` followed by the version
number, e.g. ``v2.1.0``.

6. Merge the release branch back into the main branch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the (pre-)release is tagged, it is time to merge the release branch back into `main`.
We do this for two reasons, namely, one, to mark the point up to which commits in `main`
have been considered for inclusion into the present release, and, two, to inform
setuptools-scm about the version number so that it creates the correct version number in
`main`.
However, unlike in a normal merge, we do not want to integrate any of the changes from the
release branch into main.
This is because all changes that should be in both branches, i.e. bug fixes, originate from
`main` anyway and the only other changes in the release branch relate to the release itself.
To take this into account, we perform the merge in this case on the command line using `the
ours merge strategy <https://git-scm.com/docs/merge-strategies#Documentation/merge-strategies.txt-ours-1>`__
(``git merge -s ours``), not to be confused with the ``ours`` option to the ort merge strategy
(``git merge -X ours``).
For details about merge strategies, see the above-linked page.
To execute the merge use following sequence of steps

.. code-block:: bash

   git fetch
   git checkout main
   git pull
   git merge -s ours v2.1.x
   git push

Note that the release branch remains intact and you should continue any work on the release
on that branch.

7. Create and upload the PyPI package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The package is automatically uploaded to the
`PyPI <https://pypi.org/project/ESMValTool/>`__
by a GitHub action.
If has failed for some reason, build and upload the package manually by
following the instructions below.

Follow these steps to create a new Python package:

-  Check out the tag corresponding to the release,
   e.g. ``git checkout tags/v2.1.0``
-  Make sure your current working directory is clean by checking the output
   of ``git status`` and by running ``git clean -xdf`` to remove any files
   ignored by git.
-  Install the required packages:
   ``python3 -m pip install --upgrade pep517 twine``
-  Build the package:
   ``python3 -m pep517.build --source --binary --out-dir dist/ .``
   This command should generate two files in the ``dist`` directory, e.g.
   ``ESMValTool-2.1.0-py3-none-any.whl`` and ``ESMValTool-2.1.0.tar.gz``.
-  Upload the package:
   ``python3 -m twine upload dist/*``
   You will be prompted for an API token if you have not set this up
   before, see
   `here <https://pypi.org/help/#apitoken>`__ for more information.

You can read more about this in
`Packaging Python Projects <https://packaging.python.org/tutorials/packaging-projects/>`__.

8. Create the Conda package
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``esmvaltool`` package is published on the `conda-forge conda channel
<https://anaconda.org/conda-forge>`__.
This is done via a pull request on the `esmvaltool-suite-feedstock repository
<https://github.com/conda-forge/esmvaltool-suite-feedstock>`__.

After the upload of the PyPI package, this pull request is automatically opened
by a bot.
An example pull request can be found `here
<https://github.com/conda-forge/esmvaltool-suite-feedstock/pull/5>`__.
Follow the instructions by the bot to finalize the pull request.
This step mostly contains updating dependencies that have been changed during
the last release cycle.
Once approved by the `feedstock maintainers
<https://github.com/conda-forge/esmvaltool-suite-feedstock#feedstock-maintainers>`__
they will merge the pull request, which will in turn publish the package on
conda-forge some time later.
Contact the feedstock maintainers if you want to become a maintainer yourself.

9. Check the Docker images
~~~~~~~~~~~~~~~~~~~~~~~~~~

There are three main Docker container images available for ESMValTool on
`Dockerhub <https://hub.docker.com/r/esmvalgroup/esmvaltool/tags>`_:

- ``esmvalgroup/esmvaltool:stable``, built from `docker/Dockerfile <https://github.com/ESMValGroup/ESMValTool/blob/main/docker/Dockerfile>`_,
  this is a tag that is always the same as the latest released version.
  This image is only built by Dockerhub when a new release is created.
- ``esmvalgroup/esmvaltool:development``, built from `docker/Dockerfile.dev <https://github.com/ESMValGroup/ESMValTool/blob/main/docker/Dockerfile.dev>`_,
  this is a tag that always points to the latest development version of
  ESMValTool.
  This image is built by Dockerhub every time there is a new commit to the
  ``main`` branch on Github.
- ``esmvalgroup/esmvaltool:experimental``, built from `docker/Dockerfile.exp <https://github.com/ESMValGroup/ESMValTool/blob/main/docker/Dockerfile.exp>`_,
  this is a tag that always points to the latest development version of
  ESMValTool with the latest development version of ESMValCore.
  Note that some recipes may not work as expected with this image because
  the ESMValTool development version has been designed to work with the latest
  release of ESMValCore (i.e. not with the development version).
  This image is built by Dockerhub every time there is a new commit to the
  ESMValTool ``main`` branch on Github.

In addition to the three images mentioned above, there is an image available
for every release (e.g. ``esmvalgroup/esmvaltool:v2.5.0``).
When working on the Docker images, always try to follow the
`best practices <https://docs.docker.com/develop/develop-images/dockerfile_best-practices/>`__.

After making the release, check that the Docker image for that release has been
built correctly by

1. checking that the version tag is available on `Dockerhub`_ and the ``stable``
   tag has been updated,
2. running some recipes with the ``stable`` tag Docker container, for example one
   recipe for Python, NCL, R, and Julia,
3. running a recipe with a Singularity container built from the ``stable`` tag.

If there is a problem with the automatically built container image, you can fix
the problem and build a new image locally.
For example, to
`build <https://docs.docker.com/engine/reference/commandline/build/>`__ and
`upload <https://docs.docker.com/engine/reference/commandline/push/>`__
the container image for v2.5.0 of the tool run:

.. code-block:: bash

   git checkout v2.5.0
   git clean -x
   docker build -t esmvalgroup/esmvaltool:v2.5.0 . -f docker/Dockerfile
   docker push esmvalgroup/esmvaltool:v2.5.0

and if it is the latest release that you are updating, also run

.. code-block:: bash

   docker tag esmvalgroup/esmvaltool:v2.5.0 esmvalgroup/esmvaltool:stable
   docker push esmvalgroup/esmvaltool:stable

Note that the ``docker push`` command will overwrite the existing tags on
Dockerhub.

If you would like to make a small change to an existing Docker container image,
it is also possible to do just that using the
`docker commit <https://docs.docker.com/engine/reference/commandline/commit/>`__
command.
Note that this is only recommended for very small changes, as it is not
reproducible and it will add an extra layer, increasing the size of the image.
To do this, start the container with
``docker run -it --entrypoint /bin/bash esmvalgroup/esmvaltool:v2.5.0``
and make your changes.
Exit the container by pressing `ctrl+d` and find it back by running
``docker ps -a``.
Find the `CONTAINER ID` of the image you would like to save and run
``docker commit -c 'ENTRYPOINT ["conda", "run", "--name", "esmvaltool", "esmvaltool"]' 633696a8b53a esmvalgroup/esmvaltool:v2.5.0``
where ``633696a8b53c`` is the an example of a container ID, replace it by
by the actual ID.

Changelog
---------
- 2020-09-09 Converted to rst and added to repository (future changes tracked by git)
- 2020-09-03 Update during video conference (present: Bouwe Andela, Niels Drost, Javier Vegas, Valeriu Predoi, Klaus Zimmermann)
- 2020-07-27 Update including tidying up and Glossary by Klaus Zimmermann and Bouwe Andela
- 2020-07-23 Update to timeline format by Bouwe Andela and Klaus Zimmermann
- 2020-06-08 First draft by Klaus Zimmermann and Bouwe Andela

.. _Bouwe Andela: https://github.com/bouweandela
.. _Rémi Kazeroni: https://github.com/remi-kazeroni
.. _Axel Lauer: https://github.com/axel-lauer
.. _Saskia Loosveldt Tomas: https://github.com/sloosvel
.. _Valeriu Predoi: https://github.com/valeriupredoi
.. _Manuel Schlund: https://github.com/schlunma
.. _Javier Vegas-Regidor: https://github.com/jvegasbsc
.. _Klaus Zimmermann: https://github.com/zklaus
.. _Emma Hogan: https://github.com/ehogan
.. _Chris Billows: https://github.com/chrisbillowsMO
.. _Ed Gillett: https://github.com/mo-gill
.. _Julien Lenhardt: https://github.com/jlenh
