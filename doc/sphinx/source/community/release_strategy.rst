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
:ref:`contact`.


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

- 2.4.0 (Release Manager: `Klaus Zimmermann`_)

+------------+--------------------------+
| 2021-10-04 |ESMValCore feature freeze |
+------------+--------------------------+
| 2021-10-11 |ESMValCore release        |
+------------+--------------------------+
| 2021-10-18 |ESMValTool feature freeze |
+------------+--------------------------+
| 2021-10-25 |ESMValTool release        |
+------------+--------------------------+

- 2.5.0 (Release Manager: `Manuel Schlund`_)

+------------+--------------------------+
| 2022-02-07 |ESMValCore feature freeze |
+------------+--------------------------+
| 2022-02-14 |ESMValCore release        |
+------------+--------------------------+
| 2022-02-21 |ESMValTool feature freeze |
+------------+--------------------------+
| 2022-02-28 |ESMValTool release        |
+------------+--------------------------+

- 2.6.0 (Release Manager: TBD)

+------------+--------------------------+
| 2022-06-06 |ESMValCore feature freeze |
+------------+--------------------------+
| 2022-06-13 |ESMValCore release        |
+------------+--------------------------+
| 2022-06-20 |ESMValTool feature freeze |
+------------+--------------------------+
| 2022-06-27 |ESMValTool release        |
+------------+--------------------------+

Past releases
^^^^^^^^^^^^^

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



Detailed timeline steps
~~~~~~~~~~~~~~~~~~~~~~~

These are the detailed steps to take to make a release.

1. Populate the milestone

   - The core development team will make sure it adds issues that it intends to work on as early as possible.
   - Any contributor is welcome to add issues or pull requests that they intend to work on themselves to a milestone.


2. ESMValCore feature freeze

   - A release branch is created and branch protection rules are set up so only the release manager (i.e. the person in charge of the release branch) can push commits to that branch.
   - The creation of the release branch is announced to the ESMValTool development team along with the procedures to use the branch for testing and making last-minute changes (see next step)


3. Some additional testing of ESMValCore

   - Run all the recipes (optionally with a reduced amount of data) to check that they still work
   - If a bug is discovered that needs to be fixed before the release, a pull request can be made to the main branch to fix the bug. The person making the pull request can then ask the release manager to cherry-pick that commit into the release branch.


4. ESMValCore release

   - Make the release by following the :ref:`ESMValCore release instructions <esmvalcore:how-to-make-a-release>`.
   - Ask the user engagement team to announce the release to the user mailing list, the development team mailing list, on twitter


5. ESMValTool feature freeze

   - A release branch is created and branch protection rules are set up so only the release manager (i.e. the person in charge of the release branch) can push commits to that branch.
   - The creation of the release branch is announced to the ESMValTool development team along with the procedures to use the branch for testing and making last-minute changes (see next step)


6. Some additional testing of ESMValTool

   - Run all the recipes to check that they still work and ask authors to review the plots
   - If a bug is discovered that needs to be fixed before the release, a pull request can be made to the main branch to fix the bug. The person making the pull request can then ask the release manager to cherry-pick that commit into the release branch.


7. ESMValTool release

   - Make the release by following :ref:`How to make a release`
   - Ask the user engagement team to announce the release to the user mailing list, the development team mailing list, and on twitter


8. Core development team meets to coordinate the content of next milestone

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

1. One or more issues are resolved that are deemed (by the core development team) to warrant a bugfix release.
2. A release branch is created from the last release tag and the commit that fixes the bug/commits that fix the bugs are cherry-picked into it from the main branch.
3. Some additional testing of the release branch takes place.
4. The release takes place.

Compatibility between ESMValTool and ESMValCore is ensured by the appropriate version pinning of ESMValCore by ESMValTool.

Glossary
--------

Feature freeze
~~~~~~~~~~~~~~
The date on which no new features may be submitted for the upcoming release. After this date, only critical bug fixes can still be included.

Milestone
~~~~~~~~~
A milestone is a list of issues and pull-request on GitHub. It has a due date, this date is the date of the feature freeze. Adding an issue or pull request indicates the intent to finish the work on this issue before the due date of the milestone. If the due date is missed, the issue can be included in the next milestone.

Release manager
~~~~~~~~~~~~~~~
The person in charge of making the release, both technically and organizationally. Appointed for a single release.

Release branch
~~~~~~~~~~~~~~
The release branch can be used to do some additional testing before the release, while normal development work continues in the main branch. It will be branched off from the main branch after the feature freeze and will be used to make the release on the release date. The only way to still get something included in the release after the feature freeze is to ask the release manager to cherry-pick a commit from the main branch into this branch.


.. _How to make a release:

How to make an ESMValTool release
---------------------------------

The release manager makes the release, assisted by the release manager of the
previous release, or if that person is not available, another previous release
manager. Perform the steps listed below with two persons, to reduce the risk of
error.

To make a new release of the package, follow these steps:

1. Check the tests on GitHub Actions and CircleCI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check the ``nightly``
`build on CircleCI <https://circleci.com/gh/ESMValGroup/ESMValTool/tree/main>`__
and the
`GitHub Actions run <https://github.com/ESMValGroup/ESMValTool/actions>`__.
All tests should pass before making a release (branch).

2. Increase the version number
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The version number is stored in ``esmvaltool/__init__.py``,
``package/meta.yaml``, ``CITATION.cff``. Make sure to update all files.
Also update the release date in ``CITATION.cff``.
See https://semver.org for more information on choosing a version number.
Make a pull request and get it merged into ``main``.

3. Add release notes
~~~~~~~~~~~~~~~~~~~~
Use the script :ref:`draft_release_notes.py` to create create a draft of the
release notes.
This script uses the titles and labels of merged pull requests since the
previous release.
Review the results, and if anything needs changing, change it on GitHub and
re-run the script until the changelog looks acceptable.
Copy the result to the file ``doc/sphinx/source/changelog.rst``.
Make a pull request and get it merged into ``main``.

4. Create a release branch
~~~~~~~~~~~~~~~~~~~~~~~~~~
Create a branch off the ``main`` branch and push it to GitHub.
Ask someone with administrative permissions to set up branch protection rules
for it so only you and the person helping you with the release can push to it.
Announce the name of the branch in an issue and ask the members of the
`ESMValTool development team <https://github.com/orgs/ESMValGroup/teams/esmvaltool-developmentteam>`__
to run their favourite recipe using this branch.

5. Cherry pick bugfixes into the release branch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If a bug is found and fixed (i.e. pull request merged into the
``main`` branch) during the period of testing, use the command
``git cherry-pick COMMIT_HASH``, where ``COMMIT_HASH`` is the commit hash of the
commit that needs to be cherry-picked, to include the commit for this bugfix
into the release branch.
Cherry-pick any new contributions in the order they were merged, to avoid
conflicts.
When the testing period is over, make a pull request to update
the release notes with the latest changes (do not forget to include the pull
request itself into the changelog), get it merged into ``main`` and
cherry-pick it into the release branch.

6. Make the release on GitHub
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Do a final check that all tests on CircleCI and GitHub Actions completed
successfully.
Then click the
`releases tab <https://github.com/ESMValGroup/ESMValTool/releases>`__
and create the new release from the release branch (i.e. not from ``main``).
The release tag always starts with the letter ``v`` followed by the version
number, e.g. ``v2.1.0``.

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

8. Update the conda-forge packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The upload to PyPI will automatically trigger an update PR on the
esmvaltool-suite-feedstock_. Check that it builds correctly and merge
the PR to update the conda-forge packages.

.. _esmvaltool-suite-feedstock: https://github.com/conda-forge/esmvaltool-suite-feedstock

Changelog
---------

- 2020-09-09 Converted to rst and added to repository (future changes tracked by git)
- 2020-09-03 Update during video conference (present: Bouwe Andela, Niels Drost, Javier Vegas, Valeriu Predoi, Klaus Zimmermann)
- 2020-07-27 Update including tidying up and Glossary by Klaus Zimmermann and Bouwe Andela
- 2020-07-23 Update to timeline format by Bouwe Andela and Klaus Zimmermann
- 2020-06-08 First draft by Klaus Zimmermann and Bouwe Andela

.. _Bouwe Andela: https://github.com/bouweandela
.. _Valeriu Predoi: https://github.com/valeriupredoi
.. _Manuel Schlund: https://github.com/schlunma
.. _Javier Vegas-Regidor: https://github.com/jvegasbsc
.. _Klaus Zimmermann: https://github.com/zklaus
