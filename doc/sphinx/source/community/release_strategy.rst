Release Strategy for ESMValCore and ESMValTool
================================================

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
If your project requires longer maintenance or you have other concerns about the release strategy, please contact the ESMValTool core development team.


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


Release schedule
~~~~~~~~~~~~~~~~

With the following release schedule, we strive to have three releases per year and to avoid releases too close to holidays, as well as avoiding weekends.


- 2.0.0 (Release Manager: Bouwe Andela)

+------------+--------------------------+
| 2020-07-01 |ESMValCore feature freeze |
+------------+--------------------------+
| 2020-07-20 |ESMValCore release        |
+------------+--------------------------+
| 2020-07-22 |ESMValTool feature freeze |
+------------+--------------------------+
| 2020-08-03 |ESMValTool release        |
+------------+--------------------------+

- 2.1.0 (Release Manager: Valeriu Predoi)

+------------+--------------------------+
| 2020-10-05 |ESMValCore feature freeze |
+------------+--------------------------+
| 2020-10-12 |ESMValCore release        |
+------------+--------------------------+
| 2020-10-19 |ESMValTool feature freeze |
+------------+--------------------------+
| 2020-10-26 |ESMValTool release        |
+------------+--------------------------+

- 2.2.0 (Release Manager: tbd)

+------------+--------------------------+
| 2021-02-01 |ESMValCore feature freeze |
+------------+--------------------------+
| 2021-02-07 |ESMValCore release        |
+------------+--------------------------+
| 2021-02-14 |ESMValTool feature freeze |
+------------+--------------------------+
| 2021-02-21 |ESMValTool release        |
+------------+--------------------------+

- 2.3.0 (Release Manager: tbd)

+------------+--------------------------+
| 2021-06-07 |ESMValCore feature freeze |
+------------+--------------------------+
| 2021-06-14 |ESMValCore release        |
+------------+--------------------------+
| 2021-06-21 |ESMValTool feature freeze |
+------------+--------------------------+
| 2021-06-28 |ESMValTool release        |
+------------+--------------------------+


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
   - If a bug is discovered that needs to be fixed before the release, a pull request can be made to the master branch to fix the bug. The person making the pull request can then ask the release manager to cherry-pick that commit into the release branch.


4. ESMValCore release

   - Make the release by following the `ESMValCore release instructions`_.
   - Ask the user engagement team to announce the release to the user mailing list, the development team mailing list, on twitter


5. ESMValTool feature freeze

   - A release branch is created and branch protection rules are set up so only the release manager (i.e. the person in charge of the release branch) can push commits to that branch.
   - The creation of the release branch is announced to the ESMValTool development team along with the procedures to use the branch for testing and making last-minute changes (see next step)




6. Some additional testing of ESMValTool

   - Run all the recipes to check that they still work and ask authors to review the plots
   - If a bug is discovered that needs to be fixed before the release, a pull request can be made to the master branch to fix the bug. The person making the pull request can then ask the release manager to cherry-pick that commit into the release branch.


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
2. A release branch is created from the last release tag and the commit that fixes the bug/commits that fix the bugs are cherry-picked into it from the master branch.
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
The release branch can be used to do some additional testing before the release, while normal development work continues in the master branch. It will be branched off from the master branch after the feature freeze and will be used to make the release on the release date. The only way to still get something included in the release after the feature freeze is to ask the release manager to cherry-pick a commit from the master branch into this branch.

Changelog
---------

- 2020-09-09 Converted to rst and added to repository (future changes tracked by git)
- 2020-09-03 Update during video conference (present: Bouwe Andela, Niels Drost, Javier Vegas, Valeriu Predoi, Klaus Zimmermann)
- 2020-07-27 Update including tidying up and Glossary by Klaus Zimmermann and Bouwe Andela
- 2020-07-23 Update to timeline format by Bouwe Andela and Klaus Zimmermann
- 2020-06-08 First draft by Klaus Zimmermann and Bouwe Andela


.. _ESMValCore release instructions: https://docs.esmvaltool.org/projects/esmvalcore/en/latest/contributing.html#how-to-make-a-release
