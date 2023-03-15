.. _backward-compatibility-policy:

ESMValTool policy on backward compatibility
===========================================

Motivation
----------

Development of recipes or conducting project-related work may require a
rather long period of time during which new versions of the ESMValTool
might become available. For a good user experience and a seamless
workflow, users and developers need to know before upgrading to a new
version if and how their work might be affected (backward
compatibility). This includes, for instance, information about changes
to technical features such as syntax of recipes and configuration files,
and interfaces of shared functions, but also changes that affect the
results of an ESMValTool run, e.g. modification of algorithms or changes
to the order of operators. It is therefore essential that users and
developers have the best advice on how and when to upgrade to new
versions.

While trying to minimise the impact of new code developments on users
and developers by maintaining backward compatibility where possible,
this cannot always be guaranteed. A very restrictive policy might delay
the ESMValTool development and make it more complex for developers to
contribute.

This document outlines the key principles of an updated ESMValTool policy
on backward compatibility.

Definitions
-----------

**Release:** A numbered version of ESMValCore / ESMValTool that has been
released to the community, e.g. 2.4.0. This policy relates only to
backward compatibility of releases, not to interim revisions of the main
branch. Release numbers are of the format x.y.z, where:

-  x indicates a major release
-  y indicates a minor release
-  z indicates a patch release

**Backward-incompatible change:** A change in ESMValCore or ESMValTool that causes a
recipe to no longer run successfully (a *breaking change*), or which
results in scientifically significant changes in results (a *science
change*).

**Breaking change:** A change which causes a previously working recipe
to no longer run successfully.

**Science change:** A change that alters scientific results. We do not
formally distinguish between trivial science changes (e.g. from changes
in the order of calculations) and more significant changes that would
affect interpretation, although the detail that we communicate will
share any understanding that we have regarding expected impact.

**Benign third-party dependency changes:** A change over which we have
no control, but which we believe will only have trivial technical
impacts (such as a change in font). Such changes are outside of the
scope of this policy, though we will communicate about those we are
aware of.

**Developer of backward-incompatible change:** For the purpose of this
policy, developer is the individual that is responsible for the pull
request (PR) that is not backward compatible.

**Recipe developer:** Someone who is developing a recipe that is not
(yet) integrated into the repository.

**Recipe user:** For the purpose of this policy, a *recipe user* is
anyone who runs a recipe using a *release* of ESMValTool. In this
context, someone can be both a *recipe developer* and a *recipe user*,
but they perform different activities in each capacity.

**Recipe maintainer:** First contact point for *integrated recipes* in
case of problems with that recipe (see also :ref:`Maintaining a recipe<recipe-maintainer>`).

**Integrated recipes:** Recipes that are contained within the main
branch of the ESMValTool repository, and can therefore be updated by any
developer in line with the above guidance. Note that the recipe can be
updated by someone other than the original author.

**User recipes:** Recipes developed by any developer outside of the main
branch of the repository (i.e. on a dev/feature branch or outside the
repository completely), and therefore cannot be updated by anyone else.

Scope
-----

The ESMValTool and ESMValCore policy on backward compatibility aims at balancing two
competing needs: the occasional need of improvements or maintenance to
break backward compatibility and the need for stability for existing
users and developers. The following aspects are covered by this policy:

-  Key principles and approaches to backward compatibility
-  Guidelines and requirements for *developers of backward-incompatible
   changes*
-  Communication with users and developers about *backward-incompatible
   changes*

Not within the scope of this policy are:

-  Versioning scheme of ESMValTool
-  Breakage of recipes due to changes in input data or dependencies.
   This is covered by the :ref:`broken recipe policy<broken-recipe-policy>`.

Expectations of developers, users & funders
-------------------------------------------

Stakeholders and their expectations and aims:

Projects / Funders

-  Aim to facilitate scientific discovery
-  Expect deliverables, e.g. new features/recipes
-  Expect reproducible results

*Recipe users*

-  Expect the recipe to work
-  Expect the recipe to be easy to run
-  Expect reproducible results
-  Expect easy installation of ESMValTool

*Recipe developers*

-  Develop recipes
-  Expect their recipe to keep working with every new *release* of
   ESMValCore
-  Expect ESMValCore bugfixes and new features to become available
   quickly
-  Expect reproducible results
-  Expect easy installation of ESMValTool

Core developers and *recipe maintainers*

-  Fix bugs
-  Add ESMValCore features requested by *recipe developers*
-  Try to accommodate ESMValCore features contributed to by *recipe
   developers*
-  Maintain existing recipes
-  Add new recipes
-  Try to help (other) *recipe developers* with contributing their
   recipe
-  Try to make installation as easy as possible

There is a tension between making new features available and keeping
everything as is. New features facilitate scientific discovery because
they enable *recipe developers* to do new research (e.g. analyse more
data, new data, or perform a different analysis). Ensuring that every
recipe ever made works with every new feature is technically a lot of
work, more than we have funding for. Therefore we need to make sure that
new features are added regularly, but we respect the timescale on which
*recipe developers* work when removing outdated features. Writing a
paper and getting it published may take up to a year, so this seems a
good timescale for larger changes. For changes that only affect a few
users, shorter timescales could be acceptable. It is also good to note
that we are part of a large software ecosystem (ESMValTool currently
depends on over 500 different software packages), so we may not always
be able to control at what pace changes are made to the software that we
depend upon.

Two-way communication about new and removed features is needed to make
this work. This requires active involvement from both the people
developing the new features and the *recipe developers*. ESMValTool core
developers and ESMValCore core developers need to make sure they clearly
communicate changes. In the first place, this is done by writing good
descriptions in issues and pull requests on GitHub, but some of this
material also makes it to the changelog (where the GitHub pull requests
are linked). It is highly recommended to communicate a relevant
selection (e.g. important new, scheduled for removal, and removed
features) also by other means, to ensure we reach as many people
potentially affected as possible (see :ref:`Guidance on handling
*backward-incompatible changes*<guidance-on-backward-incompatiable-changes>`
section below).
We organize :ref:`monthly community <monthly-meetings>` meetings where
*recipe developers* can learn about the latest developments and everyone is
welcome to join, ask questions, and provide feedback.

To meet the needs of users and funders, we should take reproducibility
of older results seriously, but this should not hold us back from
advancing our tools. We can support this by uploading a well tested
container image to an archive that provides a DOI and by providing clear
instructions on how to use such containers.

Helping developers to upgrade
-----------------------------

*Recipe users* of ESMValTool should be able to successfully run
*integrated recipes* using a *release*, since all
*backward-incompatible changes* introduced between *releases* will have
been fixed before the *release* is created. Please note the
:ref:`broken recipe policy<broken-recipe-policy>`.

However, *recipe developers* working on *user recipes* must be provided
with information to enable them to adapt their code to resolve issues
related to *backward-incompatible changes* when *backward-incompatible
changes* are introduced to the main branch / when a *release* of
ESMValTool is created.

.. _guidance-on-backward-incompatiable-changes:

Guidance on handling *backward-incompatible changes*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As well as helping users to handle *backward-incompatible changes*, the
policy and surrounding tools must help developers avoid making
*backward-incompatible changes*. Not many ideas are developed on this yet,
but components should include:

-  Testing; *backward-incompatible changes* should be discovered as
   early in the development process as possible. This motivates
   continued investment in automated testing.
   To discover *backward-incompatible changes* early on in the development cycle,
   every night a selection of recipes is run on
   `CircleCI <https://app.circleci.com/pipelines/github/ESMValGroup/ESMValTool?branch=main>`__.
   A recipe can be added to the test suite by adding it to the directory
   `esmvaltool/recipes/testing <https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/recipes/testing>`__.
   Only add recipes that require a small amount of data, i.e. considerably less
   than a gigabyte.
-  Guidance on how to minimise the likelihood of introducing
   *backward-incompatible changes* and how to use deprecation warnings
   when needed (see :ref:`developer guidance <esmvalcore:backward_compatibility>`).
-  :ref:`Instructions on how to provide text for the release notes <add-release-notes>`
   to assist *recipe developers* to adapt their recipe in light of the
   *backward-incompatible change*
-  General instructions for *recipe developers* working on *user
   recipes* to enable them to adapt their code related to
   *backward-incompatible changes* (see `ESMValTool_Tutorial: issue
   #263 <https://github.com/ESMValGroup/ESMValTool_Tutorial/issues/263>`__).
-  The developer or reviewer must tag the core development team to
   notify them of the *backward-incompatible change*, and give at least
   2 weeks for objections to be raised before merging to the main
   branch. If a strong objection is raised the backward-incompatible
   change should not be merged until the objection is resolved.


.. _guidance-on-releasing-backward-incompatible-changes:

Guidance on releasing *backward-incompatible changes*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

During the *release* process, the following information must be
provided:

-  **Release notes:** The *release* notes are already documented in the
   :ref:`ESMValTool Changelog <changelog>` and
   :ref:`ESMValCore Changelog <esmvalcore:changelog>`, and
   “*backward-incompatible changes*” is the first section after
   “Highlights”.

   -  **backward-incompatible changes:** This section must include
      clear instructions detailing how a *recipe developer* should adapt
      their code for each item in this section, whether the adapted code
      would introduce a *science change*, and the list of affected or
      fixed *integrated recipes* that had to be updated due to the
      *backward-incompatible changes*, if applicable (to provide
      further examples to *recipe developers* working on *user recipes*
      of how to adapt code).
   -  **Developer guidance:** *Developers* *of backward-incompatible
      changes* must:

      -  write and include the information required for the
         “*backward-incompatible changes*” section in the PR that
         introduces the *backward-incompatible change*
      -  share details of the *backward-incompatible change* at the
         next monthly ESMValTool community meeting

   -  **Communication:** The *release* notes must be shared with the
      community (for example, via the :ref:`mailing-list` and the
      `Community <https://github.com/ESMValGroup/Community>`__
      repository) at the point the first *release* candidate is made,
      highlighting the “*backward-incompatible changes*” section. The
      User Engagement Team should organise the communication of new
      *releases* together with the :ref:`release_manager`.
