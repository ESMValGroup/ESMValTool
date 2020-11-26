.. _reviewing:

Review of pull requests
=======================

In the ESMValTool community we use pull request reviews to ensure all code and
documentation contributions are of good quality.
An introduction to code reviews can be found in `The Turing Way`_, including
`why code reviews are important`_ and advice on
`how to have constructive reviews`_.

Most pull requests will need two reviews before they can be merged.
First a technical review takes place and then a scientific review.
Once both reviewers have approved a pull request, it can be merged.
These three steps are described in more detail below.
If a pull request contains only technical changes, e.g. a pull request that
corrects some spelling errors in the documentation or a pull request that
fixes some installation problem, a scientific review is not needed.

If you are a regular contributor, please try to review a bit more than two
other pull requests for every pull request you create yourself, to make sure
that each pull request gets the attention it deserves.


1. Technical review
-------------------

Technical reviews are done by the technical review team.
This team consists of regular contributors that have a strong interest and
experience in software engineering.

Technical reviewers use the technical checklist from the
`pull request template`_ to make sure the pull request follows the standards we
would like to uphold as a community.
The technical reviewer also keeps an eye on the design and checks that no major
design changes are made without the approval from the technical lead development
team.
If needed, the technical reviewer can help with programming questions, design
questions, and other technical issues.

The technical review team can be contacted by writing
`@ESMValGroup/tech-reviewers`_ in a comment on an issue or pull request on
GitHub.

2. Scientific review
--------------------

Scientific reviews are done by the scientific review team.
This team consists of contributors that have a strong interest and
experience in climate science or related domains.

Scientific reviewers use the scientific checklist from the
`pull request template`_ to make sure the pull request follows the standards we
would like to uphold as a community.

The scientific review team can be contacted by writing
`@ESMValGroup/science-reviewers`_ in a comment on an issue or pull request on
GitHub.

3. Merge
--------

Pull requests are merged by the :ref:`core-team`.
The team member who does the merge first checks that both the technical and
scientific reviewer approved the pull request and that the reviews were
conducted thoroughly.
He or she looks at the list of files that were changed
in the pull request and checks that all relevant checkboxes from the checklist
in the pull request template have been added and ticked.

The core development team can be contacted by writing `@ESMValGroup/esmvaltool-coreteam`_
in a comment on an issue or pull request on GitHub.

Frequently asked questions
--------------------------

How do I request a review of my pull request?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you know a suitable reviewer, e.g. because your pull request fixes an issue
that they opened or they are otherwise interested in the work you are
contributing, you can ask them for a review by clicking the cogwheel next to
'Reviewers' on the pull request 'Conversation' tab and clicking on that person.
When changing code, it is a good idea to ask the original authors of that code
for a review.
An easy way to find out who previously worked on a particular piece of code is
to use `git blame`_.
GitHub will also suggest reviewers based on who previously worked on the files
changed in a pull request.
Every recipe has a maintainer and authors listed in the recipe, it is a good
idea to ask these people for a review.

If there is no obvious reviewer, you can attract the attention of the relevant
team of reviewers by writing to `@ESMValGroup/tech-reviewers`_ or
`@ESMValGroup/science-reviewers`_ in a comment on your pull request.

How do I optimize for a fast review?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When authoring a pull request, please keep in mind that it is easier and
faster to review a pull request that does not contain many changes.
Try to add one new feature per pull request and change only a few files.
For the ESMValTool repository, try to limit changes to a few hundred lines of
code and new diagnostics to not much more than a thousand lines of code.
For the ESMValCore repository, a pull request should ideally change no more
than about a hundred lines of existing code, though adding more lines for unit
tests and documentation is fine.

If you are a regular contributor, make sure you regularly review other people's
pull requests, that way they will be more inclined to return the favor by
reviewing your pull request.

How do I find a pull request to review?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please pick pull requests to review yourself based on your interest or
expertise.
We try to be self organizing, so there is no central authority that will assign
you to review anything.
If someone knows you have expertise on a certain topic, they might request your
review on a pull request though.
If your review is requested, please try to respond within a few days if at all
possible.
If you do not have the time to review the pull request, notify the author and
try to find a replacement reviewer.

How do I actually do a review?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To do a review, go to the pull request on GitHub, the list of all pull requests
is available here https://github.com/ESMValGroup/ESMValCore/pulls for the ESMValCore
and here https://github.com/ESMValGroup/ESMValTool/pulls for the ESMValTool, click the
pull request you would like to review.
The top comment should contain (a selection of) the checklist available in the
`pull request template`_.
If it is not there, copy the relevant items from the `pull request template`_.
Which items from the checklist are relevant, depends on which files are changed
in the pull request. A list of items to check with brief explanations is given in
section :ref:`checklists`. The items are grouped by technical and scientific review.
To see which files have changed, click the tab 'Files changed'.
To comment on specific lines of code or documentation, click the 'plus' icon
next to a line of code and write your comment.
When you are done reviewing, use the 'Review changes' button in the top right
corner to comment on, request changes to, or approve the pull request.

What if the author and reviewer disagree?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the author and the reviewer of a pull request have difficulty agreeing
on what needs to be done before the pull request can be approved, it is usually
both more pleasant and more efficient to schedule a meeting or co-working
session, for example using `Google meet`_ or `Jitsi meet`_.

When reviewing a pull request, try to refrain from making changes to the pull
request yourself, unless the author specifically agrees to those changes, as
this could potentially be perceived as offensive.

If talking about the pull requests in a meeting still does not resolve the
disagreement, ask a member of the `@ESMValGroup/esmvaltool-coreteam`_ for
their opinion and try to find a solution.


.. _checklists:

Checklists for reviewing a pull request
---------------------------------------

Below are general checklists for doing technical and scientific reviews including brief descriptions of the tasks to do. Reviewing
CMORizer scripts consists mostly of technical tasks but differs slightly from the technical review tasks and is therefore listed
in a third table below.

Technical reviews
~~~~~~~~~~~~~~~~~

+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Item                                | Comments                                                                                         |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Documentation added to user’s guide | Check that the scientific documentation of the new diagnostic has been added to the user’s guide |
|                                     |                                                                                                  |
|                                     | * A file ./doc/sphinx/source/recipes/recipe_<diagnostic>.rst exists                              |
|                                     | * New documentation is included in ./doc/sphinx/source/recipes/index.rst                         |
|                                     | * documentation follows template (./doc/sphinx/source/recipes/recipe_template.rst.template)      |
|                                     | * configuration options                                                                          |
|                                     | * variables                                                                                      |
|                                     | * valid image files                                                                              |
|                                     | * resolution of image files (~150 dpi is usually enough; file size should be kept small)         |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| recipe                              | Check yaml syntax and that new recipe contains                                                   |
|                                     |                                                                                                  |
|                                     | * documentation: description, authors, maintainer, references, projects                          |
|                                     | * provenance tags: themes, realms                                                                |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| diagnostic script                   | Check that the new diagnostic script(s) meet(s) standards. This includes the following items:    |
|                                     |                                                                                                  |
|                                     | * In-code documentation                                                                          |
|                                     | * Code quality checks                                                                            |
|                                     |                                                                                                  |
|                                     |   (1) code quality (e.g. no hardcoded pathnames)                                                 |
|                                     |   (2) no Codacy errors reported                                                                  |
|                                     | * Re-use of existing functions whenever possible                                                 |
|                                     | * Provenance implemented                                                                         |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| run recipe                          | Make sure new diagnostic(s) is working by running the ESMValTool                                 |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Check output of diagnostic          | After successfully running the new recipe, check that                                            |
|                                     |                                                                                                  |
|                                     | * Netcdf output has been written                                                                 |
|                                     | * Output contains (some) valid values (e.g. not only nan or zeros)                               |
|                                     | * Provenance information has been written                                                        |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Check autoamted tests               | Check for errors reported by automated tests                                                     |
|                                     |                                                                                                  |
|                                     | * Codacy                                                                                         |
|                                     | * CircleCI                                                                                       |
|                                     | * documentation build                                                                            |
+-------------------------------------+--------------------------------------------------------------------------------------------------+


Scientific reviews
~~~~~~~~~~~~~~~~~~

+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Item                                | Comments                                                                                         |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Documentation added to user’s guide | Check that the scientific documentation of the new diagnostic                                    |
|                                     | ./doc/sphinx/source/recipes/recipe_<diagnostic>.rst                                              |
|                                     |                                                                                                  |
|                                     | * meets scientific documentation standard (brief description of method, references, typos,       |
|                                     |   understandable language)                                                                       |
|                                     | * references are complete                                                                        |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| recipe                              | Check that new recipe contains valid                                                             |
|                                     |                                                                                                  |
|                                     | * documentation: description, references                                                         |
|                                     | * provenance tags: themes, realms                                                                |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| diagnostic script                   | Check that the new diagnostic script(s) meet(s) scientific standards. This can include the       |
|                                     | following items:                                                                                 |
|                                     |                                                                                                  |
|                                     | * Clear and understandable in-code documentation including brief description of diagnostic       |
|                                     | * References                                                                                     |
|                                     | * Method / equations match reference(s) given                                                    |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| run recipe                          | Make sure new diagnostic(s) is working by running the ESMValTool                                 |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Check output of diagnostic          | After successfully running the new recipe, check that                                            |
|                                     |                                                                                                  |
|                                     | * Output contains (some) valid values (e.g. not only nan or zeros)                               |
|                                     | * If applicable, check plots and compare with corresponding plots in the paper(s) cited          |
+-------------------------------------+--------------------------------------------------------------------------------------------------+

CMORizer scripts
~~~~~~~~~~~~~~~~

Reviewing CMORizer scripts differs slightly from reviewing technical changes or scientific reviews of new diagnostics. A review typically
contains mostly technical aspects given in the checklist below.

+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Dataset description added to user’s | Check that new dataset has been added to the table of observations defined in the ESMValTool     |
| guide                               | user’s guide in section “Obtaining input data” (./doc/sphinx/source/input.rst).                  |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| BibTeX info file                    | Check that a BibTeX file (i.e. <dataset>.bibtex) defining the reference(s) for the new dataset   |
|                                     | has been created in ./esmvaltool/references/.                                                    |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| recipe_check_obs.yml                | Check that new dataset has been added to the testing recipe                                      |
|                                     | ./esmvaltool/recipes/examples/recipe_check_obs.yml                                               |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| CMORizer script                     | Check that the new CMORizer script (./esmvaltool/cmorizers/obs/cmorize_obs_<dataset>.py/.ncl/.r) |
|                                     | meets standards. This includes the following items:                                              |
|                                     |                                                                                                  |
|                                     | * In-code documentation (header) contains                                                        |
|                                     |                                                                                                  |
|                                     |   (1) download instructions                                                                      |
|                                     |   (2) reference(s)                                                                               |
|                                     | * Code quality checks                                                                            |
|                                     |                                                                                                  |
|                                     |   (1) code quality (e.g. no hardcoded pathnames)                                                 |
|                                     |   (2) no Codacy errors reported                                                                  |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Config file                         | If present, check config file <dataset>.yml in ./esmvaltool/cmorizers/obs/cmor_config/.          |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Run CMORizer                        | Make sure CMORizer is working by running ''cmorize_obs -c <config-file> -o <dataset>''           |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Check output of CMORizer            | After successfully running the new CMORizer, check that                                          |
|                                     |                                                                                                  |
|                                     | * Output contains (some) valid values (e.g. not only nan or zeros)                               |
|                                     | * Metadata is defined properly                                                                   |
|                                     |                                                                                                  |
|                                     | Run ./esmvaltool/recipes/examples/recipe_check_obs.yml for new dataset                           |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| RAW data                            | Contact person in charge of ESMValTool data pool and request to copy RAW data                    |
|                                     | to RAWOBS/Tier2 (Tier3)                                                                          |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| CMORized data                       | Contact person in charge of ESMValTool data pool and request to                                  |
|                                     |                                                                                                  |
|                                     | * Copy CMORized dataset to OBS/Tier2 (Tier3)                                                     |
|                                     | * Set file access rights for new dataset                                                         |
+-------------------------------------+--------------------------------------------------------------------------------------------------+


After merging a pull request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
After merging a pull request successfully, the :ref:`core-team` will:

*	Close related issue if existent
*	Delete feature branch


.. _`The Turing Way`: https://the-turing-way.netlify.app/reproducible-research/reviewing.html
.. _`why code reviews are important`: https://the-turing-way.netlify.app/reproducible-research/reviewing/reviewing-motivation.html
.. _`how to have constructive reviews`: https://the-turing-way.netlify.app/reproducible-research/reviewing/reviewing-recommend.html
.. _`@ESMValGroup/tech-reviewers`: https://github.com/orgs/ESMValGroup/teams/tech-reviewers
.. _`@ESMValGroup/science-reviewers`: https://github.com/orgs/ESMValGroup/teams/science-reviewers
.. _`@ESMValGroup/esmvaltool-coreteam`: https://github.com/orgs/ESMValGroup/teams/esmvaltool-coreteam
.. _`pull request template`: https://raw.githubusercontent.com/ESMValGroup/ESMValTool/master/.github/pull_request_template.md
.. _`Google meet`: https://meet.google.com
.. _`Jitsi meet`: https://meet.jit.si
.. _`git blame`: https://www.freecodecamp.org/news/git-blame-explained-with-examples/
