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

.. _technical_review:

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

.. _scientific_review:

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

Pull requests are merged by the `@ESMValGroup/esmvaltool-coreteam`_.
Specifically, pull requests containing a :ref:`CMORizer script<new-dataset>` can only be merged by
`@remi-kazeroni`_, who will then add the CMORized data to the OBS data pool at
DKRZ and CEDA-Jasmin.
The team member who does the merge first checks that both the technical and
scientific reviewer approved the pull request and that the reviews were
conducted thoroughly.
He or she looks at the list of files that were changed in the pull request and
checks that all relevant checkboxes from the checklist in the pull request
template have been added and ticked.
Finally, he or she checks that the :ref:`pull_request_checks` passed and
merges the pull request.
The person doing the merge commit edits the merge commit message so it
contains a concise and meaningful text.

Any issues that were solved by the pull request can be closed after merging.
It is always a good idea to check with the author of an issue and ask if it is
completely solved by the related pull request before closing the issue.

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
You can also label your pull request with one of the labels
`looking for technical reviewer <https://github.com/ESMValGroup/ESMValTool/labels/looking%20for%20technical%20reviewer>`_
or
`looking for scientific reviewer <https://github.com/ESMValGroup/ESMValTool/labels/looking%20for%20scientific%20reviewer>`_,
though asking people for a review directly is probably more effective.

.. _easy_review:

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
People may advertise that they are looking for a reviewer by applying the label
`looking for technical reviewer <https://github.com/ESMValGroup/ESMValTool/labels/looking%20for%20technical%20reviewer>`_
or `looking for scientific reviewer <https://github.com/ESMValGroup/ESMValTool/labels/looking%20for%20scientific%20reviewer>`_.
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
in the pull request.
To see which files have changed, click the tab 'Files changed'.
Please make sure you are familiar with all items from the checklist by reading
the content linked from :ref:`pull_request_checklist` and check all items
that are relevant.
Checklists with some of the items to check are available:
:ref:`recipe and diagnostic checklist <diagnostic_checklist>` and
:ref:`dataset checklist <dataset_checklist>`.

In addition to the items from the checklist, good questions to start a review
with are 'Do I understand why these changes improve the tool?' (if not, ask the
author to improve the documentation contained in the pull request and/or the
description of the pull request on GitHub) and 'What could possibly go wrong if
I run this code?'.

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


.. _`The Turing Way`: https://the-turing-way.netlify.app/reproducible-research/reviewing.html
.. _`why code reviews are important`: https://the-turing-way.netlify.app/reproducible-research/reviewing/reviewing-motivation.html
.. _`how to have constructive reviews`: https://the-turing-way.netlify.app/reproducible-research/reviewing/reviewing-recommend.html
.. _`@ESMValGroup/tech-reviewers`: https://github.com/orgs/ESMValGroup/teams/tech-reviewers
.. _`@ESMValGroup/science-reviewers`: https://github.com/orgs/ESMValGroup/teams/science-reviewers
.. _`@ESMValGroup/esmvaltool-coreteam`: https://github.com/orgs/ESMValGroup/teams/esmvaltool-coreteam
.. _`@remi-kazeroni`: https://github.com/remi-kazeroni
.. _`pull request template`: https://raw.githubusercontent.com/ESMValGroup/ESMValTool/master/.github/pull_request_template.md
.. _`Google meet`: https://meet.google.com
.. _`Jitsi meet`: https://meet.jit.si
.. _`git blame`: https://www.freecodecamp.org/news/git-blame-explained-with-examples/
