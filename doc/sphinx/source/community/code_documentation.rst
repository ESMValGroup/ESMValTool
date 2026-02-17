.. _contributing_code_docs:

Contributing code and documentation
===================================

If you would like to contribute a new diagnostic and recipe or a new feature,
please discuss your idea with the development team before getting started, to
avoid double work and/or disappointment later.
A good way to do this is to open an
`issue on GitHub <https://github.com/ESMValGroup/ESMValTool/issues>`__.
This is also a good way to get help with the implementation.

We value the time you invest in contributing and strive to make the process as
easy as possible.
If you have suggestions for improving the process of contributing, please do
not hesitate to propose them, for example by starting a discussion on our
`discussions page <https://github.com/ESMValGroup/ESMValTool/discussions>`__.

Getting started
---------------

See :ref:`install_from_source` for instructions on how to set up a development
installation.

New development should preferably be done in the
`ESMValTool <https://github.com/ESMValGroup/ESMValTool>`__
GitHub repository.
However, for scientists requiring confidentiality, private repositories are
available, see :ref:`private_repository` for more information.
The default git branch is ``main``. Use
this branch to create a new feature branch from and make a pull request
against.
This
`page <https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow>`__
offers a good introduction to git branches, but it was written for
BitBucket while we use GitHub, so replace the word BitBucket by GitHub
whenever you read it.

It is recommended that you open a `draft pull
request <https://github.blog/2019-02-14-introducing-draft-pull-requests/>`__
early, as this will cause :ref:`CircleCI to run the unit tests <tests>`,
:ref:`Codacy to analyse your code <code_quality>`, and
:ref:`readthedocs to build the documentation <documentation>`.
It’s also easier to get help from other developers if
your code is visible in a pull request.

Please review the results of the automatic checks below your pull request.
If one of the tests shows a red cross instead of a green checkmark, please click
the ``Details`` link behind the failing check and try to solve the issue.
Ask `@ESMValGroup/tech-reviewers`_ for help if you do not know how to fix the
failing check.
Note that this kind of automated checks make it easier to
:ref:`review code <reviewing>`, but they are not flawless.
Preferably Codacy code quality checks pass, however a few remaining hard to
solve Codacy issues are still acceptable.
If you suspect Codacy may be wrong, please ask by commenting on your pull
request.

.. _pull_request_checklist:

Checklist for pull requests
---------------------------

To clearly communicate up front what is expected from a pull request, we have
the following checklist.
Please try to do everything on the list before requesting a review.
If you are unsure about something on the list, please ask the
`@ESMValGroup/tech-reviewers`_ or `@ESMValGroup/science-reviewers`_ for help
by commenting on your (draft) pull request or by starting a new
`discussion <https://github.com/ESMValGroup/ESMValTool/discussions>`__.

In the ESMValTool community we use
:ref:`pull request reviews <reviewing>` to ensure all code and
documentation contributions are of good quality.
The icons indicate whether the item will be checked during the
:ref:`🛠 Technical review <technical_review>` or
:ref:`🧪 Scientific review <scientific_review>`.

All pull requests
~~~~~~~~~~~~~~~~~

- 🛠 :ref:`The pull request has a descriptive title <descriptive_pr_title>`
- 🛠 Code is written according to the :ref:`code quality guidelines <code_quality>`
- 🛠 Documentation_ is available
- 🛠 Tests_ run successfully
- 🛠 The :ref:`list of authors <authors>` is up to date
- 🛠 Changed dependencies are :ref:`added or removed correctly <dependencies>`
- 🛠 The :ref:`checks shown below the pull request <pull_request_checks>` are successful

If a pull request introduces a change that causes a recipe to
no longer run successfully (*breaking change*), or which results in scientifically
significant changes in results (*science change*), additional requirements
defined in the :ref:`backward compatibility policy<backward-compatibility-policy>` apply.
These include in particular:

- 🛠 Instructions for the release notes to assist *recipe
  developers* to adapt their recipe in light of the *backward-incompatible change*
  available.
- 🛠 If applicable, instructions for *recipe developers* working on *user
  recipes* to enable them to adapt their code related to
  *backward-incompatible changes* available (see `ESMValTool_Tutorial: issue
  #263 <https://github.com/ESMValGroup/ESMValTool_Tutorial/issues/263>`__)
  available.
- 🛠 Core development team tagged to notify them of the
  *backward-incompatible change*, and give at least
  2 weeks for objections to be raised before merging to the main
  branch. If a strong objection is raised the backward-incompatible
  change should not be merged until the objection is resolved.
- 🛠 Information required for the “*backward-incompatible changes*”
  section in the PR  that introduces the *backward-incompatible change*
  available.

New or updated recipe and/or diagnostic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See :ref:`new-diagnostic` for detailed instructions.

- 🧪 :ref:`Recipe runs successfully <testing_recipes>`
- 🧪 :ref:`recipe_documentation` is available
- 🧪 :ref:`Figure(s) and data <diagnostic_output>` look as expected from literature
- 🛠 :ref:`Provenance information <recording-provenance>` has been added

New or updated data reformatting script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See :ref:`new dataset <new-dataset>` for detailed instructions.

- 🛠 :ref:`dataset-documentation` is available
- 🛠 The dataset has been :ref:`added to the CMOR check recipe <dataset-test>`
- 🧪 Numbers and units of the data look :ref:`physically meaningful <dataset-sanity-check>`

.. _descriptive_pr_title:

Pull request title
------------------

The title of a pull request should clearly describe what the pull request changes.
If you need more text to describe what the pull request does, please add it in
the description.
The titles of pull requests are used to compile the :ref:`changelog`, therefore
it is important that they are easy to understand for people who are not
familiar with the code or people in the project.
Descriptive pull request titles also makes it easier to find back what was
changed when, which is useful in case a bug was introduced.

.. _code_quality:

Code quality
------------

To increase the readability and maintainability or the ESMValTool source
code, we aim to adhere to best practices and coding standards.
For code in all languages, it is highly recommended that you split your code up
in functions that are short enough to view without scrolling, e.g. no more than
50 lines long.

We include checks for Python, R, NCL, and yaml files, most of which are
described in more detail in the sections below.
This includes checks for invalid syntax and formatting errors.
:ref:`pre-commit` is a handy tool that can run all of these checks automatically
just before you commit your code.
It knows which tool to run for each filetype, and therefore provides
a convenient way to check your code.

To enable pre-commit checks (these will be run every time you commit something),
go to the ``ESMValTool`` directory and run:

.. code-block:: bash

    pre-commit install

just once, after you have cloned the ESMValTool repository.

If you have installed the pre-commit checks, but for some reason need to skip
them (not recommended), use the ``-n`` flag:

.. code-block:: bash

    git commit -n

To run the automatic code formatting and basic code quality checks, use the command:

.. code-block:: bash

    pre-commit run -a

Python
~~~~~~

The standard document on best practices for Python code is
`PEP8 <https://www.python.org/dev/peps/pep-0008/>`__ and there is
`PEP257 <https://www.python.org/dev/peps/pep-0257/>`__ for code documentation.
We make use of
`numpy style docstrings <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`__
to document Python functions that are visible on
`readthedocs <https://docs.esmvaltool.org>`__.

To automatically format and check if code adheres to the standard, we use
`ruff <https://docs.astral.sh/ruff/>`_. Automatic formatting and basic code
quality checks are enforced through pre-commit. In addition to these basic
quality checks, we use ruff to make recommendations on how your code could be
even better.

To automatically format your code and perform basic quality checks, run:

.. code-block:: bash

    pre-commit run -a

To perform the additional code quality checks and automatically fix many of
the reported issues, run:

.. code-block:: bash

   ruff check --fix esmvaltool/diag_scripts/your_diagnostic/your_script.py

Each code quality check performed by ruff is called a ``rule`` and has a unique
code, e.g. ``RET504``. For a detailed explanation of the purpose of these rules,
you can visit the
`rules page on the ruff website <https://docs.astral.sh/ruff/rules>`__ or
run the command:

.. code-block:: bash

    ruff rule RET504

where ``RET504`` is an arbitrary example of a rule code.

When you make a pull request, adherence to the Python development best practices
is checked in two ways:

#. As part of the unit tests, pre-commit is run by
   `CircleCI <https://app.circleci.com/pipelines/github/ESMValGroup/ESMValTool>`_,
   see the section on Tests_ for more information.
#. `Codacy <https://app.codacy.com/gh/ESMValGroup/ESMValTool/pullRequests>`_
   is a service that runs ruff (and other code quality tools) on changed
   files and reports the results.
   Click on the 'Codacy Static Code Analysis' check entry below your pull request
   and then click the 'Resolve' button or the
   'View more details on Codacy Production' link to see the results of the static
   code analysis done by Codacy_.
   If you need to log in, you can do so using your GitHub account.
   If Codacy says that there are issues, but you cannot see them on the Codacy
   webpage for some reason, please run the command ``ruff check`` on your own
   computer to see what the issues are.

A pull request should preferably not introduce new code quality issues.
However, we understand that there is a limit to how much time can be spent on
polishing code, so up to 10 new issues reported by Codacy/ruff is still an acceptable
amount.
However, a pull request where the CircleCI/pre-commit tests are failing cannot be merged.

.. note::

    The automatic code quality checks by ruff are really helpful to
    improve the quality of your code, but they are not flawless.
    If you suspect ruff or Codacy may be wrong or you just need help,
    please ask the `@ESMValGroup/tech-reviewers`_ by commenting on your pull request.

NCL
~~~

.. warning::

  Writing new NCL code is not recommended because the
  `NCL interpreter <https://github.com/NCAR/ncl>`__ is no longer maintained.

Because there is no standard best practices document for NCL, we use
`PEP8 <https://www.python.org/dev/peps/pep-0008/>`__ for NCL code as
well, with some minor adjustments to accommodate for differences in the
languages. The most important difference is that for NCL code the
indentation should be 2 spaces instead of 4.
Use the command ``nclcodestyle /path/to/file.ncl`` to check if your code
follows the style guide.
More information on the ``nclcodestyle`` command can be found
:ref:`here <nclcodestyle>`.

R
~

Best practices for R code are described in `The tidyverse style
guide <https://style.tidyverse.org/>`__. We currently don't check adherence
to this, but we strongly encourage any new R diagnostic to be checked against the
style guide by using `lintr <https://cran.r-project.org/web/packages/lintr/index.html>`__.
Please use `styler <https://styler.r-lib.org/>`__ to automatically format your code according
to this style guide. To run these tools, uncomment these lines:

.. literalinclude:: ../../../../.pre-commit-config.yaml
    :language: yaml
    :start-at: # - repo: https://github.com/lorenzwalthert/precommit/  # Checks for R
    :end-at: - id: lintr
    :caption: R configuration in ``.pre-commit-config.yaml``.

and run the commands:

.. code-block:: bash

    pre-commit run style-files --files esmvaltool/diag_scripts/your_script.R
    pre-commit run lintr --files esmvaltool/diag_scripts/your_script.R

These pre-commit hooks have been disabled by default because over the past few
years all diagnostics that have been contributed were written in Python.
However, if R contributions pick up again, we would be happy to enable these
pre-commit hooks again.
In the future we would also like to make use of
`goodpractice <https://cran.r-project.org/web/packages/goodpractice/index.html>`__
to assess the quality of R code.

YAML
~~~~

Please use `yamllint <https://yamllint.readthedocs.io>`_ to check that your
YAML files do not contain mistakes.
``yamllint`` checks for valid syntax, common mistakes like key repetition and
cosmetic problems such as line length, trailing spaces, wrong indentation, etc.
When the tool complains about the maximum line length or too many spaces, please
use your own best judgement about whether solving the issue will make your
recipe more readable.

Any text file
~~~~~~~~~~~~~

A generic tool to check for common spelling mistakes is
`codespell <https://pypi.org/project/codespell/>`__.

.. _documentation:

Documentation
-------------

The documentation lives on `docs.esmvaltool.org <https://docs.esmvaltool.org>`_
and is built using `Sphinx <https://www.sphinx-doc.org>`_.
There are two main ways of adding documentation:

#. As written text in the directory
   `doc/sphinx/source <https://github.com/ESMValGroup/ESMValTool/tree/main/doc/sphinx/source>`__.
   When writing
   `reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_
   (``.rst``) files, please try to limit the line length to 80 characters and
   always start a sentence on a new line.
   This makes it easier to review changes to documentation on GitHub.

#. As docstrings or comments in code.
   For Python code, the
   `docstrings <https://www.python.org/dev/peps/pep-0257/>`__
   of Python modules, classes, and functions
   that are mentioned in
   `doc/sphinx/source/api <https://github.com/ESMValGroup/ESMValTool/tree/main/doc/sphinx/source/api>`__
   are used to generate documentation.
   This results in the :ref:`api`.

.. _doc_howto:

What should be documented
~~~~~~~~~~~~~~~~~~~~~~~~~

See also :ref:`recipe_documentation` and :ref:`dataset-documentation`.

Any code documentation that is visible on `docs.esmvaltool.org`_
should be well written and adhere to the standards for documentation for the
respective language.
Note that there is no need to write extensive documentation for functions that
are not visible in the online documentation.
However, a short description in the docstring helps other contributors to
understand what a function is intended to do and and what its capabilities are.
For short functions, a one-line docstring is usually sufficient, but more
complex functions might require slightly more extensive documentation.

How to build and view the documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Whenever you make a pull request or push new commits to an existing pull
request, readthedocs will automatically build the documentation.
The link to the documentation will be shown in the list of checks below your
pull request, click 'Details' behind the check
``docs/readthedocs.org:esmvaltool`` to preview the documentation.
If all checks were successful, you may need to click 'Show all checks' to see
the individual checks.

To build the documentation on your own computer, go to the directory where the
repository was cloned and run

::

   sphinx-build doc/sphinx/source/ doc/sphinx/build/

or

::

   sphinx-build -Ea doc/sphinx/source/ doc/sphinx/build/

to build it from scratch.
Make sure that your newly added documentation builds without warnings or
errors and looks correctly formatted.
CircleCI_ will build the documentation with the command

.. code-block:: bash

   sphinx-build -W doc/sphinx/source/ doc/sphinx/build/

to catch mistakes that can be detected automatically.

The configuration file for Sphinx_ is
`doc/sphinx/source/conf.py <https://github.com/ESMValGroup/ESMValTool/blob/main/doc/sphinx/source/conf.py>`_
and the configuration file for ReadTheDocs is
`.readthedocs.yaml <https://github.com/ESMValGroup/ESMValTool/blob/main/.readthedocs.yaml>`_.

When reviewing a pull request, always check that the documentation checks
shown below the pull request were successful.
Successful checks have a green ✓ in front, a ❌ means the test job failed.

.. _esmvalcore-documentation-integration:

Integration with the ESMValCore documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The `ESMValCore documentation <https://docs.esmvaltool.org/projects/ESMValCore/en/latest/index.html>`_
is hosted as a
`subproject <https://docs.readthedocs.io/en/stable/subprojects.html>`_
of the ESMValTool documentation on readthedocs.
To link to a section from the ESMValCore documentation from the reStructuredText
(``.rst``) files, use the usual ``:ref:`` but prefix the reference with
``esmvalcore:``.
For example, ``:ref:`esmvalcore:recipe``` to link to
:ref:`esmvalcore:recipe`.

.. _tests:

Tests
-----

To check various aspects of the recipes and code, there tests available in the
`tests <https://github.com/ESMValGroup/ESMValTool/tree/main/tests>`__
directory.

Whenever you make a pull request or push new commits to an existing pull
request, these tests will be run automatically on CircleCI_.
The results appear at the bottom of the pull request.
Click on 'Details' for more information on a specific test job.
To see some of the results on CircleCI, you may need to log in.
You can do so using your GitHub account.

To run the tests on your own computer, go to the directory where the repository
is cloned and run the command ``pytest``.

Have a look at :ref:`testing_recipes` for information on testing recipes.

Every night, more extensive tests are run to make sure that problems with the
installation of the tool are discovered by the development team before users
encounter them.
These nightly tests have been designed to mimic the installation procedures
described in the documentation, e.g. in the :ref:`install` chapter.
The nightly tests are run using both CircleCI and GitHub Actions, the
result of the tests ran by CircleCI can be seen on the
`CircleCI project page <https://app.circleci.com/pipelines/github/ESMValGroup/ESMValTool?branch=main>`__
and the result of the tests ran by GitHub Actions can be viewed on the
`Actions tab <https://github.com/ESMValGroup/ESMValTool/actions>`__
of the repository.

The configuration of the tests run by CircleCI can be found in the directory
`.circleci <https://github.com/ESMValGroup/ESMValTool/blob/main/.circleci>`__,
while the configuration of the tests run by GitHub Actions can be found in the
directory
`.github/workflows <https://github.com/ESMValGroup/ESMValTool/blob/main/.github/workflows>`__.

When reviewing a pull request, always check that all test jobs on CircleCI_ were
successful.
Successful test jobs have a green ✓ in front, a ❌ means the test job failed.

.. _authors:

List of authors
---------------

If you make a contribution to ESMValTool and you would like to be listed as an
author (e.g. on `Zenodo <https://zenodo.org/record/4562215>`__), please add your
name to the list of authors in ``CITATION.cff`` and generate the entry for the
``.zenodo.json`` file by running the commands

::

   pip install cffconvert
   cffconvert --infile CITATION.cff --format zenodo --outfile .zenodo.json

Presently, this method unfortunately discards entries `communities`
and `grants` from that file; please restore them manually.

Note that authors of recipes and/or diagnostics also need to be added to the file
`esmvaltool/config-references.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/config-references.yml>`__,
see :ref:`recording-provenance` for more information.

.. _dependencies:

Dependencies
------------

We follow the guidelines set out in `SPEC-0 <https://scientific-python.org/specs/spec-0000/>`_
regarding the supported Python versions and dependencies.

Before considering adding a new dependency, carefully check that the
`license <https://book.the-turing-way.org/reproducible-research/licensing/>`__
of the dependency you want to add and any of its dependencies are
`compatible <https://book.the-turing-way.org/reproducible-research/licensing/licensing-compatibility/>`__
with the
`Apache 2.0 <https://github.com/ESMValGroup/ESMValTool/blob/main/LICENSE/>`_
license that applies to the ESMValTool.
Note that GPL version 2 license is considered incompatible with the Apache 2.0
license, while the compatibility of GPL version 3 license with the Apache 2.0
license is questionable.
See this `statement <https://www.apache.org/licenses/GPL-compatibility.html>`__
by the authors of the Apache 2.0 license for more information.

When adding or removing dependencies, please consider applying the changes in
the following files:

- ``environment.yml``
  contains dependencies that cannot be installed from
  `PyPI <https://pypi.org/>`__
- ``environment_osx.yml``
  contains development dependencies for MacOSX. Should be the same as ``environment.yml``,
  but currently without multi language support.
- ``pyproject.toml``
  contains all Python dependencies, regardless of their installation source

Note that packages may have a different name on
`conda-forge <https://conda-forge.org/>`__ than on PyPI or CRAN.

Several test jobs on CircleCI_ related to the installation of the tool will only
run if you change the dependencies.
These will be skipped for most pull requests.

When reviewing a pull request where dependencies are added or removed, always
check that the changes have been applied in all relevant files.

.. _pull_request_checks:

Pull request checks
-------------------

To check that a pull request is up to standard, several automatic checks are
run when you make a pull request.
Read more about it in the Tests_ and Documentation_ sections.
Successful checks have a green ✓ in front, a ❌ means the check failed.

If you need help with the checks, please ask the technical reviewer of your pull
request for help.
Ask `@ESMValGroup/tech-reviewers`_ if you do not have a technical reviewer yet.

If the checks are broken because of something unrelated to the current
pull request, please check if there is an open issue that reports the problem
and create one if there is no issue yet.
You can attract the attention of the `@ESMValGroup/esmvaltool-coreteam`_ by
mentioning them in the issue if it looks like no-one is working on solving the
problem yet.
The issue needs to be fixed in a separate pull request first.
After that has been merged into the ``main`` branch and all checks are green
again on the ``main`` branch, merge it into your own branch to get the tests
to pass.

When reviewing a pull request, always make sure that all checks were successful.
If the Codacy check keeps failing, please run prospector locally.
If necessary, ask the pull request author to do the same and to address the
reported issues.
See the section on code_quality_ for more information.
Never merge a pull request with failing CircleCI or readthedocs checks.

.. _`@ESMValGroup/esmvaltool-coreteam`: https://github.com/orgs/ESMValGroup/teams/esmvaltool-coreteam
.. _`@ESMValGroup/esmvaltool-developmentteam`: https://github.com/orgs/ESMValGroup/teams/esmvaltool-developmentteam
.. _`@ESMValGroup/tech-reviewers`: https://github.com/orgs/ESMValGroup/teams/tech-reviewers
.. _`@ESMValGroup/science-reviewers`: https://github.com/orgs/ESMValGroup/teams/science-reviewers
