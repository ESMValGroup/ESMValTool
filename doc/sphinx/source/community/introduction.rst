.. _contributing:

Contribution guidelines
=======================

Contributions are very welcome
------------------------------

We greatly value contributions of any kind.
Contributions could include, but are not limited to documentation improvements,
bug reports, new or improved diagnostic code, scientific and technical code
reviews, infrastructure improvements, mailing list and chat participation,
community help/building, education and outreach.
We value the time you invest in contributing and strive to make the process as
easy as possible.
If you have suggestions for improving the process of contributing, please do
not hesitate to propose them, for example by starting a discussion on our
`discussions page <https://github.com/ESMValGroup/ESMValTool/discussions>`__.

If you have a bug or other issue to report, please open an issue on the
`issues tab on the ESMValTool github
repository <https://github.com/ESMValGroup/ESMValTool/issues>`__.

If you would like to contribute a new diagnostic and recipe or a new
feature, please discuss your idea with the development team before
getting started, to avoid double work and/or disappointment later. A
good way to do this is to open an `issue on
GitHub <https://github.com/ESMValGroup/ESMValTool/issues>`__. This is
also a good way to get help with the implementation.

Getting started
---------------

See :ref:`install_from_source` for instructions on how to set up a development
installation.
Most contributions are either a :ref:`new recipe and diagnostic <new-diagnostic>`
or a :ref:`new dataset <new-dataset>`, have a look at those chapters for an
introduction on how to contribute those.

Checklist for pull requests
---------------------------

To clearly communicate up front what is expected from a pull request, we have
the following checklist.
Please try to do everything on the list before requesting a review.
If you are unsure about something on the list, please ask the
`@ESMValGroup/tech-reviewers`_ or `@ESMValGroup/science-reviewers`_ for help
by commenting on your (draft) pull request or by starting a new
`discussion <https://github.com/ESMValGroup/ESMValTool/discussions>`__.

The icons indicate whether the item will be checked during the
:ref:`🛠 Technical review <technical_review>` or
:ref:`🧪 Scientific review <scientific_review>`.

- 🛠 :ref:`The pull request has a descriptive title <descriptive_pr_title>`
- 🛠 Code is written according to the :ref:`code quality guidelines <code_quality>`
- 🛠 Documentation_ is available
- 🛠 Tests_ run successfully
- 🛠 New dependencies are :ref:`added or removed correctly <dependencies>`

New or updated recipe and/or diagnostic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- 🛠 :ref:`Provenance information <recording-provenance>` has been added
- 🧪 Documentation_ for the recipe/diagnostic clearly describes what it calculates from a scientific point of view
- 🧪 :ref:`Recipe runs successfully <testing_recipes>`
- 🧪 Figure(s) and data look as expected from literature
- 🧪 Code is :ref:`well documented <doc_howto>` and scientifically sound

New or updated data reformatting script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- 🛠 :ref:`Documentation <dataset-documentation>` is available
- 🛠 The dataset has been :ref:`added to the CMOR check recipe <dataset-test>`
- 🧪 Numbers and units of the data look :ref:`physically meaningful <dataset-sanity-check>`

.. _descriptive_pr_title:

Use a descriptive title for your pull request
---------------------------------------------

The title of a pull request should clearly describe what the pull request changes.
The titles of pull requests are used to compile the :ref:`changelog`, therefore
it is important that they are easy to understand for people who are not
familiar with the code or people in the project.
Using descriptive pull request titles also makes it easier to find back what was
changed when in case bugs were introduced.

.. _code_quality:

Code and documentation quality
------------------------------

To increase the readability and maintainability or the ESMValTool source
code, we aim to adhere to best practices and coding standards.
For code in all languages, it is highly recommended that you split your code up
in functions that are short enough to view without scrolling, e.g. no more than
50 lines long.

We include checks for Python, R, NCL, and yaml files, most of which are
described in more detail in the sections below.
This includes checks for invalid syntax and formatting errors.
:ref:`pre-commit` is a handy tool that can run all of these checks automatically.
It knows knows which tool to run for each filetype, and therefore provides
a simple way to check your code!

Python
~~~~~~

The standard document on best practices for Python code is
`PEP8 <https://www.python.org/dev/peps/pep-0008/>`__ and there is
`PEP257 <https://www.python.org/dev/peps/pep-0257/>`__ for
documentation. We make use of `numpy style
docstrings <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`__
to document Python functions that are visible on
`readthedocs <https://docs.esmvaltool.org>`__.

Most formatting issues in Python code can be fixed automatically by
running the commands

::

   isort some_file.py

to sort the imports in `the standard way <https://www.python.org/dev/peps/pep-0008/#imports>`__
using `isort <https://pycqa.github.io/isort/>`__ and

::

   yapf -i some_file.py

to add/remove whitespace as required by the standard using `yapf <https://github.com/google/yapf>`__,

::

   docformatter -i your_script.py

to run `docformatter <https://github.com/myint/docformatter>`__ which helps formatting the doc strings (such as line length, spaces).

To check if your code adheres to the standard, go to the directory where
the repository is cloned, e.g. ``cd ESMValTool``, and run `prospector <http://prospector.landscape.io/>`__

::

   prospector esmvaltool/diag_scripts/your_diagnostic/your_script.py

Run

::

   python setup.py lint

to see the warnings about the code style of the entire project.

We use `flake8 <https://flake8.pycqa.org/en/latest/>`__ on CircleCI to automatically check that there are
no formatting mistakes and Codacy for monitoring (Python) code quality.
Running prospector locally will give you quicker and sometimes more
accurate results.

NCL
~~~

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
guide <https://style.tidyverse.org/>`__. We check adherence to this
style guide by using
`lintr <https://cran.r-project.org/web/packages/lintr/index.html>`__ on
CircleCI. Please use `styler <https://styler.r-lib.org/>`__ to
automatically format your code according to this style guide. In the
future we would also like to make use of
`goodpractice <https://cran.r-project.org/web/packages/goodpractice/index.html>`__
to assess the quality of R code.

YAML
~~~~

Please use ``yamllint`` to check that your YAML files do not contain
mistakes.

Any text file
~~~~~~~~~~~~~

A generic tool to check for common spelling mistakes is
`codespell <https://pypi.org/project/codespell/>`__.


Documentation
-------------

.. _doc_howto:

What should be documented
~~~~~~~~~~~~~~~~~~~~~~~~~

Any code documentation that is visible on
`docs.esmvaltool.org <https://docs.esmvaltool.org>`__
should be well written and adhere to the standards for documentation for the
respective language.
Note that there is no need to write extensive documentation for functions that
are not visible in the online documentation.
However, a short description in the docstring helps other contributors to
understand what a function is intended to do and and what its capabilities are.
For short functions, a one-line docstring is usually sufficient, but more
complex functions might require slightly more extensive documentation.

Recipes should have a page in the :ref:`recipes` chapter.
This is also the place to document recipe options for the diagnostic scripts
used in those recipes.
When adding a new recipe, please start from the
`template <https://github.com/ESMValGroup/ESMValTool/blob/master/doc/sphinx/source/recipes/recipe_template.rst.template>`_
and do not forget to add your recipe to the
`index <https://github.com/ESMValGroup/ESMValTool/blob/master/doc/sphinx/source/recipes/index.rst>`_.

Functions implementing scientific formula should contain comments with
references to the paper and formula number(s).

How to view the documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

   python setup.py build_sphinx

or

::

   python setup.py build_sphinx -Ea

to build it from scratch.
Make sure that your newly added documentation builds without warnings or
errors and looks correctly formatted.
`CircleCI <https://app.circleci.com/pipelines/github/ESMValGroup/ESMValTool>`_
will build the documentation with the command

::

   python setup.py build_sphinx --warning-is-error

to catch mistakes that can be detected automatically.

When reviewing a pull request, always check that the documentation checks were
successful.
Successful checks have a green ✓ in front, a ❌ means the test job failed.

.. _tests:

Tests
-----

To check various aspects of the recipes and code, there tests available in the
`tests <https://github.com/ESMValGroup/ESMValTool/tree/master/tests>`__
directory.

Whenever you make a pull request or push new commits to an existing pull
request, these tests will be run automatically on CircleCI_.
The results appear at the bottom of the pull request.
Click on 'Details' for more information on a specific test job.

To run the tests on your own computer, go to the directory where the repository
is cloned and run the command ``pytest``.

When reviewing a pull request, always check that all test jobs on CircleCI_ were
successful.
Successful test jobs have a green ✓ in front, a ❌ means the test job failed.

.. _testing_recipes:

Testing recipes
~~~~~~~~~~~~~~~

To test a recipe, you can run it yourself on your local infrastructure or you
can ask the `@esmvalbot <https://github.com/apps/esmvalbot>`_ to run it for you.
To request a run of ``recipe_xyz.yml``, write the following comment below a pull
request:

::

   @esmvalbot Please run recipe_xyz.yml

Note that only members of the `@ESMValGroup/esmvaltool-developmentteam`_
can request runs. The memory of the `@esmvalbot`_ is limited to 16 GB and it only
has access to data available at DKRZ.

When reviewing a pull request, at the very least check that a recipes runs
without any modifications.
For a more thorough check, you might want to try out different datasets or
changing some settings if the diagnostic scripts support those.
A simple :ref:`tool <recipe_test_tool>` is available for testing recipes
with various settings.


.. _dependencies:

Adding or removing dependencies
-------------------------------

Before considering adding a new dependency, carefully check that the license of
the dependency you want to add and any of its dependencies are compatible with
the
`Apache 2.0 <https://github.com/ESMValGroup/ESMValTool/blob/master/LICENSE/>`_
license that applies to the ESMValTool.
Note that GPL version 2 license is considered incompatible with the Apache 2.0
license, while the compatibility of GPL version 3 license with the Apache 2.0
license is questionable.
See this `statement <https://www.apache.org/licenses/GPL-compatibility.html>`__
by the authors of the Apache 2.0 license for more information.

The following files contain lists of dependencies

- ``environment.yml``
  contains development dependencies that cannot be installed from
  PyPI/CRAN/Julia package repository
- ``docs/sphinx/source/requirements.txt``
  contains Python dependencies needed to build the documentation that can be
  installed from PyPI
- ``docs/sphinx/source/conf.py``
  contains a list of Python dependencies needed to build the documentation that
  cannot be installed from PyPI and need to be mocked when building the
  documentation
- ``esmvaltool/install/R/r_requirements.txt``
  contains R dependencies that can be installed from CRAN
- ``esmvaltool/install/Julia/Project.toml``
  contains Julia dependencies that can be installed from the Julia package
  repository
- ``setup.py``
  contains all Python dependencies, regardless of their installation source
- ``package/meta.yaml``
  contains dependencies for the conda package, all Python and compiled
  dependencies that can be installed from conda should be listed here, but no R
  or Julia dependencies, because this would make it impossible to solve the
  conda environment

Note that packages may have a different name on conda than on PyPI or CRAN.

Several test jobs on CircleCI_ related to the installation of the tool will only
run if you change the dependencies, these will be skipped for most pull
requests.

When reviewing a pull request where dependencies are added or removed, always
check that the changes have been applied in all relevant places.

.. _pre-commit:

Pre-commit
----------

To run ``pre-commit`` on your code, go to the ESMValTool directory
(``cd ESMValTool``) and run

::

   pre-commit run

By default, pre-commit will only run on the files that have been changed,
meaning those that have been staged in git (i.e. after
``git add your_script.py``).

To make it only check some specific files, use

::

   pre-commit run --files your_script.py

or

::

   pre-commit run --files your_script.R

Alternatively, you can configure ``pre-commit`` to run on the staged files before
every commit (i.e. ``git commit``), by installing it as a `git hook <https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks>`__ using

::

   pre-commit install

Pre-commit hooks are used to inspect the code that is about to be committed. The
commit will be aborted if files are changed or if any issues are found that
cannot be fixed automatically. Some issues cannot be fixed (easily), so to
bypass the check, run

::

   git commit --no-verify

or

::

   git commit -n

or uninstall the pre-commit hook

::

   pre-commit uninstall


Branches, pull requests and code review
---------------------------------------

New development should preferably be done in the
`ESMValTool <https://github.com/ESMValGroup/ESMValTool>`__
GitHub repository.
However, for scientists requiring confidentiality, private repositories are
available, see :ref:`private_repository` for more information.
The default git branch is ``master``. Use
this branch to create a new feature branch from and make a pull request
against.
This
`page <https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow>`__
offers a good introduction to git branches, but it was written for
BitBucket while we use GitHub, so replace the word BitBucket by GitHub
whenever you read it.

It is recommended that you open a `draft pull
request <https://github.blog/2019-02-14-introducing-draft-pull-requests/>`__
early, as this will cause CircleCI to run the unit tests, Codacy to
analyse your code, and readthedocs to build the documentation.
It’s also easier to get help from other developers if
your code is visible in a pull request.

You can view the results of the automatic checks below your pull
request by clicking on ``Details``. If one of the tests shows a red cross instead of a green
approval sign, please click the link and try to solve the issue. Note
that this kind of automated checks make it easier to review code, but
they are not flawless. Preferably Codacy code quality checks pass, however
a few remaining hard to solve Codacy issues are still acceptable.
If you suspect Codacy may be wrong, please ask by commenting.

The documentation can be seen by clicking on `Details`. Make sure the
documentation is nicely formatted, and (if necessary) add the link to the
top of the pull request.

Make sure your pull request has a descriptive title that can be used in the
`changelog <https://docs.esmvaltool.org/en/latest/changelog.html>`__.

Diagnostic script contributions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A pull request with diagnostic code should preferably not introduce new
Codacy issues. However, we understand that there is a limit to how much
time can be spend on polishing code, so up to 10 new (non-trivial)
issues is still an acceptable amount.

List of authors
~~~~~~~~~~~~~~~

If you make a (significant) contribution to ESMValTool, please add your
name to the list of authors in CITATION.cff and regenerate the file
.zenodo.json by running the command

::

   pip install cffconvert
   cffconvert --ignore-suspect-keys --outputformat zenodo --outfile .zenodo.json

.. _`@ESMValGroup/esmvaltool-developmentteam`: https://github.com/orgs/ESMValGroup/teams/esmvaltool-developmentteam
.. _`@ESMValGroup/tech-reviewers`: https://github.com/orgs/ESMValGroup/teams/tech-reviewers
.. _`@ESMValGroup/science-reviewers`: https://github.com/orgs/ESMValGroup/teams/science-reviewers
