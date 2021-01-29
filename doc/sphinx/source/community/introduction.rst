.. _contributing:

Contribution guidelines
=======================

Contributions are very welcome
------------------------------

We greatly value contributions of any kind.
Contributions could include, but are not limited to documentation improvements, bug reports, new or improved diagnostic code, scientific and technical code reviews, infrastructure improvements, mailing list and chat participation, community help/building, education and outreach.
We value the time you invest in contributing and strive to make the process as easy as possible.
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

Running tests
-------------

Go to the directory where the repository is cloned and run
``pytest``. Tests will also be run automatically by
`CircleCI <https://circleci.com/gh/ESMValGroup/ESMValTool>`__.

These automated checks are run automatically when you add new commits to your pull request.
They appear at the bottom of the pull request. Click on `Details` for more information


Code style
----------

To increase the readability and maintainability or the ESMValTool source
code, we aim to adhere to best practices and coding standards. All pull
requests are reviewed and tested by one or more members of the core
development team. For code in all languages, it is highly recommended
that you split your code up in functions that are short enough to view
without scrolling.

We include checks for Python, R, NCL, and yaml files, most of which are
described in more detail in the sections below.
This includes checks for invalid syntax and formatting errors.
`Pre-commit <https://pre-commit.com/>`__ is a handy tool that can run
all of these checks automatically.
It knows knows which tool to run for each filetype, and therefore provides
a simple way to check your code!

Pre-commit
~~~~~~~~~~

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
indentation should be 2 spaces instead of 4. Use the command
``nclcodestyle /path/to/file.ncl`` to check if your code follows the
style guide.

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

What should be documented
~~~~~~~~~~~~~~~~~~~~~~~~~

Any code documentation that is visible on
`docs.esmvaltool.org <https://docs.esmvaltool.org>`__
should be well written and adhere to the standards for documentation for the
respective language.
Recipes should have a page in the :ref:`recipes` section.
This is also the place to document recipe options for the diagnostic scripts
used in those recipes.
When adding a new recipe, please start from the
`template <https://github.com/ESMValGroup/ESMValTool/blob/master/doc/sphinx/source/recipes/recipe_template.rst.template>`_
and do not forget to add your recipe to the
`index <https://github.com/ESMValGroup/ESMValTool/blob/master/doc/sphinx/source/recipes/index.rst>`_.
Note that there is no need to write extensive documentation for functions that
are not visible in the online documentation.
However, a short description in the docstring helps other contributors to
understand what a function is intended to do and and what its capabilities are.
For short functions, a one-line docstring is usually sufficient, but more
complex functions might require slightly more extensive documentation.

How to build the documentation locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Go to the directory where the repository is cloned and run

::

   python setup.py build_sphinx -Ea

Make sure that your newly added documentation builds without warnings or
errors.

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
