# Contributions are very welcome

If you would like to contribute a new diagnostic and recipe or a new feature, please discuss your idea with the development team before getting started, to avoid double work and/or disappointment later. A good way to do this is to open an [issue on GitHub](https://github.com/ESMValGroup/ESMValTool/issues). This is also a good way to get help.

If you have a bug to report, please do so using the [issues tab on the ESMValTool github repository](https://github.com/ESMValGroup/ESMValTool/issues).

To get started developing, follow the instructions below. More detailed instructions can be found in the [manual](https://esmvaltool.readthedocs.io) under Developer's Guide.

## Getting started
To install in development mode, follow these instructions.
  - [Download and install conda](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html) (this should be done even if the system in use already has a preinstalled version of conda, as problems have been reported with NCL when using such a version)
  - To make the `conda` command availble, add `source <prefix>/etc/profile.d/conda.sh` to your `.bashrc` file and restart your shell. If using (t)csh shell, add `source <prefix>/etc/profile.d/conda.csh` to your `.cshrc`/`.tcshrc` file instead.
  - Update conda: `conda update -y conda`
  - Clone the ESMValTool public github repository: `git clone git@github.com:ESMValGroup/ESMValTool.git`, or one of the private github repositories (e.g. `git clone git@github.com:ESMValGroup/ESMValTool-private.git`)
  - Go to the esmvaltool directory: `cd ESMValTool`
  - Create the esmvaltool conda environment `conda env create --name esmvaltool --file environment.yml`
  - Activate the esmvaltool environment: `conda activate esmvaltool`
  - Install in development mode: `pip install -e '.[develop]'`. If you are installing behind a proxy that does not trust the usual pip-urls you can declare them with the option `--trusted-host`, e.g. `pip install --trusted-host=pypi.python.org --trusted-host=pypi.org --trusted-host=files.pythonhosted.org -e .[develop]`
  - If you want to use R diagnostics, run `Rscript esmvaltool/install/R/setup.R` to install the R dependences.
  - If you want to use Julia diagnostics, run `julia esmvaltool/install/Julia/setup.jl` to install the Julia dependences.
  - Test that your installation was succesful by running `esmvaltool -h`.
  - If you log into a cluster or other device via `ssh` and your origin machine sends the `locale` environment via the `ssh` connection, make sure the environment is set correctly, specifically `LANG` and `LC_ALL` are set correctly (for GB English UTF-8 encoding these variables must be set to `en_GB.UTF-8`; you can set them by adding `export LANG=en_GB.UTF-8` and `export LC_ALL=en_GB.UTF-8` in your origin or login machines' `.profile`)

## Running tests
Go to the directory where the repository is cloned and run `python setup.py test --installation`. Tests will also be run automatically by [CircleCI](https://circleci.com/gh/ESMValGroup/ESMValTool).

## Code style
To increase the readability and maintainability or the ESMValTool source code, we aim to adhere to best practices and coding standards. All pull requests are reviewed and tested by one or more members of the core development team. For code in all languages, it is highly recommended that you split your code up in functions that are short enough to view without scrolling.

### Python
The standard document on best practices for Python code is [PEP8](https://www.python.org/dev/peps/pep-0008/) and there is [PEP257](https://www.python.org/dev/peps/pep-0257/) for documentation. We make use of [numpy style docstrings](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html) to document Python functions that are visible on [readthedocs](https://esmvaltool.readthedocs.io).

Most formatting issues in Python code can be fixed automatically by running the commands
```
isort some_file.py
```
to sort the imports in the standard way and
```
yapf -i some_file.py
```
to add/remove whitespace as required by the standard.

To check if your code adheres to the standard, go to the directory where the repository is cloned, e.g. `cd ESMValTool`.
and run
```
prospector esmvaltool/diag_scripts/your_diagnostic/your_script.py
```
Run
```
python setup.py lint
```
to see the warnings about the code style of the entire project.

We use `pycodestyle` on CircleCI to automatically check that there are no formatting mistakes and Codacy for monitoring (Python) code quality. Running prospector locally will give you quicker and sometimes more accurate results.

### NCL
Because there is no standard best practices document for NCL, we use [PEP8](https://www.python.org/dev/peps/pep-0008/) for NCL code as well, with some minor adjustments to accomodate for differences in the languages. The most important difference is that for NCL code the indentation should be 2 spaces instead of 4.

### R
A document on best practices for R is [Hadley Wickham's R Style Guide](http://r-pkgs.had.co.nz/style.html). We partially check adherence to this style guide by using [lintr](https://cran.r-project.org/web/packages/lintr/index.html) on CircleCI. In the future we would also like to make use of [goodpractice](https://cran.r-project.org/web/packages/goodpractice/index.html) to assess the quality of R code.

### YAML
Please use `yamllint` to check that your YAML files do not contain mistakes.

## Documentation

### What should be documented

Any code documentation that is visible on [readthedocs](https://esmvaltool.readthedocs.io) should be well written and adhere to the standards for documentation for the respective language. Recipes should have a page in the *Recipes* section on readthedocs. This is also the place to document recipe options for the diagnostic scripts used in those recipes. Note that there is no need to write extensive documentation for functions that are not visible on readthedocs. However, adding a one line docstring describing what a function does is always a good idea.

### How to build the documentation locally
Go to the directory where the repository is cloned and run
```
python setup.py build_sphinx -Ea
```
Make sure that your newly added documentation builds without warnings or errors.

## Pull requests and code review
New development should preferably be done in a new git branch in the main ESMValTool github repository. However, for scientists requiring confidentiality, private repositories are available. It is recommended that you open a pull request early, as this will cause CircleCI to run the unit tests and Codacy to analyse your code. It's also easier to get help from other developers if your code is visible in a pull request.

You can view the results of the automatic checks below your pull request. If one of the tests shows a red cross instead of a green approval sign, please click the link and try to solve the issue. Note that this kind of automated checks make it easier to review code, but they are not flawless, so occasionally Codacy will report false positives.

### Diagnostic script contributions
A pull request with diagnostic code should preferably not introduce new Codacy issues. However, we understand that there is a limit to how much time can be spend on polishing code, so up to 10 new (non-trivial) issues is still an acceptable amount.

Never make changes to the esmvaltool core, e.g. a new preprocessor function, in diagnostic script pull requests. If you need to make this kind of change, create a separate pull request for it in the public repository.

### Contributing to the core of ESMValTool
Contributions to the core of ESMValTool should
  - Go into the public repository.
  - Preferably be covered by unit tests. Unit tests are mandatory for new preprocessor functions or modifications to existing functions. If you do not know how to start with writing unit tests, let us know in a comment on the pull request and a core development team member will try to help you get started.
 - Be accompanied by appropriate documentation.
 - Introduce no new issues on Codacy (but note that style issues reported in unit test code are not worth the effort of fixing).
