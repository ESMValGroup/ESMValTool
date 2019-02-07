Contributions are very welcome.

If you have a bug to report, please do so using the [issues tab on the ESMValTool github repository](https://github.com/ESMValGroup/ESMValTool/issues).

If you would like to contribute a new diagnostic and recipe or a new feature, please discuss your idea with the development team before getting started, to avoid disappointment later. A good way to do this is to open an [issue on GitHub]((https://github.com/ESMValGroup/ESMValTool/issues). This is also a good way to get help.

# Developing
To get started developing esmvaltool or developing/porting diagnostics, follow the instructions below. More detailed instructions can be found in the [manual](https://esmvaltool.readthedocs.io) under Developer's Guide.

## Getting started
To install in development mode, follow these instructions.
- Install gcc, g++ and gfortran if these are not available on your system. On Debian based systems, this can be done by running `apt install build-essential gfortran`, on managed systems you can often use the `module avail` command to see what compilers are available.
- [Download and install conda](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html) (this should be done even if the system in use already has a preinstalled version of conda, as problems have been reported with NCL when using such a version)
- To make the `conda` command availble, add `source <prefix>/etc/profile.d/conda.sh` to your `.bashrc` file and restart your shell. If using (t)csh shell, add `source <prefix>/etc/profile.d/conda.csh` to your `.cshrc`/`.tcshrc` file instead.
- Update conda: `conda update -y conda`
- Create a conda environment: `conda create -y -n esmvaltool python=3`
- Activate the esmvaltool environment: `conda activate esmvaltool`
- Clone the ESMValTool github repository: `git clone git@github.com:ESMValGroup/ESMValTool`
- Go to the esmvaltool directory: `cd ESMValTool`
- Update the esmvaltool conda environment `conda env update`
- Install in development mode: `pip install -e '.[develop]'`. If you are installing behind a proxy that does not trust the usual pip-urls you can declare them with the option `--trusted-host`, e.g. `pip install --trusted-host=pypi.python.org --trusted-host=pypi.org --trusted-host=files.pythonhosted.org -e .[develop]`
- If you want to use R diagnostics, run `Rscript esmvaltool/install/R/setup.R` to install the R dependences.
- Test that your installation was succesful by running `esmvaltool -h`.

## Running ESMValTool
- Review `config-user.yml`. To customize for your system, create a copy, edit and use the command line option `-c` to instruct `esmvaltool` to use your custom configuration.
- Available recipes are located in the directory `esmvaltool/recipes`.
- Run e.g. `esmvaltool -c ~/config-user.yml examples/recipe_python.yml

## Running tests
Go to the directory where the repository is cloned and run `python setup.py test --installation`. Tests will also be run automatically by [CircleCI](https://circleci.com/gh/ESMValGroup/ESMValTool).

## Code style
To increase the readability and maitainability or the ESMValTool source code, we aim to adhere to best practices and coding standards. All Pull Requests are reviewed and tested by one more members of the core development team. For code in all languages, it is highly recommended that you split your code up in functions that are short enough to view without scrolling.

### Python
The standard document on best practices for Python is [PEP8](https://www.python.org/dev/peps/pep-0008/) and [PEP257](https://www.python.org/dev/peps/pep-0257/) for documentation. We make use of [numpy style docstrings](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html) to document Python functions.

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
and run `prospector esmvaltool/diag_scripts/your_diagnostic/your_script.py`.
Run `python setup.py lint` to see the warnings about the code style of the entire project.

We use `pycodestyle` on CircleCI to automatically check that there are no formatting mistakes and Codacy for monitoring (Python) code quality. Running prospector locally will give you quicker and sometimes more accurate results.

### NCL
Because there is no standard best practices document for NCL, we use [PEP8](https://www.python.org/dev/peps/pep-0008/) for NCL code as well, with some minor adjustments to accomodate for differences in the languages. The most important difference is that for NCL code the indentation should be 2 spaces instead of 4.

### R
A document on best practices for R is [Hadley Wickham's R Style Guide](http://r-pkgs.had.co.nz/style.html). We partially check adherence to this style guide by using [lintr](https://cran.r-project.org/web/packages/lintr/index.html) on CircleCI. In the future we would also like to make use of [goodpractice](https://cran.r-project.org/web/packages/goodpractice/index.html) to assess the quality of R code.

### YAML
Please use `yamllint` to check that your YAML files do not contain mistakes.

## Documentation

### What should be documented

Any code documentation that is visible on [readthedocs](https://esmvaltool.readthedocs.io) should be well written and adhere to the standards for documentation for the respective language. Recipes should have a page in the *Recipes* section on readthedocs. This is also the place to document recipe options for the diagnostic scripts used in those recipes.

### How to build the documentation locally
Go to the directory where the repository is cloned and run `python setup.py build_sphinx -Ea`. Make sure that your newly added documentation builds without warnings or errors.
