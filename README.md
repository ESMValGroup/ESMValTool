# ESMValTool
[![DOIBadge](https://img.shields.io/badge/DOI-10.17874%2Fac8548f0315-blue.svg)](http://dx.doi.org/10.17874/ac8548f0315)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/79bf6932c2e844eea15d0fb1ed7e415c)](https://www.codacy.com/app/ESMValGroup/ESMValTool?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ESMValGroup/ESMValTool&amp;utm_campaign=Badge_Grade)
[![Codacy Coverage Badge](https://api.codacy.com/project/badge/Coverage/79bf6932c2e844eea15d0fb1ed7e415c)](https://www.codacy.com/app/ESMValGroup/ESMValTool?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ESMValGroup/ESMValTool&amp;utm_campaign=Badge_Coverage)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/ESMValGroup?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![CircleCI](https://circleci.com/gh/ESMValGroup/ESMValTool.svg?style=svg)](https://circleci.com/gh/ESMValGroup/ESMValTool)
[![Docker Build Status](https://img.shields.io/docker/build/esmvalgroup/esmvaltool.svg)](https://hub.docker.com/r/esmvalgroup/esmvaltool/)

ESMValTool: A community diagnostic and performance metrics tool for routine evaluation of Earth system models in CMIP

## Developing
This is the development branch for version 2 of ESMValTool. To get started developing esmvaltool or developing/porting diagnostics, follow the instructions below. More detailed instructions can be found in the [manual](http://esmvaltool.readthedocs.io/en/refactoring_backend/).

### Getting started
To install in development mode, follow these instructions.
- [Download and install conda](https://conda.io/docs/user-guide/install/linux.html)
- Update conda: `conda update -y conda`
- Create a conda environment: `conda create -y -n esmvaltool python=3`
- Activate the esmvaltool environment: `source activate esmvaltool`
- Clone the ESMValTool github repository: `git clone git@github.com/ESMValGroup/ESMValTool`
- Go to the esmvaltool directory: `cd ESMValTool`
- Check out the version 2 development branch: `git checkout version2_development`
- Update the esmvaltool conda environment `conda env update`
- Install in development mode: `pip install -e .[develop]`
- Test that your installation was succesful by running `esmvaltool -h`.
- Review `config-user.yml`. To customize for your system, create a copy, edit and use the command line option `-c` to instruct `esmvaltool` to use your custom configuration.
- Available namelists are located in the directory `esmvaltool/namelists`.

### Running tests
Go to the directory where the repository is cloned and run `./setup.py test`. Tests will also be run automatically by CircleCI.

### Code style
First go to the directory where the repository is cloned, e.g. `cd ESMValTool`.
- To review if your own code follows our coding standards, run `prospector esmvaltool/diag_scripts/your_diagnostic/your_script.py`.
- Run `./setup.py lint` to see the warnings about the code style of the entire project.

We use Codacy for monitoring (Python) code quality. However, running prospector locally will generally give you quicker and sometimes more accurate results. Note that Codacy does not install dependencies, so getting a warning "Unable to import 'external_library'" is probably not a real issue.   

### Building documentation
Go to the directory where the repository is cloned and run `./setup.py build_sphinx`
