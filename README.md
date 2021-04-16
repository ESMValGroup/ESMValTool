[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Documentation Status](https://readthedocs.org/projects/esmvaltool/badge/?version=latest)](https://esmvaltool.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3401363.svg)](https://doi.org/10.5281/zenodo.3401363)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/ESMValGroup?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![CircleCI](https://circleci.com/gh/ESMValGroup/ESMValTool/tree/master.svg?style=svg)](https://circleci.com/gh/ESMValGroup/ESMValTool/tree/master)
[![example branch parameter](https://github.com/github/docs/actions/workflows/main.yml/badge.svg?branch=master)](https://github.com/ESMValGroup/ESMValTool/actions)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/79bf6932c2e844eea15d0fb1ed7e415c)](https://www.codacy.com/gh/ESMValGroup/ESMValTool?utm_source=github.com&utm_medium=referral&utm_content=ESMValGroup/ESMValTool&utm_campaign=Badge_Coverage)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/79bf6932c2e844eea15d0fb1ed7e415c)](https://www.codacy.com/gh/ESMValGroup/ESMValTool?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ESMValGroup/ESMValTool&amp;utm_campaign=Badge_Grade)
[![Docker Build Status](https://img.shields.io/docker/cloud/build/esmvalgroup/esmvaltool.svg)](https://hub.docker.com/r/esmvalgroup/esmvaltool/)
[![Anaconda-Server Badge](https://anaconda.org/esmvalgroup/esmvaltool/badges/installer/conda.svg)](https://conda.anaconda.org/esmvalgroup)

![esmvaltoollogo](https://github.com/ESMValGroup/ESMValTool/blob/master/doc/sphinx/source/figures/ESMValTool-logo.png)

- [**Documentation**](https://docs.esmvaltool.org/en/latest/)
- [**ESMValTool Website**](https://www.esmvaltool.org/)
- [**ESMValTool Tutorial**](https://esmvalgroup.github.io/ESMValTool_Tutorial/index.html)
- [**ESMValGroup Project on GitHub**](https://github.com/ESMValGroup)
- [**Gallery**](https://docs.esmvaltool.org/en/latest/gallery.html)

# Introduction

ESMValTool is a community-developed climate model diagnostics and evaluation software package, driven
both by computational performance and scientific accuracy and reproducibility. ESMValTool is open to both
users and developers, encouraging open exchange of diagnostic source code and evaluation results from the
Coupled Model Intercomparison Project [CMIP](https://www.wcrp-climate.org/wgcm-cmip) ensemble. For a
comprehensive introduction to ESMValTool please visit our
[documentation](https://docs.esmvaltool.org/en/latest/introduction.html) page.

# Running esmvaltool

Diagnostics from ESMValTool are run using [recipe](https://docs.esmvaltool.org/en/latest/recipes/index.html)
files that contain pointers to the requested data types, directives for the preprocessing steps that data
will be subject to, and directives for the actual diagnostics that will be run with the now preprocessed data.
Data preprocessing is done via the [ESMValCore](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/quickstart/index.html) package, a pure Python, highly-optimized scientific library, developed by the ESMValTool core developers,
and that performs a number of common analysis tasks
such as regridding, masking, levels extraction etc. [Diagnostics](https://docs.esmvaltool.org/en/latest/develop/diagnostic.html) are written in a variety of programming languages (Python, NCL, R, Julia) and are developed by the wider
scientific community, and included after a scientific and technical review process.

# Input data

ESMValTool can run with the following types of data as input:

- CMIP5
- CMIP6
- OBS, OBS6
- obs4mips
- ana4mips
- CORDEX

# Getting started

Please see [getting started](https://docs.esmvaltool.org/en/latest/quickstart/index.html) on readthedocs as well as [ESMValTool tutorial](https://esmvalgroup.github.io/ESMValTool_Tutorial/index.html). The tutorial is a set of lessons that together teach skills needed to work with ESMValTool in climate-related domains.

## Getting help

The easiest way to get help if you cannot find the answer in the documentation on [readthedocs](https://docs.esmvaltool.org), is to open an [issue on GitHub](https://github.com/ESMValGroup/ESMValTool/issues).

## Contributing

If you would like to contribute a new diagnostic or feature, please have a look at our [contribution guidelines](https://docs.esmvaltool.org/en/latest/community/index.html).
