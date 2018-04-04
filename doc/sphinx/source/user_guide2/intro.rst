************
Introduction
************

The Earth System Model Evaluation Tool (ESMValTool) is a community-development that aims at improving diagnosing and understanding of the causes and effects of model biases and inter-model spread. The ESMValTool is open to both users and developers encouraging open exchange of diagnostic source code and evaluation results from the Coupled Model Intercomparison Project (CMIP) ensemble. This will facilitate and improve ESM evaluation beyond the state-of-the-art and aims at supporting the activities within CMIP and at individual modelling centers. We envisage running the ESMValTool routinely on the CMIP model output utilizing observations available through the Earth System Grid Federation (ESGF) in standard formats (obs4MIPs) or made available at ESGF nodes.

The goal is to develop a benchmarking and evaluation tool that produces well-established analyses as soon as model output from CMIP simulations becomes available, e.g., at one of the central repositories of the ESGF. This is realized through standard namelists that reproduce a certain set of diagnostics and performance metrics that have demonstrated its importance in benchmarking Earth System Models (ESMs) in a paper or assessment report, such as Chapter 9 of the Intergovernmental Panel on Climate Change (IPCC) Fifth Assessment Report (AR5) (Flato et al., 2013). The expectation is that in this way a routine and systematic evaluation of model results can be made more efficient, thereby enabling scientists to focus on developing more innovative methods of analysis rather than constantly having to "reinvent the wheel".

In parallel to standardization of model output, the ESGF also hosts observations for Model Intercomparison Projects (obs4MIPs) and reanalyses data (ana4MIPs). obs4MIPs provides open access data sets of satellite data that are comparable in terms of variables, temporal and spatial frequency, and periods to CMIP model output (Taylor et al., 2012). The ESMValTool utilizes these observations and reanalyses from ana4MIPs plus additionally available observations in order to evaluate the models performance. In many diagnostics and metrics, more than one observational data set or meteorological reanalysis is used to assess uncertainties in observations.

Objectives and approach
=======================

The main idea of the ESMValTool is to provide a broad suite of diagnostics which can be performed easily when new model simulations are run. The suite of diagnostics needs to be broad enough to reflect the diversity and complexity of Earth System Models, but must also be robust enough to be run routinely or semi-operationally.
In order the address these challenging objectives the ESMValTool is conceived as a framework which allows community contributions to be bound into a coherent framework.

License
=======

The ESMValTool is released under the Apache License, version 2.0 and citation
of the ESMValTool paper ("Software Documentation Paper") is kindly requested
upon use alongside with the software doi (doi:10.17874/ac8548f0315)
and version number:

* Eyring et al., ESMValTool (v1.0) - a community diagnostic and performance metrics tool for routine evaluation of Earth System Models in CMIP, Geosci. Model Dev., 9, 1747-1802, 2016.*

Besides the above citation, users are kindly asked to register any journal
articles (or other scientific documents) that use the software at the
ESMValTool webpage (http://www.esmvaltool.org/). Citing the Software
Documentation Paper and registering your paper(s) will serve to document the
scientific impact of the Software, which is of vital importance for securing
future funding. You should consider this an obligation if you have taken
advantage of the ESMValTool, which represents the end product of considerable
effort by the development team.

Architecture
============

:numref:`fig_schematic` shows a schematic of the ESMValTool architecture: the workflow manager (controlled by the Python script "main.py") runs a set of diagnostics on data provided by, for instance, a data archive. The configuration and the settings of each diagnostic are specified in namelists read and passed to the diagnostics by the workflow manager. The results which typically comprise of netCDF files and/or plots are stored in output folders along with log-files summarizing the data used, references, and technical details to ensure traceability and reproducibility of the results.

.. _fig_schematic:
.. figure::  ../figures/schematic.png
   :align:   center

   Schematic of the system architecture. The workflow manager (main.py) passes information to the diagnostics; results and log-files are written to dedicated folders.
