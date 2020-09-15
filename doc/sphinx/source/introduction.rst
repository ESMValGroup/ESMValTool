Introduction
************

About
=====

The Earth System Model Evaluation Tool (ESMValTool) is a
community-development that aims at improving diagnosing and
understanding of the causes and effects of model biases and inter-model
spread. The ESMValTool is open to both users and developers encouraging
open exchange of diagnostic source code and evaluation results from the
Coupled Model Intercomparison Project (CMIP) ensemble. This will
facilitate and improve ESM evaluation beyond the state-of-the-art and
aims at supporting the activities within CMIP and at individual
modelling centers. We envisage running the ESMValTool routinely on the
CMIP model output utilizing observations available through the Earth
System Grid Federation (ESGF) in standard formats (obs4MIPs) or made
available at ESGF nodes.

The goal is to develop a benchmarking and evaluation tool that produces
well-established analyses as soon as model output from CMIP simulations
becomes available, e.g., at one of the central repositories of the ESGF.
This is realized through standard recipes that reproduce a certain set
of diagnostics and performance metrics that have demonstrated its
importance in benchmarking Earth System Models (ESMs) in a paper or
assessment report, such as Chapter 9 of the Intergovernmental Panel on
Climate Change (IPCC) Fifth Assessment Report (AR5) (Flato et al.,
2013). The expectation is that in this way a routine and systematic
evaluation of model results can be made more efficient, thereby enabling
scientists to focus on developing more innovative methods of analysis
rather than constantly having to "reinvent the wheel".

In parallel to standardization of model output, the ESGF also hosts
observations for Model Intercomparison Projects (obs4MIPs) and
reanalyses data (ana4MIPs). obs4MIPs provides open access data sets of
satellite data that are comparable in terms of variables, temporal and
spatial frequency, and periods to CMIP model output (Taylor et al.,
2012). The ESMValTool utilizes these observations and reanalyses from
ana4MIPs plus additionally available observations in order to evaluate
the models performance. In many diagnostics and metrics, more than one
observational data set or meteorological reanalysis is used to assess
uncertainties in observations.

The main idea of the ESMValTool is to provide a broad suite of
diagnostics which can be performed easily when new model simulations are
run. The suite of diagnostics needs to be broad enough to reflect the
diversity and complexity of Earth System Models, but must also be robust
enough to be run routinely or semi-operationally. In order the address
these challenging objectives the ESMValTool is conceived as a framework
which allows community contributions to be bound into a coherent
framework.

License
=======

The ESMValTool is released under the Apache License, version 2.0.
Citation of the ESMValTool paper ("Software Documentation Paper") is
kindly requested upon use, alongside with the software DOI for
ESMValTool
(`doi:10.5281/zenodo.3401363 <https://doi.org/10.5281/zenodo.3401363>`__)
and ESMValCore
(`doi:10.5281/zenodo.3387139 <https://doi.org/10.5281/zenodo.3387139>`__)
and version number:

-  Righi, M., Andela, B., Eyring, V., Lauer, A., Predoi, V., Schlund,
   M., Vegas-Regidor, J., Bock, L., Brötz, B., de Mora, L., Diblen, F.,
   Dreyer, L., Drost, N., Earnshaw, P., Hassler, B., Koldunov, N.,
   Little, B., Loosveldt Tomas, S., and Zimmermann, K.: Earth System
   Model Evaluation Tool (ESMValTool) v2.0 – technical overview, Geosci.
   Model Dev., 13, 1179–1199, https://doi.org/10.5194/gmd-13-1179-2020,
   2020. 

Besides the above citation, users are kindly asked to register any
journal articles (or other scientific documents) that use the software
at the ESMValTool webpage (http://www.esmvaltool.org/). Citing the
Software Documentation Paper and registering your paper(s) will serve to
document the scientific impact of the Software, which is of vital
importance for securing future funding. You should consider this an
obligation if you have taken advantage of the ESMValTool, which
represents the end product of considerable effort by the development
team.

What ESMValTool can do for you
==============================

The ESMValTool applies a great variety of standard diagnostics and
metrics, and produces a collection of netCDF and graphical files
(plots). Thus, the tool needs a certain amount of input from the user so
that it can:

-  establish the correct input and output parameters and the structured
   workflow;
-  acquire the correct data;
-  execute the workflow; and
-  output the desired collective data and media.

To facilitate these four steps, the user has control over the tool via
two main input files: the :ref:`user configuration file <config-user>`
and the :ref:`recipe <esmvalcore:recipe>`. The configuration file sets
user and site-specific parameters (like input and output paths, desired
output graphical formats, logging level, etc.), whereas the recipe file
sets data, preprocessing and diagnostic-specific parameters (data
parameters grouped in the datasets sections, preprocessing steps for
various preprocessors sections, variables' parameters and
diagnostic-specific instructions grouped in the diagnostics sections).
The configuration file may be used for a very large number of runs with
very minimal changes since most of the parameters it sets are
recyclable; the recipe file can be used for a large number of
applications, since it may include as many datasets, preprocessors and
diagnostics sections as the user deems useful.

Once the user configuration files and the recipe are at hand, the user
can start the tool. A schematic overview of the ESMValTool workflow is
depited in the figure below.

.. container::
   :name: figarch

   .. figure:: figures/schematic.png
      :alt: Schematic of the system architecture.
      :figclass: align-center

      Schematic of the system architecture.

For a generalized run scenario, the tool will perform the following
ordered procedures.

Data finding
------------

-  read the data requirements from the :ref:`datasets section
   <esmvalcore:Datasets>` of the recipe and assemble the data request to
   locate the data;
-  find the data using the specified root paths and DRS types in the
   configuration file (note the flexibility allowed by the
   :ref:`data finder
   <esmvalcore:findingdata>`);

Data selection
--------------

-  data selection is performed using the parameters specified in the
   :ref:`datasets section <esmvalcore:Datasets>` (including e.g. type of
   experiment, type of ensemble, time boundaries etc); data will be
   retrieved and selected for each variable that is specified in the
   :ref:`diagnostics <esmvalcore:Diagnostics>` section of the recipe;

Data fixing
-----------

-  the ESMValTool requires data to be in CMOR format; since errors in
   the data formatting are not uncommon, the ESMValTool performs
   :ref:`checks against the
   CMOR library and fixes small irregularities <esmvalcore:CMOR check and
   dataset-specific fixes>` (note that the degree of leniency is not
   very high).

Variable derivation
-------------------

-  :ref:`variable derivation <esmvalcore:Variable derivation>` (in the
   case of non CMOR-standard variables, most likely associated with
   observational datasets) is performed automatically before running the
   preprocessor;
-  if the variable definitions are already in the database then the user
   will just have to specify the variableto be derived in the
   :ref:`diagnostics
   <esmvalcore:Diagnostics>` section (as any other standard variable,
   just setting ``derive: true``).

Run the preprocessor
--------------------

-  if any :ref:`preprocessor section <esmvalcore:preprocessor>` is
   specified in the recipe file, then data will be loaded in memory as
   iris cubes and passed through the preprocessing steps required by the
   user and specified in the preprocessor section, using the specific
   preprocessing step parameters provided by the user as keys (for the
   parameter name) and values (for the paramater value); the
   preprocessing order is very imprtant since a number of steps depend
   on prior excution of other steps (e.g. :ref:`multimodel
   statistics <esmvalcore:Multi-model statistics>` can not be computed
   unless all models are on a common grid, hence a prior
   :ref:`regridding
   <esmvalcore:Horizontal regridding>` on a common grid is necessary);
   the preprocessor steps order can be set by the user as custom or the
   default order can be used;
-  once preprocessing has finished, the tool writes the data output to
   disk as netCDF files so that the diagnostics can pick it up and use
   it; the user will also be provided with a metadata file containing a
   summary of the preprocessing and pointers to its output. Note that
   writing data to disk between the preprocessing and the diagnostic
   phase is required to ensure multi-language support for the latter.

Run the diagnostics
-------------------

-  the last and most important phase can now be run: using output files
   from the preprocessor, the diagnostic scripts are executed using the
   provided diagnostics parameters.
