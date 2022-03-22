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

.. _Support-and-Contact:

Support
=======

Support for ESMValTool can be found in `ESMValTool Discussions page
<https://github.com/ESMValGroup/ESMValTool/discussions>`__
where users can open an issue and a member of the `User Engagement Team
<mailto:esmvaltool_user_engagement_team@listserv.dfn.de>`_ of ESMValTool
will reply as soon as possible. 
This is open for all general and technical questions on the ESMValTool:
installation, application, development, or any other question or comment
you may have.

.. _mailing-list:

User mailing list
-----------------

Subscribe to the ESMValTool announcements mailing list
`esmvaltool@listserv.dfn.de <mailto:esmvaltool@listserv.dfn.de>`__
to stay up to date about new releases, monthly online meetings, upcoming workshops, and trainings.

To subscribe, send an email to
`sympa@listserv.dfn.de <mailto:sympa@listserv.dfn.de?subject=subscribe%20esmvaltool>`_
with the following subject line:

-  *subscribe esmvaltool*

or

-  *subscribe esmvaltool YOUR_FIRSTNAME YOUR_LASTNAME*

The mailing list also has a `public archive <https://www.listserv.dfn.de/sympa/arc/esmvaltool>`_ online.




Monthly meetings
----------------

We have monthly online meetings using `zoom <https://zoom.us/>`__, anyone with
an interest in the ESMValTool is welcome to join these meetings to connect with
the community.
These meetings are always announced in an issue
on the `ESMValTool <https://github.com/ESMValGroup/ESMValTool/issues>`_
repository and on the mailing-list_.

.. _core-team:

Core development team
---------------------

-  Deutsches Zentrum für Luft- und Raumfahrt (DLR), Institut für Physik
   der Atmosphäre, Germany (Co-PI)

   - ESMValTool Core Co-PI and Developer: contact for requests to use the ESMValTool and for collaboration with the development team, access to the PRIVATE GitHub repository.

-  Met Office, United Kingdom (Co-PI)
-  Alfred Wegener institute (AWI) Bremerhaven, Germany
-  Barcelona Supercomputing Center (BSC), Spain
-  Netherlands eScience Center (NLeSC), The Netherlands
-  Ludwig Maximilian University of Munich, Germany
-  Plymouth Marine Laboratory (PML), United Kingdom
-  Swedish Meteorological and Hydrological Institute (SMHI), Sweden
-  University of Bremen, Germany
-  University of Reading, United Kingdom

Recipes and diagnostics
-----------------------

Contacts for specific diagnostic sets are the respective authors, as
listed in the corresponding :ref:`recipe and diagnostic documentation<recipes>`
and in the source code.


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
