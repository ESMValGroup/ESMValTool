.. _recipe:

*****************
ESMValTool recipe
*****************

Recipes are the instructions telling ESMValTool
about the user who wrote the recipe,
the datasets which need to be run,
the preprocessors that need to be applied,
and the diagnostics which need to be run over the preprocessed data.
This information is provided to ESMValTool in the recipe sections:
`Documentation`_, `Datasets`_, `Preprocessors`_ and `Diagnostics`_, respectively.


Documentation
=============

The documentation section includes:

- The recipe's author's user name
- A description of the recipe
- The user name of the maintainer
- A list of scientific references
- the project or projects associated with the recipe.

For example, please see the documentation section from the recipe:
recipe_ocean_amoc.yml.

.. code-block:: yaml

    documentation:
      description: |
        Recipe to produce time series figures of the derived variable, the
        Atlantic meriodinal overturning circulation (AMOC).
        This recipe also produces transect figures of the stream functions for
        the years 2001-2004.

      authors:
        - demo_le

      maintainer:
        - demo_le

      references:
        - demora2018gmd

      projects:
        - ukesm

Note that the author, project, and references will need to be included in the
config-references.yml file. The author name uses the format:
last name (clipped after 4 letters) underscore first name (clipped to two
letters) and all letters are lowercase. For instance, Mickey Mouse would be:
mous_mi. Also note that this username is unlikely to be the same as the github
user name.



Datasets
========

The datasets section includes:

- dataset name
- Project (CMIP5 or 6, observations...)
- experiment (historical/ RCP8.5 etc...)
- Ensemble member
- The time range
- The model grid, gn or gr, (CMIP6 only).

For example, a datasets section could be:

.. code-block:: yaml

    datasets:
      - {dataset: CanESM2, project: CMIP5, exp: historical, ensemble: r1i1p1, start_year: 2001, end_year: 2004}
      - {dataset: UKESM1-0-LL, project: CMIP6, exp: historical, ensemble: r1i1p1f2, start_year: 2001, end_year: 2004, grid: gn}


Note that this section is not required, as datasets can also be provided in the
`Diagnostics`_ section.


Preprocessors
=============

The preprocessor section of the recipe includes one or more preprocesors, each
of which may call one or several preprocessor functions.

Each preprocessor section includes:

- A preprocessor name.
- A list of preprocesor functions to apply
- Any Arguments given to the preprocessor functions.
- The order that the preprocesor functions are applied can also be specified.

The following preprocessor is an example of a preprocessor that contains
multiple preprocessor functions:

.. code-block:: yaml

    preprocessors:
      prep_timeseries: # For regional volume averaged (at the surface)
        custom_order: true
        extract_region:
          start_longitude: -80.
          end_longitude: 30.
          start_latitude: -80.
          end_latitude: 80.
        extract_volume:
          z_min: 0.
          z_max: 100.
        volume_statistics:
          operator: mean
        multi_model_statistics:
          span: overlap
          statistics: [mean ]

If only the default preprocessor is needed, then this section can be omitted.


Diagnostics
===========

The diagnostics section includes one or more
diagnostics. Each diagnostics will have:

- A list of which variables to load
- A description of the variables (optional)
- Which preprocessor to apply to each variable
- The script to run
- The diagnostics can also include an additional_datasets section, if specific
  datasets are linked to specific diagnostics.

The following example, taken from recipe_ocean_amoc.yml, shows a diagnostic
named `diag_timeseries_amoc`, which calculates a derived variable and
sends it to the ocean time series scipt.

.. code-block:: yaml

    diagnostics:
      diag_timeseries_amoc:
        description: atlantic_meridional_overturning_circulation
        variables:
          amoc:
            mip: Omon
            derive: true
            force_derivation: false
        additional_datasets:
         - {dataset: CanESM2, project: CMIP5, exp: historical, ensemble: r1i1p1, start_year: 1860, end_year: 2004}
        scripts:
          AMOC_timeseries:
            script: ocean/diagnostic_timeseries.py
            moving_average: 6 years


As mentionned above, the datasets are provided in the `Diagnostics`_ section
in this section. However, they could also be included in the `Datasets`_
section.



Brief introduction to YAML
==========================

While .yaml is a relatively common format, maybe users may not have
encountered this language before. The key information about this format is:

- Yaml is a human friendly markup language.
- Yaml is commonly used for configuration files.
- the syntax is relatively straightforward
- Indentation matters a lot (like python)!
- yaml is case sensitive
- A yml tutorial is available here: https://learnxinyminutes.com/docs/yaml/
- A yml quick reference card is available here: https://yaml.org/refcard.html
- ESMValTool uses the yamllint linter tool: http://www.yamllint.com
