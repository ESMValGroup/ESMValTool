.. _recipe:

*****************
ESMValTool recipe
*****************

Recipes are the instructions telling ESMValTool about the user who wrote the
recipe, the datasets which need to be run, the preprocessors that need to be
applied, and the diagnostics which need to be run over the preprocessed data.
This information is provided to ESMValTool in the recipe sections:
`Documentation`_, `Datasets`_, `Preprocessors`_ and `Diagnostics`_,
respectively.


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

Note that the authors, projects, and references will need to be included in the
``config-references.yml`` file. The author name uses the format:
`surname_name`. For instance, Mickey Mouse would be: `mouse_mickey`.
Also note that this username is unlikely to be the same as the github
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
- The order that the preprocesor functions are applied can also be specified using the ``custom_order`` preprocesor function.

The following preprocessor is an example of a preprocessor that contains
multiple preprocessor functions:

.. code-block:: yaml

    preprocessors:
      prep_map:
        regrid:
          target_grid: 1x1
          scheme: linear
        time_average:
        multi_model_statistics:
          span: overlap
          statistics: [mean ]

If only the default preprocessor is needed, then this section can be omitted.


Diagnostics
===========

The diagnostics section includes one or more diagnostics. Each diagnostics will
have:

- A list of which variables to load
- A description of the variables (optional)
- Which preprocessor to apply to each variable
- The script to run
- The diagnostics can also include an optional ``additional_datasets`` section.

The ``additional_datasets`` can add datasets beyond those listed in the the
`Datasets`_ section. This is useful if specific datasets need to be linked with
a specific diagnostics. The addition datasets can be used to add variable
specific datasets. This is also a good way to add observational datasets can be
added to the diagnostic.

The following example, taken from recipe_ocean_example.yml, shows a diagnostic
named `diag_map`, which loads the temperature at the ocean surface between
the years 2001 and 2003 and then passes it to the prep_map preprocessor.
The result of this process is then passed to the ocean diagnostic map scipt,
``ocean/diagnostic_maps.py``.

.. code-block:: yaml

    diagnostics:

    diag_map:
      description: Global Ocean Surface regridded temperature map
      variables:
        tos: # Temperature at the ocean surface
          preprocessor: prep_map
          start_year: 2001
          end_year: 2003
      scripts:
        Global_Ocean_Surface_regrid_map:
          script: ocean/diagnostic_maps.py

To define a variable/dataset combination, the keys in the diagnostic section
are combined with the keys from datasets section. If two versions of the same
key are provided, then the key in the datasets section will take precedence
over the keys in variables section. For many recipes it makes more sense to
define the ``start_year`` and ``end_year`` items in the variable section, because the
diagnostic script assumes that all the data has the same time range.

Note that the path to the script provided in the `script` option should be
either:

1. the absolute path to the script.
2. the path relative to the ``esmvaltool/diag_scripts`` directory.


As mentioned above, the datasets are provided in the `Diagnostics`_ section
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
