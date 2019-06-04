.. _recipe:

*****************
ESMValTool recipe
*****************

Recipes are the instructions telling ESMValTool
about the user who wrote the recipe,
the datasets which need to be run,
the preprocessors that need to be applied,
and the diagnostics which need to be run over the preprocessed data.

Documentation:User information:
- The recipe's author's user name
- A description of the recipe
- the name of the maintainer
- a list of scientific references
- the project or projects associated with the recipe.

Note that the author, project, and references will need to be included in the
config-references.yml file.

Datasets:
- Which Models to run
- Project (CMIP5 or 6, observations...)
- Scenario (historicalk/ RCP8.5 etc...)
- Ensemble member
- The Time resolution
- The time range

Preprocessors:
- Which preprocesors to apply
- The order to apply them
- Any Arguments given to the preprocessors
- A name is required for each group (or chain) of preprocessors.

Diagnostics:
- Which variables to load
- A description of the variables
- Which preprocessor to apply to each variable
- Which diagnostics to run
- The diagnostics can also include n addition_datasets section, if specific datasets are linked to specific diagnostics.






Brief introduction to YAML
==========================

While .yaml is a relatively common format, maybe users may not have 
encountered this language before. The key information about this format is:
- Yaml is a human friendly markup language.
- yaml is commonly used for configuration files.
- the syntax is relatively straightforward
- Indentation matters a lot (like python)!
- yaml is case sensitive
- More details available here: https://learn.getgrav.org/16/advanced/yaml


