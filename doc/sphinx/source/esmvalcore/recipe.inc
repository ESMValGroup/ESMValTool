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

Note that the author, project, and references will need to be included in the
config-references.yml file. The author name uses the format:
last name (clipped after 4 letters) underscore first name (clipped to two letters)
and all letters are lowercase. 
For instance, Mickey Mouse would be: mous_mi.
This is unlikely to be the same as the github user name.

Datasets
========

The datasets section includes:

- dataset name
- Project (CMIP5 or 6, observations...)
- experiment (historical/ RCP8.5 etc...)
- Ensemble member
- The time range

Note that this section is not required, as datasets 
can also be provided in the `Diagnostics`_ section.


Preprocessors
=============

The preprocessor section of the recipe includes 
one or more preprocesors, each of which
may call one or several preprocessor functions.

Each preprocessor section includes.

- A preprocessor name.
- Which preprocesor functions to apply
- Any Arguments given to the preprocessors
- The order that the preprocesor functions are applied can also be specified here.

If only the default preprocessor is needed, then this section can be omitted.


Diagnostics
===========

The diagnostics section includes one or more
diagnostics. Each diagnostics will have: 

- A list of which variables to load
- A description of the variables (optional)
- Which preprocessor to apply to each variable
- The script to run
- The diagnostics can also include an additional_datasets section, if specific datasets are linked to specific diagnostics.






Brief introduction to YAML
==========================

While .yaml is a relatively common format, maybe users may not have 
encountered this language before. The key information about this format is:

- Yaml is a human friendly markup language.
- Yaml is commonly used for configuration files.
- the syntax is relatively straightforward
- Indentation matters a lot (like python)!
- yaml is case sensitive
- More details available here: https://learn.getgrav.org/16/advanced/yaml


