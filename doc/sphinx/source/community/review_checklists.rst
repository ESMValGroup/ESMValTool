.. _checklists:

Checklists for reviewing a pull request
=======================================

General tasks
-------------
After merging a pull request successfully:

*	Close related issue if existent
*	Delete feature branch

ESMValTool diagnostics
----------------------

+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Item                                | Comments                                                                                         |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Documentation added to user’s guide | Check that the scientific documentation of the new diagnostic has been added to the user’s guide |
|                                     |                                                                                                  |
|                                     | * ./doc/sphinx/source/recipes/recipe_<diagnostic>.rst                                            |
|                                     | * include documentation in ./doc/sphinx/source/recipes/index.rst                                 |
|                                     | * documentation follows template (./doc/sphinx/source/recipes/recipe_template.rst.template)      |
|                                     | * typos                                                                                          |
|                                     | * references                                                                                     |
|                                     | * configuration options                                                                          |
|                                     | * variables                                                                                      |
|                                     | * observations                                                                                   |
|                                     | * valid image files                                                                              |
|                                     | * resolution of image files (~150 dpi is usually enough; file size should be kept small)         |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| recipe                              | Check that new recipe has been added to the testing recipe                                       |
|                                     | ./esmvaltool/recipes/examples/recipe_check_obs.yml                                               |
|                                     |                                                                                                  |
|                                     | * documentation: description, authors, maintainer, references, projects                          |
|                                     | * provenance: themes, realms                                                                     |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| diagnostic script                   | Check that the new diagnostic script(s) meet(s) standards. This includes the following items:    |
|                                     |                                                                                                  |
|                                     | * In-code documentation                                                                          |
|                                     | * Code quality checks                                                                            |
|                                     |                                                                                                  |
|                                     |   (1) code quality (e.g. no hardcoded pathnames)                                                 |
|                                     |   (2) no Codacy errors reported                                                                  |
|                                     | * Re-use of existing functions whenever possible                                                 |
|                                     | * Provenance implemented                                                                         |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| run recipe                          | Make sure new diagnostic(s) is working by running the ESMValTool                                 |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Check output of diagnostic          | After successfully running the new recipe, check that                                            |
|                                     |                                                                                                  |
|                                     | * Netcdf output has been written                                                                 |
|                                     | * Output contains (some) valid values (e.g. not only nan or zeros)                               |
|                                     | * Provenance information has been written                                                        |
|                                     | * If applicable, check plots and compare with corresponding plots in the paper(s) cited          |
+-------------------------------------+--------------------------------------------------------------------------------------------------+

ESMValTool CMORizer scripts
---------------------------

+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Dataset description added to user’s | Check that new dataset has been added to the table of observations defined in the ESMValTool     |
| guide                               | user’s guide in section “Obtaining input data” (./doc/sphinx/source/input.rst).                  |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| BibTeX info file                    | Check that a BibTeX file (i.e. <dataset>.bibtex) defining the reference(s) for the new dataset   |
|                                     | has been created in ./esmvaltool/references/.                                                    |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| recipe_check_obs.yml                | Check that new dataset has been added to the testing recipe                                      |
|                                     | ./esmvaltool/recipes/examples/recipe_check_obs.yml                                               |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| CMORizer script                     | Check that the new CMORizer script (./esmvaltool/cmorizers/obs/cmorize_obs_<dataset>.py/.ncl/.r) |
|                                     | meets standards. This includes the following items:                                              |
|                                     |                                                                                                  |
|                                     | * In-code documentation (header) contains                                                        |
|                                     |                                                                                                  |
|                                     |   (1) download instructions                                                                      |
|                                     |   (2) reference(s)                                                                               |
|                                     | * Code quality checks                                                                            |
|                                     |                                                                                                  |
|                                     |   (1) code quality (e.g. no hardcoded pathnames)                                                 |
|                                     |   (2) no Codacy errors reported                                                                  |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Config file                         | If present, check config file <dataset>.yml in ./esmvaltool/cmorizers/obs/cmor_config/.          |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Run CMORizer                        | Make sure CMORizer is working by running ''cmorize_obs -c <config-file> -o <dataset>''           |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| Check output of CMORizer            | After successfully running the new CMORizer, check that                                          |
|                                     |                                                                                                  |
|                                     | * Output contains (some) valid values (e.g. not only nan or zeros)                               |
|                                     | * Metadata is defined properly                                                                   |
|                                     | Run ./esmvaltool/recipes/examples/recipe_check_obs.yml for new dataset                           |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| RAW data                            | Contact person in charge of ESMValTool data pool (on Mistral: Axel) and request to copy RAW data |
|                                     | to RAWOBS/Tier2 (Tier3) (on Mistral: /work/bd0854/DATA/ESMValTool2/RAWOBS)                       |
+-------------------------------------+--------------------------------------------------------------------------------------------------+
| CMORized data                       | Contact person in charge of ESMValTool data pool (on Mistral: Axel) and request to               |
|                                     |                                                                                                  |
|                                     | * Copy CMORized dataset to OBS/Tier2 (Tier3) (on Mistral: /work/bd0854/DATA/ESMValTool2/OBS)     |
|                                     | * Set file access rights for new dataset (on Mistral: /work/bd0854/DATA/set_rights.csh)          |
+-------------------------------------+--------------------------------------------------------------------------------------------------+

ESMValCore pull requests
------------------------

+---------------+----------------------------------------+
| Item          | Comments                               |
+---------------+----------------------------------------+
| Documentation	| In-code documentation                  |
+---------------+----------------------------------------+
| user’s guide  |                                        |
+---------------+----------------------------------------+
| Code quality  |                                        |
+---------------+----------------------------------------+
| Unit tests    | Check that unit test has been provided |
+---------------+----------------------------------------+
| Tests	Codacy  |                                        |
+---------------+----------------------------------------+
| CircleCI      |                                        |
+---------------+----------------------------------------+
| VM            |                                        |
+---------------+----------------------------------------+
| Installation  |                                        |
+---------------+----------------------------------------+
