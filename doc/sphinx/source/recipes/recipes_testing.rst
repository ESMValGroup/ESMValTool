.. _recipes_testing:

Short test versions of scientific recipes to check for backward compatibility.
==========================================================================================================================

Overview
--------


These recipes are created to cover typical functionalities in the ESMValTool and allow to test them quickly. 
Each recipe should run less than 5 minutes to facilitate fast tests.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/testing/

   * recipe_deangelis15nat_fig1_fast.yml
   
Diagnostics are stored in diag_scripts/

   * deangelis15nat/deangelisf1b.py


User settings in recipes
-----------------------

The recipe recipe_deangelis15nat_fig1_fast.yml calls the first diagnostic (deangelisf1b.py) from the original recipe recipe_deangelis15nat.yml.
It can be run with CMIP5 and CMIP6 models for any duration.
Several flux variables (W m\ :sup:`-2`\) and up to 6 different model experiments can be handled.
Each variable needs to be given for each model experiment. The same experiments must
be given for all models. For testing purpose it was reduce to two models, 3 experiments and one year.
For a more detailed documentation see :ref: recipes_deangelis15nat 
