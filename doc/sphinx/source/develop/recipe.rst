Recipe
******

Writing a basic recipe
======================
The user will need to write a basic recipe to be able to run their own personal diagnostic.
An example of such a recipe is found in `esmvaltool/recipes/recipe_my_personal_diagnostic.yml`.
For general guidelines with regards to ESMValTool recipes please consult the User Guide;
the specific parameters needed by a recipe that runs a personal diagnostic are:

.. code-block:: yaml

  scripts:
    my_diagnostic:
    script: /path/to/your/my_little_diagnostic.py

i.e. the full path to the personal diagnostic that the user needs to run.
