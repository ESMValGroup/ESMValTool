.. _broken-recipe-policy:

Broken/outdated recipe policy
==============================

Recipes might stop working for different reasons. Among those are, for instance, withdrawal of datasets
used by the recipe (i.e. the recipe contains data that are no longer publicly available), backward incompatible development
of the ESMValTool including new features or retiring old ones as well as
changes to Python or used dependencies such as library functions. Recipes might also be considered outdated by the :ref:`recipe maintainer<recipe-maintainer>`.

In case a recipe stopped working, the :ref:`recipe maintainer<recipe-maintainer>` is contacted by the technical lead development team (`@ESMValGroup/technical-lead-development-team`_) to find
a solution, fixing the affected recipe and checking the scientific output after applying the fix. If no recipe maintainer is
available, such recipes will be flagged by the release manager during the
:ref:`Release schedule and procedure<preparation-new-release>` as "broken".
For this, the affected recipe will be added to the :ref:`list of broken recipes <broken-recipe-list>`, together with the version
number of the last known release in which the recipe was still working.
If a recipe continues to be broken for three releases of the ESMValTool (about one year) and no recipe maintainer could be found
during this time, or the affected recipe and diagnostics will be retired. This means the recipe and diagnostic code are
moved from the folders ``esmvaltool/recipes`` and ``esmvaltool/diag_scripts`` to ``archive/legacy_recipes`` and ``archive/legacy_diag_scripts`` by the release manager. In case of recipes considered
outdated by the :ref:`recipe maintainer<recipe-maintainer>`, a pull request can be opened to move recipe and diagnostic code to the ``archive`` folders at any time.

In both cases, the scientific documentation of the recipe (and diagnostics) will be kept in the user and developer guide with an
additional note and link to the last release in which the recipe was still fully functional and a reference to the ``archive`` folder.
The documentation of the retired recipe in the user and developer guide is moved into the category "legacy recipes" to distinguish
between currently available and functional recipes and such legacy recipes.

Special care needs to be taken when retiring diagnostic scripts to avoid unintentionally affecting other recipes or diagnostics.

.. _`@ESMValGroup/technical-lead-development-team`: https://github.com/orgs/ESMValGroup/teams/technical-lead-development-team
