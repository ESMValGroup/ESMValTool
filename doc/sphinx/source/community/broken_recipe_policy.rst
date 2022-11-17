.. _broken-recipe-policy:

Broken recipe policy
====================

Recipes might stop working for different reasons. Among those are, for instance, withdrawal of datasets
used by the recipe (i.e. the recipe contains data that are no longer publicly available), backward incompatible development
of the ESMValTool including new features or retiring old ones as well as
changes to Python or used dependencies such as library functions.
In such cases, the :ref:`Maintaining a recipe<recipe-maintainer>` is contacted by the technical lead development team (`@ESMValGroup/technical-lead-development-team`_) to find
a solution, fixing the affected recipe and checking the scientific output after applying the fix. If no recipe maintainer is
available, such recipes will be flagged by the release manager during the
:ref:`Release schedule and procedure<preparation-new-release>` as "broken".
For this, the affected recipe will be listed under "broken recipes" in the :ref:`Changelog`, together with the version
number of the last known release in which the recipe was still working.
If a recipe continues to be broken for three releases of the ESMValTool (about one year) and no recipe maintainer could be found
during this time, the affected recipe and diagnostics will be retired. This means the recipe and diagnostic code are
removed from the ESMValTool main branch by the release manager and thus will not be included in future releases.
Only the scientific documentation of the recipe (and diagnostics) will be kept in the user and developer guide with an
additional note and link to the last release in which the recipe was still fully functional.

.. _`@ESMValGroup/technical-lead-development-team`: https://github.com/orgs/ESMValGroup/teams/technical-lead-development-team
