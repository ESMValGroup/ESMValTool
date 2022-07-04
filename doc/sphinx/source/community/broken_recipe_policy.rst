Broken recipe policy
====================

Recipes might stop to work properly because of different reasons. Among those is the continued development of the ESMValTool
introducing new features or retiring old ones as well as e.g. changes of Python or used dependencies such as library functions.
In such cases, the :ref:`recipe maintainer<Maintaining a recipe>` is contacted by the technical core development team to find
a solution for fixing the affected recipe and checking the scientific validity of the fix. If no recipe maintainer is available,
such recipes will be flagged by the release manager during the
:ref:`Release schedule and procedure for ESMValCore and ESMValTool release procedure` as "broken".
If a recipe continues to be broken for three releases of the ESMValTool (about one year) and no recipe maintainer could be found
during this time, the affected recipe and diagnostics will be retired. This means the recipe and diagnostic code are
removed from the ESMValTool main branch by the release manager and thus will not be included in future releases.
Only the scientific documentation of the recipe (and diagnostics) will be kept in the user and developer guide with an
additional note and link to the last release in which the recipe was still fully functional.
