How to add a recipe to the |RTW|
===================================================

.. include:: common.txt

**Please note**: Before you follow these steps to add your recipe, you must be
able to successfully run the recipe with the latest version of ESMValTool on
the compute server you use at your site, as detailed by the ``platform`` option
in the ``[[COMPUTE]]`` section in the site-specific ``.cylc`` file in the
``esmvaltool/utils/recipe_test_workflow/recipe_test_workflow/site/`` directory.

#. Obtain the duration and memory usage of the recipe from the messages printed
   to screen, or at the end of the ``run/main_log.txt`` file in the output
   directory after running your recipe on the compute cluster you use at your
   site; these messages will look something like::

    YYYY-MM-DD HH:MM:SS:sss UTC [12345] INFO    Time for running the recipe was: 0:02:13.334742
    YYYY-MM-DD HH:MM:SS:sss UTC [12345] INFO    Maximum memory used (estimate): 2.4 GB
    [...]
    YYYY-MM-DD HH:MM:SS:sss UTC [12345] INFO    Run was successful

#. Add the recipe to the ``[task parameters]`` section in the
   ``esmvaltool/utils/recipe_test_workflow/recipe_test_workflow/flow.cylc``
   file. If the recipe takes less than 10 minutes to run then it should be
   added to the ``fast`` option. Recipes that take longer than ten minutes
   should be added to the ``medium`` option.

#. If the duration of the recipe is larger than the value specified by the
   ``execution time limit`` option in the ``[[COMPUTE]]`` section in the
   aforementioned site-specific ``.cylc`` file, and / or the memory usage of
   the recipe is larger than the value specified by the ``--mem`` option in the
   ``[[[directives]]]`` section in the ``[[COMPUTE]]`` section, add a section
   (in alphabetical order) to this file as shown below (round the duration to
   the nearest second).

   Supplementary to the paragraph above, we advise execution with a duration of
   less than ~1m45s should use an execution time limit of PT2M, and round up to
   the nearest integer for the memory usage.::

    [[process<fast=recipe_albedolandcover>]]
    # Actual: 0m31s, 2.5 GB on 2024-04-08.
    execution time limit = PT2M
    [[[directives]]]
        --mem = 3G


   **Please note**: When you add the your modified version of the above example
   to your site-specific ``.cylc`` file, the key *fast* should match the
   category you placed your recipe within the ``[task parameters]`` section
   of the flow.cylc

#. Stop any running ``recipe_test_workflow`` workflows::

    cylc stop recipe_test_workflow

#. Run the recipe test workflow (it is expected that the ``compare`` task will
   fail).

#. Update the KGO:

        #. Locate the workflow run folder of the workflow you just completed.

        #. Recursively copy the recipe output directory (i.e.
           ``recipe_<recipe>_<date>_<time>``) from the
           ``${HOME}/cylc-run/recipe_test_workflow/run1/share/cycle/<cycle>``
           directory to your site-specific KGO directory, as detailed by the
           ``KGO_ROOT_PATH`` option in the site-specific ``.conf`` file in the
           ``esmvaltool/utils/recipe_test_workflow/recipe_test_workflow/opt/``
           directory::

            cp -r <directory_of_recipe_output_cycle_folder> <KGO_rootpath_folder>

        #. Change directory to the rootpath KGO directory::

            cd <KGO_rootpath_directory>

        #. Allow write permissions for all users on the directory and it's
           subdirectories of the recipe you've added to the KGO folder.::

            chmod a+w <the_directory_of_the_recipe_you_have_copied_into_the_KGO_folder>

#. Stop any running ``recipe_test_workflow`` workflows::

    cylc stop recipe_test_workflow

#. Run the RTW again. The workflow should now succeed.::

    cylc vip -O <site>

#. Add the recipe to the documentation; add a link to the recipe to the list of
   "Currently tested recipes" in ``doc/source/tested_recipes.rst`` in
   alphabetical order.

#. Commit and push your changes.
