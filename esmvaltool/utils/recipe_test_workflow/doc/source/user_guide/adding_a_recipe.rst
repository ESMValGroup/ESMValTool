Adding a recipe to the Recipe Test Workflow (|RTW|)
===================================================

.. include:: ../common.txt

**Please note**: Before you follow these steps to add it, you must be
able to successfully run it with the latest version of ESMValTool on the
compute server you use at your site, as detailed by the ``platform`` option in
the ``[[COMPUTE]]`` section in the site-specific ``.cylc`` file in the
``esmvaltool/utils/recipe_test_workflow/recipe_test_workflow/site/`` directory.

#. Stop any running recipe_test_workflow workflows::

    cylc stop "a_running_recipe_test_workflow"

#. Open it in your preferred code editor.

#. Locate the ``[[COMPUTE]]`` section, it should look something like this::

    [[COMPUTE]]
    platform = <your_platform_here>
    execution time limit = PT2M
    [[[directives]]]
    --wckey = RTW
    --ntasks = {{ MAX_PARALLEL_TASKS }}
    --mem = 2G

#. Run your recipe with ESMValTool on your compute server for your site.

#. If the recipe takes less than 10 minutes to run then it should be added as a
   "fast" recipe in the ``flow.cylc`` file within the ``[task parameters]``
   section.If it takes longer than ten minutes it should be included in "medium
   ".

#. If either of the memory readings from your run are larger than the
   values specified in the ``[[COMPUTE]]`` section, you need to add your recipe as
   another `process` similar to::

    [[process<fast=recipe_albedolandcover>]]
    # Actual: 0m31s, 2.5 GB on 2024-04-08.
    execution time limit = PT2M
    [[[directives]]]
        --mem = 3G

#. Run the recipe test workflow.

#. The ``compare`` task will fail. These next steps should fix it,
   and even if the ``compare`` task has passed, they must still be completed.

#. Locate the workflow run folder of the workflow you just completed.

#. Copy the output files from the "cycle" folder (``run/share/cycle``) of the
   workflow you just ran, to your site specific KGO rootpath folder
   (this folder should be set as the value for the variable ``KGO_ROOT_PATH=``
   in your ``rose-suite<your_site>.conf`` file found in the ``/data/users/esmva
   l/KGO`` directory)::

    cp -r <directory_of_recipe_output_cycle_folder> <KGO_rootpath_folder>

#. Change directory to the rootpath KGO directory::

    cd <KGO_rootpath_directory>

#. Allow write permissions for all users on the directory and it's
   subdirectories of the recipe you've added to the KGO folder::

    chmod a+w <the_directory_of_the_recipe_you_have_copied_into_the_KGO_folder>

#. Stop any running recipe_test_workflow workflows::

    cylc stop "a_running_recipe_test_workflow"

#. Run the RTW again::

    cylc vip -O <your_site_name>

#. The workflow should now succeed.

#. Add your recipe to the list of "Currently tested recipes" in ``tested_recipes.rst``
   (this file is located within ``recipe_test_workflow/doc/source``) in
   it's correct position (the list is ordered alphabetically), this will add
   it to the documentation.

#. If the workflow succeeds then your recipe has successfully been added to the
   workflow. You can now commit your changes and push them onto GitHub.
