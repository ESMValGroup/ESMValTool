Adding a recipe to the Recipe Test Workflow (|RTW|)
===================================================

.. include:: ../common.txt

**Please note:** Before you follow these steps to add your recipe, you must be
able to successfully run it with ESMValtool without error.

* Run ESMValTool locally with your recipe, make sure to take note of the memory
  and time expenditure provided by the terminal output::

    esmvaltool run <your_recipe_name.yml>

* Run the same recipe on JASMIN, again taking note of the memory and time
  expenditure printed to the terminal.

* Stop any running recipe_test_workflow workflows::

    cylc stop "a_running_recipe_test_workflow"

* Add the recipe to the ``[task parameters]`` section of the |RTW| ``flow.cylc``
  workflow file, make sure your recipe name fits the formatting of ``recipe_new_
  recipe \\``.

* Run the workflow. The process task for the new recipe should succeed, but its
  compare task should fail with an error that the reference data (KGO) does
  not exist::

    cylc vip -O <your_site_name>

* Locate the workflow run folder of the workflow you just completed.

* Copy the output files from the "cycle" folder (run/share/bin/cycle) of the
  workflow you just ran, to your site specific KGO rootpath folder
  (this folder should be set as the value for the variable "KGO_ROOT_PATH="
  in your rose-suite<your_site>.conf file)::

    cp -r <directory_of_recipe_output_cycle_folder> <KGO_rootpath_folder>

* Change directory to the rootpath KGO directory::

    cd <KGO_rootpath_directory>

* Allow write permissions for all users on the directory and it's
  subdirectories of the recipe you've added to the KGO folder::

    chmod a+w <the_directory_of_the_recipe_you_have_copied_into_the_KGO_folder>

* Run the RTW again::

    cylc vip -O <your_site_name>

* The workflow should now succeed.

* Take note of how long the run took  to complete on cylc. This can be found in the `job.time` section of the task listed as `process_<your_recipe>`.

* Locate your local <site>.cylc file, found in (recipe_test_workflow/site).

* Open it in your preferred code editor.

* Locate the `COMPUTE` section, it should look something like this::

    [[COMPUTE]]
    platform = <your_platform_here>
    execution time limit = PT2M
    [[[directives]]]
    --wckey = RTW
    --ntasks = {{ MAX_PARALLEL_TASKS }}
    --mem = 2G

* Compare the `execution time limit` and --mem (memory) units here with the
  readings you took locally. If your local readings do not exceed these, then
  you have successfully added the recipe to the workflow, and can now commit
  and push your changes.


* If either of the time/memory readings from your local run are larger than the
  values specified in the `COMPUTE` section, you need to add your recipe as
  another `process` similar to::

    [[process<fast=recipe_albedolandcover>]]
    # Actual: 0m31s, 2.5 GB on 2024-04-08.
    execution time limit = PT2M
    [[[directives]]]
        --mem = 3G

* Variable (fast, medium) must be consistent with flow.cylc.

* The commented "Actual" reading should be the time/memory
  reading you recorded from cylc review.

* Adjust the values for `execution_time` and `--mem` to be larger than the
  values you recorded from cylc review.

* Stop any running workflows.

* Run the recipe test workflow again.

* If the workflow succeeds then your recipe has successfully been added to the
  workflow. You can now commit your changes and push them onto GitHub.
