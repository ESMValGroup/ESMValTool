.. include:: ../common.txt

Adding a Recipe to the RTW
==========================

Can your recipe be added to the RTW?
------------------------------------
The |RTW| will run a recipe and compare the output to a |KGO|.
The compare task will fail if the output doesn not exactly match the |KGO|,
therefore, if the recipe output is not reproducible it cannot be compared to a |KGO|, and the compare task will always fail.

It is recommended that recipes are added to the |RTW| as this will make the process of a release easier,
and should a breaking change be introduced ESMValTool, the author will be notified when the recipe on the |RTW| fails.

Prerequisites
-------------

To add a recipe to the |RTW|, the first thing that must be done is run the recipe.

Take note of the resources required to run the recipe and the time it takes.

Inspect the results of your recipe and verify they are as expected.

Add the verified output of the recipe to the KGO_ROOT_PATH of each site. These can be found in in the
`recipe_test_workflow/opt/rose-suite-<site>.conf` files. Please ensure the |KGO| for the recipe to be
added is placed in the root path for every site.

If for any reason you are unable to access the KGO_ROOT_PATH directory, please contact the maintainer of the |RTW|.

Add your recipe to the workflow
-------------------------------
* Checkout the |RTW|::

    git clone git@github.com:ESMValGroup/ESMValTool.git

* Make a branch of the |RTW|::

    git checkout recipe_test_workflow_prototype
    git checkout -b add_recipe_branch

* Every recipe in the |RTW| is sorted into one of three groups; fast, medium, and slow, depending on what resources
  they require and how long they take to run.

  * fast - Anything under 2 minutes to run and under 2G memory required.
  * medium - Anything between 2 and 6 minutes to run and under 2G memory required.
  * slow - Anything above this.


* Open the flow.cylc file and edit the task parameters of the group the added recipe belongs to::

    [task parameters]
       fast = radiation_budget, <your_recipe_name>
       medium = ensclus, heatwaves_coldwaves

  The name of the task parameter must match the file name of the recipe to be added

* Run the branch of |RTW| on all sites and ensure the added recipe passes.

* Make a pull request to merge the branch with the added recipe with `recipe_test_workflow_prototype`
