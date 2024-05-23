Adding a recipe to the |RTW|
============================

.. include:: ../common.txt

* Run ESMValTool locally with your recipe, make sure to take note of the memory
  and time expenditure provided by the terminal output::

    esmvaltool run <your_recipe_name.yml>

* Run the same recipe on JASMIN, again taking note of the memory and time
  expenditure printed to the terminal.

* Stop any running workflows::

    cylc stop "*"

* Add the recipe to the ``[task parameters]`` section of the |RTW| ``flow.cylc``
  workflow file (make sure to place it within the "fast" or "medium" category
  depending on how long it took to run):

* Run the workflow (it should fail)::

    cylc vip -O <your_site_name>

* Locate the workflow run folder of the workflow you just completed.

* Copy the output files from the "cycle" folder (run/share/bin/cycle) of the
  workflow you just ran, to the KGO rootpath folder (/data/users/esmval/KGO)::

    cp -r <directory_of_recipe_output_cycle_folder> <KGO_rootpath_folder>

* Change directory to the rootpath KGO directory::

    cd <KGO_rootpath_directory>

* Allow write permissions for all users on the directory and it's
  subdirectories of the recipe you've added to the KGO folder::

    chmod a+w <the_directory_of_the_recipe_you_have_copied_into_the_KGO_folder>

* Run the RTW again::

    cylc vip -O <your_site_name>

* If the workflow succeeds then your recipe has successfully been added to the
  workflow. You can now commit your changes and push them onto GitHub.
