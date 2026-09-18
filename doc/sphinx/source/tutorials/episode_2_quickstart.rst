.. _episode_2_quickstart:

Quickstart guide
================

.. admonition:: Overview
   :class: note

   .. grid:: 2
      :gutter: 1
      :margin: 3 3 0 5

      .. grid-item-card:: Timings
         :columns: 12

         * Teaching: 2 min
         * Exercises: 8 min

      .. grid-item-card:: Questions

         * What is the purpose of the quickstart guide?
         * How do I load and check the ESMValTool environment?
         * How do I configure ESMValTool?
         * How do I run a recipe?

      .. grid-item-card:: Learning outcomes

         * Understand the purpose of the quickstart guide
         * Load and check the ESMValTool environment
         * Configure ESMValTool
         * Run a recipe

      .. grid-item-card:: Compatibility
         :columns: 12

         ESMValTool v2.15.0


What is the purpose of the quickstart guide?
--------------------------------------------

The purpose of the quickstart guide is to enable a user of ESMValTool to run
ESMValTool as quickly as possible by making the minimum number of changes.

.. admonition:: Purpose of the quickstart guide
   :class: discussion

   The purpose of this quickstart guide is to get you running with ESMValTool
   as quickly as possible without requiring a full tutorial walkthrough.


How do I load and check the ESMValTool environment?
---------------------------------------------------

For this quickstart guide, it is assumed that ESMValTool has already been
installed at the site where it will be run. If this is not the case, see the
:doc:`Installation <episode_3_installation>` episode in this tutorial.

Load the ESMValTool environment by following the instructions in the
:ref:`ESMValTool installation and environment activation guide <install_on_hpc>`.
This will typically involving loading a module like so:

.. code-block:: bash

   module load esmvaltool

Check the ESMValTool environment by accessing the help for ESMValTool:

.. code-block:: bash

   esmvaltool --help

.. admonition:: Check the ESMValTool environment
   :class: admonition-todo

   Run the command above to confirm that your ESMValTool environment is loaded
   correctly and that the command-line interface is available.


How do I configure ESMValTool?
------------------------------

Create the ESMValTool user configuration file. By default, this file is written
to ``~/.config/esmvaltool/config-user.yml``:

.. code-block:: bash

   esmvaltool config copy defaults/config-user.yml

Edit the ESMValTool user configuration file using your favourite text editor to
uncomment the lines relating to the site where ESMValTool will be run.

For more details about the ESMValTool user configuration file, see the
:doc:`Configuration <episode_4_configuration>` episode in this tutorial.

.. admonition:: Configure the tool
   :class: admonition-todo

   Create the user configuration file and adjust the settings for the machine on
   which you are running ESMValTool.


How do I run a recipe?
----------------------

Run the example Python recipe:

.. code-block:: bash

   esmvaltool run examples/recipe_python.yml

Wait for the recipe to complete. If the recipe completes successfully, the last
line printed to the screen at the end of the log will look something like:

.. code-block:: bash

   YYYY-MM-DD HH:mm:SS, NNN UTC [NNNNN] INFO    Run was successful

View the output of the recipe by opening the HTML file produced by ESMValTool.
The location of this file is printed to the screen near the end of the log:

.. code-block:: bash

   YYYY-MM-DD HH:mm:SS, NNN UTC [NNNNN] INFO    Wrote recipe output to:
   file:///$HOME/esmvaltool_output/recipe_python_``<date>_<time>``/index.html

For more details about running recipes, see the
:doc:`episode_5_recipe` episode in this tutorial.

.. admonition:: Run an example recipe
   :class: admonition-todo

   Try running the example recipe above and confirm that the run completes
   successfully before moving on to other recipes.


Conclusion
----------

This quickstart guide introduces the minimum steps needed to get ESMValTool
running on your system. Once the environment is loaded, configured, and a recipe
has been run successfully, you are ready to explore the wider tutorial and more
advanced workflows.

.. admonition:: Key points
   :class: important

   * The purpose of the quickstart guide is to enable a user to run ESMValTool
     as quickly as possible without having to go through the whole tutorial.
   * Use ``module load`` to load the ESMValTool environment.
   * Use ``esmvaltool --help`` to check that the environment is working.
   * Use ``esmvaltool config copy defaults/config-user.yml`` to create the user
     configuration file.
   * Use ``esmvaltool run ``<recipe>``.yml`` to run a recipe.
