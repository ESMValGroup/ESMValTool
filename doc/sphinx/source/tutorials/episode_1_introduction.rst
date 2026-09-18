.. _epsiode_1:

Introduction
=====================

.. admonition:: Overview
   :class: note

   .. grid:: 2
      :gutter: 1
      :margin: 3 3 0 5

      .. grid-item-card:: Timings
         :columns: 12

         * Teaching: 5 min
         * Exercises: 10 min

      .. grid-item-card:: Questions

         * What is ESMValTool?
         * Who are the people behind ESMValTool?

      .. grid-item-card:: Learning outcomes

         * Familiarization with ESMValTool
         * Synchronize expectations

      .. grid-item-card:: Compatibility
         :columns: 12

         ESMValTool v2.15.0


What is ESMValTool?
-------------------

This tutorial is a first introduction to ESMValTool. Before diving into the
technical steps, let’s talk about what ESMValTool is all about.

.. admonition:: What is ESMValTool?
   :class: admonition-todo

   What do you already know about or expect from ESMValTool?

   .. dropdown:: ESMValTool is...
      :color: secondary
      :icon: eye

      ESMValTool is many things, but in this tutorial we will focus on the
      following traits:

      * **A tool to analyse climate data**
      * **A collection of diagnostics for reproducible climate science**
      * **A community effort**

A tool to analyse climate data
------------------------------

ESMValTool takes care of finding, opening, checking, fixing, concatenating,
and preprocessing CMIP data and several other supported datasets.

The central component of ESMValTool that we will see in this tutorial is the
**recipe**. Any ESMValTool recipe is basically a set of instructions to
reproduce a certain result. The basic structure of a recipe is as follows:

* **Documentation** with relevant (citation) information
* **Datasets** that should be analysed
* **Preprocessor** steps that must be applied
* **Diagnostic** scripts performing more specific evaluation steps

An example recipe could look like this:

.. code-block:: bash

   documentation:
   title: This is an example recipe.
   description: Example recipe
   authors:
      - lastname_firstname

   datasets:
   - {dataset: UKESM1-0-LL, project: CMIP6, exp: historical, mip: Amon,
      ensemble: r1i1p1f2, start_year: 1960, end_year: 2005}

   preprocessors:
   global_mean:
      area_statistics:
         operator: mean

   diagnostics:
   hockeystick_plot:
      description: plot of global mean temperature change
      variables:
         temperature:
         short_name: tas
         preprocessor: global_mean
      scripts: hockeystick.py

.. admonition:: Understanding the different sections of the recipe
   :class: admonition-todo

   Try to figure out the meaning of the different dataset keys. Hint: they can
   be found in the documentation of ESMValTool.

   .. dropdown:: Solution
      :color: secondary
      :icon: eye

      The keys are explained in the ESMValCore documentation, in **The
      recipe format** section under
      :ref:`Overview <esmvalcore:recipe_overview>`.


A collection of diagnostics for reproducible climate science
------------------------------------------------------------

More than a tool, ESMValTool is a collection of publicly available recipes and
diagnostic scripts. This makes it possible to easily reproduce important
results.

.. admonition:: Explore the available recipes
   :class: admonition-todo

   Go to the :doc:`ESMValTool documentation </index>`
   and explore the **Recipes** section in the sidebar. Which recipe(s) would
   you like to try?

A community effort
------------------

ESMValTool is built and maintained by an active community of scientists and
software engineers. It is an open source project to which anyone can
contribute. Many of the interactions take place on GitHub. Here, we briefly
introduce you to some of the most important pages.

.. admonition:: Meet the ESMValGroup
   :class: admonition-todo

   Go to https://github.com/ESMValGroup. This is the
   GitHub page of our 'organization'. Have a look around. How many
   collaborators are there? Do you know any of them?

   Near the top of the page there are 3 pinned repositories: ESMValTool,
   ESMValCore and Community. Visit each of the repositories. How many people
   have contributed to each of them? Can you also find out how many people
   have contributed to this tutorial?

.. admonition:: Issues and pull requests
   :class: admonition-todo

   Go back to the repository pages of `ESMValTool <https://github.com/
   ESMValGroup/ESMValTool>`_ or `ESMValCore <https://github.com/ESMValGroup/
   ESMValCore>`_. There are tabs for ‘issues’ and ‘pull requests’. You can use
   the labels to navigate them a bit more. How many open issues are about
   enhancements of ESMValTool? And how many bugs have been fixed in
   ESMValCore? There is also an ‘insights’ tab, where you can see a summary
   of recent activity. How many issues have been opened and closed in the
   past month?

Conclusion
----------

This concludes the introduction of the tutorial. You now have a basic
knowledge of ESMValTool and its community. The following episodes will walk
you through the installation, configuration and running your first recipes.

.. admonition:: Key points
   :class: important

   * ESMValTool provides a reliable interface to analyse and evaluate climate
     data
   * A large collection of recipes and diagnostic scripts is already available
   * ESMValTool is built and maintained by an active community of scientists
     and developers

.. _`PyData Theme documentation: Admonitions`: https://pydata-sphinx-theme.readthedocs.io/en/stable/examples/kitchen-sink/admonitions.html
.. _`Font Awesome`: https://fontawesome.com/search?ic=free-collection
