How to add a recipe to the |RTW|
================================

.. include:: common.txt

.. note::
   Before you follow these steps to add your recipe, you must be able to
   successfully run the recipe with the latest version of ESMValTool on the
   compute server you use at your site, as detailed by the ``platform`` option
   in the ``[[COMPUTE]]`` section in the site-specific ``.cylc`` file in the
   ``esmvaltool/utils/recipe_test_workflow/site/`` directory.

#. Open a `new ESMValTool issue`_ on GitHub, assign yourself to the issue, and
   add the ``Recipe Test Workflow (RTW)`` label to the issue, see
   `ESMValTool issue #3663`_ for an example.

#. Create a branch.

#. Obtain the duration and memory usage of the recipe from the messages printed
   to screen, or at the end of the ``run/main_log.txt`` file in the recipe
   output directory after running your recipe on the compute cluster you use at
   your site; these messages will look something like::

     YYYY-MM-DD HH:MM:SS:sss UTC [12345] INFO    Time for running the recipe was: 0:02:13.334742
     YYYY-MM-DD HH:MM:SS:sss UTC [12345] INFO    Maximum memory used (estimate): 2.4 GB
     [...]
     YYYY-MM-DD HH:MM:SS:sss UTC [12345] INFO    Run was successful

.. _run_a_recipe_in_the_rtw:

Run a Recipe in the |RTW|
-------------------------

#. The recipe should now be added to your ``<site>-recipes.cylc`` file. (Where
   ``<site>`` is your site). Find this in the
   ``esmvaltool/utils/recipe_test_workflow/site/`` directory.

#. ``<site>-recipes.cylc`` contains **two lists of dictionaries**. The lists
   are ``FAST_RECIPES`` and ``MEDIUM_RECIPES``.

   .. hint::
      ``FAST_RECIPES`` take *less* than 10 minutes to run at your site.
      ``MEDIUM_RECIPES`` take *more* than 10 minutes.

#. Add the recipe to one of lists as as a dictionary of ``key, value`` pairs.
   E.g.::

      {
         'recipe_path': 'recipe_a_fast_recipe',
         'actual': '1m30s, 1.5 GB on 2025-01-01',
         'max_time': 'PT2M',
         'max_memory': '2G',
      }

   Add the following information for each key:

   .. list-table::
      :widths: 25 75
      :header-rows: 1

      * - Key
        - Value
      * - ``recipe_path``
        - The path to the recipe. Recipe paths are specified relative to
          ``esmvaltool/recipes``. For recipes in subdirectories, ``--`` stands
          for ``/`` since the latter is an illegal char.
      * - ``actual``
        - A note of  the recipe's actual resource usage. From the successful run
          of the recipe on the compute server at your site.
      * - ``max_time``
        - Set the maximum amount of time the recipe has to complete. Use the
          ISO8601 duration format (see `Cylc ISO8601 Durations`_ for more).
      * - ``max_memory``
        - Set the memory to allocate to running the recipe. Default units are
          megabytes. Different units can be specified using the
          suffix [K|M|G|T] (see `Slurm sbatch --mem`_ for more).

   .. hint::
      Set the ``max_time`` to 10-20% more than the actual duration. For actual
      durations of up to ``1m45s``, set ``max_time`` to ``PT2M`` (2 minutes).

   .. hint::
      Try not to regularly waste more than 500 MiB in memory usage. Typically,
      rounding the actual memory usage up to the nearest integer is acceptable.

#. Once the recipe dictionary is added to either ``FAST_RECIPES`` or
   ``MEDIUM_RECIPES``, the recipe will run as part of the |RTW| at your site.
   Next, the recipe will need |KGOs| to run successfully in the RTW.

   .. _jinja2_templating_language:

   .. note::
      The ``<site>-recipes.cylc`` file is actually written in the `Jinja2`_
      templating language. Jinja2 gives |Cylc| many powerful features (see
      `Cylc Jinja2`_). This is beyond the scope of this guide. Follow the links
      for more information.

Update the KGOs
---------------

#. Stop any running ``recipe_test_workflow`` workflows::

    cylc stop recipe_test_workflow/*

#. Run the |RTW|, as detailed in the :ref:`quick_start_guide`; it is expected
   that the ``compare`` task will fail.

#. Update the Known Good Outputs (|KGOs|):

   * Recursively copy the recipe output directory (i.e.
     ``recipe_<recipe>_<date>_<time>/``) from the
     ``${HOME}/cylc-run/recipe_test_workflow/runN/share/cycle/<cycle>/``
     directory to your site-specific |KGO| directory, as detailed by the
     ``KGO_ROOT_PATH`` option in the site-specific ``.conf`` file in the
     ``esmvaltool/utils/recipe_test_workflow/opt/`` directory::

       cp -r ${HOME}/cylc-run/recipe_test_workflow/runN/share/cycle/<cycle>/recipe_<recipe>_<date>_<time> <KGO_ROOT_PATH>

     .. note::
        Cylc is typically configured such that
        ``${HOME}/cylc-run/recipe_test_workflow/runN/share`` is a symbolic link
        to a share directory located on a scratch disk; the recipe output
        directory is not stored in ``${HOME}``.

   * Enable write permissions for all users on the recipe output directory in
     your site-specific |KGO| directory::

       chmod -R a+w <KGO_ROOT_PATH>/recipe_<recipe>_<date>_<time>

#. Stop any running ``recipe_test_workflow`` workflows::

     cylc stop recipe_test_workflow/*

#. Run the |RTW| again, as detailed in the :ref:`quick_start_guide`; the
   ``compare`` task should now succeed.

Add the recipe to the documentation
------------------------------------

#. Add the recipe to the documentation:

   * Add a URL for the recipe in ``doc/sphinx/source/utils/RTW/common.txt``
     under the ``.. Links`` section in alphabetical order (follow the format
     ``.. _<name>: <URL>``).

   * Add the recipe to the list of :ref:`currently_tested_recipes` in
     alphabetical order (follow the format ``* `<name>`_``).

#. Commit and push your changes, create a PR, assign yourself to the PR, and
   add the ``Recipe Test Workflow (RTW)`` label to the PR, see
   `ESMValTool PR #3664`_ for an example.
