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

#. Add the recipe to the ``[task parameters]`` section in the
   ``esmvaltool/utils/recipe_test_workflow/flow.cylc`` file.

   .. hint::
      If the recipe takes less than 10 minutes to run then it should be added
      to the ``fast`` option. Recipes that take longer than ten minutes should
      be added to the ``medium`` option.

   .. hint::
      The line added should follow the format of ``recipe_new_recipe, \``,
      unless the line is the last one in the list, in which case the line added
      should follow the format of ``recipe_new_recipe``.

#. If the duration of the recipe is larger than the value specified by the
   ``execution time limit`` option in the ``[[COMPUTE]]`` section in the
   aforementioned site-specific ``.cylc`` file, and / or the memory usage of
   the recipe is larger than the value specified by the ``--mem`` option in the
   ``[[[directives]]]`` section in the ``[[COMPUTE]]`` section, add a section
   (in alphabetical order) to this file as shown below (round the duration to
   the nearest second)::

     [[process<fast=recipe_albedolandcover>]]
     # Actual: 0m31s, 2.5 GB on 2024-04-08.
     execution time limit = PT2M
     [[[directives]]]
         --mem = 3G

   .. hint::
      The ``fast`` key in the example task definition above
      (``[[process<fast=recipe_albedolandcover>]]``) should match name of the
      option the recipe was added to in the ``[task parameters]`` section in
      the ``esmvaltool/utils/recipe_test_workflow/flow.cylc`` file

   .. hint::
      Set the ``execution time limit`` to 10-20% more than the actual duration.
      For actual durations of up to ``1m45s``, set the ``execution time limit``
      to ``PT2M`` (2 minutes).

   .. hint::
      Try not to regularly waste more than 500 MiB in memory usage. Typically,
      rounding the actual memory usage up to the nearest integer is acceptable.

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

#. Add the recipe to the documentation:

   * Add a URL for the recipe in ``doc/sphinx/source/utils/RTW/common.txt``
     under the ``.. Links`` section in alphabetical order (follow the format
     ``.. _<name>: <URL>``).

   * Add the recipe to the list of :ref:`currently_tested_recipes` in
     alphabetical order (follow the format ``* `<name>`_``).

#. Commit and push your changes, create a PR, assign yourself to the PR, and
   add the ``Recipe Test Workflow (RTW)`` label to the PR, see
   `ESMValTool PR #3664`_ for an example.
