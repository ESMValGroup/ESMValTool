How to add a recipe to the |RTW|
================================

.. include:: common.txt

Overview
--------

To add a recipe to the |RTW| you will:

* Run the recipe at your site
* Note the actual duration and memory usage
* Edit your site's recipe file
* Create the recipe's KGOs
* Request a review

The recipe will then run at your site whenever the |RTW| is run.

Preparation
-----------

#. Open a `new ESMValTool issue`_ on GitHub.
   Assign yourself to the issue
   and add the ``Recipe Test Workflow (RTW)`` label.
   `ESMValTool issue #3663`_ provides an example.

#. Create a branch.

#. Run the recipe:

   * with the latest version of ESMValTool
   * on the compute server you use at your site

   .. hint::
      Your compute server is defined in the ``site/<site>/runtime.cylc`` file
      as follows::

         [[COMPUTE]]
            platform = <your compute server>

#. Obtain the actual duration and memory usage of the recipe.
   This can be found either in the message printed to screen,
   or at the end of the ``run/main_log.txt`` file
   in the recipe output directory.
   The relevant lines will look something like::

      YYYY-MM-DD HH:MM:SS:sss UTC [12345] INFO    Time for running the recipe was: 0:02:13.334742
      YYYY-MM-DD HH:MM:SS:sss UTC [12345] INFO    Maximum memory used (estimate): 2.4 GB
      [...]
      YYYY-MM-DD HH:MM:SS:sss UTC [12345] INFO    Run was successful

Adding the recipe
-----------------

#. Add the recipe in alphabetical order to either ``FAST_RECIPES`` or
   ``MEDIUM_RECIPES`` in the ``site/<site>/recipes.jinja`` file.
   It should look something like::

      {
         'recipe_path': 'recipe_a_fast_recipe',
         'actual': '2m13s, 2.4 GB on YYYY-MM-DD',
         'max_time': 'PT3M',
         'max_memory': '3G',
      }

   .. important::
      Add the recipe to ``FAST_RECIPES`` if it takes *less* than 10 mins to
      run at your site. Add the recipe to ``MEDIUM_RECIPES`` if it takes *more*
      than 10 mins.

   .. hint::
      The :ref:`site_recipes_file` file provides more information.

   .. hint::
      Set the ``max_time`` to 10-20% more than the actual duration. For actual
      durations of up to ``1m45s``, set ``max_time`` to ``PT2M`` (2 minutes).

   .. hint::
      Try not to regularly waste more than 500 MiB in memory usage. Typically,
      rounding the actual memory usage up to the nearest integer is acceptable.

Create the |KGOs|
-----------------

#. Run the |RTW|, as detailed in the :ref:`quick_start_guide`; it is expected
   that the ``compare`` task will fail.

   .. important::
      The ``compare`` task fails because the |KGOs| for the recipe do not yet
      exist. This run of the |RTW| will generate the outputs that will be
      used as |KGOs|.

#. Recursively copy the recipe output directory
   ``recipe_<recipe>_<date>_<time>/`` to your site-specific |KGO| directory::

      cp -r ${HOME}/cylc-run/recipe_test_workflow/runN/share/cycle/<cycle>/recipe_<recipe>_<date>_<time> <KGO_ROOT_PATH>

   .. hint::
      ``<cycle>`` will look something like: ``20250101T0900Z``.
      The recipe output directory
      will look something like: ``recipe_python_20250101_090000``

   .. hint::
      Find your site-specific ``KGO_ROOT_PATH``
      in the ``opt/rose-suite-<site>.conf`` file.

   .. note::
      Cylc is typically configured such that
      ``${HOME}/cylc-run/recipe_test_workflow/runN/share`` is a symbolic link
      to a share directory located on a scratch disk; the recipe output
      directory is not stored in ``${HOME}``.

#. Enable write permissions for all users on the recipe output directory in
   your site-specific |KGO| directory::

      chmod -R a+w <KGO_ROOT_PATH>/recipe_<recipe>_<date>_<time>

#. Stop any running ``recipe_test_workflow`` workflows::

      cylc stop recipe_test_workflow/*

#. Run the |RTW| again, as detailed in the :ref:`quick_start_guide`; the
   ``compare`` task should now succeed.

   .. important::
      Unless you want the |RTW| to continue running every night, repeat the
      previous step to stop the workflow.

Request a review
----------------

#. Commit and push your changes.

#. Create a PR.
   Assign yourself to the PR,
   and add the ``Recipe Test Workflow (RTW)`` label to the PR.
   `ESMValTool PR #3664`_ provides an example.

   .. note::
      Reviewers will automatically be assigned to your PR.

Congratulations!
----------------

The recipe will now run at your site whenever the RTW is run.
