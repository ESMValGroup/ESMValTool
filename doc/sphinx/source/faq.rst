.. _faq:

Frequently Asked Questions
**************************

Is there a mailing list?
========================

Yes, you can subscribe to the ESMValTool user mailing list and join the discussion on general topics (installation, configuration, etc). See :ref:`mailing-list`.

What is YAML?
=============

While ``.yaml`` or ``.yml`` is a relatively common format, users may not have
encountered this language before. The key information about this format is:

- yaml is a human friendly markup language;
- yaml is commonly used for configuration files (gradually replacing the
  venerable ``.ini``);
- the syntax is relatively straightforward;
- indentation matters a lot (like ``Python``)!
- yaml is case sensitive;

More information can be found in the `yaml tutorial
<https://learnxinyminutes.com/docs/yaml/>`_ and `yaml quick reference card
<https://yaml.org/refcard.html>`_. ESMValTool uses the `yamllint
<http://www.yamllint.com>`_ linter tool to check recipe syntax.


.. _rerunning:

Re-running diagnostics
======================

If a diagnostic fails, you will get the message

.. code:: bash

   INFO    To re-run this diagnostic script, run:

If you run the command in the stdout you will be able to re-run the
diagnostic without having to re-run the whole preprocessor. If you add the ``-f``
argument (available only for Python diagnostics, check your options with ``--help``)
that will force an overwrite, and it will delete not just the failed diagnostic,
but the contents of its ``work_dir`` and ``plot_dir`` directories - this is useful when needing to
redo the whole work. Adding ``-i`` or ``--ignore-existing`` will not delete any existing files,
and it can be used to skip work that was already done successfully, provided
that the diagnostic script supports this.


Enter interactive mode with iPython
===================================

Sometimes it is useful to enter an interactive session to have a look what's going on.
Insert a single line in the code where you want to enter IPython:
``import IPython; IPython.embed()``

This is a useful functionality because it allows the user to `fix` things on-the-fly and after
quitting the Ipython console, code execution continues as per normal.


Using multiple configuration directories
========================================

By default, ESMValTool will read YAML configuration files from the user
configuration directory ``~/.config/esmvaltool``, which can be changed with the
``ESMVALTOOL_CONFIG_DIR`` environment variable.
If required, users can specify the command line option ``--config_dir`` to
select another configuration directory, which is read **in addition** to the
user configuration directory
See the section on configuration :ref:`config_yaml_files` for details on this.


Create a symbolic link to the latest output directory
=====================================================

When running multiple times the same recipe, the tool creates separate output directories
sorted by the time tag that they were created at; sometimes, when running quite a few times,
it is not straightforward to detect which one is the `latest` output directory, so a symbolic
link attached to it would make things more clear e.g.:

.. code:: bash

   recipe_example_20190905_163431
   recipe_example_20190905_163519
   recipe_example_latest -> recipe_example_20190905_163519


You can do that by running the tool using the latest output as basis and creating
a symbolic link to it so it gets picked up at every re-run iteration:

.. code:: bash

   esmvaltool run recipe_example.yml; \
   ln -sfT $(ls -1d ~/esmvaltool_output/recipe_example_* | tail -1) ~/esmvaltool_output/recipe_example_latest


.. uncomment when feature plopped in main
.. # Running a dry run
.. =================

.. You can run in dry-run mode with

.. .. code:: bash

..   esmvaltool run recipe_xxx.yml --dry-run


.. This mode activated will run through the data finding and CMOR checks and fixes
.. and will highlight on screen and in `run/main_log.txt` every time certain data is
.. missing or there are issues with the CMOR checks; note that no data is written
.. to disk and no diagnostics are run; you don't have to modify your recipe in any
.. way to have this mode run. The information provided will help you obtain any data
.. that is missing and/or create fixes for the datasets and variables that failed the
.. CMOR checks and could not be fixed on the fly.


Can ESMValTool plot arbitrary model output?
===========================================

:ref:`recipe_model_evaluation` provides a set of recipes that can be used for a
basic climate model evaluation with observational data.
This is especially useful to get an overview of the general performance of a
simulation.

Furthermore, recipe :ref:`recipe_monitor` allows for the plotting of any
preprocessed model.
The plotting parameters are set through a yaml configuration file, and the type
of plots to be generated are determined in the recipe.

Similarly, recipe :ref:`recipes_seaborn_diag` and the corresponding diagnostic
:ref:`seaborn_diag.py <api.esmvaltool.diag_scripts.seaborn_diag>` provide a
high-level interface to the `Seaborn <https://seaborn.pydata.org>`__ package
which can also be used to create a large variety of different plots.

See also :ref:`general_purpose_diags`.
