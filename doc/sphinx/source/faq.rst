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


Use multiple config-user.yml files
==================================

The user selects the configuration yaml file at run time. It's possible to
have several configurations files. For instance, it may be practical to have one
config file for debugging runs and another for production runs.

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

Recipe :ref:`recipe_monitor` allows for the plotting of any preprocessed model.
The plotting parameters are set through a yaml configuration file, and the
type of plots to be generated are determined in the recipe.

Debugging FAQ
=============

Something's going wrong in preprocessing, how can I find out where?
-------------------------------------------------------------------

Sometimes there is a problem in the preprocessing, often because of faulty data.

To locate the problem, the first step is to edit your `config-user.yml` to set

.. code-block:: yaml

   log_level: debug

Then run `esmvaltool` again and try to understand the error message. Don't
forget that you can open the error log in the run directory instead of working
only with the terminal output.  Usually, the error log ends in a python
traceback.  See for example `here <https://realpython.com/python-traceback/>`_
for an introduction to understanding tracebacks.

It can also be helpful to look in the error log directly above the traceback for
more information on what caused the problem.

If the error message is not enough to diagnose the problem, you can "spread out"
the effects of the different preprocessors. For that, edit again
`config-user.yml` to set

.. code-block:: yaml

   save_intermediary_cubes: true
   remove_preproc_dir: false

This way, you can inspect the output of every preprocessing step individually by
looking at the individual outputs in the `preproc` directory. Usually, the last
existing file can give you a clue on what went wrong.

.. warning::

   This can be costly both in storage space and computing time, so be careful to
   use as small a dataset as possible to diagnose your problem and to deactivate
   the setting when you are done debugging!

I get a `iris.exceptions.ConcatenateError`. What's the problem?
---------------------------------------------------------------

Broadly speaking, the problem is that data that was spread over different files
was supposed to be part of one dataset, but doesn't follow the same conventions.
This can be due to different coordinates, different attributes, or a number of
other differences. Unfortunately, it can be difficult to understand what exactly
is the difference between the data.

To diagnose the problem, first, figure out which netcdf files cause the
issue. If it happens in the preprocessing, refer to `this FAQ
<somethings-going-wrong-in-preprocessing-how-can-I-find-out-where>`_ for some
pointers.

Next, try to understand the differences between the two files. Perhaps the error
message already pointed to the area of difference, for example some attributes
or particular coordinates. In that case you can focus on those first.

How can I compare two (or more) netcdf files?
---------------------------------------------

There are three popular tools/methods to understand the differences between
netcdf files. The first two, `looking at the metadata of the individual files
<inspecting-the-netcdf-file-headers>`_ and `visualizing their data
<visualizing-the-data>`_, are almost always possible and can give you a quick
idea of what's going on. They are, however, often not suited to detailed
comparisons, for example when the actual data contains subtle problems or when
the problem is in the values of coordinates. In those cases, turn to `a more
in-depth look <comparing-data-with-iris>`_ at the data. This last method also
works better when you have a larger collection of files and are not sure where
exactly the problem lies.

Inspecting the netcdf file headers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can use the `ncdump` utility to look at the metadata of netcdf files. More
specifically, `ncdump -h my_file.nc` will output the so-called header of the
file that contains information about all the variables present in the file,
along with their units, etc. This includes the coordinate variables. It also
will give you all the global, ie file level, and local, ie variable level
attributes.  It will, however, not give you the coordinate values or the data
values itself.  It will also not help you with comparison, so you can only
compare the output for two files manually.

Nevertheless, when you want to get a quick overview, or when you have some idea
what to look for, or when the other methods don't work, this is a good starting
point.

Visualizing the data
~~~~~~~~~~~~~~~~~~~~

Often a quick look at the data can reveal the most glaring problems, like
missing or wrong masks, or wrong units.  There are a number of tools that can be
used for this, but one of the most basic and almost universally available ones
is `ncview`. This can also be installed with conda from conda-forge.

Just do a quick `ncview my_file.nc` to open the gui. Then select the variable
and use the navigation buttons to flip through the timesteps and, in the case of
3d data, the levels.

Again, this tool doesn't help you with comparison, but you can open two
instances to compare two files in side-by-side.

Comparing data with iris
~~~~~~~~~~~~~~~~~~~~~~~~

`Iris <https://scitools.org.uk/iris/docs/latest/>`_ is a powerful python
framework for working with `CF <http://cfconventions.org/>`_ compliant netcdf
files. It also includes some support for other formats like grib, and pp files.
Iris is the foundation for esmvaltool, but it can also be used on its own. Two
convenient ways to do that are `ipython <https://ipython.org/>`_ and `jupyter
notebooks <https://jupyter.org/>`_. This FAQ cannot cover any of these tools in
great detail, so please refer to their respective documentation for an overview
and in-depth questions. Instead, here we will give just a quick way to
understand differences in a set of netcdf files with iris.

Iris loads data with one of several functions whose names start with
`load`. These functions take a single string or a list of strings that are
interpreted as filenames. They may also contain directories and globs, like
`./my_data/tas_*.nc` to load all surface temperature data from the `my_data`
subdirectory of the current working directory. The iris function will return
either a single `Cube
<https://scitools.org.uk/iris/docs/latest/iris/iris/cube.html#iris.cube.Cube>`_,
or a list of cubes in the form of a `CubeList
<https://scitools.org.uk/iris/docs/latest/iris/iris/cube.html#iris.cube.CubeList>`_.

Often, we want to combine different parts of dataset to form the whole dataset
again. This is accomplished in iris using either a `merge` or a
`concatenate`. Refer to `the iris user guide
<https://scitools.org.uk/iris/docs/latest/userguide/merge_and_concat.html>`_ to
understand the difference.

In fact, it is this `concatenate` that has failed you when you are reading this
FAQ entry.

To find out what went wrong, let's assume we have boiled down the problem to the
two files `tas_1969.nc` and `tas_1970.nc`. Then one way to learn more about the
issue is the following.

.. code-block:: python

   import iris
   cube_list = iris.load('tas_19*.nc')
   cube_list.concatenate_cube()

This will likely produce the exact same error that we started with. While we do
want to do a concatenation, since both files contain several timesteps each, the
output of the corresponding merge function is often more illuminating. So let's
give that a try.

.. code-block:: python

   cube_list.merge_cube()

Hopefully, that gives you a better idea of the problem.

If you are dealing with more than two files it can be difficult to figure out
which files are actually the culprits. In that case try

.. code-block:: python

   import iris
   cube_list = iris.load('tas_*.nc')
   concatenated_cube_list = cube_list.concatenate()
   print(concatenated_cube_list)

This will concatenate all the consecutive bits it can manage and keep only those
parts separate, that cannot be merged. This way, you can often boil down the
problem to only two or a few cubes. Then you can diagnose those in the same way
as discussed above, eg

.. code-block:: python

   concatenated_cube_list.concatenate_cube()  # or
   concatenated_cube_list.merge_cube()

You can also inspect the attributes and coordinates of the resulting cubes with
a particular focus on the areas of problems pointed to by the output of
`merge_cube`.
