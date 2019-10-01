.. _running:

******************
Running ESMValTool
******************

ESMValTool is mostly used as a command line tool. Whenever your
conda environment for ESMValTool is active, you can just run the command
``esmvaltool <options> <recipe>``. One of the options that must be specified
is the user configuration file, which is specified using the
option ``-c /path/to/config-user.yaml``. An example recipe is available in the
ESMValTool installation folder as ``examples/recipe_python.yml``:

.. code-block:: yaml

	# ESMValTool
	# recipe_python.yml
	---
	documentation:
	  description: |
	    Example recipe that plots the mean precipitation and temperature.

	  authors:
	    - andela_bouwe
	    - righi_mattia

	  maintainer:
	    - schlund_manuel

	  references:
	    - acknow_project

	  projects:
	    - esmval
	    - c3s-magic

	datasets:
	  - {dataset: CanESM2,  project: CMIP5,  mip: Amon,  exp: historical,  ensemble: r1i1p1,  start_year: 2000,  end_year: 2002}
	  - {dataset: MPI-ESM-LR,  project: CMIP5,  mip: Amon,  exp: historical,  ensemble: r1i1p1,  start_year: 2000,  end_year: 2002}

	preprocessors:

	  preprocessor1:
	    extract_levels:
	      levels: 85000
	      scheme: nearest
	    regrid:
	      target_grid: 1x1
	      scheme: linear
	    multi_model_statistics:
	      span: overlap
	      statistics: [mean, median]

	diagnostics:

	  diagnostic1:
	    description: Air temperature and precipitation Python tutorial diagnostic.
	    themes:
	      - phys
	    realms:
	      - atmos
	    variables:
	      ta:
	        preprocessor: preprocessor1
	      pr:
	    scripts:
	      script1:
	        script: examples/diagnostic.py
	        quickplot:
	          plot_type: pcolormesh

This recipe finds data from CanESM2 and MPI-ESM-LR for 2000 - 2002,
extracts a single level (850 hPa), regrids it to a 1x1 degree mesh and runs
a diagnostic script that creates some plots of Air temperature and
precipitation flux. You can copy the recipe above and save it in your project
directory as (e.g.) ``example_recipe_python.yml`` and then run ESMValTool with

.. code:: bash

	esmvaltool -c /path/to/config-user.yml example recipe_python.yml -s

The ``-s`` option tells ESMValTool to use Synda to search for and download
the necessary datasets.

ESMValTool will also find recipes that are stored in its install directory.
A copy of the example recipe is available in ``examples/recipe_python``. Thus,
the following also works:

.. code:: bash

	esmvaltool -c /path/to/config-user.yml examples/recipe_python.yml

Note that this command does not call Synda. The required data should thus be
located in the directories specified in your user config file.
Recall that the chapter :ref:`User configuration file <config-user>`
provides an explanation of how to create your own config-user.yml file.

To get help on additional commands, please use

.. code:: bash

	esmvaltool --help


Available diagnostics and metrics
=================================

See Section :doc:`Recipes <../recipes/index>` for a description of all
available recipes.
