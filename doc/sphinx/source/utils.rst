.. _utils:

Utilities
*********

This section provides information on tools that are useful when developing
ESMValTool.
Tools that are specific to ESMValTool live in the
`esmvaltool/utils <https://github.com/ESMValGroup/ESMValTool/tree/master/esmvaltool/utils>`_
directory, while others can be installed using the usual package managers.

.. _pre-commit:

Pre-commit
==========

`pre-commit <https://pre-commit.com/>`__ is a handy tool that can run many
tools for checking code quality with a single command.
Usually it is used just before committing, to avoid accidentally committing
mistakes.
It knows knows which tool to run for each filetype, and therefore provides
a convenient way to check your code!


To run ``pre-commit`` on your code, go to the ESMValTool directory
(``cd ESMValTool``) and run

::

   pre-commit run

By default, pre-commit will only run on the files that have been changed,
meaning those that have been staged in git (i.e. after
``git add your_script.py``).

To make it only check some specific files, use

::

   pre-commit run --files your_script.py

or

::

   pre-commit run --files your_script.R

Alternatively, you can configure ``pre-commit`` to run on the staged files before
every commit (i.e. ``git commit``), by installing it as a `git hook <https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks>`__ using

::

   pre-commit install

Pre-commit hooks are used to inspect the code that is about to be committed. The
commit will be aborted if files are changed or if any issues are found that
cannot be fixed automatically. Some issues cannot be fixed (easily), so to
bypass the check, run

::

   git commit --no-verify

or

::

   git commit -n

or uninstall the pre-commit hook

::

   pre-commit uninstall


Note that the configuration of pre-commit lives in
`.pre-commit-config.yaml <https://github.com/ESMValGroup/ESMValTool/blob/master/.pre-commit-config.yaml>`_.

.. _nclcodestyle:

nclcodestyle
============

A tool for checking the style of NCL code, based on pycodestyle.
Install ESMValTool in development mode (``pip install -e '.[develop]'``) to make it available.
To use it, run

.. code-block:: bash

    nclcodestyle /path/to/file.ncl

.. _recipe_test_tool:

Colormap samples
================
Tool to generate colormap samples for ESMValTool's default Python and NCL colormaps.

Run

.. code-block:: bash

    esmvaltool colortables python

or

.. code-block:: bash

    esmvaltool colortables ncl

to generate the samples.

Testing recipes
===============

Tools for testing recipes.

Test recipe settings
--------------------

A tool for generating recipes with various diagnostic settings, to test of those work.
Install ESMValTool in development mode (``pip install -e '.[develop]'``) to make it available.
To use it, run

.. code-block:: bash

    test_recipe --help


.. _draft_release_notes.py:

draft_release_notes.py
======================

`draft_release_notes.py <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/utils/draft_release_notes.py>`__
is a script for drafting release notes based on the titles and labels of
the GitHub pull requests that have been merged since the previous release.

To use it, install the package pygithub_:

.. code-block:: bash

   pip install pygithub

Create a `GitHub access token`_ (leave all boxes for additional
permissions unchecked) and store it in the file ``~/.github_api_key``.

Edit the script and update the date and time of the previous release and run
the script:

.. code-block:: bash

   python esmvaltool/utils/draft_release_notes.py ${REPOSITORY}

``REPOSITORY`` can be either ``esmvalcore`` or ``esmvaltool`` depending on the
release notes you want to create.

Review the resulting output (in ``.rst`` format) and if anything needs changing,
change it on GitHub and re-run the script until the changelog looks acceptable.
In particular, make sure that pull requests have the correct label, so they are
listed in the correct category.
Finally, copy and paste the generated content at the top of the changelog.

Converting Version 1 Namelists to Version 2 Recipes
===================================================

The
`xml2yml <https://github.com/ESMValGroup/ESMValTool/tree/master/esmvaltool/utils/xml2yml>`_
converter can turn the old xml namelists into new-style yml
recipes. It is implemented as a xslt stylesheet that needs a processor
that is xslt 2.0 capable. With this, you simply process your old
namelist with the stylesheet xml2yml.xsl to produce a new yml recipe.

After the conversion you need to manually check the mip information in
the variables! Also, check the caveats below!

Howto
-----

One freely available processor is the Java based
`saxon <http://saxon.sourceforge.net/>`__. You can download the free he
edition
`here <https://sourceforge.net/projects/saxon/files/latest/download?source=files>`__.
Unpack the zip file into a new directory. Then, provided you have Java
installed, you can convert your namelist simply with:

::

   java -jar $SAXONDIR/saxon9he.jar -xsl:xml2yml.xsl -s:namelist.xml -o:recipe.yml

Caveats/Known Limitations
-------------------------

-  At the moment, not all model schemes (OBS, CMIP5, CMIP5_ETHZ…) are
   supported. They are, however, relatively easy to add, so if you need
   help adding a new one, please let me know!
-  The documentation section (namelist_summary in the old file) is not
   automatically converted.
-  In version 1, one could name an exclude, similar to the reference
   model. This is no longer possible and the way to do it is to include
   the models with another ``additional_models`` tag in the variable
   section. That conversion is not performed by this tool.

Authored by **Klaus Zimmermann**, direct questions and comments to
klaus.zimmermann@smhi.se

.. _GitHub access token: https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line
.. _pygithub: https://pygithub.readthedocs.io/en/latest/introduction.html


Recipe filler
=============

If you need to fill in a blank recipe with additional datasets, you can do that with
the command `recipe_filler`. This runs a tool to obtain a set of additional datasets when
given a blank recipe, and you can give an arbitrary number of data parameters. The blank recipe
should contain, to the very least, a list of diagnostics, each with their variable(s).
Example of running the tool:

.. code-block:: bash

    recipe_filler recipe.yml

where `recipe.yml` is the recipe that needs to be filled with additional datasets; a minimal
example of this recipe could be:

.. code-block:: yaml

    diagnostics:
      diagnostic:
        variables:
          ta:
            mip: Amon  # required
            start_year: 1850  # required
            end_year: 1900  # required


Key features
------------

- you can add as many variable parameters as are needed; if not added, the
  tool will use the ``"*"`` wildcard and find all available combinations;
- you can restrict the number of datasets to be looked for with the ``dataset:``
  key for each variable, pass a list of datasets as value, e.g.
  ``dataset: [MPI-ESM1-2-LR, MPI-ESM-LR]``;
- you can specify a pair of experiments, e.g. ``exp: [historical, rcp85]``
  for each variable; this will look for each available dataset per experiment
  and assemble an aggregated data stretch from each experiment to complete
  for the total data length specified by ``start_year`` and ``end_year``; equivalent to
  ESMValTool's syntax on multiple experiments; this option needs an ensemble
  to be declared explicitly; it will return no entry if there are gaps in data;
- ``start_year`` and ``end_year`` are required and are used to filter out the
  datasets that don't have data in the interval; as noted above, the tool will not
  return datasets with partial coverage from ``start_year`` to ``end_year``;
  if you want all possible years hence no filtering on years just use ``"*"``
  for start and end years;
- ``config-user: rootpath: CMIPX`` may be a list, rootpath lists are supported;
- all major DRS paths (including ``default``, ``BADC``, ``ETHZ`` etc) are supported;
- speedup is achieved through CMIP mip tables lookup, so ``mip`` is required in recipe;

Caveats
-------

- the tool doesn't yet work with derived variables; it will not return any available datasets;
- operation restricted to CMIP data only, OBS lookup is not available yet.
