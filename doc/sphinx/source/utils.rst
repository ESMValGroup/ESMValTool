.. _utils:

Utilities
*********

This section provides information on small tools that are available in the
`esmvaltool/utils <https://github.com/ESMValGroup/ESMValTool/tree/master/esmvaltool/utils>`_
directory.

.. _draft_release_notes.py:

draft_release_notes.py
======================

This is a script for drafting release notes based on the titles and labels of
the GitHub pull requests that have been merged since the previous release.

To use the tool, install the package pygithub_:

.. code-block:: bash

   pip install pygithub

Create a `GitHub access token`_ (leave all boxes for additional
permissions unchecked) and store it in the file ``~/.github_api_key``.

Edit the script and update the date and time of the previous release. If needed,
change it so it uses the correct repository.

Run the script:

.. code-block:: bash

   python esmvaltool/utils/draft_release_notes.py

Review the resulting output (in ``.rst`` format) and if anything needs changing,
change it on GitHub and re-run the script until the changelog looks acceptable.
In particular, make sure that pull requests have the correct label, so they are
listed in the correct category.
Finally, copy and paste the generated content at the top of the changelog.

nclcodestyle
============

A tool for checking the style of NCL code, based on pycodestyle.
Install ESMValTool in development mode (``pip install -e '.[develop]'``) to make it available.
To use it, run

.. code-block:: bash

    nclcodestyle /path/to/file.ncl


xml2yml
=======

A tool for converting version 1 recipes to version 2 recipes.
See the README.md file in the directory esmvaltool/utils/xml2yml for detailed usage instructions.


testing
=======

Tools for testing recipes.

test recipe settings
--------------------

A tool for generating recipes with various diagnostic settings, to test of those work.
Install ESMValTool in development mode (``pip install -e '.[develop]'``) to make it available.
To use it, run

.. code-block:: bash

    test_recipe --help

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
  return datasets with partial coverage from `start_year` to `end_year`;
- `config-user: rootpath: CMIPX` may be a list, rootpath lists are supported;
- all major DRS paths (including `default`, `BADC`, `ETHZ` etc) are supported;
- speedup is achieved through CMIP mip tables lookup, so `mip` is required in recipe;

Caveats
-------

- the tool doesn't yet work with derived variables; it will not return any available datasets;
- operation restricted to CMIP data only, OBS lookup is not available yet.
