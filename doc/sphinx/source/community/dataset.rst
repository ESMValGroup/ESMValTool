.. _new-dataset:

Making a new dataset
********************

If you are contributing a new dataset, please have a look at :ref:`new-cmorizer` for how to do so.
If you need the new dataset for a new recipe, please make a separate pull
for the CMORizer script.

.. _dataset-documentation:

Dataset documentation
=====================

The documentation required for a CMORizer script is the following:

- Make sure that the new dataset is added to the list of
  :ref:`supported_datasets`
- The in code documentation should contain clear instructions on how to obtain
  the data
- A BibTeX file named ``<dataset>.bibtex`` defining the reference for the new
  dataset should be placed in the directory ``esmvaltool/references/``, see
  :ref:`adding_references` for detailed instructions.

For more general information on writing documentation, see :ref:`documentation`.

.. _dataset-test:

Testing
=======

When contributing a new script, add an entry for the CMORized data to
`recipes/examples/recipe_check_obs.yml <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/recipes/examples/recipe_check_obs.yml>`__
and run the recipe, to make sure the CMOR checks pass without warnings or errors.

To test a pull request for a new CMORizer script:

#. Download the data following the instructions included in the script and place
   it in the ``RAWOBS`` path specified in your ``config-user.yml``
#. Run the CMORizer script by running ``cmorize_obs -c <config-file> -o <dataset>``
#. Copy the resulting data to the ``OBS`` (for CMIP5 compliant data) or ``OBS6``
   (for CMIP6 compliant data) path specified in your
   ``config-user.yml``
#. Run ``recipes/examples/recipe_check_obs.yml`` with the new dataset to check that
   the data can be used

.. _dataset-sanity-check:

Scientific sanity check
=======================

When contributing a new dataset, we expect that the numbers and units of the dataset look physically meaningful.
The scientific reviewer needs to check this.

Data availability
=================

Once your pull request has been approved by the reviewers, ask
`@remi-kazeroni <https://github.com/remi-kazeroni>`_
to add the new dataset to the data pool at DKRZ and CEDA-Jasmin.
He is also the person in charge of merging CMORizer pull requests.

.. _dataset_checklist:

Detailed checklist for reviews
==============================

This (non-exhaustive) checklist provides ideas for things to check when reviewing
pull requests for new or updated CMORizer scripts.

Dataset description
-------------------

Check that new dataset has been added to the table of observations defined in
the ESMValTool guide user’s guide in section :ref:`inputdata`
(generated from ``doc/sphinx/source/input.rst``).

BibTeX info file
----------------

Check that a BibTeX file, i.e. ``<dataset>.bibtex`` defining the reference for
the new dataset has been created in ``esmvaltool/references/``.

recipe_check_obs.yml
--------------------

Check that new dataset has been added to the testing recipe
``esmvaltool/recipes/examples/recipe_check_obs.yml``

CMORizer script
---------------

Check that the new CMORizer script
``esmvaltool/cmorizers/obs/cmorize_obs_<dataset>.{py,ncl}``
meets standards.
This includes the following items:

* In-code documentation (header) contains

  1. Download instructions
  2. Reference(s)

* Code quality checks

  1. Code quality (e.g. no hardcoded pathnames)
  2. No Codacy errors reported


Config file
-----------

If present, check config file ``<dataset>.yml`` in
``esmvaltool/cmorizers/obs/cmor_config/`` for correctness.
Use ``yamllint`` to check for syntax errors and common mistakes.

Run CMORizer
------------

Make sure CMORizer is working by running ``cmorize_obs -c <config-file> -o <dataset>``

Check output of CMORizer
------------------------

After successfully running the new CMORizer, check that:

* Output contains (some) valid values (e.g. not only nan or zeros)
* Metadata is defined properly

Run ``esmvaltool/recipes/examples/recipe_check_obs.yml`` for new dataset.


RAW data
--------

Contact person in charge of ESMValTool data pool (`@remi-kazeroni`_) and
request to copy RAW data to RAWOBS/Tier2 (Tier3).


CMORized data
-------------

Contact person in charge of ESMValTool data pool (`@remi-kazeroni`_) and
request to

* Merge the pull request
* Copy CMORized dataset to OBS/Tier2 (Tier3)
* Set file access rights for new dataset
