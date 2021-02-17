.. _new-dataset:

********************
Making a new dataset
********************

If you are contributing a new dataset, please have a look at :ref:`new-cmorizer` for how to do so.


.. _dataset-documentation:

Documentation
=============

Make sure that the new dataset is added to :ref:`supported_datasets`.

The documentation should have clear instructions on how to obtain the data.

Testing
=======

Add a test for the CMORized data to ``recipes/example/recipe_check_obs.yml`` and run the recipe, to make sure the CMOR checks pass without errors.

CMORizer output
===============

When contributing a new dataset, we expect that the numbers and units of the dataset look physically meaningful.


Adding your dataset to the OBS data pool
========================================

Once your contribution has been accepted, tag `@remi-kazeroni <https://github.com/remi-kazeroni>`_ in your pull request, so that the new dataset can be added to the OBS data pool at DKRZ and synchronized with CEDA-Jasmin
