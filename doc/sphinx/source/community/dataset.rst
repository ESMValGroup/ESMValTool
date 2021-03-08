.. _new-dataset:

********************
Making a new dataset
********************

If you are contributing a new dataset, please have a look at :ref:`new-cmorizer` for how to do so.
If you need the new dataset for a new recipe, please make a separate pull
for the CMORizer script.

.. _dataset-documentation:

Dataset documentation
=====================

Make sure that the new dataset is added to the list of
:ref:`supported_datasets`.
The documentation should contain clear instructions on how to obtain the data.

For more information on writing documentation, see :ref:`documentation`.

.. _dataset-test:

Testing
=======

Add a test for the CMORized data to
`recipes/examples/recipe_check_obs.yml <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/recipes/examples/recipe_check_obs.yml>`__
and run the recipe, to make sure the CMOR checks pass without warnings or errors.

.. _dataset-sanity-check:

Scientific sanity check
=======================

When contributing a new dataset, we expect that the numbers and units of the dataset look physically meaningful.
The scientific reviewer needs to check this.

Data availability
=================

Once your contribution has been accepted, ask
`@remi-kazeroni <https://github.com/remi-kazeroni>`_
to add the new dataset to the data pool at DKRZ and CEDA-Jasmin.
He is also the person in charge of merging CMORizer pull requests.
