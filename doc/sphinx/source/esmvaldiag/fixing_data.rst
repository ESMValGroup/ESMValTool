.. fixing_data:

**********
Data fixes
**********

One of the most common sources of frustation for diagnostic developers is data check errors.
They appear often and don´t let you test you run.

What can you do to get rid of them? ESMValTool has support for automatic on the fly fixes for
those datasets that had known errors that can be fixed

Fix structure
=============

Fixes are clases stored in ``esmvaltool/cmor/_fixes/[PROJECT]/[DATASET].py``
that derive from :class:`esmvalcore.cmor._fixes.fix.Fix` and
named after the short name of the variable they fix. You can use the name
``allvars`` if you want the fix to be applied to the whole dataset

.. warning::
    Be careful to replace any ``-`` with ``_`` in your dataset name.
    We need this replacement to have proper python module names

They are automatically loaded and applied when the dataset is preprocessed.

Fixing a dataset
================

To illustrate the process of creating a fix we are going to construct a new
one from scratch for a fictional dataset. We need to fix a CMIP7 model
called PERFECT-MODEL that is reporting a missing latitude coordinate for
variable tas.

Check the output
----------------

Next to the error message, you should see some info about the iris cube: size,
coordinates. In our example it looks like this:

.. code-block:: python

   cube output with latitude called altitude


So now the mistake is clear: the latitude coordinate is badly named and the
fix should just rename it

Create the fix
--------------

We start by creating the module file. In our example the path will be
``esmvaltool/cmor/_fixes/CMIP7/PERFECT_MODEL.py``. If it already exists
just add the class to the file, there is no limit in the number of fixes
we can have in any given file

Then we have to create the class for the fix deriving from
:class:`esmvalcore.cmor._fixes.Fix`

.. code-block:: python

    """Fixes for PERFECT-MODEL"""
    from esmvalcore.cmor.fix import Fix

    class tas(Fix):
         """Fixes for tas class""""

And now we must choose the method to use between the ones offered by the
Fix class:

- ``fix_file`` : should be used only to fix errors that prevent data loading.
  As a rule of thumb, you should only use it you are not even
  reaching the checks

- ``fix_metadata`` : you want to change something in the cube that is not the data.

- ``fix_data``: you need to fix the data. Beware: coordinates data values are part
  of the metadata

In our case we need to rename the coordinate ``altitude`` to ``latitude``,
so we will implement the ``fix_metadata`` method

.. code-block:: python

    """Fixes for PERFECT-MODEL."""
    from esmvalcore.cmor.fix import Fix

    class tas(Fix):
        """Fixes for tas variable.""""

        def fix_metadata(self, cubes):
            """
            Fix metadata for tas

            Fix the name of the latitude coordinate, which is called altitude
            in the original file.
            """"
            # Sometimes Iris will interpret the data as multiple cubes.
            # Good CMOR datasets will only show one but we support the
            # multiple cubes case to be able to fix the errors that are
            # leading to that extra cubes.
            # In our case this means that we can safely assume that the
            # tas cube is the first one
            tas_cube = cubes[0]
            latitude = tas_cube.coord('altitude')

            # Fix the names. Latitude values, units and
            latitude.short_name = 'lat'
            latitude.standard_name = 'latitude'
            latitude.long_name = 'latitude'

And that's all. The next time you run ESMValTool you will find that the error is fixed on
the fly and, hopefully, your recipe will run free of errors.

Sometimes other errors can appear after you fix the first one because they were hidden by it.
In our case,  the latitude coordinate could have bad units or values outside the
valid range for example. Just extend your fix to fix those errors and keep going

Finishing
---------

Chances are that you are not the only one that wants to use that dataset and variable. Other users
will be very grateful to have your fixes available as soon as possible. Please, create a separated
pull request for the fix and submit it.

It will also be very helpful if you just scan a couple of other variables from the same dataset and
check if they share this error. In case that you find that it is a general one, you can change the fix
name to ``allvars`` so it gets executed for the full dataset. If you find that this is shared only by
a handful of similar vars you can just make the fix for those new vars derive from the one you just created:

.. code-block:: python

    """Fixes for PERFECT-MODEL."""
    from esmvalcore.cmor.fix import Fix

    class tas(Fix):
        """Fixes for tas variable.""""

        def fix_metadata(self, cubes):
            """
            Fix metadata for tas

            Fix the name of the latitude coordinate, which is called altitude
            in the original file.
            """"
            # Sometimes Iris will interpret the data as multiple cubes.
            # Good CMOR datasets will only show one but we support the
            # multiple cubes case to be able to fix the errors that are
            # leading to that extra cubes.
            # In our case this means that we can safely assume that the
            # tas cube is the first one
            tas_cube = cubes[0]
            latitude = tas_cube.coord('altitude')

            # Fix the names. Latitude values, units and
            latitude.short_name = 'lat'
            latitude.standard_name = 'latitude'
            latitude.long_name = 'latitude'


    class ps(tas):
        """Fixes for ps variable."""


Common errors
=============

Our example covered one of the most common cases: variables / coordinates that have names that
do not match the expected. But there are some others that use to appear frequently. This section
will describe them

Bad units declared
------------------

Is quite common that a variable declares to be using some units but the data
is stored in another. This can be solved ovwerwriting the units attribute
with the real data units.

.. code-block:: python

    ...
        def fix_metadata(self, cubes):
            ...
            cube.units = 'real_units'
            ...

Detecting this error can be tricky if the units are similar enough. It also
has a good chance of going undetected until you notice strange results in
your diagnostic


Coordinates missing
-------------------

Another common error is to have missing coordinates. Usually it just means
that the file does not follow the CF-conventions and Iris can not interpret it.

If that is the case, you should see a warning about ESMValTool discarding some
cubes in the fix metadata step. Above that warning you should see the full
list of cubes as readed by Iris. If that list contains your missing coordinate
you can create a fix for it following this model:

.. code-block:: bash
    ...
        def fix_metadata(self, cubes):
            coord_cube = cubes.extract_strict('COORDINATE_NAME')
            # Usually this will correspond to an auxiliary coordinate
            # because the most common error is to forget adding it to the
            # coordinates attribute
            coord = iris.coords.AuxCoord(
                coord_cube.data,
                var_name = coord_cube.var_name,
                standard_name = coord_cube.standard_name,
                long_name = coord_cube.long_name,
                units = coord_cube.units,
                attributes =
            }

            # It may also have bounds as another cube
            coord.bounds = cubes.extract_strict('BOUNDS_NAME').data
            ...
