'''
Module to provide an Iris constraint mechanism for Cell Methods
'''

from iris import Constraint
from itertools import izip_longest


class MethodConstraint(Constraint):
    """Provides a simple Cube-method based :class:`Constraint`."""
    def __init__(self, coord, **kwargs):
        """
        Accept any cube that has a class:`CellMethod` with the specified
        coordinate that optionally matches the given method or interval.

        Args:

        * coord:
            The coordinate name for which we want any method to include

        Kwargs:

        * method:
            A specific method that we want for the coordinate

        * interval:
            A description of the sampling interval for the method, can be
            callable function (default units is hours as per lbtim parameter
            in pp metadata)

        Example usage::

            iris.MethodConstraint('time', method='mean')

            iris.MethodConstraint('time',
                interval=lambda interval: interval<24)

        """
        self._coord = coord
        self._kwargs = kwargs
        Constraint.__init__(self, cube_func=self._cube_func)

    def _cube_func(self, cube):
        match = False

        # coord must be a coordinate in the cube!
        if not cube.coords(self._coord):
            raise Exception('Coord %s is not in cube' % coord)

        # Cube must have cell method for given coordinate
        nomethod = True
        for cm in cube.cell_methods:
            if self._coord in cm.coord_names:
                nomethod = False
                match = True

                # Now check method condition if it exists
                if 'method' in self._kwargs:
                    if self._kwargs['method'] != cm.method:
                        match = False

                # Now check interval condition if it exists
                # This assumes that interval is given in hours as per lbtim
                #  parameter in pp metadata, but it needs to be extended
                #  for use beyond this
                if 'interval' in self._kwargs:
                    cell_intervals = dict(([x, y] for x, y in
                                           izip_longest(cm.coord_names,
                                                        cm.intervals)))
                    # Interval is a string combining integer value of interval
                    # and units. How to extend to take account of units?
                    # What if cell_intervals[coord] is None?
                    interval = int(cell_intervals[self._coord].split(' ')[0])
                    if callable(self._kwargs['interval']):
                        if not self._kwargs['interval'](interval):
                            match = False
                    else:
                        if interval != self._kwargs['interval']:
                            match = False

                if match:
                    break

        # Now deal with instantaneous edge case, i.e. no cell_method will
        # exist for coord. For this case the user must specify "method=None"
        # in constraint. Essentially the above code will find no matches for
        # the coord and return that the cube will fail, however setting
        # "method=None" will reverse this decision as we demand no matches.
        if nomethod:
            if 'method' in self._kwargs:
                if self._kwargs['method'] is None:
                    match = True

        # Look away from the screen!!!
        # Naughty fudge for specific stash codes like landseamask
        if str(cube.attributes.get('STASH', '')) in ('m01s00i030', ):
            match = True
        # Okay, can look back now!!!

        return match

    def __repr__(self):
        return 'MethodConstraint({})'.format(self._kwargs)
