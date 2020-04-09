"""Common utilities for KCS scripts

Note that the name of this subpackage is "double" namespaced: kcsutils
inside a kcs directory. The reason is to avoid problems that importing
`utils` may lead to importing a different module: `import kcsutils (as
utils)` should be abundantly clear.

"""

from .io import read_data, concat_cubes
from .matching import match
from .date import make_year_constraint_all_calendars
from .constraints import EqualConstraint, RangeConstraint
