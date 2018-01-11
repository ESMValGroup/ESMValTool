#########################################################################
# FILE OPERATIONS
#########################################################################
import logging

import iris

from ..preprocessor._io import concatenate_callback

logger = logging.getLogger(__name__)


# merge multiple files assigned to a same diagnostic and variable
def glob(file_list, varname):
    """
    Function that takes a list of nc files and globs them into a single one
    """
    # there may be the case where the nc file contains multiple cubes
    # and/or the time calendars are crooked
    # -> these exceptional cases are nicely solved by applying Javier's nice
    # iris fixing (merge_callback and get_attr_from_field_coord are also
    # in preprocess.py but keeping them here in case we will have to
    # change things)

    var_name = varname

    def cube_var_name(raw_cube):
        return raw_cube.var_name == var_name

    var_cons = iris.Constraint(cube_func=cube_var_name)
    # force single cube; this function defaults a list of cubes
    cl = [
        iris.load(a, var_cons, callback=concatenate_callback)[0]
        for a in file_list
    ]

    c = iris.cube.CubeList(cl)

    try:
        concatenated = c.concatenate()
        try:
            logger.info("Successfully concatenated cubes")
            return concatenated[0]
        except (OSError, iris.exceptions.IrisError) as exc:
            logger.warning("Could not save concatenated cube, keeping a "
                           "list of files - %s ", exc)
            return 0
    except iris.exceptions.ConcatenateError as exc:
        error_message = "Problem trying to concatenate the following cubes:\n"
        for cube in cl:
            error_message += cube.summary(shorten=True) + '\n'
        logger.warning(
            "Could not concatenate cubes, keeping a list of files - %s",
            error_message)
        return 0
