"""
Module for custom Jinja2 Tests.

In Jinja2 "tests" are used to test a variable against an expression.
Refer to Jinja2 tests:
https://jinja.palletsprojects.com/en/stable/templates/#tests

Jinja2 has a number of Builtin tests:
https://jinja.palletsprojects.com/en/stable/templates/#list-of-builtin-tests

Jinja2 also allows for writing custom tests in Python:
https://jinja.palletsprojects.com/en/stable/api/#custom-tests

Cylc supports custom filters, tests and globals as described here:
https://cylc.github.io/cylc-doc/stable/html/user-guide/writing-workflows/jinja2.html#custom-jinja2-filters-tests-and-globals

"""

from pathlib import Path


def file_exists(file_path):
    """
    Test if a file exists.

    Parameters
    ----------
    file_path : str
        A file path.

    Returns
    -------
    bool
        Returns True if the file exists.
    """
    run_directory = Path.cwd()
    site_recipes_file_path = run_directory / file_path
    return site_recipes_file_path.is_file()
