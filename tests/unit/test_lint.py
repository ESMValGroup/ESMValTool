""" Lint tests """
from __future__ import print_function

import os
import textwrap

import pycodestyle  # formerly known as pep8

from esmvaltool.utils.nclcodestyle import nclcodestyle


def test_pep8_conformance():
    """Test that we conform to PEP-8."""
    check_paths = [
        'esmvaltool',
        'tests',
    ]
    exclude_paths = [
        'esmvaltool/doc',
    ]

    print("PEP8 check of directories: {}\n".format(', '.join(check_paths)))

    # Get paths wrt package root
    package_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    for paths in (check_paths, exclude_paths):
        for i, path in enumerate(paths):
            paths[i] = os.path.join(package_root, path)

    style = pycodestyle.StyleGuide()
    style.options.exclude.extend(exclude_paths)

    success = style.check_files(check_paths).total_errors == 0

    if not success:
        print(textwrap.dedent("""
            Your Python code does not conform to the official Python style
            guide (PEP8), see https://www.python.org/dev/peps/pep-0008

            A list of warning and error messages can be found above,
            prefixed with filename:line number:column number.

            Run `yapf -i yourfile.py` to automatically fix most errors.
            Run `yapf -d yourfile.py` to preview what would be changed.
            Run `pip install --upgrade yapf` to install the latest version
            of yapf.
        """))

    assert success, "Your code does not conform to PEP8"


def test_nclcodestyle():
    """Test that NCL code is formatted according to our standards."""
    check_paths = [
        'esmvaltool',
        'tests',
    ]

    print("Formatting check of NCL code in directories: {}\n".format(
        ', '.join(check_paths)))

    style = nclcodestyle.StyleGuide()
    package_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    check_paths = [os.path.join(package_root, path) for path in check_paths]
    success = style.check_files(check_paths).total_errors == 0

    if not success:
        print(textwrap.dedent("""
            Your NCL code does not follow our formatting standards.

            A list of warning and error messages can be found above,
            prefixed with filename:line number:column number.

            Please fix the mentioned issues.
        """))

    assert success, "Your NCL code does not follow our formatting standards."


def test_r_lint():
    """Test R lint"""
    from rpy2.robjects.packages import importr

    lintr = importr('lintr')

    check_paths = [
        'esmvaltool',
        'tests',
    ]
    root_folder = os.path.abspath(os.path.join(__file__, '..', '..', '..'))

    has_errors = False

    linters = lintr.with_defaults()

    for path in check_paths:
        for dirpath, _, filenames in os.walk(os.path.join(root_folder, path)):
            for filename in filenames:
                if os.path.splitext(filename)[1].lower() == '.r':
                    errors = lintr.lint(os.path.join(dirpath, filename),
                                        linters)
                    for error in errors:
                        print(str(error)[0:-1])
                        has_errors = True
    assert not has_errors, 'Your R code does not follow our R standards'

