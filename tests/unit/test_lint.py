"""Lint tests."""
import os
import subprocess
import textwrap

from esmvaltool.utils.nclcodestyle import nclcodestyle


def test_nclcodestyle():
    """Test that NCL code is formatted according to our standards."""
    check_paths = [
        'esmvaltool',
        'tests',
    ]

    exclude_paths = [
        'esmvaltool/diag_scripts/cvdp/cvdp',
    ]

    print("Formatting check of NCL code in directories: {}\n".format(
        ', '.join(check_paths)))

    # Get paths wrt package root
    package_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    for paths in (check_paths, exclude_paths):
        for i, path in enumerate(paths):
            paths[i] = os.path.join(package_root, path)

    style = nclcodestyle.StyleGuide()
    style.options.exclude.extend(exclude_paths)

    success = style.check_files(check_paths).total_errors == 0

    if not success:
        print(
            textwrap.dedent("""
            Your NCL code does not follow our formatting standards.

            A list of warning and error messages can be found above,
            prefixed with filename:line number:column number.

            Please fix the mentioned issues.
        """))

    assert success, "Your NCL code does not follow our formatting standards."


def test_r_lint(monkeypatch):
    """Test R lint."""
    monkeypatch.setenv("LINTR_COMMENT_BOT", "FALSE")
    package_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    checker = os.path.join(package_root, 'tests', 'unit', 'check_r_code.R')
    try:
        output = subprocess.check_output(('Rscript', checker, package_root),
                                         stderr=subprocess.STDOUT,
                                         universal_newlines=True)
        print(output)
        return
    except subprocess.CalledProcessError as ex:
        print(
            textwrap.dedent("""
            Your R code does not follow our formatting standards.

            Please fix the following issues:
        """))
        print(ex.output)

    assert False,\
        'Your R code does not follow our formatting standards.'
