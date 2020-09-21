"""Lint tests."""
import os
import pytest
import subprocess
import textwrap
from pathlib import Path

import esmvaltool
from esmvaltool.utils.nclcodestyle import nclcodestyle


def test_nclcodestyle():
    """Test that NCL code is formatted according to our standards."""
    package_root = Path(esmvaltool.__file__).absolute().parent
    check_paths = [
        package_root,
    ]

    print("Formatting check of NCL code in directories: {}\n".format(
        ', '.join(str(p) for p in check_paths)))

    exclude_paths = [
        package_root / 'diag_scripts' / 'cvdp' / 'cvdp',
    ]

    style = nclcodestyle.StyleGuide()
    style.options.exclude.extend(str(p) for p in exclude_paths)

    success = style.check_files(str(p) for p in check_paths).total_errors == 0

    if not success:
        print(
            textwrap.dedent("""
            Your NCL code does not follow our formatting standards.

            A list of warning and error messages can be found above,
            prefixed with filename:line number:column number.

            Please fix the mentioned issues.
        """))

    assert success, "Your NCL code does not follow our formatting standards."


@pytest.mark.installation
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
