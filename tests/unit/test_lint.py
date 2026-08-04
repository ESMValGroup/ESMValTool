"""Lint tests."""

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

    print(
        "Formatting check of NCL code in directories: {}\n".format(
            ", ".join(str(p) for p in check_paths),
        ),
    )

    exclude_paths = [
        package_root / "diag_scripts" / "cvdp" / "cvdp",
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
        """),
        )

    assert success, "Your NCL code does not follow our formatting standards."
