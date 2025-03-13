#!/usr/bin/env python
"""ESMValTool installation script."""

import json
import os
import re
import sys
from pathlib import Path

from setuptools import Command, setup

PACKAGES = [
    "esmvaltool",
]

REQUIREMENTS = {
    # Installation script (this file) dependencies
    "setup": [
        "setuptools_scm",
    ],
    # Installation dependencies
    # Use with pip install . to install from source
    "install": [
        "aiohttp",
        "cartopy",
        "cdo",
        "cdsapi",
        "cf-units",
        "cfgrib",
        "cftime",
        "cmocean",
        "dask!=2024.8.0",  # https://github.com/dask/dask/issues/11296
        "distributed",
        "ecmwf-api-client",
        "eofs",
        "ESMPy",  # not on PyPI
        "esmvalcore",
        "esmf-regrid>=0.10.0",  # iris-esmf-regrid #342
        "fiona",
        "fire",
        "fsspec",
        "GDAL",
        "ipython<9.0",  # github.com/ESMValGroup/ESMValCore/issues/2680
        "jinja2",
        "joblib",
        "lime",
        "mapgenerator>=1.0.5",
        "matplotlib",
        "natsort",
        "nc-time-axis",
        "netCDF4",
        "numba",
        "numpy!=1.24.3",  # severe masking bug
        "openpyxl",
        "packaging",
        "pandas",
        "progressbar2",
        "psyplot>=1.5.0",  # psy*<1.5.0 are not py312 compat
        "psy-maps>=1.5.0",
        "psy-reg>=1.5.0",
        "psy-simple>=1.5.0",
        "pyproj>=2.1",
        "pys2index",
        "python-dateutil",
        "pyyaml",
        "rasterio>=1.3.10",
        "requests",
        "ruamel.yaml",
        "scikit-image",
        "scikit-learn>=1.4.0",  # github.com/ESMValGroup/ESMValTool/issues/3504
        "scipy",
        "scitools-iris>=3.11",
        "seaborn",
        "seawater",
        "shapely>=2",
        "xarray>=0.12.0",
        "xesmf>=0.7.1",
        "xgboost>1.6.1",  # github.com/ESMValGroup/ESMValTool/issues/2779
        "xlsxwriter",
        "zarr",
    ],
    # Test dependencies (unit tests)
    # Execute `pip install .[test]` once and then use `pytest` to run tests
    "test": [
        "pre-commit",
        "pytest>=3.9,!=6.0.0rc1,!=6.0.0",
        "pytest-cov>=2.10.1",
        "pytest-env",
        "pytest-html!=2.1.0",
        "pytest-metadata>=1.5.1",
        "pytest-mock",
        "pytest-xdist",
    ],
    # Documentation dependencies
    "doc": [
        "autodocsumm>=0.2.2",
        "nbsphinx",
        "sphinx>=6.1.3",
        "pydata-sphinx-theme",
    ],
    # Development dependencies
    # Use pip install -e .[develop] to install in development mode
    "develop": [
        "codespell",
        "docformatter",
        "imagehash",
        "pre-commit",
        "prospector[with_pyroma]>=1.12",
        "vprof",
        "yamllint",
    ],
}


def discover_python_files(paths, ignore):
    """Discover Python files."""

    def _ignore(path):
        """Return True if `path` should be ignored, False otherwise."""
        return any(re.match(pattern, path) for pattern in ignore)

    for path in sorted(set(paths)):
        for root, _, files in os.walk(path):
            if _ignore(path):
                continue
            for filename in files:
                filename = os.path.join(root, filename)
                if filename.lower().endswith(".py") and not _ignore(filename):
                    yield filename


class RunLinter(Command):
    """Class to run a linter and generate reports."""

    user_options = []

    def initialize_options(self):
        """Do nothing."""

    def finalize_options(self):
        """Do nothing."""

    def install_deps_temp(self):
        """Try to temporarily install packages needed to run the command."""
        if self.distribution.install_requires:
            self.distribution.fetch_build_eggs(
                self.distribution.install_requires
            )
        if self.distribution.tests_require:
            self.distribution.fetch_build_eggs(self.distribution.tests_require)

    def run(self):
        """Run prospector and generate a report."""
        check_paths = PACKAGES + [
            "setup.py",
            "tests",
            "util",
        ]
        ignore = [
            "doc/",
        ]

        # try to install missing dependencies and import prospector
        try:
            from prospector.run import main
        except ImportError:
            # try to install and then import
            self.distribution.fetch_build_eggs(["prospector[with_pyroma]"])
            from prospector.run import main

        self.install_deps_temp()

        # run linter

        # change working directory to package root
        package_root = os.path.abspath(os.path.dirname(__file__))
        os.chdir(package_root)

        # write command line
        files = discover_python_files(check_paths, ignore)
        sys.argv = ["prospector"]
        sys.argv.extend(files)

        # run prospector
        errno = main()

        sys.exit(errno)


def read_authors(filename):
    """Read the list of authors from .zenodo.json file."""
    with Path(filename).open() as file:
        info = json.load(file)
        authors = []
        for author in info["creators"]:
            name = " ".join(author["name"].split(",")[::-1]).strip()
            authors.append(name)
        return ", ".join(authors)


def read_description(filename):
    """Read the description from .zenodo.json file."""
    with Path(filename).open() as file:
        info = json.load(file)
        return info["description"]


setup(
    name="ESMValTool",
    author=read_authors(".zenodo.json"),
    description=read_description(".zenodo.json"),
    long_description=Path("README.md").read_text(),
    long_description_content_type="text/markdown",
    url="https://www.esmvaltool.org",
    download_url="https://github.com/ESMValGroup/ESMValTool",
    license="Apache License, Version 2.0",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "Topic :: Scientific/Engineering :: GIS",
        "Topic :: Scientific/Engineering :: Hydrology",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    packages=PACKAGES,
    # Include all version controlled files
    include_package_data=True,
    setup_requires=REQUIREMENTS["setup"],
    install_requires=REQUIREMENTS["install"],
    tests_require=REQUIREMENTS["test"],
    extras_require={
        "develop": REQUIREMENTS["develop"]
        + REQUIREMENTS["test"]
        + REQUIREMENTS["doc"],
        "doc": REQUIREMENTS["doc"],
        "test": REQUIREMENTS["test"],
    },
    entry_points={
        "console_scripts": [
            "nclcodestyle = esmvaltool.utils.nclcodestyle.nclcodestyle:_main",
            "test_recipe = "
            "esmvaltool.utils.testing.recipe_settings.install_expand_run:main",
        ],
        "esmvaltool_commands": [
            "colortables = "
            "esmvaltool.utils.color_tables.show_color_tables:ColorTables",
            "install = esmvaltool.install:Install",
            "data = esmvaltool.cmorizers.data.cmorizer:DataCommand",
        ],
    },
    cmdclass={
        "lint": RunLinter,
    },
    zip_safe=False,
)
