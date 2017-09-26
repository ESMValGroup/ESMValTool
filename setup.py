#!/usr/bin/env python
"""ESMValTool installation script"""
# This script installs dependencies available on PyPI
#
# Dependencies that need to be installed some other way (e.g. conda):
# - ncl
# - iris
# - basemap

import sys
import unittest

from setuptools import Command, setup

PACKAGES = [
    'esmvaltool',
]


class RunTests(Command):
    """Class to run tests and generate a coverage report htmlcov/index.html"""
    user_options = []

    def initialize_options(self):
        """Do nothing"""

    def finalize_options(self):
        """Do nothing"""

    def run(self):
        """Run tests and generate a coverage report."""
        import coverage

        cov = coverage.Coverage(source=PACKAGES)
        cov.set_option("run:branch", True)
        cov.set_option("html:title", 'Coverage report for ESMValTool')
        cov.start()

        loader = unittest.TestLoader()
        tests = loader.discover(
            'tests', pattern='test_*.py', top_level_dir='.')
        runner = unittest.TextTestRunner(verbosity=2)
        results = runner.run(tests)

        cov.stop()
        cov.save()
        cov.report()
        cov.html_report()
        cov.xml_report()

        sys.exit(0 if results.wasSuccessful() else 1)


with open('README.md') as readme:
    setup(
        name='ESMValTool',
        version="2.0.0",
        description=readme.read(),
        packages=PACKAGES,
        # Include all version controlled files
        include_package_data=True,
        use_scm_version=True,
        setup_requires=[
            'setuptools_scm',
        ],
        install_requires=[
            'cdo',
            'cf_units',
            'coverage',
            'cython',
            'esgf-pyclient',
            'numpy',
            'netCDF4',
            'matplotlib',
            'pyyaml',
            'shapely',
        ],
        tests_require=[
            'coverage',
            'easytest',
            'mock',
            'nose',
            'pytest',
        ],
        entry_points={
            'console_scripts': [
                'esmvaltool = esmvaltool.main:run',
            ],
        },
        cmdclass={
            'test': RunTests,
        }, )
