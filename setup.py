#!/usr/bin/env python
"""ESMValTool installation script"""

import sys
import unittest

import coverage
from setuptools import Command, setup

import esmvaltool


class RunTests(Command):
    """Run tests and generate a coverage report htmlcov/index.html"""

    cov = coverage.Coverage(source=['esmvaltool'])
    cov.set_option("run:branch", True)
    cov.set_option("html:title", 'Coverage report for ESMValTool')
    cov.start()

    loader = unittest.TestLoader()
    tests = loader.discover('tests', pattern='test_*.py', top_level_dir='.')
    runner = unittest.TextTestRunner(verbosity=2)
    results = runner.run(tests)

    cov.stop()
    cov.save()
    cov.report()
    cov.html_report()

    sys.exit(0 if results.wasSuccessful() else 1)


with open('README.md') as readme:
    setup(
        name='ESMValTool',
        version=esmvaltool.__version__,
        description=readme.read(),
        packages=[
            'esmvaltool',
        ],
        # Include all version controlled files
        include_package_data=True,
        use_scm_version=True,
        setup_requires=['setuptools_scm'],
        install_requires=[
            'numpy',
            'netCDF4',
            'matplotlib'
            'basemap',
            'iris',
            'esgf-pyclient',
            'scientificpython',
        ],
        tests_requires=[
            'nose',
            'easytest',
            'mock',
            'pytest',
            'coverage',
        ],
        entry_points={
            'console_scripts': ['esmvaltool = esmvaltool.main:main'],
        },
        cmdclass={
            'test': RunTests,
        },
    )
