#!/usr/bin/env python
"""ESMValTool installation script."""
# This script only installs dependencies available on PyPI
#
# Dependencies that need to be installed some other way (e.g. conda):
# - ncl
# - iris
# - python-stratify

import os
import re
import sys

from setuptools import Command, setup

from esmvaltool._version import __version__

PACKAGES = [
    'esmvaltool',
]

REQUIREMENTS = {
    # Installation script (this file) dependencies
    'setup': [
        'setuptools_scm',
    ],
    # Installation dependencies
    # Use with pip install . to install from source
    'install': [
        'cartopy',
        'cdo',
        'cf_units',
        'cython',
        # 'scitools-iris',  # Only iris 2 is on PyPI
        'matplotlib<3',
        'netCDF4',
        'numba',
        'numpy',
        'pillow',
        'prov[dot]',
        'psutil',
        'pyyaml',
        'shapely',
        'six',
        'stratify',
        'vmprof',
        'yamale',
    ],
    # Test dependencies
    # Execute 'python setup.py test' to run tests
    'test': [
        'easytest',
        # TODO: add dummydata package, see environment.yml
        'mock',
        'nose',
        'pycodestyle',
        'pytest>=3.9',
        'pytest-cov',
        'pytest-env',
        'pytest-html',
        'pytest-metadata>=1.5.1',
    ],
    # Development dependencies
    # Use pip install -e .[develop] to install in development mode
    'develop': [
        'isort',
        'prospector[with_pyroma]',
        'pycodestyle',
        'pydocstyle',
        'pylint',
        'sphinx',
        'sphinx_rtd_theme',
        'yamllint',
        'yapf',
    ],
}

if sys.version_info.major == 2:
    REQUIREMENTS['test'].append('more-itertools<6')
    for i, req in enumerate(REQUIREMENTS['install']):
        if req.startswith('cdo'):
            REQUIREMENTS['install'][i] = 'cdo!=1.5.*'


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
                if (filename.lower().endswith('.py')
                        and not _ignore(filename)):
                    yield filename


class CustomCommand(Command):
    """Custom Command class."""

    def install_deps_temp(self):
        """Try to temporarily install packages needed to run the command."""
        if self.distribution.install_requires:
            self.distribution.fetch_build_eggs(
                self.distribution.install_requires)
        if self.distribution.tests_require:
            self.distribution.fetch_build_eggs(self.distribution.tests_require)


class RunTests(CustomCommand):
    """Class to run tests and generate reports."""

    user_options = [('installation', None,
                     'Run tests that require installation.')]

    def initialize_options(self):
        """Initialize custom options."""
        self.installation = False

    def finalize_options(self):
        """Do nothing."""

    def run(self):
        """Run tests and generate a coverage report."""
        self.install_deps_temp()

        import pytest

        version = sys.version_info[0]
        report_dir = 'test-reports/python{}'.format(version)
        args = [
            'tests',
            'esmvaltool',  # for doctests
            '--ignore=esmvaltool/cmor/tables/',
            '--doctest-modules',
            '--cov=esmvaltool',
            '--cov-report=term',
            '--cov-report=html:{}/coverage_html'.format(report_dir),
            '--cov-report=xml:{}/coverage.xml'.format(report_dir),
            '--junit-xml={}/report.xml'.format(report_dir),
            '--html={}/report.html'.format(report_dir),
        ]
        if self.installation:
            args.append('--installation')
        errno = pytest.main(args)

        sys.exit(errno)


class RunLinter(CustomCommand):
    """Class to run a linter and generate reports."""

    user_options = []

    def initialize_options(self):
        """Do nothing."""

    def finalize_options(self):
        """Do nothing."""

    def run(self):
        """Run prospector and generate a report."""
        check_paths = PACKAGES + [
            'setup.py',
            'tests',
            'util',
        ]
        ignore = [
            'doc/',
        ]

        # try to install missing dependencies and import prospector
        try:
            from prospector.run import main
        except ImportError:
            # try to install and then import
            self.distribution.fetch_build_eggs(['prospector[with_pyroma]'])
            from prospector.run import main

        self.install_deps_temp()

        # run linter

        # change working directory to package root
        package_root = os.path.abspath(os.path.dirname(__file__))
        os.chdir(package_root)

        # write command line
        files = discover_python_files(check_paths, ignore)
        sys.argv = ['prospector']
        sys.argv.extend(files)

        # run prospector
        errno = main()

        sys.exit(errno)


with open('README.md') as readme:
    setup(
        name='ESMValTool',
        version=__version__,
        description='Earth System Models eValuation Tool',
        long_description=readme.read(),
        url='https://www.esmvaltool.org',
        download_url='https://github.com/ESMValGroup/ESMValTool',
        license='Apache License, Version 2.0',
        classifiers=[
            'Environment :: Console',
            'License :: OSI Approved :: Apache Software License',
            'Programming Language :: Python',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.6',
        ],
        packages=PACKAGES,
        # Include all version controlled files
        include_package_data=True,
        setup_requires=REQUIREMENTS['setup'],
        install_requires=REQUIREMENTS['install'],
        tests_require=REQUIREMENTS['test'],
        extras_require={
            'develop': REQUIREMENTS['develop'] + REQUIREMENTS['test']
        },
        entry_points={
            'console_scripts': [
                'esmvaltool = esmvaltool._main:run',
                'cmorize_obs = esmvaltool.'
                'utils.cmorizers.obs.cmorize_obs:execute_cmorize',
                'nclcodestyle = esmvaltool.'
                'utils.nclcodestyle.nclcodestyle:_main',
                'mip_convert_setup = esmvaltool.'
                'utils.cmorizers.mip_convert.esmvt_mipconv_setup:main'
            ],
        },
        cmdclass={
            'test': RunTests,
            'lint': RunLinter,
        },
        zip_safe=False,
    )
