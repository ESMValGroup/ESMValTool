#!/usr/bin/env python
"""ESMValTool installation script."""
import os
import re
import sys
from pathlib import Path

from setuptools import Command, setup

from esmvaltool import __version__

PACKAGES = [
    'esmvaltool',
]

REQUIREMENTS = {
    # Installation script (this file) dependencies
    'setup': [
        'pytest-runner',
        'setuptools_scm',
    ],
    # Installation dependencies
    # Use with pip install . to install from source
    'install': [
        'cartopy',
        'cdo',
        'cf-units',
        'cython',
        'ecmwf-api-client',
        'eofs',
        'esmvalcore>=2.0.0b3,<2.1',
        'fiona',
        'jinja2',
        'matplotlib<3',
        'nc-time-axis',  # needed by iris.plot
        'netCDF4',
        'numpy',
        'pandas',
        'pyyaml',
        'scitools-iris>=2.2',
        'scikit-learn',
        'seaborn',
        'shapely',
        'stratify',
        'xarray>=0.12',
        'xesmf',
        'xlrd',
        'xlsxwriter',
    ],
    # Test dependencies
    # Execute 'python setup.py test' to run tests
    'test': [
        'easytest',
        'mock',
        'nose',
        'pycodestyle',
        'pytest>=3.9',
        'pytest-cov',
        'pytest-env',
        'pytest-flake8',
        'pytest-html',
        'pytest-metadata>=1.5.1',
    ],
    # Development dependencies
    # Use pip install -e .[develop] to install in development mode
    'develop': [
        'isort',
        'prospector[with_pyroma]!=1.1.6.3,!=1.1.6.4',
        'sphinx',
        'sphinx_rtd_theme',
        'vmprof',
        'yamllint',
        'yapf',
    ],
}


def read_authors(citation_file):
    """Read the list of authors from .cff file."""
    authors = re.findall(
        r'family-names: (.*)$\s*given-names: (.*)',
        Path(citation_file).read_text(),
        re.MULTILINE,
    )
    return ', '.join(' '.join(author[::-1]) for author in authors)


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

        report_dir = 'test-reports'
        args = [
            'tests',
            'esmvaltool',  # for doctests
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


setup(
    name='ESMValTool',
    version=__version__,
    author=read_authors('CITATION.cff'),
    description='Earth System Models eValuation Tool',
    long_description=Path('README.md').read_text(),
    url='https://www.esmvaltool.org',
    download_url='https://github.com/ESMValGroup/ESMValTool',
    license='Apache License, Version 2.0',
    classifiers=[
        'Development Status :: 4 - Beta', 'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English', 'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: GIS',
        'Topic :: Scientific/Engineering :: Hydrology',
        'Topic :: Scientific/Engineering :: Physics'
    ],
    packages=PACKAGES,
    # Include all version controlled files
    include_package_data=True,
    setup_requires=REQUIREMENTS['setup'],
    install_requires=REQUIREMENTS['install'],
    tests_require=REQUIREMENTS['test'],
    extras_require={
        'develop': (set(REQUIREMENTS['develop'] + REQUIREMENTS['test']) -
                    {'pycodestyle'}),
    },
    entry_points={
        'console_scripts': [
            'cmorize_obs = esmvaltool.cmorizers.obs.cmorize_obs:main',
            'mip_convert_setup = '
            'esmvaltool.cmorizers.mip_convert.esmvt_mipconv_setup:main',
            'nclcodestyle = esmvaltool.utils.nclcodestyle.nclcodestyle:_main',
            'showcolortables = '
            'esmvaltool.utils.color_tables.show_color_tables:run',
            'test_recipe = '
            'esmvaltool.utils.testing.recipe_settings.install_expand_run:main'
        ],
    },
    cmdclass={
        'test': RunTests,
        'lint': RunLinter,
    },
    zip_safe=False,
)
