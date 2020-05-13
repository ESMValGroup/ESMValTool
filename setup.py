#!/usr/bin/env python
"""ESMValTool installation script."""
import json
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
        'cdsapi',
        'cf-units',
        'cftime',
        'cmocean',
        'dask>=2.12',
        'ecmwf-api-client',
        'eofs',
        'ESMPy',
        'esmvalcore>=2.0.0,<2.1',
        'fiona',
        'GDAL',
        'jinja2',
        'joblib',
        'matplotlib',
        'nc-time-axis',  # needed by iris.plot
        'netCDF4',
        'numpy',
        'pandas',
        'pynio',
        'pyproj>=2.1'
        'pyyaml',
        'progressbar2',
        'scikit-learn',
        'scipy',
        'scitools-iris>=2.2.1',
        'seaborn',
        'seawater',
        'shapely',
        'xarray>=0.12',
        'xesmf',
        'xlsxwriter',
    ],
    # Test dependencies
    # Execute 'python setup.py test' to run tests
    'test': [
        'pytest>=3.9,!=6.0.0rc1,!=6.0.0',
        'pytest-cov',
        'pytest-env',
        'pytest-flake8',
        'pytest-html!=2.1.0',
        'pytest-metadata>=1.5.1',
    ],
    # Development dependencies
    # Use pip install -e .[develop] to install in development mode
    'develop': [
        'autodocsumm',
        'codespell',
        'isort',
        'prospector[with_pyroma]!=1.1.6.3,!=1.1.6.4',
        'sphinx>2',
        'sphinx_rtd_theme',
        'vmprof',
        'yamllint',
        'yapf',
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
                if (filename.lower().endswith('.py')
                        and not _ignore(filename)):
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
                self.distribution.install_requires)
        if self.distribution.tests_require:
            self.distribution.fetch_build_eggs(self.distribution.tests_require)

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


def read_authors(filename):
    """Read the list of authors from .zenodo.json file."""
    with Path(filename).open() as file:
        info = json.load(file)
        authors = []
        for author in info['creators']:
            name = ' '.join(author['name'].split(',')[::-1]).strip()
            authors.append(name)
        return ', '.join(authors)


def read_description(filename):
    """Read the description from .zenodo.json file."""
    with Path(filename).open() as file:
        info = json.load(file)
        return info['description']


setup(
    name='ESMValTool',
    version=__version__,
    author=read_authors('.zenodo.json'),
    description=read_description('.zenodo.json'),
    long_description=Path('README.md').read_text(),
    long_description_content_type='text/markdown',
    url='https://www.esmvaltool.org',
    download_url='https://github.com/ESMValGroup/ESMValTool',
    license='Apache License, Version 2.0',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: GIS',
        'Topic :: Scientific/Engineering :: Hydrology',
        'Topic :: Scientific/Engineering :: Physics',
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
        'test': REQUIREMENTS['test'],
    },
    entry_points={
        'console_scripts': [
            'cmorize_obs = esmvaltool.cmorizers.obs.cmorize_obs:main',
            'mip_convert_setup = '
            'esmvaltool.cmorizers.mip_convert.esmvt_mipconv_setup:main',
            'nclcodestyle = esmvaltool.utils.nclcodestyle.nclcodestyle:_main',
            'test_recipe = '
            'esmvaltool.utils.testing.recipe_settings.install_expand_run:main'
        ],
        'esmvaltool_commands': [
            'colortables = '
            'esmvaltool.utils.color_tables.show_color_tables:ColorTables',
            'install = esmvaltool.install:Install',
        ]
    },
    cmdclass={
        'lint': RunLinter,
    },
    zip_safe=False,
)
