"""
Run the CMORization module as a utility executable.

This utility allows the user to call and execute CMOR reformatting
scripts (support for NCL and Python at the moment), that will use
two I/O variables passed by this utility: an input directory as
specified in config-user.yml by the RAWOBS key, and an output dir
created in the form of output_dir/CMOR_DATE_TIME/TierTIER/DATASET.
The user can specify a list of DATASETS that the CMOR reformatting
can by run on by using -o (--obs-list-cmorize) command line argument.
The CMOR reformatting scripts are to be found in:
esmvalcore.cmor/cmorizers/obs
"""
import datetime
import importlib
import logging
import os
import subprocess
import shutil
from pathlib import Path

import esmvalcore
from esmvalcore._config import configure_logging, read_config_user_file
from esmvalcore._task import write_ncl_settings

from esmvaltool.cmorizers.obs.utilities import read_cmor_config

logger = logging.getLogger(__name__)


class Formatter():

    def __init__(self,):
        self.datasets = []
        self.config = None

    def start(self, command, datasets, config_file, options):
        """
        Read configuration and set up formatter for data processing.

        Parameters
        ----------
        command: str
            Name of the command to execute
        datasets: str
            List of datasets to process, comma separated
        config_file: str
            Config file to use
        options: dict()
            Extra options to overwrite config user file
        """
        if datasets:
            self.datasets = datasets.split(',')

        self.config = read_config_user_file(
            config_file, f'data_{command}', options)

        if not os.path.isdir(self.run_dir):
            os.makedirs(self.run_dir)

        # configure logging
        log_files = configure_logging(
            output_dir=self.run_dir, console_log_level=self.log_level)
        logger.info("Writing program log files to:\n%s", "\n".join(log_files))

        # run
        timestamp1 = datetime.datetime.utcnow()
        timestamp_format = "%Y-%m-%d %H:%M:%S"

        logger.info("Starting the CMORization Tool at time: %s UTC",
                    timestamp1.strftime(timestamp_format))

        logger.info(70 * "-")
        logger.info("input_dir  = %s", self.rawobs)
        # check if the inputdir actually exists
        if not os.path.isdir(self.rawobs):
            logger.error("Directory %s does not exist", self.rawobs)
            raise ValueError
        logger.info("output_dir = %s", self.output_dir)
        logger.info(70 * "-")

    @property
    def rawobs(self):
        """Raw obs folder path"""
        return self.config["rootpath"]["RAWOBS"][0]

    @property
    def output_dir(self):
        """Output folder path"""
        return self.config['output_dir']

    @property
    def run_dir(self):
        """Run dir folder path"""
        return os.path.join(self.config['output_dir'], 'run')

    @property
    def log_level(self):
        """"Console log level"""
        return self.config['log_level']

    @staticmethod
    def _dataset_to_module(dataset):
        return dataset.lower().replace('-', '_')

    def download(self, start_date, end_date, overwrite):
        """
        Download all datasets.

        Parameters
        ----------
        start_date: datetime
            First date to download
        end_date: datetime
            Last date to download
        overwrite: boolean
            If True, download again existing files
        """
        if not self.datasets:
            logger.error('Missing datasets to download')
        logger.info("Downloading original data...")
        # master directory
        failed_datasets = []
        for dataset in self.datasets:
            try:
                self.download_dataset(
                    dataset, start_date, end_date, overwrite)
            except Exception:
                logger.exception('Failed to download %s', dataset)
                failed_datasets.append(dataset)
        if failed_datasets:
            logger.error('Download failed for datasets %s', self.datasets)

    def download_dataset(self, dataset, start_date, end_date, overwrite):
        """
        Download a single dataset.

        Parameters
        ----------
        dataset: str
            Dataset name
        start_date: datetime
            First date to download
        end_date: datetime
            Last date to download
        overwrite: boolean
            If True, download again existing files
        """
        dataset_module = self._dataset_to_module(dataset)
        logger.info('Downloading %s', dataset)
        logger.debug("Download module: %s", dataset_module)
        try:
            downloader = importlib.import_module(
                f'.{dataset_module}',
                package='esmvaltool.cmorizers.obs.downloaders.datasets'
            )
        except ImportError:
            logger.exception('Could not find cmorizer for %s', dataset)
            raise
        downloader.download_dataset(
            self.config, dataset, start_date, end_date, overwrite)
        logger.info('%s downloaded', dataset)

    def format(self, start, end, install):
        """Format all available datasets."""
        logger.info("Running the CMORization scripts.")
        # datsets dictionary of Tier keys
        datasets = self._assemble_datasets()
        if not datasets:
            logger.warning("Check input: could not find required %s in %s",
                           self.datasets, self.rawobs)
        logger.info("Processing datasets %s", datasets)

        # loop through tier/datasets to be cmorized
        failed_datasets = []
        for tier in datasets:
            for dataset in datasets:
                if not self.format_dataset(dataset, start, end, install):
                    failed_datasets.append(dataset)

            if failed_datasets:
                raise Exception(
                    f'Format failed for datasets {" ".join(failed_datasets)}'
                )

    def _assemble_datasets(self):
        """Get my datasets as dictionary keyed on Tier."""
        # check for desired datasets only (if any)
        # if not, walk all over rawobs dir
        # assume a RAWOBS/TierX/DATASET input structure

        # get all available tiers in source dir
        tiers = ['Tier{}'.format(i) for i in [2, 3]]
        tiers = [
            tier for tier in tiers
            if os.path.exists(os.path.join(self.rawobs, tier))
        ]
        datasets = []
        if self.datasets:
            return self.datasets
        else:
            for tier in tiers:
                for dataset in os.listdir(os.path.join(self.rawobs, tier)):
                    datasets.append(dataset)

        return datasets

    def format_dataset(self, dataset, start, end, install):
        """
        Format a single dataset.

        Parameters
        ----------
        dataset: str
            Dataset name
        """
        reformat_script_root = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'formatters',
            'datasets',
            self._dataset_to_module(dataset)
        )
        tier = self._get_dataset_tier(dataset)
        if tier is None:
            logger.error('Data for %s not found', dataset)
            return False

        # in-data dir; build out-dir tree
        in_data_dir = os.path.join(self.rawobs, tier, dataset)
        logger.info("Input data from: %s", in_data_dir)
        out_data_dir = os.path.join(self.output_dir, tier, dataset)
        logger.info("Output will be written to: %s", out_data_dir)
        if not os.path.isdir(out_data_dir):
            os.makedirs(out_data_dir)

        # all operations are done in the working dir now
        os.chdir(out_data_dir)
        # figure out what language the script is in
        logger.info("Reformat script: %s", reformat_script_root)
        if os.path.isfile(reformat_script_root + '.ncl'):
            reformat_script = reformat_script_root + '.ncl'
            success = self._run_ncl_script(
                in_data_dir,
                out_data_dir,
                dataset,
                reformat_script,
                start,
                end
            )
        elif os.path.isfile(reformat_script_root + '.py'):
            success = self._run_pyt_script(
                in_data_dir, out_data_dir, dataset, start, end)
        else:
            logger.error('Could not find formatter for %s', dataset)
            return False
        if not success:
            logger.error('Foramtting failed for datset %s', dataset)
            return False
        if install:
            rootpath = self.config['rootpath']
            target_dir = rootpath.get('OBS', rootpath['default'])[0]
            target_dir = os.path.join(target_dir, tier, dataset)
            if os.path.isdir(target_dir):
                logger.info(
                    'Automatic installation of dataset %s skipped: '
                    'target folder %s already exists',
                    dataset, target_dir)
            else:
                logger.info(
                    'Installing dataset %s in folder %s', dataset, target_dir)
                shutil.copytree(out_data_dir, target_dir)
        return True

    def _get_dataset_tier(self, dataset):
        for tier in [2, 3]:
            if os.path.isdir(
                    os.path.join(self.rawobs, f"Tier{tier}", dataset)):
                return f"Tier{tier}"
        return None

    def _write_ncl_settings(self, project_info, dataset, run_dir,
                            reformat_script, start_year, end_year):
        """Write the information needed by the ncl reformat script."""
        if start_year is None:
            start_year = 0
        else:
            start_year = start_year.year
        if end_year is None:
            end_year = 0
        else:
            end_year = end_year.year
        settings = {
            'cmorization_script': reformat_script,
            'input_dir_path': project_info[dataset]['indir'],
            'output_dir_path': project_info[dataset]['outdir'],
            'config_user_info': {
                'log_level': self.config['log_level'],
            },
            'start_year': start_year,
            'end_year': end_year,
        }
        settings_filename = os.path.join(run_dir, dataset, 'settings.ncl')
        if not os.path.isdir(os.path.join(run_dir, dataset)):
            os.makedirs(os.path.join(run_dir, dataset))
        # write the settings file
        write_ncl_settings(settings, settings_filename)
        return settings_filename

    def _run_ncl_script(self, in_dir, out_dir, dataset, script, start, end):
        """Run the NCL cmorization mechanism."""
        logger.info("CMORizing dataset %s using NCL script %s",
                    dataset, script)
        project = {}
        project[dataset] = {}
        project[dataset]['indir'] = in_dir
        project[dataset]['outdir'] = out_dir
        settings_file = self._write_ncl_settings(
            project, dataset, self.run_dir, script, start, end)
        esmvaltool_root = os.path.dirname(
            os.path.dirname(os.path.dirname(os.path.dirname(script))))

        # put settings in environment
        env = dict(os.environ)
        env['settings'] = settings_file
        env['esmvaltool_root'] = esmvaltool_root
        env['cmor_tables'] = str(
            Path(esmvalcore.cmor.__file__).parent / 'tables')
        logger.info("Using CMOR tables at %s", env['cmor_tables'])
        # call NCL
        ncl_call = ['ncl', script]
        logger.info("Executing cmd: %s", ' '.join(ncl_call))
        process = subprocess.Popen(ncl_call,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT,
                                   env=env)
        output, err = process.communicate()
        for oline in str(output.decode('utf-8')).split('\n'):
            logger.info('[NCL] %s', oline)
        if err:
            logger.error('[NCL][subprocess.Popen ERROR] %s', err)
            return False
        return True

    def _run_pyt_script(self, in_dir, out_dir, dataset, start, end):
        """Run the Python cmorization mechanism."""
        module_name = 'esmvaltool.cmorizers.obs.formatters.datasets.{}'.format(
            dataset.lower().replace("-", "_"))
        module = importlib.import_module(module_name)
        logger.info("CMORizing dataset %s using Python script %s",
                    dataset, module.__file__)
        cmor_cfg = read_cmor_config(dataset)
        module.cmorization(in_dir, out_dir, cmor_cfg, self.config, start, end)
        return True


class DataCommand():
    """Download and format data to use with ESMValTool"""

    def download(self, datasets, config_file=None, start=None, end=None,
                 overwrite=False, **kwargs):
        """
        Download datasets.

        Parameters
        ----------
        datasets

        config_file
        """
        start = self._parse_date(start)
        end = self._parse_date(end)

        formatter = Formatter()
        formatter.start('download', datasets, config_file, kwargs)
        formatter.download(start, end, overwrite)

    def format(self, datasets=None, config_file=None, start=None, end=None,
               install=False, **kwargs):
        start = self._parse_date(start)
        end = self._parse_date(end)

        formatter = Formatter()
        formatter.start('download', datasets, config_file, kwargs)
        formatter.format(start, end, install)

    def prepare(self, datasets=None, config_file=None,
                start=None, end=None,
                overwrite=False, install=False, **kwargs):
        """
        Download a format a set of datasets.

        Parameters
        ----------
        datasets
        """
        start = self._parse_date(start)
        end = self._parse_date(end)

        formatter = Formatter()
        formatter.start('download', datasets, config_file, kwargs)
        formatter.download(start, end, overwrite)
        formatter.format(start, end, install)

    @staticmethod
    def _parse_date(date):
        if date is None:
            date = ''
        date = str(date)
        if len(date) == 8:
            return datetime.datetime.strptime(date, "%Y%m%d")
        if len(date) == 6:
            return datetime.datetime.strptime(date, "%Y%m")
        if len(date) == 4:
            return datetime.datetime.strptime(date, "%Y")
