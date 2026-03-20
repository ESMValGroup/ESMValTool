"""Download and formatting of non-ESGF datasets.

This module adds new commands to the ESMValTool to allow the user to get
and reformat to the ESMValTool's data format a set of observations and
reanalysis.
"""

from __future__ import annotations

import datetime
import importlib
import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import esmvalcore
import yaml
from esmvalcore._task import write_ncl_settings
from esmvalcore.config import CFG
from esmvalcore.config._dask import get_distributed_client
from esmvalcore.config._logging import configure_logging

from esmvaltool.cmorizers.data.utilities import read_cmor_config

if TYPE_CHECKING:
    from esmvalcore.config import Session

    from esmvaltool.cmorizers.data.typing import DatasetInfo

logger = logging.getLogger(__name__)


DATASETS_FILE = Path(__file__).parent / "datasets.yml"


class _Formatter:
    """
    Class to manage the download and formatting of datasets.

    Parameters
    ----------
    info : dict
        Datasets information
    """

    def __init__(
        self,
        info: dict[Literal["datasets"], dict[str, DatasetInfo]],
    ) -> None:
        self.datasets: list[str] = []
        self.datasets_info = info
        self.config: Session | None = None

    def start(
        self,
        command: str,
        datasets: str | list[str],
        original_data_dir: Path | str | None,
        config_dir: Path | None,
        options: dict,
    ) -> None:
        """Read configuration and set up formatter for data processing.

        Parameters
        ----------
        command:
            Name of the command to execute.
        datasets:
            List of datasets to process, comma separated.
        original_data_dir:
            Directory containing the original data.
        config_dir:
            Config directory to use.
        options:
            Extra options to overwrite configuration.

        """
        if isinstance(datasets, str):
            self.datasets = datasets.split(",")
        else:
            self.datasets = datasets

        if config_dir is not None:
            config_dir = (
                Path(os.path.expandvars(config_dir)).expanduser().absolute()
            )
            if not config_dir.is_dir():
                msg = (
                    f"Invalid --config_dir given: {config_dir} is not an "
                    f"existing directory",
                )
                raise NotADirectoryError(msg)
            CFG.update_from_dirs([config_dir])
        CFG.nested_update(options)
        self.config = CFG.start_session(f"data_{command}")
        self.run_dir.mkdir(parents=True, exist_ok=True)

        # configure logging
        log_files = configure_logging(
            output_dir=self.run_dir,
            console_log_level=self.log_level,
        )
        logger.info("Writing program log files to:\n%s", "\n".join(log_files))

        # Locate the input data.
        if original_data_dir is None:
            original_data_dir = Path.cwd()
            # TODO: remove the lines below in ESMValTool v2.16.0.
            rawobs = self.config.get("rootpath", {}).get("RAWOBS", None)
            if rawobs is not None:
                logger.warning(
                    (
                        "Using the 'rootpath: RAWOBS' setting to specify the "
                        "input data directory is deprecated and this will stop "
                        "working in ESMValTool v2.16.0. Please use the "
                        "'--original-data-dir' argument instead."
                    ),
                )
                original_data_dir = Path(rawobs[0])

        self.original_data_dir = Path(
            os.path.expandvars(original_data_dir),
        ).expanduser()

        # run
        timestamp1 = datetime.datetime.now(datetime.UTC)
        timestamp_format = "%Y-%m-%d %H:%M:%S"

        logger.info(
            "Starting the CMORization Tool at time: %s UTC",
            timestamp1.strftime(timestamp_format),
        )

        logger.info(70 * "-")
        logger.info("input_dir  = %s", self.original_data_dir)
        # check if the inputdir actually exists
        if not self.original_data_dir.is_dir():
            msg = (
                f"The path '{self.original_data_dir}' is not a directory. "
                "Please specify the correct path to the input data using the "
                "--original-data-dir flag."
            )
            raise NotADirectoryError(msg)
        logger.info("output_dir = %s", self.output_dir)
        logger.info(70 * "-")

    @property
    def output_dir(self):
        """Output folder path."""
        return self.config.session_dir

    @property
    def run_dir(self):
        """Run dir folder path."""
        return self.config.run_dir

    @property
    def log_level(self):
        """Console log level."""
        return self.config["log_level"]

    @staticmethod
    def _dataset_to_module(dataset: str) -> str:
        return dataset.lower().replace("-", "_")

    def download(
        self,
        start_date: datetime.datetime | None,
        end_date: datetime.datetime | None,
        *,
        overwrite: bool,
    ) -> bool:
        """Download all datasets.

        Parameters
        ----------
        start_date:
            First date to download
        end_date:
            Last date to download
        overwrite:
            If True, download again existing files

        Returns
        -------
        :
            :obj:`True` if all datasets were downloaded, :obj:`False` otherwise.

        """
        if not self.datasets:
            logger.error("Missing datasets to download")
        logger.info("Downloading original data...")
        # master directory
        failed_datasets = []
        for dataset in self.datasets:
            try:
                self.download_dataset(
                    dataset,
                    start_date,
                    end_date,
                    overwrite=overwrite,
                )
            except ValueError:
                logger.exception("Failed to download %s", dataset)
                failed_datasets.append(dataset)
        if failed_datasets:
            logger.error("Download failed for datasets %s", failed_datasets)
            return False
        return True

    def download_dataset(
        self,
        dataset: str,
        start_date: datetime.datetime | None,
        end_date: datetime.datetime | None,
        *,
        overwrite: bool,
    ) -> None:
        """Download a single dataset.

        Parameters
        ----------
        dataset:
            Dataset name
        start_date:
            First date to download
        end_date:
            Last date to download
        overwrite:
            If True, download again existing files
        """
        if not self.has_downloader(dataset):
            msg = f"Dataset {dataset} does not have an automatic downloader"
            raise ValueError(msg)
        dataset_module = self._dataset_to_module(dataset)
        logger.info("Downloading %s", dataset)
        logger.debug("Download module: %s", dataset_module)
        try:
            downloader = importlib.import_module(
                f".{dataset_module}",
                package="esmvaltool.cmorizers.data.downloaders.datasets",
            )
        except ImportError:
            logger.exception("Could not find cmorizer for %s", dataset)
            raise

        downloader.download_dataset(
            original_data_dir=self.original_data_dir,
            dataset=dataset,
            dataset_info=self.datasets_info["datasets"][dataset],
            start_date=start_date,
            end_date=end_date,
            overwrite=overwrite,
        )
        logger.info("%s downloaded", dataset)

    def format(
        self,
        start: datetime.datetime | None,
        end: datetime.datetime | None,
        *,
        install: bool,
    ) -> None:
        """Format all available datasets.

        Parameters
        ----------
        start:
            Start of the period to format
        end:
            End of the period to format
        install:
            If True, automatically moves the data to the final location if
            there is no
        """
        logger.info("Running the CMORization scripts.")
        # datasets dictionary of Tier keys
        datasets = self._assemble_datasets()
        if not datasets:
            logger.warning(
                "Check input: could not find required %s in %s",
                self.datasets,
                self.original_data_dir,
            )
        logger.info("Processing datasets %s", datasets)

        with get_distributed_client():
            # loop through tier/datasets to be cmorized
            failed_datasets = [
                dataset
                for dataset in datasets
                if not self.format_dataset(
                    dataset,
                    start,
                    end,
                    install=install,
                )
            ]

        if failed_datasets:
            msg = f"Format failed for datasets {' '.join(failed_datasets)}"
            raise RuntimeError(msg)

    @staticmethod
    def has_downloader(dataset: str) -> bool:
        """Check if a given datasets has an automatic downloader.

        Parameters
        ----------
        dataset : str
            Name of the dataset to check

        Returns
        -------
        str
            'Yes' if the downloader exists, 'No' otherwise
        """
        try:
            importlib.import_module(
                f".{dataset.lower().replace('-', '_')}",
                package="esmvaltool.cmorizers.data.downloaders.datasets",
            )
        except ImportError:
            return False
        return True

    def _assemble_datasets(self) -> list[str]:
        """Get the datasets to CMORize."""
        # check for desired datasets only (if any)
        if self.datasets:
            return self.datasets
        # if not, look in configuration file and at the available data.
        return [
            dataset
            for dataset, info in self.datasets_info["datasets"].items()
            if (
                self.original_data_dir / f"Tier{info['tier']}" / dataset
            ).is_dir()
        ]

    def format_dataset(
        self,
        dataset: str,
        start: datetime.datetime | None,
        end: datetime.datetime | None,
        *,
        install: bool,
    ) -> bool:
        """Format a single dataset.

        Parameters
        ----------
        dataset:
            Dataset name
        start:
            Start of the period to format
        end:
            End of the period to format
        install:
            If True, automatically moves the data to the final location if
            there is no data there.
        """
        reformat_script_root = Path(
            Path(__file__).parent.resolve(),
            "formatters",
            "datasets",
            self._dataset_to_module(dataset),
        )
        tier = self._get_dataset_tier(dataset)

        # in-data dir; build out-dir tree
        in_data_dir = self.original_data_dir / tier / dataset
        logger.info("Input data from: %s", in_data_dir)
        if not in_data_dir.is_dir():
            msg = (
                f"Data for dataset '{dataset}' not found. "
                f"Path to original data '{in_data_dir}' is not a directory'"
            )
            raise NotADirectoryError(msg)
        out_data_dir = self.output_dir / tier / dataset
        logger.info("Output will be written to: %s", out_data_dir)
        if not out_data_dir.is_dir():
            out_data_dir.mkdir(parents=True)

        # all operations are done in the working dir now
        os.chdir(out_data_dir)
        # figure out what language the script is in
        logger.info("Reformat script: %s", reformat_script_root)
        if reformat_script_root.with_suffix(".ncl").is_file():
            reformat_script = reformat_script_root.with_suffix(".ncl")
            success = self._run_ncl_script(
                in_data_dir,
                out_data_dir,
                dataset,
                reformat_script,
                start,
                end,
            )
        elif reformat_script_root.with_suffix(".py").is_file():
            success = self._run_pyt_script(
                in_data_dir,
                out_data_dir,
                dataset,
                start,
                end,
            )
        else:
            logger.error("Could not find formatter for %s", dataset)
            return False
        if success:
            logger.info("Formatting successful for dataset %s", dataset)
        else:
            logger.error("Formatting failed for dataset %s", dataset)
            return False
        if install:
            target_dir = self._get_install_dir(dataset)
            if target_dir is None:
                logger.warning(
                    "Unable determine install path for dataset '%s', please "
                    "move files from %s to the desired location manually.",
                    dataset,
                    out_data_dir,
                )
            elif target_dir.is_dir():
                logger.info(
                    "Automatic installation of dataset %s skipped: "
                    "target folder %s already exists",
                    dataset,
                    target_dir,
                )
            else:
                logger.info(
                    "Installing dataset %s in folder %s",
                    dataset,
                    target_dir,
                )
                shutil.move(out_data_dir, target_dir)
        return True

    def _get_install_dir(self, dataset) -> Path | None:
        """Get the installation directory for a dataset.

        Parameters
        ----------
        dataset:
            Dataset name.

        Returns
        -------
        :
            Path to the installation directory, or None if it cannot be
            determined.
        """
        tier = self._get_dataset_tier(dataset)
        if "rootpath" in self.config:
            # May not be available since ESMValCore v2.14, will be removed in
            # ESMValCore v2.16.
            rootpath = self.config["rootpath"]
            target_dir = rootpath.get("OBS", rootpath["default"])[0]
            return Path(target_dir, tier, dataset)

        try:
            # Only available since ESMValCore v2.14
            import esmvalcore.io
            import esmvalcore.io.local
        except ImportError:
            return None

        target_dir = None
        # Normalize the attributes so they can be used as facets in the
        # directory name template of a LocalDataSource.
        try:
            attributes = read_cmor_config(dataset)["attributes"]
        except FileNotFoundError:
            # NCL scripts do not use a cmor config file.
            return None
        for attr in ["dataset", "project"]:
            if attr not in attributes:
                attributes[attr] = attributes[f"{attr}_id"]
        # Load data sources from configuration.
        data_sources = esmvalcore.io.load_data_sources(
            self.config,
            project=attributes["project"],
        )
        # Loop over potential target directories and try if the right attributes
        # are available to format the directory name template. Use the first
        # one that works.
        for data_source in sorted(
            data_sources,
            key=lambda ds: ds.priority,
        ):
            if isinstance(
                data_source,
                esmvalcore.io.local.LocalDataSource,
            ):
                try:
                    target_dir = (
                        data_source.rootpath
                        / data_source.dirname_template.format(
                            **attributes,
                        )
                    )
                except KeyError:
                    pass
                else:
                    break
        return target_dir

    def _get_dataset_tier(self, dataset: str) -> str:
        return f"Tier{self.datasets_info['datasets'][dataset]['tier']}"

    def _write_ncl_settings(
        self,
        project_info: dict[str, dict[str, str]],
        dataset: str,
        run_dir: Path,
        reformat_script: Path,
        start: datetime.datetime | None,
        end: datetime.datetime | None,
    ) -> Path:
        """Write the information needed by the ncl reformat script."""
        start_year = 0 if start is None else start.year
        end_year = 0 if end is None else end.year
        settings = {
            "cmorization_script": reformat_script,
            "input_dir_path": project_info[dataset]["indir"],
            "output_dir_path": project_info[dataset]["outdir"],
            "config_user_info": {
                "log_level": self.config["log_level"],
            },
            "start_year": start_year,
            "end_year": end_year,
        }
        settings_filename = run_dir / dataset / "settings.ncl"
        (run_dir / dataset).mkdir(parents=True, exist_ok=True)
        # write the settings file
        write_ncl_settings(settings, settings_filename)
        return settings_filename

    def _run_ncl_script(
        self,
        in_dir: Path,
        out_dir: Path,
        dataset: str,
        script: Path,
        start: datetime.datetime | None,
        end: datetime.datetime | None,
    ) -> bool:
        """Run the NCL cmorization mechanism."""
        logger.info(
            "CMORizing dataset %s using NCL script %s",
            dataset,
            script,
        )
        project = {}
        project[dataset] = {}
        project[dataset]["indir"] = str(in_dir)
        project[dataset]["outdir"] = str(out_dir)
        settings_file = self._write_ncl_settings(
            project_info=project,
            dataset=dataset,
            run_dir=self.run_dir,
            reformat_script=script,
            start=start,
            end=end,
        )

        # put settings in environment
        env = dict(os.environ)
        env["settings"] = str(settings_file)
        env["esmvaltool_root"] = str(script.parents[3])
        env["cmor_tables"] = str(
            Path(esmvalcore.cmor.__file__).parent / "tables",
        )
        logger.info("Using CMOR tables at %s", env["cmor_tables"])
        # call NCL
        ncl_call = ["ncl", str(script)]
        logger.info("Executing cmd: %s", " ".join(ncl_call))
        with subprocess.Popen(
            ncl_call,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env=env,
        ) as process:
            output, err = process.communicate()
        for oline in str(output.decode("utf-8")).split("\n"):
            logger.info("[NCL] %s", oline)
        if err:
            logger.error("[NCL][subprocess.Popen ERROR] %s", err)
            return False
        return True

    def _run_pyt_script(
        self,
        in_dir: str,
        out_dir: str,
        dataset: str,
        start: datetime.datetime | None,
        end: datetime.datetime | None,
    ) -> Literal[True]:
        """Run the Python cmorization mechanism."""
        module_name = (
            "esmvaltool.cmorizers.data.formatters.datasets."
            + dataset.lower().replace("-", "_")
        )
        module = importlib.import_module(module_name)
        logger.info(
            "CMORizing dataset %s using Python script %s",
            dataset,
            module.__file__,
        )
        cmor_cfg = read_cmor_config(dataset)
        module.cmorization(in_dir, out_dir, cmor_cfg, self.config, start, end)
        logger.info("CMORization of dataset %s finished!", dataset)
        return True


class DataCommand:
    """Download and format data to use with ESMValTool."""

    def __init__(self) -> None:
        self._info = yaml.safe_load(DATASETS_FILE.read_text(encoding="utf8"))
        self.formatter = _Formatter(self._info)

    def _has_downloader(self, dataset: str) -> Literal["Yes", "No"]:
        return "Yes" if self.formatter.has_downloader(dataset) else "No"

    def list(self) -> None:
        """List all supported datasets."""
        print()
        print(f"| {'Dataset name':30} | Tier | Auto-download | Last access |")
        print("-" * 71)
        for dataset, dataset_info in self._info["datasets"].items():
            date = datetime.datetime.strptime(
                str(dataset_info["last_access"]),
                "%Y-%m-%d",
            )
            print(
                f"| {dataset:30} | {dataset_info['tier']:4} "
                f"| {self._has_downloader(dataset):13} "
                f"|  {date.strftime('%Y-%m-%d')} |",
            )
        print("-" * 71)

    def info(self, dataset: str) -> None:
        """Show detailed info about a specific dataset.

        Parameters
        ----------
        dataset:
            dataset to show
        """
        dataset_info = self._info["datasets"][dataset]
        print(dataset)
        print()
        print(f"Tier: {dataset_info['tier']}")
        print(f"Source: {dataset_info['source']}")
        print(f"Automatic download: {self._has_downloader(dataset)}")
        print()
        print(dataset_info["info"])

    def download(
        self,
        datasets: str | list[str],
        *,
        original_data_dir: Path | str | None = None,
        start: str | None = None,
        end: str | None = None,
        overwrite: bool = False,
        config_dir: Path | None = None,
        **kwargs,
    ) -> None:
        """Download datasets.

        Parameters
        ----------
        datasets: list(str)
            List of datasets to format
        original_data_dir:
            Directory where original data will be stored.
        start:
            Start of the interval to process, by default None. Valid formats
            are YYYY, YYYYMM and YYYYMMDD.
        end:
            End of the interval to process, by default None. Valid formats
            are YYYY, YYYYMM and YYYYMMDD.
        overwrite:
            If true, download already present data again
        config_dir:
            Path to additional ESMValTool configuration directory. See
            https://docs.esmvaltool.org/projects/ESMValCore/en/latest/quickstart/configure.html#yaml-files
            for details.

        """
        start_date = self._parse_date(start)
        end_date = self._parse_date(end)

        self.formatter.start(
            "download",
            datasets=datasets,
            original_data_dir=original_data_dir,
            config_dir=config_dir,
            options=kwargs,
        )
        self.formatter.download(start_date, end_date, overwrite=overwrite)

    def format(
        self,
        datasets: str | list[str],
        *,
        original_data_dir: Path | str | None = None,
        start: str | None = None,
        end: str | None = None,
        install: bool = False,
        config_dir: Path | None = None,
        **kwargs,
    ) -> None:
        """Format datasets.

        Parameters
        ----------
        datasets:
            List of datasets to format
        original_data_dir:
            Directory where original data is stored.
        start:
            Start of the interval to process, by default None. Valid formats
            are YYYY, YYYYMM and YYYYMMDD.
        end:
            End of the interval to process, by default None. Valid formats
            are YYYY, YYYYMM and YYYYMMDD.
        install:
            If true, move processed data to the folder, by default False
        config_dir:
            Path to additional ESMValTool configuration directory. See
            https://docs.esmvaltool.org/projects/ESMValCore/en/latest/quickstart/configure.html#yaml-files
            for details.

        """
        start_date = self._parse_date(start)
        end_date = self._parse_date(end)

        self.formatter.start(
            "formatting",
            datasets=datasets,
            original_data_dir=original_data_dir,
            config_dir=config_dir,
            options=kwargs,
        )
        self.formatter.format(start_date, end_date, install=install)

    def prepare(
        self,
        datasets: str | list[str],
        *,
        original_data_dir: Path | str | None = None,
        start: str | None = None,
        end: str | None = None,
        overwrite: bool = False,
        install: bool = False,
        config_dir: Path | None = None,
        **kwargs,
    ) -> None:
        """Download and format a set of datasets.

        Parameters
        ----------
        datasets:
            List of datasets to format
        original_data_dir:
            Directory where original data is stored or will be stored after download.
        start:
            Start of the interval to process, by default None. Valid formats
            are YYYY, YYYYMM and YYYYMMDD.
        end:
            End of the interval to process, by default None. Valid formats
            are YYYY, YYYYMM and YYYYMMDD.
        install:
            If true, move processed data to the folder, by default False
        overwrite:
            If true, download already present data again
        config_dir:
            Path to additional ESMValTool configuration directory. See
            https://docs.esmvaltool.org/projects/ESMValCore/en/latest/quickstart/configure.html#yaml-files
            for details.

        """
        start_date = self._parse_date(start)
        end_date = self._parse_date(end)

        self.formatter.start(
            "preparation",
            datasets=datasets,
            original_data_dir=original_data_dir,
            config_dir=config_dir,
            options=kwargs,
        )
        if self.formatter.download(start_date, end_date, overwrite=overwrite):
            self.formatter.format(start_date, end_date, install=install)
        else:
            logger.warning("Download failed, skipping format step")

    @staticmethod
    def _parse_date(date: str | None) -> datetime.datetime | None:
        if date is None:
            return None
        date_string = str(date)
        date_formats = {
            4: "%Y",
            6: "%Y%m",
            8: "%Y%m%d",
        }
        format_string = date_formats.get(len(date_string), None)
        if format_string is None:
            msg = (
                f"Unsupported date format for {date}. "
                'Supported formats for "start" and "end" are: '
                '"None", "YYYY", "YYYYMM", "YYYYMMDD"',
            )
            raise ValueError(msg)
        return datetime.datetime.strptime(date_string, format_string)
