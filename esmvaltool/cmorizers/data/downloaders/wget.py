"""wget based downloader."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING

from .downloader import BaseDownloader

if TYPE_CHECKING:
    from esmvaltool.cmorizers.data.typing import DatasetInfo

logger = logging.getLogger(__name__)


class WGetDownloader(BaseDownloader):
    """Data downloader based on wget."""

    def download_folder(self, server_path, wget_options):
        """Download folder.

        Parameters
        ----------
        server_path: str
            Path to remote folder
        wget_options: list(str)
            Extra options for wget
        """
        if self.overwrite:
            raise ValueError(
                "Overwrite does not work with downloading directories through "
                "wget. Please, remove the unwanted data manually",
            )
        command = (
            ["wget"]
            + wget_options
            + [
                "--no-clobber",
                f"--directory-prefix={self.local_folder}",
                "--recursive",
                "--no-directories",
                f"{server_path}",
            ]
        )
        logger.debug(command)
        subprocess.check_output(command)

    def download_file(self, server_path, wget_options, output_filename=None):
        """Download file.

        Parameters
        ----------
        server_path: str
            Path to remote file
        wget_options: list(str)
            Extra options for wget
        output_filename: str, optional
            Name of the downloaded file. If not given, use the one given by
            ``server_path``.

        """
        output_dir = Path(self.local_folder)
        output_dir.mkdir(parents=True, exist_ok=True)
        if output_filename is None:
            output_path = output_dir / Path(server_path).name
        else:
            output_path = output_dir / output_filename
        output_options = []

        # If no specific output filename is desired (i.e., the option -O can be
        # omitted), wget can be used with the --no-clobber and
        # --directory-prefix options to avoid overwriting data. Otherwise, we
        # will need to check file existence manually here (-O and --no-clobber
        # do not work well together).
        if not self.overwrite and output_filename is None:
            output_options.append(f"--directory-prefix={output_dir!s}")
            output_options.append("--no-clobber")
        else:
            if (
                not self.overwrite
                and output_filename is not None
                and output_path.exists()
            ):
                logger.info("File %s exists, skipping download", output_path)
                return
            output_options.append("-O")
            output_options.append(str(output_path))

        command = (
            ["wget"]
            + wget_options
            + output_options
            + ["--no-directories", server_path]
        )
        logger.debug(command)
        subprocess.check_output(command)

    def login(self, server_path, wget_options):
        """Login.

        Parameters
        ----------
        server_path: str
            Path to remote file
        wget_options: list(str)
            Extra options for wget
        """
        command = ["wget"] + wget_options + [server_path]
        logger.debug(command)
        subprocess.check_output(command)


class NASADownloader(WGetDownloader):
    """Downloader for the NASA repository."""

    def __init__(
        self,
        original_data_dir: Path,
        dataset: str,
        dataset_info: DatasetInfo,
        *,
        overwrite: bool,
    ) -> None:
        super().__init__(
            original_data_dir=original_data_dir,
            dataset=dataset,
            dataset_info=dataset_info,
            overwrite=overwrite,
        )

        self._wget_common_options = [
            "--load-cookies=~/.urs_cookies",
            "--save-cookies=~/.urs_cookies",
            "--auth-no-challenge=on",
            "--keep-session-cookies",
            "--no-check-certificate",
        ]

    def download_folder(self, server_path, wget_options=None):
        """Download folder.

        Parameters
        ----------
        server_path: str
            Path to remote folder
        wget_options: list(str)
            Extra options for wget, by default None
        """
        if wget_options is None:
            wget_options = []
        wget_options = (
            self._wget_common_options
            + ["-np", "--accept=nc,nc4,hdf"]
            + wget_options
        )
        super().download_folder(server_path, wget_options)

    def download_file(self, server_path, wget_options=None):
        """Download file.

        Parameters
        ----------
        server_path: str
            Path to remote folder
        wget_options: list(str)
            Extra options for wget, by default None
        """
        if wget_options is None:
            wget_options = []
        super().download_file(
            server_path,
            self._wget_common_options + wget_options,
        )
