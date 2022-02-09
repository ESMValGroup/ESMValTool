"""wget based downloader."""

import logging
import os
import subprocess

from .downloader import BaseDownloader

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
                'Overwrite does not work with downloading directories through '
                'wget. Please, remove the unwanted data manually')
        command = ['wget'] + wget_options + self.overwrite_options + [
            f'--directory-prefix={self.local_folder}',
            '--recursive',
            '--no-directories',
            f'{server_path}',
        ]
        logger.debug(command)
        subprocess.check_output(command)

    def download_file(self, server_path, wget_options):
        """Download file.

        Parameters
        ----------
        server_path: str
            Path to remote file
        wget_options: list(str)
            Extra options for wget
        """
        command = ['wget'] + wget_options + self.overwrite_options + [
            f'--directory-prefix={self.local_folder}',
            '--no-directories',
            server_path,
        ]
        if self.overwrite:
            command.append(f'-O {os.path.basename(server_path)}')
        logger.debug(command)
        subprocess.check_output(command)

    @property
    def overwrite_options(self):
        """Get overwrite options as configured in downloader."""
        if not self.overwrite:
            return [
                '--no-clobber',
            ]
        return []


class NASADownloader(WGetDownloader):
    """Downloader for the NASA repository."""
    def __init__(self, config, dataset, dataset_info, overwrite):
        super().__init__(config, dataset, dataset_info, overwrite)

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
        wget_options = self._wget_common_options + ["-np", "--accept=nc,nc4"
                                                    ] + wget_options
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
        super().download_file(server_path,
                              self._wget_common_options + wget_options)
