"""Downloader for FTP repositories."""

import ftplib
import logging
import os
import re

from progressbar import (
    ETA,
    Bar,
    DataSize,
    FileTransferSpeed,
    Percentage,
    ProgressBar,
)

from .downloader import BaseDownloader

logger = logging.getLogger(__name__)


class FTPDownloader(BaseDownloader):
    """Downloader for FTP repositories.

    Parameters
    ----------
    config : dict
        ESMValTool's user configuration
    server : str
        FTP server URL
    dataset : str
        Dataset to download
    dataset_info : dict
        Dataset information from the datasets.yml file
    overwrite : bool
        Overwrite already downloaded files
    """
    def __init__(self, config, server, dataset, dataset_info, overwrite):
        super().__init__(config, dataset, dataset_info, overwrite)
        self._client = None
        self.server = server

    def connect(self):
        """Connect to the FTP server."""
        self._client = ftplib.FTP(self.server)
        logger.info(self._client.getwelcome())
        self._client.login()

    def set_cwd(self, path):
        """Set current working directory in the remote.

        Parameters
        ----------
        path : str
            Remote path to set as current working directory.
        """
        logger.debug('Current working directory: %s', self._client.pwd())
        logger.debug('Setting working directory to %s', path)
        self._client.cwd(path)
        logger.debug('New working directory: %s', self._client.pwd())

    def list_folders(self, server_path='.'):
        """List folder in the remote.

        Parameters
        ----------
        server_path : str, optional
            Folder to list, by default '.'

        Returns
        -------
        list(str)
            List of folder names
        """
        filenames = self._client.mlsd(server_path, facts=['type'])
        return [
            filename for filename, facts in filenames if facts['type'] == 'dir'
        ]

    def exists(self, server_path):
        """Check if a given path exists in the server.

        Parameters
        ----------
        server_path : str
            Path to check for existence.
        """
        return server_path in self._client.nlst()

    def download_folder(self, server_path, sub_folder='', filter_files=None):
        """Download files from a given folder.

        Parameters
        ----------
        server_path : str
            Folder to download
        sub_folder : str, optional
            Name of the local subfolder to store the results in, by default ''
        filter_files : str, optional
            If set, only download files that match this regular expression,
            by default None
        """
        # get filenames within the directory
        filenames = self._client.nlst(server_path)
        logger.info('Downloading files in %s', server_path)
        if filter_files:
            expression = re.compile(filter_files)
            filenames = [
                filename for filename in filenames
                if expression.match(os.path.basename(filename))
            ]
        for filename in filenames:
            self.download_file(filename, sub_folder)

    def download_file(self, server_path, sub_folder=''):
        """Download a file from the server.

        Parameters
        ----------
        server_path : str
            Path to the file
        sub_folder : str, optional
            Name of the local subfolder to store the results in, by default ''
        """
        os.makedirs(os.path.join(self.local_folder, sub_folder), exist_ok=True)
        local_path = os.path.join(self.local_folder, sub_folder,
                                  os.path.basename(server_path))
        if not self.overwrite and os.path.isfile(local_path):
            logger.info('File %s already downloaded. Skipping...', server_path)
            return
        logger.info('Downloading %s', server_path)
        logger.debug('Downloading to %s', local_path)

        self._client.sendcmd("TYPE i")
        size = self._client.size(server_path)

        widgets = [
            DataSize(),
            Bar(),
            Percentage(), ' ',
            FileTransferSpeed(), ' (',
            ETA(), ')'
        ]

        progress = ProgressBar(max_value=size, widgets=widgets)
        progress.start()

        with open(local_path, 'wb') as file_handler:

            def _file_write(data):
                file_handler.write(data)
                nonlocal progress
                progress += len(data)

            try:
                self._client.retrbinary(f'RETR {server_path}', _file_write)
            except Exception:
                file_handler.close()
                if os.path.exists(local_path):
                    os.remove(local_path)
                raise

        progress.finish()


class CCIDownloader(FTPDownloader):
    """Downloader for the CDA ESA-CCI repository.

    Parameters
    ----------
    config : dict
        ESMValTool's user configuration
    dataset : str
        Dataset to download
    dataset_info : dict
        Dataset information from the datasets.yml file
    overwrite : bool
        Overwrite already downloaded files
    """
    def __init__(self, config, dataset, dataset_info, overwrite):
        super().__init__(config, 'anon-ftp.ceda.ac.uk', dataset, dataset_info,
                         overwrite)
        self.ftp_name = self.dataset_name[7:]

    def set_cwd(self, path):
        """Set current work directory.

        Relative to the dataset root folder.

        Parameters
        ----------
        path : str
            Remote path to set as current working directory.
        """
        cwd = f'/neodc/esacci/{self.ftp_name}/data/{path}'
        super().set_cwd(cwd)

    @property
    def dataset_name(self):
        """Name of the dataset in the repository.

        Returns
        -------
        str
            Name of the dataset in the repository.
        """
        return self.dataset.lower().replace('-', '_')

    def download_year(self, year):
        """Download a specific year.

        Parameters
        ----------
        year : int
            Year to download
        """
        self.download_folder(str(year))
