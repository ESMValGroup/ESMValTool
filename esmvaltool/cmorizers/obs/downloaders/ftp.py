"""Downloader for the Climate Data Store"""

import os
import ftplib
import logging
import re

from progressbar import ProgressBar, Bar, DataSize, ETA,\
    FileTransferSpeed, Percentage

from .downloader import BaseDownloader


logger = logging.getLogger(__name__)


class FTPDownloader(BaseDownloader):

    def __init__(self, config, server, dataset):
        super().__init__(config, dataset)
        self._client = None
        self.server = server
        self.tier = 2

    def connect(self):
        self._client = ftplib.FTP(self.server)
        logger.info(self._client.getwelcome())
        self._client.login()

    def set_cwd(self, path):
        logger.debug('Current working directory: %s', self._client.pwd())
        logger.debug('Setting working directory to %s', path)
        self._client.cwd(path)
        logger.debug('New working directory: %s', self._client.pwd())

    def list_folders(self, server_path='.'):
        filenames = self._client.mlsd(server_path, facts=['type'])
        return [
            filename for filename, facts in filenames
            if facts['type'] == 'dir'
        ]

    def download_folder(self, server_path, sub_folder='', filter_files=None):
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
        os.makedirs(os.path.join(self.local_folder, sub_folder), exist_ok=True)
        local_path = os.path.join(
            self.local_folder, sub_folder, os.path.basename(server_path))
        if not self.overwrite and os.path.isfile(local_path):
            logger.info('File %s already downloaded. Skipping...', server_path)
            return
        logger.info('Downloading %s', server_path)
        logger.debug('Downloading to %s', local_path)

        self._client.sendcmd("TYPE i")
        size = self._client.size(server_path)

        widgets = [
            DataSize(), Bar(), Percentage(), ' ', FileTransferSpeed(),
            ' (', ETA(), ')'
        ]

        progress = ProgressBar(max_value=size, widgets=widgets)
        progress.start()

        with open(local_path, 'wb') as file_handler:
            def _file_write(data):
                file_handler.write(data)
                nonlocal progress
                progress += len(data)

            try:
                self._client.retrbinary(
                    f'RETR {server_path}', _file_write)
            except Exception:
                file_handler.close()
                if os.path.exists(local_path):
                    os.remove(local_path)
                raise

        progress.finish()


class CCIDownloader(FTPDownloader):

    def __init__(self, config, dataset):
        super().__init__(config, 'anon-ftp.ceda.ac.uk', dataset)
        self.ftp_name = self.dataset_name

    def set_cwd(self, path):
        cwd = f'/neodc/esacci/{self.ftp_name}/data/{path}'
        super().set_cwd(cwd)

    @ property
    def dataset_name(self):
        return self.dataset.lower().replace('-', '_')

    def download_year(self, year):
        self.download_folder(str(year))
