"""Downloader for the Climate Data Store"""

import os
import ftplib
import logging
import subprocess

from progressbar import ProgressBar, Bar, DataSize, ETA,\
    FileTransferSpeed, Percentage

logger = logging.getLogger(__name__)


class WGetDownloader():

    def __init__(self, config, dataset):
        self.config = config
        self.dataset = dataset
        self.tier = 2

    @property
    def local_folder(self):
        return os.path.join(
            self.config['rootpath']['RAWOBS'][0],
            f'Tier{self.tier}',
            self.dataset
        )

    def download_folder(self, server_path, wget_options):
        # get filenames within the directory
        subprocess.check_output(
            [
                'wget',
                f'--directory-prefix={self.local_folder}',
                '--recursive',
                '--no-directories',
                '--no-clobber',
                f'{server_path}',
            ] + wget_options
        )


class NASADownloader(WGetDownloader):

    def __init__(self, config, dataset):
        super().__init__(config, dataset)
        self.tier = 2

    def download_folder(self, folder_path):
        wget_options = [
            "--load-cookies",
            "~/.urs_cookies",
            "--save-cookies",
            "~/.urs_cookies",
            "--auth-no-challenge=on",
            "--keep-session-cookies",
            "-np",
            "--accept=nc,nc4"
        ]
        super().download_folder(folder_path, wget_options)
