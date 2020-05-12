"""Downloader for the Climate Data Store"""

import os
import ftplib
import logging

logger = logging.getLogger(__name__)


class FTPDownloader():

    def __init__(self, config, server, dataset):
        self._client = None
        self.config = config
        self.server = server
        self.dataset = dataset
        self.overwrite = False

    def connect(self):
        self._client = ftplib.FTP(self.server)
        logger.info(self._client.getwelcome())
        self._client.login()

    def list_folders(self, server_path='.'):
        filenames = self._client.mlsd(server_path, facts=['type'])
        return [
            filename for filename, facts in filenames if facts['type'] == 'dir']

    @property
    def local_folder(self):
        return os.path.join(
            self.config['rootpath']['RAWOBS'][0],
            f'Tier{self.tier}',
            self.dataset
        )

    def download_folder(self, server_path):
        # get filenames within the directory
        filenames = self._client.nlst(server_path)
        logger.info('Downloading files in %s', server_path)
        for filename in filenames:
            self.download_file(filename)

    def download_file(self, server_path):
        os.makedirs(self.local_folder, exist_ok=True)
        local_path = os.path.join(
            self.local_folder, os.path.basename(server_path))
        if os.path.isfile(local_path):
            logger.info('File %s already downloaded. Skipping...', server_path)
            return
        logger.info('Downloading %s', server_path)
        logger.debug('Downloading to %s', local_path)

        with open(local_path, 'wb') as file_handler:
            try:
                self._client.retrbinary(
                    f'RETR {server_path}', file_handler.write)
            except:
                file_handler.close()
                if os.path.exists(local_path):
                    os.remove(local_path)
                raise


class CCIDownloader(FTPDownloader):

    def __init__(self, config, dataset):
        super().__init__(config, 'anon-ftp.ceda.ac.uk', dataset)
        self.tier = 2

    def set_cwd(self, path):
        cwd = f'/neodc/esacci/{self.dataset_name[7:]}/data/{path}'
        logger.debug('Current working directory: %s', self._client.pwd())
        logger.debug('Setting working directory to %s', cwd)
        self._client.cwd(cwd)
        logger.debug('New working directory: %s', self._client.pwd())

    @ property
    def dataset_name(self):
        return self.dataset.lower().replace('-', '_')

    def download_year(self, year):
        self.download_folder(str(year))
