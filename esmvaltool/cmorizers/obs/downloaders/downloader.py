"""Downloader for the Climate Data Store"""

import os


class BaseDownloader():

    def __init__(self, config, dataset):
        self._config = config
        self.tier = 3
        self.dataset = dataset
        self.overwrite = False

    @property
    def local_folder(self):
        return os.path.join(
            self.rawobs_folder, f'Tier{self.tier}', self.dataset)

    @property
    def rawobs_folder(self):
        return self._config['rootpath']['RAWOBS'][0]
