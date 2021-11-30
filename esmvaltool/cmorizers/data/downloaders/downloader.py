"""Downloader base class."""

import os


class BaseDownloader():
    """Base class for all downloaders.

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
        self._config = config
        self.tier = dataset_info['tier']
        self.dataset = dataset
        self.dataset_info = dataset_info
        self.overwrite = overwrite

    @property
    def local_folder(self):
        """Folder to store the downloader date.

        Returns
        -------
        str
            Path to the download folder
        """
        return os.path.join(self.rawobs_folder, f'Tier{self.tier}',
                            self.dataset)

    @property
    def rawobs_folder(self):
        """RAWOBS base path.

        Returns
        -------
        str
            Path to the RAWOBS folder
        """
        return self._config['rootpath']['RAWOBS'][0]
