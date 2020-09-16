"""Downloader for the Climate Data Store."""

import os
import logging
from collections.abc import Iterable
import cdsapi

from .downloader import BaseDownloader

logger = logging.getLogger(__name__)


class CDSDownloader(BaseDownloader):
    """
    Downloader class for the climate data store.

    Parameters
    ----------
    product_name : str
        Name of the product in the CDS
    config : dict
        ESMValTool's user configuration
    request_dictionary : dict
        Common CDS request parameters
    dataset : str
        Name of the dataset
    overwrite : bool
        Overwrite already downloaded files
    extra_name : str, optional
        Some products have a subfix appended to their name for certain
        variables. This parameter is to specify it, by default ''
    """

    def __init__(self, product_name, config, request_dictionary, dataset,
                 overwrite, extra_name=''):
        super().__init__(config, dataset, overwrite)
        self._client = cdsapi.Client()
        self._product_name = product_name
        self._request_dict = request_dictionary
        self.extra_name = extra_name

    def download(self, year, month, day=None):
        """
        Download a specfic month from the CDS.

        Parameters
        ----------
        year : int
            Year to download
        month : int
            Month to download
        day : int, list(int), optional
            Day or days to download, by default None
        """
        request_dict = self._request_dict.copy()
        request_dict['year'] = f'{year}'
        request_dict['month'] = f"{month:02d}"
        if day:
            if isinstance(day, Iterable):
                request_dict['day'] = day
            else:
                request_dict['day'] = f"{day:02d}"

        date_str = f"{self.dataset}{self.extra_name}_{year}{month:02d}"
        if day:
            if not isinstance(day, Iterable):
                date_str += f"{day:02d}"

        os.makedirs(self.local_folder, exist_ok=True)
        file_path = f"{self.dataset}_{date_str}.tar"
        self.download_request(file_path, request_dict)

    def download_request(self, filename, request=None):
        """
        Download a specific request.

        Parameters
        ----------
        filename : str
            Name of the file to download
        request : dict, optional
            Request dictionary for the CDS, by default None
        """
        if request is None:
            request = self._request_dict.copy()

        os.makedirs(self.local_folder, exist_ok=True)
        filename = os.path.join(self.local_folder, filename)
        if os.path.exists(filename):
            if self.overwrite:
                os.remove(filename)
            else:
                logger.info(
                    'File %s already downloaded. Skipping...', filename)
                return
        self._client.retrieve(
            self._product_name,
            request,
            filename,
        )
