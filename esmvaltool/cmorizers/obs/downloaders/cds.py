"""Downloader for the Climate Data Store"""

import os
import cdsapi

from .downloader import BaseDownloader


class CDSDownloader(BaseDownloader):

    def __init__(self, product_name, config, request_dictionary, dataset):
        super().__init__(config, dataset)
        self._client = cdsapi.Client()
        self._product_name = product_name
        self._request_dict = request_dictionary

    def download(self, year, month, day=None):
        request_dict = self._request_dict.copy()
        request_dict['year'] = f'{year}'
        request_dict['month'] = f"{month:02d}"
        if day:
            request_dict['day'] = f"{day:02d}"

        date_str = f"{self.dataset}_{year}{month:02d}"
        if day:
            date_str += f"{day:02d}"

        os.makedirs(self.local_folder, exist_ok=True)
        file_path = os.path.join(
            self.local_folder,
            f"{self.dataset}_{date_str}.tar"
        )
        if os.path.exists(file_path):
            if self.overwrite:
                os.remove(file_path)
            else:
                return
        self._client.retrieve(
            self._product_name,
            request_dict,
            file_path,
        )
