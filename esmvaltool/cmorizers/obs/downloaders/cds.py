"""Downloader for the Climate Data Store"""

import os
import cdsapi


class CDSDownloader():

    def __init__(self, product_name, config, request_dictionary, dataset):
        self._client = cdsapi.Client()
        self._product_name = product_name
        self._request_dict = request_dictionary
        self._config = config
        self.tier = 3
        self.dataset = dataset
        self.overwrite = True

    def download(self, year, month, day=None):
        request_dict = self._request_dict.copy()
        request_dict['year'] = f'{year}'
        request_dict['month'] = f"{month:02d}"
        if day:
            request_dict['day'] = f"{day:02d}"

        date_str = f"{self.dataset}_{year}{month:02d}"
        if day:
            date_str += f"{day:02d}"

        os.makedirs(self.download_folder, exist_ok=True)
        file_path = os.path.join(
            self.download_folder,
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

    @property
    def download_folder(self):
        return os.path.join(self._rawobs, f'Tier{self.tier}', self.dataset)

    @property
    def _rawobs(self):
        return self._config['rootpath']['RAWOBS'][0]
