"""Downloader for the Copernicus Marine Service data store."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import copernicusmarine

from .downloader import BaseDownloader

if TYPE_CHECKING:
    from pathlib import Path

    from esmvaltool.cmorizers.data.typing import DatasetInfo

logger = logging.getLogger(__name__)


class CMSDownloader(BaseDownloader):
    """Downloader class for the COpernicus Marine Service.

    Parameters
    ----------
    product_name:
        Name of the product in the Copernicus Marine Service data store.
    product_filename:
        Product filename in the Copernicus Marine Service data store.
    original_data_dir:
        Directory where original data will be stored.
    dataset:
        Name of the dataset
    dataset_info:
        Dataset information from the datasets.yml file
    overwrite:
        Overwrite already downloaded files
    """

    def __init__(
        self,
        original_data_dir: Path,
        product_name: str,
        dataset: str,
        dataset_info: DatasetInfo,
        *,
        overwrite: bool,
        no_directories: bool,
        product_filename: str | None = None,
    ) -> None:
        super().__init__(
            original_data_dir=original_data_dir,
            dataset=dataset,
            dataset_info=dataset_info,
            overwrite=overwrite,
        )
        self._product_name = product_name
        self._product_filename = product_filename
        self._no_directories = no_directories
        if product_filename is None:
            self._product_filename = ""

    def download(self):
        """Download a whole dataset from the Copernicus Marine Service."""
        _ = copernicusmarine.get(
            dataset_id=self._product_name,
            output_directory=self.local_folder,
            no_directories=self._no_directories,
            filter=f"{self._product_filename}*.nc",
        )

    def download_year(self, year):
        """Download a specific year from the Copernicus Marine Service.

        Parameters
        ----------
        year : int
            Year to download
        """
        _ = copernicusmarine.get(
            dataset_id=self._product_name,
            output_directory=self.local_folder,
            no_directories=self._no_directories,
            filter=f"{self._product_filename}_{year}*.nc",
        )
