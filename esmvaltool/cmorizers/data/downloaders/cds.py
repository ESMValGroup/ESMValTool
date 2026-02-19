"""Downloader for the Climate Data Store."""

from __future__ import annotations

import logging
import os
from collections.abc import Iterable
from typing import TYPE_CHECKING

import cdsapi

from .downloader import BaseDownloader

if TYPE_CHECKING:
    from pathlib import Path

    from esmvaltool.cmorizers.data.typing import DatasetInfo

logger = logging.getLogger(__name__)


class CDSDownloader(BaseDownloader):
    """Downloader class for the climate data store.

    Parameters
    ----------
    product_name:
        Name of the product in the CDS
    original_data_dir:
        Directory where original data will be stored.
    request_dictionary:
        Common CDS request parameters
    dataset:
        Name of the dataset
    dataset_info:
        Dataset information from the datasets.yml file
    overwrite:
        Overwrite already downloaded files
    extra_name:
        Some products have a subfix appended to their name for certain
        variables. This parameter is to specify it, by default ''
    """

    def __init__(
        self,
        original_data_dir: Path,
        product_name: str,
        request_dictionary: dict,
        dataset: str,
        dataset_info: DatasetInfo,
        *,
        overwrite: bool,
        extra_name: str = "",
        cds_url: str = "https://cds.climate.copernicus.eu/api",
    ) -> None:
        super().__init__(
            original_data_dir=original_data_dir,
            dataset=dataset,
            dataset_info=dataset_info,
            overwrite=overwrite,
        )
        try:
            self._client = cdsapi.Client(url=cds_url)
        except Exception as ex:
            if str(ex).endswith(".cdsapirc"):
                logger.exception(
                    "Could not connect to the CDS due to issues with your "
                    '".cdsapirc" file. More info in '
                    "https://cds.climate.copernicus.eu/api-how-to.",
                )
            raise
        self._product_name = product_name
        self._request_dict = request_dictionary
        self.extra_name = extra_name

    def download(
        self,
        year,
        month,
        day=None,
        file_pattern=None,
        file_format="tar",
    ):
        """Download a specific month from the CDS.

        Parameters
        ----------
        year : int
            Year to download
        month : int
            Month to download
        day : int, list(int), optional
            Day or days to download, by default None
        file_pattern : str, optional
            Filename pattern, by default None
        file_format : str, optional
            File format, by default tar
        """
        request_dict = self._request_dict.copy()
        request_dict["year"] = f"{year}"
        request_dict["month"] = f"{month:02d}"
        if day:
            if isinstance(day, Iterable):
                request_dict["day"] = day
            else:
                request_dict["day"] = f"{day:02d}"

        date_str = f"{year}{month:02d}"
        if day and not isinstance(day, Iterable):
            date_str += f"{day:02d}"

        os.makedirs(self.local_folder, exist_ok=True)
        if file_pattern is None:
            file_pattern = f"{self._product_name}"
        file_path = f"{file_pattern}_{date_str}.{file_format}"
        self.download_request(file_path, request_dict)

    def download_year(self, year, file_pattern=None, file_format="zip"):
        """Download a specific year from the CDS.

        Parameters
        ----------
        year : int
            Year to download
        file_pattern : str, optional
            Filename pattern, by default None
        file_format : str, optional
            File format, by default tar
        """
        request_dict = self._request_dict.copy()
        request_dict["year"] = f"{year}"
        request_dict["month"] = [f"{m:02d}" for m in range(1, 13)]

        os.makedirs(self.local_folder, exist_ok=True)
        if file_pattern is None:
            file_pattern = f"{self._product_name}"
        file_path = f"{file_pattern}_{year}.{file_format}"
        self.download_request(file_path, request_dict)

    def download_request(self, filename, request=None):
        """Download a specific request.

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
                    "File %s already downloaded. Skipping...",
                    filename,
                )
                return
        try:
            self._client.retrieve(
                self._product_name,
                request,
                filename,
            )
        except Exception:
            logger.exception("Failed request: %s", request)
            raise
