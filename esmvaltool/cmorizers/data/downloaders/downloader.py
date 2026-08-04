"""Downloader base class."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    from esmvaltool.cmorizers.data.typing import DatasetInfo


class BaseDownloader:
    """Base class for all downloaders.

    Parameters
    ----------
    original_data_dir:
        Directory where original data will be stored.
    dataset:
        Dataset to download
    dataset_info:
        Dataset information from the datasets.yml file
    overwrite:
        Overwrite already downloaded files
    """

    def __init__(
        self,
        original_data_dir: Path,
        dataset: str,
        dataset_info: DatasetInfo,
        *,
        overwrite: bool,
    ) -> None:
        self.original_data_dir = original_data_dir
        self.tier = dataset_info["tier"]
        self.dataset = dataset
        self.dataset_info = dataset_info
        self.overwrite = overwrite

    @property
    def local_folder(self) -> str:
        """Folder to store the downloader date.

        Returns
        -------
        :
            Path to the download folder
        """
        return str(self.original_data_dir / f"Tier{self.tier}" / self.dataset)
