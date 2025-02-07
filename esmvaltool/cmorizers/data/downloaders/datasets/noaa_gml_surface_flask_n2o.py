"""Script to download NOAA Global Monitoring Lab surface flask data
for N2O from NOAA's archive."""
import logging

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(config, dataset, dataset_info, overwrite):
    """Download dataset.

    Parameters
    ----------
    config : dict
        ESMValTool's user configuration
    dataset : str
        Name of the dataset
    dataset_info : dict
         Dataset information from the datasets.yml file
    overwrite : bool
        Overwrite already downloaded files
    """
    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    path = "https://gml.noaa.gov/aftp/data/trace_gases/n2o/flask/surface/"
    file = "n2o_surface-flask_ccgg_text.tar.gz"
    downloader.download_file(
        path + file,
        wget_options=[],
    )
