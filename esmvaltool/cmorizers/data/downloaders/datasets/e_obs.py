"""Script to download E-OBS from its webpage."""
import logging

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(config, dataset, dataset_info, start_date, end_date,
                     overwrite):
    """Download dataset.

    Parameters
    ----------
    config : dict
        ESMValTool's user configuration
    dataset : str
        Name of the dataset
    dataset_info : dict
         Dataset information from the datasets.yml file
    start_date : datetime
        Start of the interval to download
    end_date : datetime
        End of the interval to download
    overwrite : bool
        Overwrite already downloaded files
    """
    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    for var in ['TG', 'TN', 'TX', 'RR', 'PP']:
        for grid in ('0.1deg', '0.25deg'):
            for version in ('20.0e', ):
                downloader.download_file(
                    "https://knmi-ecad-assets-prd.s3.amazonaws.com/ensembles/"
                    f"data/Grid_{grid}_reg_ensemble/"
                    f"{var.lower()}_ens_mean_{grid}_reg_v{version}.nc",
                    wget_options=[])
