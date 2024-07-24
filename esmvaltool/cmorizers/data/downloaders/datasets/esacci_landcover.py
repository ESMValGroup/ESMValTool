"""Script to download ESACCI-LANDCOVER pft data from the CEDA."""

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.ftp import CCIDownloader


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
    if not start_date:
        start_date = datetime(1992, 1, 1)
    if not end_date:
        end_date = datetime(2020, 12, 31)

    loop_date = start_date

    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.ftp_name = 'land_cover'
    downloader.connect()
    downloader.set_cwd('/pft/v2.0.8/')
    downloader.download_folder('.')
    while loop_date <= end_date:
        year = loop_date.year
        downloader.download_year(f'{year}')
        loop_date += relativedelta.relativedelta(years=1)
