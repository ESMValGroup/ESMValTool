"""Script to download ESACCI-OC from the Climate Data Store(CDS)."""

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
        start_date = datetime(1997, 9, 1)
    if not end_date:
        end_date = datetime(2020, 12, 1)

    loop_date = start_date

    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.ftp_name = 'ocean_colour'
    downloader.connect()

    downloader.set_cwd('v5.0-release/geographic/netcdf/chlor_a/monthly/v5.0/')
    while loop_date <= end_date:
        year = loop_date.year
        downloader.download_year(f'{year}')
        loop_date += relativedelta.relativedelta(years=1)
