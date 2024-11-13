"""Script to download ESACCI-OZONE from the Climate Data Store(CDS)."""

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
    if start_date is None:
        start_date = datetime(1997, 1, 1)
    if end_date is None:
        end_date = datetime(2010, 1, 1)

    loop_date = start_date

    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.ftp_name = 'ozone'
    downloader.connect()
    downloader.set_cwd(
        'limb_profiles/l3/merged/merged_monthly_zonal_mean/v0002')
    downloader.download_folder('.')

    downloader.set_cwd('total_columns/l3/merged/v0100/')
    while loop_date <= end_date:
        year = loop_date.year
        downloader.set_cwd('total_columns/l3/merged/v0100/')
        downloader.download_year(f'{year}')
        loop_date += relativedelta.relativedelta(years=1)
