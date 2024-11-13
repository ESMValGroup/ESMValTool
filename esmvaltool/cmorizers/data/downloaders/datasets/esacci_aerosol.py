"""Script to download ESACCI-AEROSOL from CCI CEDA ftp."""

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
        end_date = datetime(2011, 1, 1)
    loop_date = start_date

    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.connect()

    while loop_date <= end_date:
        year = loop_date.year
        if year < 2003:
            downloader.set_cwd('ATSR2_SU/L3/v4.21/MONTHLY')
        else:
            downloader.set_cwd('AATSR_SU/L3/v4.21/MONTHLY')

        downloader.download_year(loop_date.year)
        loop_date += relativedelta.relativedelta(years=1)
