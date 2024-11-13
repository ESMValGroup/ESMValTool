"""Script to download CT2019 from NOAA's webpage."""

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.ftp import FTPDownloader


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
        start_date = datetime(2000, 1, 1)
    if end_date is None:
        end_date = datetime(2018, 12, 31)
    loop_date = start_date

    downloader = FTPDownloader(config=config,
                               dataset=dataset,
                               dataset_info=dataset_info,
                               overwrite=overwrite,
                               server='aftp.cmdl.noaa.gov')
    downloader.connect()
    downloader.set_cwd(
        'products/carbontracker/co2/CT2019/molefractions/co2_total_monthly/')

    while loop_date <= end_date:
        year = loop_date.year
        month = loop_date.month

        downloader.download_file(
            f'CT2019.molefrac_glb3x2_{year}-{month:02}.nc')
        loop_date += relativedelta.relativedelta(months=1)
