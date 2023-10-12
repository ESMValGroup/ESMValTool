"""Script to download CFSR from ESGF."""
import logging
import os

from datetime import datetime

from dateutil import relativedelta

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
    if not start_date:
        start_date = datetime(1979, 1, 1)
    if not end_date:
        end_date = datetime(2019, 1, 1)
    loop_date = start_date

    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    os.makedirs(downloader.local_folder, exist_ok=True)

    url = ("https://esgf.nccs.nasa.gov/thredds/fileServer/CREATE-IP/"
           "reanalysis/NOAA-NCEP/CFSR/CFSR/mon/atmos/")

    while loop_date <= end_date:
        year = loop_date.year
        downloader.download_file(url +
                                 "hus/"
                                 "hus_Amon_reanalysis_CFSR_"
                                 f"{year}01-{year}12.nc",
                                 wget_options=[])
        downloader.download_file(url +
                                 "hur/"
                                 "hur_Amon_reanalysis_CFSR_"
                                 f"{year}01-{year}12.nc",
                                 wget_options=[])
        downloader.download_file(url +
                                 "wap/"
                                 "wap_Amon_reanalysis_CFSR_"
                                 f"{year}01-{year}12.nc",
                                 wget_options=[])
        loop_date += relativedelta.relativedelta(years=1)

    downloader.download_file(url +
                             "clt/clt_Amon_reanalysis_CFSR_197901-201912.nc",
                             wget_options=[])
    downloader.download_file(url +
                             "prw/prw_Amon_reanalysis_CFSR_197901-201912.nc",
                             wget_options=[])
    downloader.download_file(url +
                             "rlut/"
                             "rlut_Amon_reanalysis_CFSR_197901-201912.nc",
                             wget_options=[])
    downloader.download_file(url +
                             "rlutcs/"
                             "rlutcs_Amon_reanalysis_CFSR_197901-201912.nc",
                             wget_options=[])
    downloader.download_file(url +
                             "rsut/"
                             "rsut_Amon_reanalysis_CFSR_197901-201912.nc",
                             wget_options=[])
    downloader.download_file(url +
                             "rsutcs/"
                             "rsutcs_Amon_reanalysis_CFSR_197901-201912.nc",
                             wget_options=[])
    downloader.download_file(url +
                             "ts/"
                             "ts_Amon_reanalysis_CFSR_197901-201912.nc",
                             wget_options=[])
