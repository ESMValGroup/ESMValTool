"""Script to download daily and monthly ESACCI-CLOUD data."""

import logging
from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader
from esmvaltool.cmorizers.data.utilities import read_cmor_config

logger = logging.getLogger(__name__)


def download_dataset(
    original_data_dir,
    dataset,
    dataset_info,
    start_date,
    end_date,
    overwrite,
):
    """Download dataset.

    Parameters
    ----------
    original_data_dir : Path
        Directory where original data will be stored.
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
    start_date_day = False
    if start_date is None:
        start_date = datetime(1982, 1, 1)
        start_date_day = datetime(2003, 1, 1)
    if end_date is None:
        end_date = datetime(2016, 12, 31)
    loop_date = start_date

    downloader = WGetDownloader(
        original_data_dir=original_data_dir,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    # check if daily data needs to be downloaded
    cmor_config = read_cmor_config(dataset)
    daily_data = cmor_config["daily_data"]
    if not daily_data:
        logger.info(
            'If daily data needs to be downloaded change "daily_data" in the '
            'cmor_config file to "True" '
            "(esmvaltool/cmorizers/data/cmor_config/ESACCI-CLOUD.yml)",
        )

    # Base paths for L3U (daily data) and L3C (monthly data)
    base_path_l3u = (
        "https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/"
        "CLD_PRODUCTS/v3.0/L3U/"
    )
    base_path_l3c = (
        "https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/"
        "CLD_PRODUCTS/v3.0/L3C/"
    )

    wget_options = [
        "-r",
        "-e robots=off",  # Ignore robots.txt
        "--no-parent",  # Don't ascend to the parent directory
        '--reject="index.html"',  # Reject any HTML files
    ]

    while loop_date <= end_date:
        year = loop_date.year
        month = loop_date.month
        date = f"{year}{month:02}"

        if datetime(1982, 1, 1) <= loop_date < datetime(1985, 2, 1):
            sat_am = ""
            sat_pm = "AVHRR-PM/AVHRR_NOAA-7/"
        elif datetime(1985, 2, 1) <= loop_date < datetime(1988, 11, 1):
            sat_am = ""
            sat_pm = "AVHRR-PM/AVHRR_NOAA-9/"
        elif datetime(1988, 11, 1) <= loop_date < datetime(1991, 9, 1):
            sat_am = ""
            sat_pm = "AVHRR-PM/AVHRR_NOAA-11/"
        elif datetime(1991, 9, 1) <= loop_date < datetime(1994, 9, 1):
            sat_am = "AVHRR-AM/AVHRR_NOAA-12/"
            sat_pm = "AVHRR-PM/AVHRR_NOAA-11/"
        elif datetime(1994, 9, 1) <= loop_date < datetime(1995, 2, 1):
            sat_am = "AVHRR-AM/AVHRR_NOAA-12/"
            sat_pm = ""
        elif datetime(1995, 2, 1) <= loop_date < datetime(1999, 1, 1):
            sat_am = "AVHRR-AM/AVHRR_NOAA-12/"
            sat_pm = "AVHRR-PM/AVHRR_NOAA-14/"
        elif datetime(1999, 1, 1) <= loop_date < datetime(2001, 4, 1):
            sat_am = "AVHRR-AM/AVHRR_NOAA-15/"
            sat_pm = "AVHRR-PM/AVHRR_NOAA-14/"
        elif datetime(2001, 4, 1) <= loop_date < datetime(2002, 11, 1):
            sat_am = "AVHRR-AM/AVHRR_NOAA-15/"
            sat_pm = "AVHRR-PM/AVHRR_NOAA-16/"
        elif datetime(2002, 11, 1) <= loop_date < datetime(2005, 9, 1):
            sat_am = "AVHRR-AM/AVHRR_NOAA-17/"
            sat_pm = "AVHRR-PM/AVHRR_NOAA-16/"
        elif datetime(2005, 9, 1) <= loop_date < datetime(2007, 7, 1):
            sat_am = "AVHRR-AM/AVHRR_NOAA-17/"
            sat_pm = "AVHRR-PM/AVHRR_NOAA-18/"
        elif datetime(2007, 7, 1) <= loop_date < datetime(2009, 6, 1):
            sat_am = "AVHRR-AM/AVHRR_METOPA/"
            sat_pm = "AVHRR-PM/AVHRR_NOAA-18/"
        elif datetime(2009, 6, 1) <= loop_date < datetime(2017, 1, 1):
            sat_am = "AVHRR-AM/AVHRR_METOPA/"
            sat_pm = "AVHRR-PM/AVHRR_NOAA-19/"
        else:
            msg = f"Data for this date {date} is not available"
            raise ValueError(msg)

        # Download monthly data from L3C
        for sat in (sat_am, sat_pm):
            if sat != "":
                # monthly data
                logger.info("Downloading monthly data (L3C) for sat = %s", sat)
                folder_l3c = base_path_l3c + sat + f"{year}/"
                wget_options_l3c = wget_options.copy()
                wget_options_l3c.append(f"--accept={date}*.nc")
                logger.info(
                    "Download folder for monthly data (L3C): %s",
                    folder_l3c,
                )
                downloader.download_file(folder_l3c, wget_options_l3c)

                # daily data
                if daily_data:
                    if not start_date_day or (
                        start_date_day
                        and datetime(2003, 1, 1)
                        <= loop_date
                        <= datetime(2007, 2, 1)
                    ):
                        logger.info(
                            "Downloading daily data (L3U) for sat = %s",
                            sat,
                        )
                        folder_l3u = base_path_l3u + sat + f"{year}/{month:02}"
                        wget_options_l3u = wget_options.copy()
                        wget_options_l3u.append(
                            f"--accept={date}*CLD_MASKTYPE*.nc,"
                            f"{date}*CLD_PRODUCTS*.nc",
                        )
                        logger.info(
                            "Download folder for daily data (L3U): %s",
                            folder_l3u,
                        )
                        downloader.download_file(folder_l3u, wget_options_l3u)

        # Increment the loop_date by one month
        loop_date += relativedelta.relativedelta(months=1)
