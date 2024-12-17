"""Script to download daily and monthly ESACCI-CLOUD data."""
import logging
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
    if start_date is None:
        start_date = datetime(2000, 1, 1)
    if end_date is None:
        end_date = datetime(2000, 2, 1)
    loop_date = start_date

    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    # Base paths for L3U (daily data) and L3C (monthly data)
    base_path_l3u = ('https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/'
                     'CLD_PRODUCTS/v3.0/L3U/')
    base_path_l3c = ('https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/'
                     'CLD_PRODUCTS/v3.0/L3C/')

    # File patterns for daily (L3U) and monthly (L3C) data
    files_l3u = [
        "*-ESACCI-L3U_CLOUD-CLD_MASKTYPE-AVHRR_*-fv3.0.nc",
        "*-ESACCI-L3U_CLOUD-CLD_PRODUCTS-AVHRR_*-fv3.0.nc"
    ]
    files_l3c = ["*-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_*-fv3.0.nc"]

    wget_options = [
        '-r',
        '-nH',  # Disable the creation of directory structure
        '-e',
        'robots=off',  # Ignore robots.txt
        '--cut-dirs=9',
        '--no-parent',  # Don't ascend to the parent directory
        '--reject="index.html"',  # Reject any HTML files
        #'--accept=*.nc'  # Accept only .nc files
    ]

    #end_year = end_date.year
    #loop_year = start_date.year
    while loop_date <= end_date:
    #while loop_year <= end_year:
        year = loop_date.year
        month = loop_date.month
        date = f'{year}{month:02}'

        if int(date) in range(198201, 198502):
            sat_am = ''
            sat_pm = 'AVHRR-PM/AVHRR_NOAA-7/'
        elif int(date) in range(198502, 198811):
            sat_am = ''
            sat_pm = 'AVHRR-PM/AVHRR_NOAA-9/'
        elif int(date) in range(198811, 199109):
            sat_am = ''
            sat_pm = 'AVHRR-PM/AVHRR_NOAA-11/'
        elif int(date) in range(199109, 199409):
            sat_am = 'AVHRR-AM/AVHRR_NOAA-12/'
            sat_pm = 'AVHRR-PM/AVHRR_NOAA-11/'
        elif int(date) in range(199409, 199502):
            sat_am = 'AVHRR-AM/AVHRR_NOAA-12/'
            sat_pm = ''
        elif int(date) in range(199502, 199901):
            sat_am = 'AVHRR-AM/AVHRR_NOAA-12/'
            sat_pm = 'AVHRR-PM/AVHRR_NOAA-14/'
        elif int(date) in range(199901, 200104):
            sat_am = 'AVHRR-AM/AVHRR_NOAA-15/'
            sat_pm = 'AVHRR-PM/AVHRR_NOAA-14/'
        elif int(date) in range(200104, 200211):
            sat_am = 'AVHRR-AM/AVHRR_NOAA-15/'
            sat_pm = 'AVHRR-PM/AVHRR_NOAA-16/'
        elif int(date) in range(200211, 200509):
            sat_am = 'AVHRR-AM/AVHRR_NOAA-17/'
            sat_pm = 'AVHRR-PM/AVHRR_NOAA-16/'
        elif int(date) in range(200509, 200707):
            sat_am = 'AVHRR-AM/AVHRR_NOAA-17/'
            sat_pm = 'AVHRR-PM/AVHRR_NOAA-18/'
        elif int(date) in range(200707, 200906):
            sat_am = 'AVHRR-AM/AVHRR_METOPA/'
            sat_pm = 'AVHRR-PM/AVHRR_NOAA-18/'
        elif int(date) in range(200906, 201701):
            sat_am = 'AVHRR-AM/AVHRR_METOPA/'
            sat_pm = 'AVHRR-PM/AVHRR_NOAA-19/'
        else:
            logger.error("Data for this date %s is not available",
                         date)

        ## Download daily data from L3U
        #for sat in (sat_am, sat_pm):
        #    logger.info("Downloading daily data (L3U) for sat = %s", sat)
        #    if sat != '':
        #        folder_l3u = base_path_l3u + sat + f'{year}/{month:02}'
        #        logger.info("Download folder for daily data (L3U): %s",
        #                    folder_l3u)
        #        try:
        #            downloader.download_file(folder_l3u, wget_options)
        #        except Exception as e:
        #            logger.error("Failed to download daily data from %s: %s",
        #                         folder_l3u, str(e))

        # Download monthly data from L3C
        for sat in (sat_am, sat_pm):
            logger.info("Downloading monthly data (L3C) for sat = %s", sat)
            if sat != '':
                folder_l3c = base_path_l3c + sat + f'{year}/'
                wget_options_new = wget_options.copy()
                wget_options_new.append(f'--accept={date}*.nc')
                logger.info("Download folder for monthly data (L3C): %s",
                            folder_l3c)
                try:
                    downloader.download_file(folder_l3c, wget_options_new)
                except Exception as e:
                    logger.error("Failed to download monthly data from %s: %s",
                                 folder_l3c, str(e))

        # Increment the loop_date by one month
        loop_date += relativedelta.relativedelta(months=1)
