"""Script to download ESACCI-LST from CCI CEDA ftp."""

import logging

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.ftp import CCIDownloader

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
        start_date = datetime(2003, 1, 1)
    if end_date is None:
        end_date = datetime(2018, 12, 31)
    loop_date = start_date

    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
        dataset_name='land_surface_temperature',
    )
    downloader.connect()

    # download daily data

#    loop_date = start_date
#    rel_base_dir = 'AQUA_MODIS/L3C/0.01/v3.00/daily'
#    while loop_date <= end_date:
#        downloader.set_cwd(rel_base_dir)
#        if downloader.exists(f'{loop_date.year}'):
#            downloader.set_cwd(f'{rel_base_dir}/{loop_date.year}')
#            if downloader.exists(f'{loop_date.month:02}'):
#                downloader.set_cwd(f'{rel_base_dir}/{loop_date.year}/'
#                                   f'{loop_date.month:02}')
#                if downloader.exists(f'{loop_date.day:02}'):
#                    downloader.download_folder(f'{loop_date.day:02}',
#                        sub_folder=f'{loop_date.year}_daily')
#                else:
#                    logger.info('%d/%d/%d: no data found', loop_date.year,
#                                loop_date.month, loop_date.day)
#                loop_date += relativedelta.relativedelta(days=1)
#            else:
#                logger.info('%d/%d: no data found', loop_date.year,
#                            loop_date.month)
#                loop_date += relativedelta.relativedelta(months=1)
#        else:
#            logger.info('%d: no data found', loop_date.year)
#            loop_date += relativedelta.relativedelta(years=1)

    # download monthly data

    loop_date = start_date
    rel_base_dir = 'AQUA_MODIS/L3C/0.01/v3.00/monthly'
    while loop_date <= end_date:
        downloader.set_cwd(rel_base_dir)
        if downloader.exists(f'{loop_date.year}'):
            downloader.set_cwd(f'{rel_base_dir}/{loop_date.year}')
            if downloader.exists(f'{loop_date.month:02}'):
                downloader.download_folder(f'{loop_date.month:02}')
            else:
                logger.info('%d/%d: no data found', loop_date.year,
                            loop_date.month)
            loop_date += relativedelta.relativedelta(months=1)
        else:
            logger.info('%d: no data found', loop_date.year)
            loop_date += relativedelta.relativedelta(years=1)
