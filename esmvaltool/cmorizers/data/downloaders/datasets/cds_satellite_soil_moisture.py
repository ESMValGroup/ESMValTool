"""Script to download CDS-SATELLITE-SOIL-MOISTURE from the CDS."""

import calendar
import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.cds import CDSDownloader
from esmvaltool.cmorizers.data.utilities import unpack_files_in_folder


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
        start_date = datetime.datetime(1991, 9, 1)
    if not end_date:
        end_date = datetime.datetime(2020, 6, 30)

    loop_date = start_date
    downloader = CDSDownloader(
        product_name='satellite-soil-moisture',
        request_dictionary={
            'format': 'tgz',
            'variable': 'volumetric_surface_soil_moisture',
            'type_of_sensor': 'combined_passive_and_active',
            'type_of_record': 'cdr',
            'version': 'v201912.0.0',
            'time_aggregation': 'month_average',
            'day': ['01']
        },
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    monthly_downloaders = {}
    daily_downloaders = {}

    for sensor in ['combined_passive_and_active', 'passive', 'active']:
        monthly_downloaders[sensor] = get_downloader(config, dataset,
                                                     dataset_info, overwrite,
                                                     sensor, 'month')
        daily_downloaders[sensor] = get_downloader(config, dataset,
                                                   dataset_info, overwrite,
                                                   sensor, 'day')
    while loop_date <= end_date:
        for sensor, downloader in monthly_downloaders.items():
            pattern = f'cds-satellite-soil-moisture_cdr_{sensor}_monthly'
            downloader.download(loop_date.year,
                                loop_date.month,
                                file_pattern=pattern)
        loop_date += relativedelta.relativedelta(months=1)

    loop_date = start_date
    while loop_date <= end_date:
        for sensor, downloader in daily_downloaders.items():
            downloader.download(
                loop_date.year, loop_date.month, [
                    f'{i+1:02d}' for i in range(
                        calendar.monthrange(loop_date.year, loop_date.month)
                        [1])
                ], f'cds-satellite-soil-moisture_cdr_{sensor}_daily')
        loop_date += relativedelta.relativedelta(months=1)
    unpack_files_in_folder(downloader.local_folder)


def get_downloader(config, dataset, dataset_info, overwrite, sensor,
                   frequency):
    """Create download request.

    Parameters
    ----------
    config : dict
        ESMValTool's user configuration
    dataset : str
        Name of the dataset
    dataset_info : dict
         Dataset information from the datasets.yml file
    overwrite : bool
        Overwrite already downloaded files
    sensor : str
        Type of sensor
    frequency : str
        Time aggregation
    """
    if sensor == 'active':
        variable = 'surface_soil_moisture'
    else:
        variable = 'volumetric_surface_soil_moisture'
    downloader = CDSDownloader(
        product_name='satellite-soil-moisture',
        request_dictionary={
            'format': 'tgz',
            'variable': variable,
            'type_of_sensor': sensor,
            'day': '01',
            'type_of_record': 'cdr',
            'version': 'v201912.0.0',
            'time_aggregation': f'{frequency}_average',
        },
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    return downloader
