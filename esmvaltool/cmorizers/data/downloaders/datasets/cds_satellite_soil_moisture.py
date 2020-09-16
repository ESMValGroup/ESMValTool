"""Script to download cds-satellite-soil-moisture from the CDS"""

from dateutil import relativedelta
import datetime
import calendar

from esmvaltool.cmorizers.data.downloaders.cds import CDSDownloader
from esmvaltool.cmorizers.data.utilities import unpack_files_in_folder


def download_dataset(config, dataset, start_date, end_date, overwrite):
    """
    Download dataset.

    Parameters
    ----------
    config : dict
        ESMValTool's user configuration
    dataset : str
        Name of the dataset
    start_date : datetime
        Start of the interval to download
    end_date : datetime
        End of the interval to download
    overwrite : bool
        Overwrite already downloaded files
    """
    if not start_date:
        start_date = datetime.datetime(1978, 11, 1)
    if not end_date:
        end_date = datetime.datetime(2020, 6, 30)

    # loop_date = start_date
    # downloader = CDSDownloader(
    #     product_name='satellite-soil-moisture',
    #     request_dictionary={
    #         'format': 'tgz',
    #         'variable': 'volumetric_surface_soil_moisture',
    #         'type_of_sensor': 'combined_passive_and_active',
    #         'type_of_record': 'cdr',
    #         'version': 'v201912.0.0',
    #         'time_aggregation': 'month_average',
    #         'day': ['01']
    #     },
    #     config=config,
    #     dataset=dataset,
    #     overwrite=overwrite,
    # )

    # while loop_date <= end_date:
    #     downloader.download(loop_date.year, loop_date.month)
    #     loop_date += relativedelta.relativedelta(months=1)

    downloader = CDSDownloader(
        product_name='satellite-soil-moisture',
        request_dictionary={
            'format': 'tgz',
            'variable': 'volumetric_surface_soil_moisture',
            'type_of_sensor': 'combined_passive_and_active',
            'day': '01',
            'type_of_record': 'cdr',
            'version': 'v201912.0.0',
            'time_aggregation': 'day_average',
        },
        config=config,
        dataset=dataset,
        overwrite=overwrite,
        extra_name='_monthly',
    )
    loop_date = start_date
    while loop_date <= end_date:
        downloader.download(
            loop_date.year, loop_date.month,
            [f'{i+1:02d}' for i in range(
                calendar.monthrange(loop_date.year, loop_date.month)[1])])
        loop_date += relativedelta.relativedelta(months=1)

    unpack_files_in_folder(downloader.local_folder)
