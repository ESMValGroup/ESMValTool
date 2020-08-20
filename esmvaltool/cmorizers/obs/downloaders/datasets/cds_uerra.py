"""Script to download cds-satellite-albedo from the Climate Data Store(CDS)"""

from dateutil import relativedelta
import datetime
import calendar

from esmvaltool.cmorizers.obs.downloaders.cds import CDSDownloader
from esmvaltool.cmorizers.obs.utilities import unpack_files_in_folder


def download_dataset(config, dataset, start_date, end_date, overwrite):
    """Download dataset cds-satellite-albedo."""
    if not start_date:
        start_date = datetime.datetime(1961, 1, 1)
    if not end_date:
        end_date = datetime.datetime(2019, 7, 1)

    loop_date = start_date
    downloader = CDSDownloader(
        product_name='reanalysis-uerra-europe-soil-levels',
        request_dictionary={
            'format': 'netcdf',
            'origin': 'uerra_harmonie',
            'variable': 'volumetric_soil_moisture',
            'soil_level': ['1', '2', '3',],
            'time': ['00:00', '06:00', '12:00', '18:00'],
        },
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )

    loop_date = start_date
    while loop_date <= end_date:
        downloader.download(
            loop_date.year, loop_date.month,
            [f'{i+1:02d}' for i in range(calendar.monthrange(loop_date.year, loop_date.month)[1])])
        loop_date += relativedelta.relativedelta(months=1)

    unpack_files_in_folder(downloader.local_folder)
