"""Script to download CDS-UERRA from the CDS."""

import calendar
import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.cds import CDSDownloader


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
            'soil_level': [
                '1',
                '2',
                '3',
            ],
            'time': ['00:00', '06:00', '12:00', '18:00'],
        },
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    loop_date = start_date
    while loop_date <= end_date:
        downloader.download(
            loop_date.year,
            loop_date.month, [
                f'{i+1:02d}' for i in range(
                    calendar.monthrange(loop_date.year, loop_date.month)[1])
            ],
            file_format='nc')
        loop_date += relativedelta.relativedelta(months=1)
