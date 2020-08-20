"""Script to download cds-satellite-albedo from the Climate Data Store(CDS)"""

from dateutil import relativedelta
import datetime

from esmvaltool.cmorizers.obs.downloaders.cds import CDSDownloader
from esmvaltool.cmorizers.obs.utilities import unpack_files_in_folder


def download_dataset(config, dataset, start_date, end_date, overwrite):
    """Download dataset cds-satellite-albedo."""
    if not start_date:
        start_date = datetime.datetime(1998, 4, 1)
    if not end_date:
        end_date = datetime.datetime(2014, 5, 1)

    loop_date = start_date
    downloader = CDSDownloader(
        product_name='satellite-lai-fapar',
        request_dictionary={
            'variable': [
                'fapar', 'lai',
            ],
            'satellite': 'spot',
            'sensor': 'vgt',
            'horizontal_resolution': '1km',
            'product_version': 'V1',
            'nominal_day': '20',
            'format': 'tgz',
        },
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )

    while loop_date <= end_date:
        downloader.download(loop_date.year, loop_date.month)
        loop_date += relativedelta.relativedelta(months=1)

    unpack_files_in_folder(downloader.local_folder)
