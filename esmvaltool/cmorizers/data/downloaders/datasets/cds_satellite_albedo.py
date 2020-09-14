"""Script to download cds-satellite-albedo from the Climate Data Store(CDS)"""

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.cds import CDSDownloader
from esmvaltool.cmorizers.data.utilities import unpack_files_in_folder


def download_dataset(config, dataset, start_date, end_date, overwrite):
    """Download dataset cds-satellite-albedo."""
    loop_date = start_date

    downloader = CDSDownloader(
        product_name='satellite-albedo',
        request_dictionary={
            'format': 'tgz',
            'satellite': 'spot',
            'sensor': 'vgt',
            'product_version': 'V1',
            'horizontal_resolution': '1km',
            'variable': [
                'albb_bh',
                'albb_dh',
            ],
            'nominal_day': '20',
        },
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )

    while loop_date <= end_date:
        downloader.download(loop_date.year, loop_date.month)
        loop_date += relativedelta.relativedelta(months=1)

    unpack_files_in_folder(downloader.local_folder)
