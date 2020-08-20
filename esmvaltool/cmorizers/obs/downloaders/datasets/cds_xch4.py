"""Script to download cds-satellite-albedo from the Climate Data Store(CDS)"""

from dateutil import relativedelta
import datetime
import calendar

from esmvaltool.cmorizers.obs.downloaders.cds import CDSDownloader
from esmvaltool.cmorizers.obs.utilities import unpack_files_in_folder


def download_dataset(config, dataset, start_date, end_date, overwrite):
    """Download dataset cds-satellite-albedo."""
    downloader = CDSDownloader(
        product_name='satellite-methane',
        request_dictionary={
            'format': 'tgz',
            'processing_level': 'level_3',
            'variable': 'xch4',
            'sensor_and_algorithm': 'merged_obs4mips',
            'version': '4.1',
        },
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )
    downloader.download_request("CDS-XCH4.tar")
    unpack_files_in_folder(downloader.local_folder)
