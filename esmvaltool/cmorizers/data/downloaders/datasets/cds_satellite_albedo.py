"""Script to download CDS-SATELLITE-ALBEDO from the Climate Data Store."""

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
    if start_date is None:
        start_date = datetime.datetime(1998, 4, 1)
    if end_date is None:
        end_date = datetime.datetime(2013, 5, 1)

    downloader = CDSDownloader(
        product_name='satellite-albedo',
        request_dictionary={
            'format': 'tgz',
            'satellite': 'spot',
            'sensor': 'vgt',
            'product_version': 'v1',
            'horizontal_resolution': '1km',
            'variable': [
                'albb_bh',
                'albb_dh',
            ],
            'nominal_day': '20',
        },
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    loop_date = start_date
    while loop_date <= end_date:
        downloader.download(loop_date.year, loop_date.month)
        loop_date += relativedelta.relativedelta(months=1)

    unpack_files_in_folder(downloader.local_folder)
