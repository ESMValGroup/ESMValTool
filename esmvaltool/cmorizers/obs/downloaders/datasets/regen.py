"""Script to download ESACCI-AEROSOL from CCI CEDA ftp."""

from datetime import datetime
from dateutil import relativedelta

from esmvaltool.cmorizers.obs.downloaders.wget import WGetDownloader
from esmvaltool.cmorizers.obs.utilities import read_cmor_config


def download_dataset(config, dataset, start_date, end_date, overwrite):
    """Download dataset ESACCI-AEROSOL."""
    if start_date is None:
        start_date = datetime(1950, 1, 1)
    if end_date is None:
        end_date = datetime(2016, 1, 1)
    loop_date = start_date

    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )
    base_path = (
        "http://dapds00.nci.org.au/thredds/fileServer/ks32/CLEX_Data/"
        "REGEN_AllStns/v1-2019/REGEN_AllStns_{version}_{year}.nc"
    )
    version=read_cmor_config(dataset)['attributes']['version']
    while loop_date <= end_date:
        downloader.download_folder(
            base_path.format(year=loop_date.year, version=version),
            []
        )
        loop_date += relativedelta.relativedelta(years=1)
