"""Script to download ESACCI-AEROSOL from CCI CEDA ftp."""

from datetime import datetime
from dateutil import relativedelta

from esmvaltool.cmorizers.obs.downloaders.wget import WGetDownloader
from esmvaltool.cmorizers.obs.utilities import unpack_files_in_folder


def download_dataset(config, dataset, start_date, end_date, overwrite):
    """Download dataset ESACCI-AEROSOL."""
    if start_date is None:
        start_date = datetime(1984, 1, 1)
    if end_date is None:
        end_date = datetime(2016, 1, 1)
    loop_date = start_date

    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )

    while loop_date <= end_date:
        downloader.download_file(
            "https://isccp.giss.nasa.gov/pub/flux-fh/tar-nc4_MPF/"
            f"ISCCP-FH_nc4_MPF_v.0.0_{loop_date.year}.tar.gz",
            wget_options=[]
        )
        loop_date += relativedelta.relativedelta(years=1)
    unpack_files_in_folder(downloader.local_folder)
