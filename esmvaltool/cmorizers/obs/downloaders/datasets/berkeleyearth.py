"""Script to download BerkeleyEarth from its webpage."""

from esmvaltool.cmorizers.obs.downloaders.wget import WGetDownloader


def download_dataset(config, dataset, start_date, end_date, overwrite):
    """Download dataset BerkeleyEarth."""
    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )

    downloader.download_file(
        "http://berkeleyearth.lbl.gov/auto/Global/Gridded/"
        "Land_and_Ocean_LatLong1.nc",
        wget_options=[]
    )
