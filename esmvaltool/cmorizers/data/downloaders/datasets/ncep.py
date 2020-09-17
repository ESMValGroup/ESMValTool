"""Script to download NCEP."""

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.ftp import FTPDownloader


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
    downloader = FTPDownloader(
        config=config,
        server='ftp.cdc.noaa.gov',
        dataset=dataset,
        overwrite=overwrite,
    )
    downloader.connect()

    loop_date = start_date

    downloader.set_cwd("/Datasets/ncep.reanalysis.derived/pressure/")
    downloader.download_file("air.mon.mean.nc", sub_folder='pressure')
    downloader.download_file("hgt.mon.mean.nc", sub_folder='pressure')
    downloader.download_file("rhum.mon.mean.nc", sub_folder='pressure')
    downloader.download_file("shum.mon.mean.nc", sub_folder='pressure')
    downloader.download_file("uwnd.mon.mean.nc", sub_folder='pressure')
    downloader.download_file("vwnd.mon.mean.nc", sub_folder='pressure')
    downloader.download_file("omega.mon.mean.nc", sub_folder='pressure')

    downloader.set_cwd("/Datasets/ncep.reanalysis.derived/surface/")
    downloader.download_file("air.mon.mean.nc", sub_folder='surface')
    downloader.set_cwd("/Datasets/ncep.reanalysis.derived/surface_gauss/")
    downloader.download_file("prate.mon.mean.nc", sub_folder='surface')

    while loop_date <= end_date:
        year = loop_date.year
        downloader.set_cwd("/Datasets/ncep.reanalysis.dailyavgs/pressure/")
        downloader.download_file(f"uwnd.{year}.nc", sub_folder='pressure')
        downloader.download_file(f"vwnd.{year}.nc", sub_folder='pressure')
        downloader.set_cwd("/Datasets/ncep.reanalysis.dailyavgs/surface_gauss")
        downloader.download_file(
            f"prate.sfc.gauss.{year}.nc", sub_folder='surface')
        downloader.set_cwd("/Datasets/ncep.reanalysis.dailyavgs/other_gauss")
        downloader.download_file(
            f"ulwrf.ntat.gauss.{year}.nc", sub_folder='surface')

        loop_date += relativedelta.relativedelta(years=1)
