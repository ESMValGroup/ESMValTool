"""Script to download NCEP-NCAR-R1."""

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.ftp import FTPDownloader


def download_dataset(
    original_data_dir,
    dataset,
    dataset_info,
    start_date,
    end_date,
    overwrite,
):
    """Download dataset.

    Parameters
    ----------
    original_data_dir : Path
        Directory where original data will be stored.
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
        start_date = datetime(1948, 1, 1)
    if end_date is None:
        end_date = datetime(2021, 1, 1)
    downloader = FTPDownloader(
        original_data_dir=original_data_dir,
        server="ftp.cdc.noaa.gov",
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.connect()

    loop_date = start_date

    downloader.set_cwd("/Datasets/ncep.reanalysis.derived/pressure/")
    downloader.download_file("air.mon.mean.nc", sub_folder="pressure")
    downloader.download_file("hgt.mon.mean.nc", sub_folder="pressure")
    downloader.download_file("rhum.mon.mean.nc", sub_folder="pressure")
    downloader.download_file("shum.mon.mean.nc", sub_folder="pressure")
    downloader.download_file("uwnd.mon.mean.nc", sub_folder="pressure")
    downloader.download_file("vwnd.mon.mean.nc", sub_folder="pressure")
    downloader.download_file("omega.mon.mean.nc", sub_folder="pressure")

    downloader.set_cwd("/Datasets/ncep.reanalysis.derived/surface/")
    downloader.download_file("air.mon.mean.nc", sub_folder="surface")
    downloader.download_file("pr_wtr.mon.mean.nc", sub_folder="surface")
    downloader.download_file("slp.mon.mean.nc", sub_folder="surface")
    downloader.download_file("wspd.mon.mean.nc", sub_folder="surface")
    downloader.download_file("rhum.mon.mean.nc", sub_folder="surface")

    downloader.set_cwd("/Datasets/ncep.reanalysis.derived/surface_gauss/")
    downloader.download_file("air.2m.mon.mean.nc", sub_folder="surface")
    downloader.download_file("prate.mon.mean.nc", sub_folder="surface")
    downloader.download_file("tmax.2m.mon.mean.nc", sub_folder="surface")
    downloader.download_file("tmin.2m.mon.mean.nc", sub_folder="surface")

    downloader.set_cwd("/Datasets/ncep.reanalysis.derived/other_gauss/")
    downloader.download_file("tcdc.eatm.mon.mean.nc", sub_folder="surface")
    downloader.download_file("ulwrf.ntat.mon.mean.nc", sub_folder="surface")
    downloader.download_file("csulf.ntat.mon.mean.nc", sub_folder="surface")
    downloader.download_file("uswrf.ntat.mon.mean.nc", sub_folder="surface")
    downloader.download_file("csusf.ntat.mon.mean.nc", sub_folder="surface")

    while loop_date <= end_date:
        year = loop_date.year
        downloader.set_cwd("/Datasets/ncep.reanalysis.dailyavgs/pressure/")
        downloader.download_file(f"uwnd.{year}.nc", sub_folder="pressure")
        downloader.download_file(f"vwnd.{year}.nc", sub_folder="pressure")
        downloader.set_cwd("/Datasets/ncep.reanalysis.dailyavgs/surface_gauss")
        downloader.download_file(
            f"prate.sfc.gauss.{year}.nc",
            sub_folder="surface",
        )
        downloader.set_cwd("/Datasets/ncep.reanalysis.dailyavgs/other_gauss")
        downloader.download_file(
            f"ulwrf.ntat.gauss.{year}.nc",
            sub_folder="surface",
        )

        loop_date += relativedelta.relativedelta(years=1)
