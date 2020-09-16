"""Script to download CALIPSO-GOCCP from IPSL's ftp server."""

from dateutil import relativedelta
import datetime

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
        server='ftp.climserv.ipsl.polytechnique.fr',
        dataset=dataset,
        overwrite=overwrite,
    )
    downloader.connect()
    if not start_date:
        start_date = datetime.datetime(2007, 1, 1)
    if not end_date:
        end_date = datetime.datetime(2015, 1, 1)

    loop_date = start_date
    two_digits = r"\d{2}"
    while loop_date <= end_date:
        year = loop_date.year
        downloader.set_cwd(
            f"/cfmip/GOCCP_v3/3D_CloudFraction/grid_2x2xL40/{year}/avg/")
        downloader.download_folder(
            ".",
            filter_files=(
                f"3D_CloudFraction330m_{year}{two_digits}"
                "_avg_CFMIP2_sat_3.1.2.nc"
            )
        )
        loop_date += relativedelta.relativedelta(years=1)
