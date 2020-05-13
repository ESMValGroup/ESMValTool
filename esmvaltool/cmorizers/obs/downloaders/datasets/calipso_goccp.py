"""Script to download cds-satellite-albedo from the Climate Data Store(CDS)"""

from dateutil import relativedelta

from esmvaltool.cmorizers.obs.downloaders.ftp import FTPDownloader


def download_dataset(config, dataset, start_date, end_date):
    """Download dataset cds-satellite-albedo."""
    downloader = FTPDownloader(
        config=config,
        server='ftp.climserv.ipsl.polytechnique.fr',
        dataset=dataset
    )
    downloader.connect()

    loop_date = start_date
    two_digits = r"\d{2}"
    while loop_date <= end_date:
        year = loop_date.year
        downloader.set_cwd(
            f"/cfmip/GOCCP_v3/3D_CloudFraction/grid_2x2xL40/{year}/avg/")
        downloader.download_folder(
            ".",
            filter_files=f"3D_CloudFraction330m_{year}{two_digits}_avg_CFMIP2_sat_3.1.2.nc"
        )
        loop_date += relativedelta.relativedelta(years=1)
