"""Script to download ESACCI-AEROSOL from CCI CEDA ftp (AATSR)
   and via Copernicus Climate Data Store (SLSTR)."""

import logging
from datetime import datetime

import cdsapi
from dateutil import relativedelta

from pathlib import Path
from esmvaltool.cmorizers.data.downloaders.ftp import CCIDownloader

logger = logging.getLogger(__name__)


def download_dataset(
    config, dataset, dataset_info, start_date, end_date, overwrite
):
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
        aatsr_year1 = 1997
        slstr_year1 = 2017
    else:
        aatsr_year1 = start_date.year
        slstr_year1 = start_date.year
    if end_date is None:
        aatsr_year2 = 2011
        slstr_year2 = 2023
    else:
        aatsr_year2 = end_date.year
        slstr_year2 = end_date.year

    # =============================
    # Download AATSR data from CEDA
    # =============================

    if start_date is None:
        start_date = datetime(aatsr_year1, 1, 1)
    if end_date is None:
        end_date = datetime(aatsr_year2, 12, 31)

    loop_date = start_date

    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.connect()

    version = 'v4.3'
    algorithm = 'SU'

    # download monthly AATSR data

    loop_date = start_date
    while loop_date <= end_date:
        if loop_date.year < 2003:
            instrument = 'ATSR2'
        else:
            instrument = 'AATSR'
        rel_base_dir = f'{instrument}_{algorithm}/L3/{version}/MONTHLY'
        downloader.set_cwd(rel_base_dir)
        if downloader.exists(f'{loop_date.year}'):
            downloader.download_folder(f'{loop_date.year}',
                                       f'{algorithm}-{version}-monthly')
        else:
            logger.info('%d: no data found', loop_date.year)
        loop_date += relativedelta.relativedelta(years=1)

    # download daily AATSR data

    loop_date = start_date
    while loop_date <= end_date:
        if loop_date.year < 2003:
            instrument = 'ATSR2'
        else:
            instrument = 'AATSR'
        rel_base_dir = f'{instrument}_{algorithm}/L3/{version}/DAILY'
        downloader.set_cwd(rel_base_dir)
        if downloader.exists(f'{loop_date.year}'):
            downloader.set_cwd(f'{rel_base_dir}/{loop_date.year}')
            if downloader.exists(f"{loop_date.month:02}"):
                downloader.download_folder(f'{loop_date.month:02}',
                                           f'{algorithm}-{version}-daily')
            else:
                logger.info('%d/%d: no data found', loop_date.year,
                            loop_date.month)
        else:
            logger.info('%d: no data found', loop_date.year)
        loop_date += relativedelta.relativedelta(months=1)

    # ================================================
    # Download SLSTR data from CDS (daily and monthly)
    # ================================================

    raw_obs_dir = Path(config["rootpath"]["RAWOBS"][0])
    output_folder = raw_obs_dir / f"Tier{dataset_info['tier']}" / dataset
    output_folder.mkdir(parents=True, exist_ok=True)

    cds_url = "https://cds.climate.copernicus.eu/api"

    requests = {
        "aod_slstr_daily": {
            "time_aggregation": "daily_average",
            "variable": "aerosol_optical_depth",
            "sensor_on_satellite": [
                "slstr_on_sentinel_3a",
                "slstr_on_sentinel_3b"
            ],
            "algorithm": ["swansea"],
            "year": [str(y) for y in range(slstr_year1, slstr_year2)],
            "month": [f"{m:02d}" for m in range(1, 13)],
            "day": [f"{m:02d}" for m in range(1, 32)],
            "version": ["v1_12"]
        },
        "aod_slstr_monthly": {
            "time_aggregation": "monthly_average",
            "variable": "aerosol_optical_depth",
            "sensor_on_satellite": [
                "slstr_on_sentinel_3a",
                "slstr_on_sentinel_3b"
            ],
            "algorithm": ["swansea"],
            "year": [str(y) for y in range(slstr_year1, slstr_year2)],
            "month": [f"{m:02d}" for m in range(1, 13)],
            "version": ["v1_12"]
        },
        "fine_aod_slstr_daily": {
            "time_aggregation": "daily_average",
            "variable": "fine_mode_aerosol_optical_depth",
            "sensor_on_satellite": [
                "slstr_on_sentinel_3a",
                "slstr_on_sentinel_3b"
            ],
            "algorithm": ["swansea"],
            "year": [str(y) for y in range(slstr_year1, slstr_year2)],
            "month": [f"{m:02d}" for m in range(1, 13)],
            "day": [f"{m:02d}" for m in range(1, 32)],
            "version": ["v1_12"]
        },
        "fine_aod_slstr_monthly": {
            "time_aggregation": "monthly_average",
            "variable": "fine_mode_aerosol_optical_depth",
            "sensor_on_satellite": [
                "slstr_on_sentinel_3a",
                "slstr_on_sentinel_3b"
            ],
            "algorithm": ["swansea"],
            "year": [str(y) for y in range(slstr_year1, slstr_year2)],
            "month": [f"{m:02d}" for m in range(1, 13)],
            "version": ["v1_12"]
        },
    }

    cds_client = cdsapi.Client(cds_url)

    for var_name, request in requests.items():
        outdir = output_folder / f"{var_name}/"
        logger.info("Downloading %s data to %s", var_name, outdir)

        file_path = outdir / f"{var_name}.gz"

        if file_path.exists() and not overwrite:
            logger.info(
                "File %s already exists. Skipping download.",
                file_path,
            )
            continue

        cds_client.retrieve(
            "satellite-aerosol-properties",
            request,
            file_path.as_posix(),
        )

        # Handle both .gz and .zip files
        with open(file_path, "rb") as file:
            magic = file.read(2)

        if magic == b"PK":  # ZIP file signature
            logger.info("Detected ZIP file: %s", file_path)
            with zipfile.ZipFile(file_path, "r") as zip_ref:
                zip_ref.extractall(outdir)
        else:
            logger.info("Detected GZIP file: %s", file_path)
            with gzip.open(file_path, "rb") as f_in:
                with open(outdir / file_path.stem, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
