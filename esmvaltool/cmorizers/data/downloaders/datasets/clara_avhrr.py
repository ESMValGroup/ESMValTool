"""Script to download CM SAF CLARA-AHRR data from the CDS."""

import gzip
import logging
import shutil
import zipfile
from datetime import datetime
from pathlib import Path

import cdsapi

logger = logging.getLogger(__name__)


def download_dataset(
    config, dataset, dataset_info, start_date, end_date, overwrite
):
    """Download CLARA-AVHRR dataset using CDS API.

    - An ECMWF account is needed to download the datasets from
      https://cds.climate.copernicus.eu/datasets/satellite-cloud-properties
    - The file named .cdspirc containing the key associated to
      the ECMWF account needs to be saved in user's ${HOME} directory.
    - All the files will be saved in ${RAWOBS}/Tier2/CLARA-AVHRR.
    """
    cds_url = "https://cds.climate.copernicus.eu/api"

    raw_obs_dir = Path(config["rootpath"]["RAWOBS"][0])
    output_folder = raw_obs_dir / f"Tier{dataset_info['tier']}" / dataset
    output_folder.mkdir(parents=True, exist_ok=True)

    if start_date is None:
        start_date = datetime(1979, 1, 1)
    if end_date is None:
        end_date = datetime(2020, 12, 31)

    requests = {}

    # The CDS requests for daily values are done for each month separately
    # to avoid the error "cost limits exceeded".

    for year in range(start_date.year, end_date.year + 1):
        requests.update(
            {
                "clivi_monthly_" + str(year): {
                    "product_family": "clara_a3",
                    "origin": "eumetsat",
                    "variable": "cloud_physical_properties_of_the_ice_phase",
                    "climate_data_record_type": "thematic_climate_data_record",
                    "time_aggregation": "monthly_mean",
                    "year": str(year),
                    "month": [f"{m:02d}" for m in range(1, 13)],
                },
                "clt_monthly_" + str(year): {
                    "product_family": "clara_a3",
                    "origin": "eumetsat",
                    "variable": "cloud_fraction",
                    "climate_data_record_type": "thematic_climate_data_record",
                    "time_aggregation": "monthly_mean",
                    "year": str(year),
                    "month": [f"{m:02d}" for m in range(1, 13)],
                },
                "lwp_monthly_" + str(year): {
                    "product_family": "clara_a3",
                    "origin": "eumetsat",
                    "variable": "cloud_physical_properties_of_the_liquid_phase",
                    "climate_data_record_type": "thematic_climate_data_record",
                    "time_aggregation": "monthly_mean",
                    "year": str(year),
                    "month": [f"{m:02d}" for m in range(1, 13)],
                },
            }
        )
        for month in range(1, 13):
            requests.update(
                {
                    "clivi_daily_" + str(year) + f"{month:02d}": {
                        "product_family": "clara_a3",
                        "origin": "eumetsat",
                        "variable": "cloud_physical_properties_of_the_ice_phase",
                        "climate_data_record_type": "thematic_climate_data_record",
                        "time_aggregation": "daily_mean",
                        "year": str(year),
                        "month": f"{month:02d}",
                        "day": [f"{m:02d}" for m in range(1, 32)],
                    },
                    "clt_daily_" + str(year) + f"{month:02d}": {
                        "product_family": "clara_a3",
                        "origin": "eumetsat",
                        "variable": "cloud_fraction",
                        "climate_data_record_type": "thematic_climate_data_record",
                        "time_aggregation": "daily_mean",
                        "year": str(year),
                        "month": f"{month:02d}",
                        "day": [f"{m:02d}" for m in range(1, 32)],
                    },
                    "lwp_daily_" + str(year) + f"{month:02d}": {
                        "product_family": "clara_a3",
                        "origin": "eumetsat",
                        "variable": "cloud_physical_properties_of_the_liquid_phase",
                        "climate_data_record_type": "thematic_climate_data_record",
                        "time_aggregation": "daily_mean",
                        "year": str(year),
                        "month": f"{month:02d}",
                        "day": [f"{m:02d}" for m in range(1, 32)],
                    },
                }
            )

    cds_client = cdsapi.Client(cds_url)

    for var_name, request in requests.items():
        datestr = var_name.split("_")[2]
        if "daily" in var_name:
            outdir = output_folder / f"daily/{datestr}/"
        else:
            outdir = output_folder / f"monthly/{datestr}/"
        outdir.mkdir(parents=True, exist_ok=True)

        logger.info("Downloading %s data to %s", var_name, outdir)

        file_path = outdir / f"{var_name}.gz"

        if file_path.exists() and not overwrite:
            logger.info(
                "File %s already exists. Skipping download.",
                file_path,
            )
            continue

        try:
            cds_client.retrieve(
                "satellite-cloud-properties",
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
        except Exception as ex:
            logger.info("%s: no data downloaded for %s", type(ex), var_name)
