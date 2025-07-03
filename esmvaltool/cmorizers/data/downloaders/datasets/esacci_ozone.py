"""Script to download ESACCI-OZONE from the CDS and BIRA WebDAV."""

import gzip
import logging
import os
import shutil
import zipfile
from datetime import datetime
from pathlib import Path

import cdsapi
import webdav.client as wc
from dateutil import relativedelta

logger = logging.getLogger(__name__)


def download_dataset(
    config, dataset, dataset_info, start_date, end_date, overwrite,
):
    """Download ESACCI-OZONE dataset using CDS API.

    - An ECMWF account is needed to download the datasets from
      https://cds.climate.copernicus.eu/datasets/satellite-ozone-v1.
    - The file named .cdspirc containing the key associated to
      the ECMWF account needs to be saved in user's ${HOME} directory.
    - All the files will be saved in ${RAWOBS}/Tier2/ESACCI-OZONE.
    """
    if dataset == "ESACCI-OZONE":
        raw_obs_dir = Path(config["rootpath"]["RAWOBS"][0])
        output_folder = raw_obs_dir / f"Tier{dataset_info['tier']}" / dataset
        output_folder.mkdir(parents=True, exist_ok=True)

        cds_url = "https://cds.climate.copernicus.eu/api"

        requests = {
            "toz_gto_ecv": {
                "processing_level": "level_3",
                "variable": "atmosphere_mole_content_of_ozone",
                "vertical_aggregation": "total_column",
                "sensor": ["merged_uv"],
                "year": [str(y) for y in range(1995, 2024)],
                "month": [f"{m:02d}" for m in range(1, 13)],
                "version": ["v2000"],
            },
            "o3_sage_omps": {
                "processing_level": "level_3",
                "variable": "mole_concentration_of_ozone_in_air",
                "vertical_aggregation": "vertical_profiles_from_limb_sensors",
                "sensor": ["cmzm"],
                "year": [str(y) for y in range(1984, 2023)],
                "month": [f"{m:02d}" for m in range(1, 13)],
                "version": ["v0008"],
            },
            "o3_sage_megridop": {
                "processing_level": "level_3",
                "variable": "mole_concentration_of_ozone_in_air",
                "vertical_aggregation": "vertical_profiles_from_limb_sensors",
                "sensor": ["cllg"],
                "year": [str(y) for y in range(2001, 2025)],
                "month": [f"{m:02d}" for m in range(1, 13)],
                "version": ["v0005"],
            },
        }

#        client = cdsapi.Client(cds_url)
#
#        for var_name, request in requests.items():
#            logger.info("Downloading %s data to %s", var_name, output_folder)
#
#            file_path = output_folder / f"{var_name}.gz"
#
#            if file_path.exists() and not overwrite:
#                logger.info(
#                    "File %s already exists. Skipping download.", file_path,
#                )
#                continue
#
#            client.retrieve(
#                "satellite-ozone-v1", request, file_path.as_posix(),
#            )
#
#            # Handle both .gz and .zip files
#            with open(file_path, "rb") as file:
#                magic = file.read(2)
#
#            if magic == b"PK":  # ZIP file signature
#                logger.info("Detected ZIP file: %s", file_path)
#                with zipfile.ZipFile(file_path, "r") as zip_ref:
#                    zip_ref.extractall(output_folder)
#            else:
#                logger.info("Detected GZIP file: %s", file_path)
#                with gzip.open(file_path, "rb") as f_in:
#                    with open(output_folder / file_path.stem, "wb") as f_out:
#                        shutil.copyfileobj(f_in, f_out)

        # download IASI data from BIRA WebDAV (IASI data not available on CDS)
        # all the files will be saved by year (yyyy) in
        # ${RAWOBS}/Tier2/ESACCI-OZONE/IASI_yyyy

        if start_date is None:
            start_date = datetime(2008, 1, 1)
        if end_date is None:
            end_date = datetime(2023, 12, 31)

        options = {
            "webdav_hostname": "https://webdav.aeronomie.be",
            "webdav_login": "o3_cci_public",
            "webdav_password": "",
        }

        client = wc.Client(options)

        basepath = "/guest/o3_cci/webdata/Nadir_Profiles/L3/IASI_MG_FORLI/"

        loop_date = start_date
        while loop_date <= end_date:
            year = loop_date.year

            # if needed, create local output directory
            outdir = output_folder / f"IASI_{year}"
            os.makedirs(outdir, exist_ok=True)

            # directory on WebDAV server to download
            remotepath = f"{basepath}/{year}"
            files = client.list(remotepath)
            info = client.info(remotepath + "/" + files[0])
            numfiles = len(files)
            # calculate approx. download volume in Gbytes
            size = int(info["size"]) * numfiles // 1073741824
            del files

            loginfo = (
                f"downloading {numfiles} files for year {year}"
                f" (approx. {size} Gbytes)"
            )
            logger.info(loginfo)

            # synchronize local (output) directory and WebDAV server directory
            client.pull(remote_directory=remotepath, local_directory=outdir)

            loop_date += relativedelta.relativedelta(years=1)

    else:
        errmsg = f"Unknown dataset: {dataset}"
        raise ValueError(errmsg)
