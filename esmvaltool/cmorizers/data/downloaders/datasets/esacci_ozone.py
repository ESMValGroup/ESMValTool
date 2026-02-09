"""Script to download ESACCI-OZONE from the CDS."""

import gzip
import logging
import shutil
import zipfile

import cdsapi

logger = logging.getLogger(__name__)


def download_dataset(
    original_data_dir,
    dataset,
    dataset_info,
    start_date,
    end_date,
    overwrite,
):
    """Download ESACCI-OZONE dataset using CDS API.

    - An ECMWF account is needed to download the datasets from
      https://cds.climate.copernicus.eu/datasets/satellite-ozone-v1.
    - The file named .cdspirc containing the key associated to
      the ECMWF account needs to be saved in user's ${HOME} directory.
    - All the files will be saved in Tier2/ESACCI-OZONE.
    """
    cds_url = "https://cds.climate.copernicus.eu/api"

    if dataset == "ESACCI-OZONE":
        requests = {
            "toz": {
                "processing_level": "level_3",
                "variable": "atmosphere_mole_content_of_ozone",
                "vertical_aggregation": "total_column",
                "sensor": ["merged_uv"],
                "year": [str(y) for y in range(1995, 2024)],
                "month": [f"{m:02d}" for m in range(1, 13)],
                "version": ["v2000"],
            },
            "o3": {
                "processing_level": "level_3",
                "variable": "mole_concentration_of_ozone_in_air",
                "vertical_aggregation": "vertical_profiles_from_limb_sensors",
                "sensor": ["cmzm"],
                "year": [str(y) for y in range(1984, 2023)],
                "month": [f"{m:02d}" for m in range(1, 13)],
                "version": ["v0008"],
            },
        }

        client = cdsapi.Client(cds_url)
        output_folder = (
            original_data_dir / f"Tier{dataset_info['tier']}" / dataset
        )
        output_folder.mkdir(parents=True, exist_ok=True)

        for var_name, request in requests.items():
            logger.info("Downloading %s data to %s", var_name, output_folder)

            file_path = output_folder / f"{var_name}.gz"

            if file_path.exists() and not overwrite:
                logger.info(
                    "File %s already exists. Skipping download.",
                    file_path,
                )
                continue

            client.retrieve(
                "satellite-ozone-v1",
                request,
                file_path.as_posix(),
            )

            # Handle both .gz and .zip files
            with open(file_path, "rb") as file:
                magic = file.read(2)

            if magic == b"PK":  # ZIP file signature
                logger.info("Detected ZIP file: %s", file_path)
                with zipfile.ZipFile(file_path, "r") as zip_ref:
                    zip_ref.extractall(output_folder)
            else:
                logger.info("Detected GZIP file: %s", file_path)
                with gzip.open(file_path, "rb") as f_in:
                    with open(output_folder / file_path.stem, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)

    else:
        raise ValueError(f"Unknown dataset: {dataset}")
