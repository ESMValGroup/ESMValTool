# # """Script to download ESACCI-OZONE from the CDS."""

from pathlib import Path
import shutil
import gzip
import zipfile
import cdsapi


def download_dataset(config, dataset, dataset_info, start_date, end_date,
                     overwrite):
    """Download ESACCI-OZONE dataset using CDS API

        An ECMWF account is needed to download the datasets from:
        https://cds.climate.copernicus.eu/datasets/satellite-ozone-v1

        - The file named .cdspirc containing the key associated to
        the ECMWF account needs to be saved in user's ${HOME} directory.

        - All the files will be saved in ${RAWOBS}/Tier2/ESACCI-OZONE
        """

    cds_url = "https://cds.climate.copernicus.eu/api"

    if dataset == "ESACCI-OZONE":
        requests = {
            "toz": {
                "processing_level": "level_3",
                "variable": "atmosphere_mole_content_of_ozone",
                "vertical_aggregation": "total_column",
                "sensor": ["merged_uv"],
                "year": [str(y) for y in range(1995, 2023)],
                "month": [f"{m:02d}" for m in range(1, 13)],
                "version": ["v2000"]
            },
            "o3": {
                "processing_level": "level_3",
                "variable": "mole_concentration_of_ozone_in_air",
                "vertical_aggregation": "vertical_profiles_from_limb_sensors",
                "sensor": ["cmzm"],
                "year": [str(y) for y in range(1984, 2022)],
                "month": [f"{m:02d}" for m in range(1, 13)],
                "version": ["v0008"]
            }
        }

        client = cdsapi.Client(cds_url)
        raw_obs_dir = Path(config['rootpath']['RAWOBS'][0])
        output_folder = raw_obs_dir / f"Tier{dataset_info['tier']}" / dataset
        output_folder.mkdir(parents=True, exist_ok=True)

        for var_name, request in requests.items():
            print(f"Downloading {var_name} data to {output_folder}...")

            file_path = output_folder / f"{var_name}.gz"

            if file_path.exists() and not overwrite:
                print(f"File {file_path} already exists. Skipping download.")
                continue

            client.retrieve("satellite-ozone-v1", request,
                            file_path.as_posix())

            # Handle both .gz and .zip files
            with open(file_path, "rb") as f:
                magic = f.read(2)

            if magic == b'PK':  # ZIP file signature
                print(f"Detected ZIP file: {file_path}")
                with zipfile.ZipFile(file_path, 'r') as zip_ref:
                    zip_ref.extractall(output_folder)
            else:
                print(f"Detected GZIP file: {file_path}")
                with gzip.open(file_path, 'rb') as f_in:
                    with open(output_folder / file_path.stem, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

    else:
        raise ValueError(f"Unknown dataset: {dataset}")
