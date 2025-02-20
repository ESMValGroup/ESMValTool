# """Script to download ESACCI-OZONE from the CDS."""

import cdsapi
from pathlib import Path
from esmvaltool.cmorizers.data.utilities import unpack_files_in_folder


def download_dataset(config, dataset, dataset_info, start_date, end_date,
                     overwrite):
    """Download ESACCI-OZONE dataset using CDS API and create tar
     file for each variable."""

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

        client = cdsapi.Client()
        output_folder = Path(config["rootpath"][dataset][0])
        output_folder.mkdir(parents=True, exist_ok=True)

        for var_name, request in requests.items():
            print(f"Downloading {var_name} data to {output_folder}...")

            # Get the expected filename from the CDS API request
            file_path = output_folder / f"{var_name}.tar.gz"

            if file_path.exists() and not overwrite:
                print(f"File {file_path} already exists. Skipping download.")
                continue

            # Download file to the specified path
            client.retrieve("satellite-ozone-v1", request,
                            file_path.as_posix())

        unpack_files_in_folder(output_folder)

    else:
        raise ValueError(f"Unknown dataset: {dataset}")
