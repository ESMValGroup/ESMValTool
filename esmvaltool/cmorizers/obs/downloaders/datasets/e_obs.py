"""Script to download E-OBS from its webpage."""
import logging

from esmvaltool.cmorizers.obs.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(config, dataset, start_date, end_date, overwrite):
    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )
    for var in ['TG', 'TN', 'TX', 'RR', 'PP']:
        for grid in ('0.1deg', '0.25deg'):
            for version in ('20.0e', ):
                downloader.download_file(
                    "https://knmi-ecad-assets-prd.s3.amazonaws.com/ensembles/"
                    f"data/Grid_{grid}_reg_ensemble/"
                    f"{var.lower()}_ens_mean_{grid}_reg_v{version}.nc",
                    wget_options=[]
                )
