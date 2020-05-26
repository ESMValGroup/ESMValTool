"""Script to download cds-satellite-albedo from the Climate Data Store(CDS)"""

from dateutil import relativedelta

from esmvaltool.cmorizers.obs.downloaders.ftp import CCIDownloader


def download_dataset(config, dataset, start_date, end_date, overwrite):
    """Download dataset cds-satellite-albedo."""
    loop_date = start_date

    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )
    downloader.connect()
    end_of_file = 'ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA-12-fv3.0.nc'
    filler_data = {
        1994: [
            f'AVHRR_NOAA_12/1994/199409-{end_of_file}',
            f'AVHRR_NOAA_12/1994/199410-{end_of_file}',
            f'AVHRR_NOAA_12/1994/199411-{end_of_file}',
            f'AVHRR_NOAA_12/1994/199412-{end_of_file}',
        ],
        1995:  [
            f'AVHRR_NOAA_12/1995/199501-{end_of_file}',
        ],
    }

    while loop_date <= end_date:
        year = loop_date.year
        downloader.set_cwd('version3/L3C/AVHRR-PM/v3.0')
        for folder in downloader.list_folders():
            for year_folder in downloader.list_folders(folder):
                if int(year_folder) == year:
                    downloader.download_year(f'{folder}/{year_folder}')
        downloader.set_cwd('version3/L3C/AVHRR-AM/v3.0')
        for extra_file in filler_data.get(year, []):
            downloader.download_file(extra_file)
        loop_date += relativedelta.relativedelta(years=1)
