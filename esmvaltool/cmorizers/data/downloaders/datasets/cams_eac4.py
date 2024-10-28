"""Script to download CAMS data from the Climate Data Store."""

import datetime

from esmvaltool.cmorizers.data.downloaders.cds import CDSDownloader


def download_dataset(config, dataset, dataset_info, start_date, end_date,
                     overwrite):
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
        start_date = datetime.datetime(2003, 1, 1)
    if end_date is None:
        end_date = datetime.datetime(2023, 12, 31)
    #    #end_date = datetime.datetime(2022, 12, 1)

    variables = [
        'geopotential', 'hydroxyl_radical', 'methane_chemistry', 'ozone',
        'specific_humidity', 'vertical_velocity', 'nitrogen_monoxide'
    ]

    for var in variables:

        downloader = CDSDownloader(
            product_name='cams-global-reanalysis-eac4-monthly',
            request_dictionary={
                'variable': [var],
                'pressure_level': [
                    "1", "2", "3", "5", "7", "10", "20", "30", "50", "70",
                    "100", "150", "200", "250", "300", "400", "500", "600",
                    "700", "800", "850", "900", "925", "950", "1000"
                ],
                "product_type": ["monthly_mean"],
                "year": [
                    "2003", "2004", "2005", "2006", "2007", "2008", "2009",
                    "2010", "2011", "2012", "2013", "2014", "2015", "2016",
                    "2017", "2018", "2019", "2020", "2021", "2022", "2023"
                ],
                "month": [
                    "01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
                    "11", "12"
                ],
                "data_format":
                'grib'
            },
            config=config,
            dataset=dataset,
            dataset_info=dataset_info,
            overwrite=overwrite,
        )

        downloader.download_request(f"CAMS_EAC4_{var}_2003_2023.grib")
