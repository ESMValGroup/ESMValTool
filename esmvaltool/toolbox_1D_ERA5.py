import cdstoolbox as ct

@ct.application(title='Retrieve Data')
@ct.output.dataarray()
def retrieve_sample_data():
    """
    Retrieve a variable from sample dataset.
    """

    data = ct.catalogue.retrieve(
    'reanalysis-era5-single-levels',
    {
        'variable':'sea_surface_temperature',
        'product_type':'reanalysis',
        'year':'2017',
        'month':'04',
        'day':'10',
        'time':'08:00'
    }
)
    datas = ct.climate.daily_mean(data)
    return datas
