
@ct.application(title='Retrieve Data')
@ct.output.dataarray()
def retrieve_sample_data():
    """
    Application main steps:
    
    - retrieve a variable from CDS Catalogue
    - produce a link to download it in original data format (usually 'grib' if not otherwise specified).
    """

    data = ct.catalogue.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'variable': 'temperature',
            'product_type': 'reanalysis',
            'pressure_level':['50','100','500','1000'],
            'year': '2010',
            'month': '08',
            'day': ['01', '02', '03', '04', '05', '06',
                '07', '08', '09', '10', '11', '12',
                '13', '14', '15', '16', '17', '18',
                '19', '20', '21', '22', '23', '24',
                '25', '26', '27', '28', '29', '30',
                '31'],
            'time': ['00:00','12:00']
       	        }
    )

    datas = ct.climate.daily_mean(data)
return datas
