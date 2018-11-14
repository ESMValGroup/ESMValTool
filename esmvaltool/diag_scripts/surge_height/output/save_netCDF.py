import datetime
import os

import numpy as np
import xarray as xr


def save_netCDF(dates, stat, srg_est_full, cfg, dataset):
    # define coordinates, i.e. time and stations
    coord = {
        'northcor': [1.25, 61.3328],
        'f3': [4.75, 54.91639],
        'lauwerso': [6.375, 53.49978],
        'helgolan': [8.0, 54.24975],
        'helgeroa': [9.875, 59.08289],
        'immingha': [0.25, 53.74977],
        'goeree': [3.75, 51.99984],
        'bremerha': [8.625, 53.66644],
        'scheveni': [4.375, 52.24983],
        'os15': [3.625, 51.66652],
        'devonpor': [-4.0, 50.41657],
        'oscarsbo': [10.75, 59.58287],
        'leith': [-3.0, 56.08301],
        'os11': [3.625, 51.74985],
        'newlyn': [-5.375, 50.16658],
        'harlinge': [5.5, 53.33312],
        'vlaktevd': [3.25, 51.58319],
        'europlat': [3.375, 52.08317],
        'lowestof': [1.875, 52.58315],
        'felixsto': [1.5, 51.99984],
        'aberdeen': [-1.875, 57.24963],
        'portsmou': [-1.0, 50.91655],
        'newhaven': [0.0, 50.83322],
        'hoekvanh': [4.25, 52.08317],
        'holyhead': [-4.5, 53.49978],
        'westters': [5.25, 53.41645],
        'esbjerg': [8.375, 55.58303],
        'lerwick': [-0.875, 60.58283],
        'bg2': [3.75, 51.83318],
        'husum': [9.0, 54.49974],
        'terschel': [5.375, 53.49978],
        'stornowa': [-5.625, 58.16626],
        'wick': [-2.875, 58.49958],
        'torsmind': [8.125, 56.41633],
        'scarboro': [-0.25, 54.41641],
        'ilfracom': [-4.0, 51.3332],
        'roompotb': [3.75, 51.74985],
        'northshi': [-1.25, 55.08305],
        'delfzijl': [7.125, 53.41645],
        'k13a': [3.375, 53.33312],
        'denhelde': [5.25, 53.41645],
        'cadzand': [3.375, 51.49986],
        'ekofisk': [3.375, 56.66632],
        'borkums': [6.75, 53.66644],
        'vlissing': [3.625, 51.49986],
        'cromer': [1.5, 52.9998],
        'stmarys': [-6.125, 49.99992],
        'oostende': [3.0, 51.3332],
        'aukalpha': [2.25, 56.49966],
        'southend': [0.875, 51.58319],
        'huibertg': [6.5, 53.66644],
        'denoever': [5.125, 53.08313],
        'texelnoo': [4.875, 53.24979],
        'sheernes': [1.0, 51.49986],
        'innerdow': [0.625, 53.41645],
        'ijmuiden': [4.625, 52.49982],
        'kornwerd': [5.25, 53.08313],
        'meetpost': [4.375, 52.41649],
        'cuxhaven': [8.875, 53.99976],
        'tregde': [7.75, 58.08293],
        'vidaa': [8.75, 55.08305],
        'westkape': [3.5, 51.58319],
        'duinkerk': [2.5, 51.16654],
        'offharwi': [1.625, 51.91651],
        'weymouth': [-2.375, 50.66656],
        'stavange': [5.625, 59.08289],
        'zeebrugg': [3.25, 51.41653],
        'scillyis': [-6.375, 49.91659],
        'dover': [1.5, 51.24987]
    }

    srg = []
    stats_lon = []
    stats_lat = []
    stat_out = []
    for s in srg_est_full.keys():
        srg.append(srg_est_full[s])
        stats_lat.append(coord[s][0])
        stats_lon.append(coord[s][1])
    #
    now = datetime.datetime.utcnow().strftime("%a %Y-%m-%d %H:%M")
    srg_out = xr.DataArray(
        srg,
        coords=[np.asarray(stat).astype('<U8'), dates],
        dims=['station', 'time'])
    srg_out.attrs['units'] = 'meters'
    lons = xr.DataArray(
        stats_lon, coords=[np.asarray(stat).astype('<U8')], dims=['station'])
    lons.attrs['units'] = 'Degrees East'
    lats = xr.DataArray(
        stats_lat, coords=[np.asarray(stat).astype('<U8')], dims=['station'])
    lats.attrs['units'] = 'Degrees West'
    out = xr.Dataset({
        'surge': srg_out,
        'longitude': lons,
        'latitude': lats
    },
                     attrs={
                         'description':
                         'Estimate of surge height along the North Sea coast.',
                         'history':
                         'created ' + now,
                         'featureType':
                         'timeSeries'
                     })
    #
    path = os.path.join(
        cfg["work_dir"],
        dataset + '_' + cfg['fout_name'] + '_' + dates[0].strftime("%Y-%m-%d")
        + '-' + dates[-1].strftime("%Y-%m-%d") + '.nc',
    )

    out.to_netcdf(path)
