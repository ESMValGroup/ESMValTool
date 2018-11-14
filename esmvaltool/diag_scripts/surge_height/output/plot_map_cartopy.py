import os
from datetime import datetime

import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np


def plot_map_cartopy(date, srg, tidx, cfg, dataset):
    # coordinates of stations
    coords = {
        'aberdeen': [81, 111],
        'aukalpha': [114, 102],
        'bg2': [126, 46],
        'borkums': [150, 68],
        'bremerha': [165, 68],
        'cadzand': [123, 42],
        'cromer': [108, 60],
        'cuxhaven': [167, 72],
        'delfzijl': [153, 65],
        'denhelde': [138, 65],
        'denoever': [137, 61],
        'devonpor': [64, 29],
        'dover': [108, 39],
        'duinkerk': [116, 38],
        'ekofisk': [123, 104],
        'esbjerg': [163, 91],
        'europlat': [123, 49],
        'f3': [134, 83],
        'felixsto': [108, 48],
        'goeree': [126, 48],
        'harlinge': [140, 64],
        'helgeroa': [175, 133],
        'helgolan': [160, 75],
        'hoekvanh': [130, 49],
        'holyhead': [60, 66],
        'huibertg': [148, 68],
        'husum': [168, 78],
        'ijmuiden': [133, 54],
        'ilfracom': [64, 40],
        'immingha': [98, 69],
        'innerdow': [101, 65],
        'k13a': [123, 64],
        'kornwerd': [138, 61],
        'lauwerso': [147, 66],
        'leith': [72, 97],
        'lerwick': [89, 151],
        'lowestof': [111, 55],
        'meetpost': [131, 53],
        'newhaven': [96, 34],
        'newlyn': [53, 26],
        'northcor': [106, 160],
        'northshi': [86, 85],
        'offharwi': [109, 47],
        'oostende': [120, 40],
        'os11': [125, 45],
        'os15': [125, 44],
        'oscarsbo': [182, 139],
        'portsmou': [88, 35],
        'roompotb': [126, 45],
        'scarboro': [94, 77],
        'scheveni': [131, 51],
        'scillyis': [45, 23],
        'sheernes': [104, 42],
        'southend': [103, 43],
        'stavange': [141, 133],
        'stmarys': [47, 24],
        'stornowa': [51, 122],
        'terschel': [139, 66],
        'texelnoo': [135, 63],
        'torsmind': [161, 101],
        'tregde': [158, 121],
        'vidaa': [166, 85],
        'vlaktevd': [122, 43],
        'vlissing': [125, 42],
        'westkape': [124, 43],
        'westters': [138, 65],
        'weymouth': [77, 32],
        'wick': [73, 126],
        'zeebrugg': [122, 41]
    }
    # create map
    fig = plt.figure()
    #ax  = plt.axes(projection=ccrs.cfg['map_projection'])
    ax1 = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax1.coastlines('50m')
    ax1.set_extent([-7, 11.5, 49.4, 61.75], ccrs.PlateCarree())
    ax1.add_feature(cartopy.feature.OCEAN, zorder=0)
    ax1.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
    ax1.gridlines()
    lons = np.arange(-12.00, -12.00 + (0.125 * 201), 0.125).tolist()
    lats = np.arange(48.00, 48.00 + (0.08333 * 173), 0.08333).tolist()

    X, Y = np.meshgrid(lons, lats)

    #
    for stat in srg.keys():
        scat = plt.scatter(
            lons[coords[stat][0]],
            lats[coords[stat][1]],
            c=srg[stat][tidx],
            edgecolors='k',
            cmap=plt.get_cmap(cfg["colormap"], 20),
            vmin=-3,
            vmax=3,
            zorder=50)
    #
    cbar_scat = plt.colorbar(scat)
    cbar_scat.set_label('surge height (m)')
    plt.title('North Sea coastal surge ' +
              datetime.strftime(date, '(%d-%m-%Y)'))
    #
    fdates = datetime.strftime(date, '%Y-%m-%d')
    fsave = os.path.join(
        cfg["plot_dir"],
        dataset + '_' + cfg["fout_name"] + '_map_' + fdates + '.pdf')
    plt.savefig(fsave, dpi=100, format='pdf')
