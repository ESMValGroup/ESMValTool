"""Diagnostic to select grid points within a shapefile."""
import logging
import os
import sys
from copy import deepcopy

import fiona
# import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap
import iris
# import shapefile
import numpy as np
import xlsxwriter
from matplotlib import pyplot as p
from netCDF4 import Dataset, num2date, stringtochar
# from shapely.geometry import shape, mapping, Point, Polygon, MultiPolygon
from shapely.geometry import MultiPoint, shape
from shapely.ops import nearest_points

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Select grid points within shapefiles."""
    for filename, attributes in cfg['input_data'].items():
        logger.info("Processing variable %s from dataset %s",
                    attributes['standard_name'], attributes['dataset'])
        logger.debug("Loading %s", filename)
        cube = iris.load_cube(filename)
        ncts, nclon, nclat = shapeselect(cfg, cube, filename)
        name = os.path.splitext(os.path.basename(filename))[0] + '_polygon'
        if cfg['write_xlsx']:
            writexls(cfg, filename, ncts, nclon, nclat)
        if cfg['write_netcdf']:
            path = os.path.join(
                cfg['work_dir'],
                name + '.nc',
            )
        write_netcdf(path, ncts, nclon, nclat, cube, cfg)


def write_keyvalue_toxlsx(worksheet, row, key, value):
    if type(value).__name__ == 'dict':
        worksheet.write(row, 0, key)
        row += 1
        for dictkey, dictvalue in value.items():
            row = write_keyvalue_toxlsx(worksheet, row, dictkey, dictvalue)
    elif type(value).__name__ == 'list':
        for listvalue in value:
            row = write_keyvalue_toxlsx(worksheet, row, key, listvalue)
    else:
        worksheet.write(row, 0, key)
        if type(value).__name__ == 'bool':
            worksheet.write(row, 1, str(int(value)))
        else:
            worksheet.write(row, 1, value)
        row += 1

    return row


def writexls(cfg, filename, ncts, nclon, nclat):
    """Write the content of a netcdffile as .xls"""
    name = os.path.splitext(os.path.basename(filename))[0] + '_polygon'
    ncfile = Dataset(filename, 'r')
    nclon = ncfile.variables['lon']
    nclat = ncfile.variables['lat']
    otime = ncfile.variables['time']
    dtime = num2date(otime[:], otime.units, otime.calendar)
    wtime = []
    for tim in range(len(dtime)):
        wtime.append(str(dtime[tim]))
    pathx = os.path.join(cfg['work_dir'], name + '.xlsx')
    pathc = os.path.join(cfg['work_dir'], name + '.csv')
    workbook = xlsxwriter.Workbook(pathx)
    worksheet = workbook.add_worksheet('Data')
    worksheet.write(0, 0, 'Date')
    worksheet.write(0, 1, 'Lon/Lat')
    worksheet.write_column(2, 0, wtime)
    for row in range(ncts.shape[1]):
        worksheet.write(
            1, row + 1,
            str("%#.3f" % round(float(nclon[row]), 3)) + '_' + str(
                "%#.3f" % round(float(nclat[row]), 3)))
        worksheet.write_column(2, row + 1,
                               np.around(np.squeeze(ncts[:, row]), decimals=8))
        worksheet.set_column(0, row + 1, 20)
    worksheet = workbook.add_worksheet('NetCDFheader')
    worksheet.set_column(0, 0, 20)
    for row, attr in enumerate(ncfile.ncattrs()):
        worksheet.write(row, 0, attr)
        worksheet.write(row, 1, getattr(ncfile, attr))
    worksheet = workbook.add_worksheet('ESMValTool')
    worksheet.set_column(0, 0, 20)
    row = 0
    for key, value in cfg.items():
        print(key, value)
        row = write_keyvalue_toxlsx(worksheet, row, key, value)
    workbook.close()


def dates(otime):
    """ Converts iris time to human readable format """
    calendar = otime.units.calendar
    origin = otime.units.origin


def shapeselect(cfg, cube, filename):
    """ New solution with fiona! """
    shppath = cfg['shppath']
    wgtmet = cfg['wgtmet']
    if ((cube.coord('latitude').ndim == 1
         and cube.coord('longitude').ndim == 1)):
        coordpoints = [(x, y) for x in cube.coord('longitude').points
                       for y in cube.coord('latitude').points]
    else:
        logger.info("Support for 2-d coords not yet implemented!")
        sys.exit(1)
    points = MultiPoint(coordpoints)
    if cfg['evalplot']:
        # Set limits for map (This can definitely be improved!)
        shap = fiona.open(shppath)
        llcrnrlon = shap.bounds[0] - 1
        llcrnrlat = max((shap.bounds[1] - 1, -90))
        urcrnrlon = shap.bounds[2] + 1
        urcrnrlat = min((shap.bounds[3] + 1, 90))
        alons = []
        alats = []
        for lon, lat in coordpoints:
            if (lat > llcrnrlat and lat < urcrnrlat and lon > llcrnrlon
                    and lon < urcrnrlon):
                alons.append(lon)
                alats.append(lat)
        cnt = 0
        # rngx = np.arange(np.min(alons), np.max(alons), 1)
        # rngy = np.arange(np.min(alats), np.max(alats), 1)
        p.gcf().clear()
        # ax.set_xticks(rngx)
        # ax.set_yticks(rngy)
        # for shapp in shap.shapeRecords(): #shape.next()
        for shapp in shap:
            xm = [i[0] for i in shapp['geometry']['coordinates'][0][:]]
            ym = [i[1] for i in shapp['geometry']['coordinates'][0][:]]
            p.plot(xm, ym)
            p.plot(alons, alats, 'ro', markersize=2)
            cnt += 1
        p.xlabel('Longitude')
        p.ylabel('Latitude')
    with fiona.open(shppath) as shp:
        gpx = []
        gpy = []
        cnt = -1
        ncts = np.zeros((cube.coord('time').shape[0], len(shp)))
        nclon = np.zeros((len(shp)))  # Takes representative point
        nclat = np.zeros((len(shp)))
        for ishp, multipol in enumerate(shp):
            cnt += 1
            multi = shape(multipol['geometry'])
            if wgtmet == 'mean_inside':
                gpx, gpy = mean_inside(gpx, gpy, points, multi, cube, cfg)
                if len(gpx) == 0:
                    gpx, gpy = representative(gpx, gpy, points, multi, cube,
                                              cfg)
                    pth = 'g+'
                else:
                    pth = 'b+'
            elif wgtmet == 'representative':
                gpx, gpy = representative(gpx, gpy, points, multi, cube, cfg)
            if cfg['evalplot']:
                p.plot(
                    cube.coord('longitude').points[gpx],
                    cube.coord('latitude').points[gpy],
                    pth,
                    markersize=10)
            if len(gpx) == 1:
                ncts[:, ishp] = np.reshape(cube.data[:, gpy, gpx],
                                           (cube.data.shape[0], ))
            else:
                ncts[:, ishp] = np.mean(cube.data[:, gpy, gpx], axis=1)
            gx, gy = representative([], [], points, multi, cube, cfg)
            nclon[ishp] = cube.coord('longitude').points[gx]
            nclat[ishp] = cube.coord('latitude').points[gy]
    if cfg['evalplot']:
        name = os.path.splitext(os.path.basename(filename))[0]
        path = os.path.join(cfg['plot_dir'], name + '.png')
        p.savefig(path)
    return ncts, nclon, nclat


def mean_inside(gpx, gpy, points, multi, cube, cfg):
    for point in points:
        if point.within(multi):
            if point.x < 0 and np.min(cube.coord('longitude').points) >= 0:
                addx = 360.
            else:
                addx = 0.
            x, y = best_match(
                cube.coord('longitude').points,
                cube.coord('latitude').points, point.x + addx, point.y)
            gpx.append(x)
            gpy.append(y)
    return gpx, gpy


def representative(gpx, gpy, points, multi, cube, cfg):
    reprpoint = multi.representative_point()
    nearest = nearest_points(reprpoint, points)
    npx = nearest[1].coords[0][0]
    npy = nearest[1].coords[0][1]
    if npx < 0 and np.min(cube.coord('longitude').points) >= 0:
        npx = npx + 360.
    x, y = best_match(
        cube.coord('longitude').points,
        cube.coord('latitude').points, npx, npy)
    gpx.append(x)
    gpy.append(y)
    return gpx, gpy


def best_match(i, j, px, py):
    """ Identifies the grid points in 2-d with minimum distance. """
    if i.shape != 2 or j.shape != 2:
        x = deepcopy(i)
        y = deepcopy(j)
        xx = np.zeros((len(x), len(y)))
        yy = np.zeros((len(x), len(y)))
        for i in range(0, len(y)):
            xx[:, i] = x
        for j in range(0, len(x)):
            yy[j, :] = y
    else:
        xx = deepcopy(i)
        yy = deepcopy(j)
    distance = ((xx - px)**2 + (yy - py)**2)**0.5
    ind = np.unravel_index(np.argmin(distance, axis=None), distance.shape)
    return ind[0], ind[1]


def write_netcdf(path, var, plon, plat, cube, cfg):
    """Write results to a netcdf file."""
    shppath = cfg['shppath']
    # wgtmet = cfg['wgtmet']
    polyid = []
    for row in range(var.shape[1]):
        polyid.append(
            str("%#.3f" % round(plon[row], 3)) + '_' +
            str("%#.3f" % round(plat[row], 3)))
    # polyid = stringtochar(np.array(polyid, 'S4'))
    # print(polyid)
    ncout = Dataset(path, mode='w')
    ncout.createDimension('time', None)
    ncout.createDimension('polygon', len(polyid))
    times = ncout.createVariable('time', 'f8', ('time'), zlib=True)
    times.setncattr_string('standard_name', cube.coord('time').standard_name)
    times.setncattr_string('long_name', cube.coord('time').long_name)
    times.setncattr_string('calendar', cube.coord('time').units.calendar)
    times.setncattr_string('units', cube.coord('time').units.origin)
    polys = ncout.createVariable('polygon', 'S1', ('polygon'), zlib=True)
    polys.setncattr_string('standard_name', 'polygon')
    polys.setncattr_string('long_name', 'polygon')
    polys.setncattr_string('shapefile', shppath)
    lon = ncout.createVariable(
        cube.coord('longitude').var_name, 'f8', 'polygon', zlib=True)
    lon.setncattr_string('standard_name',
                         cube.coord('longitude').standard_name)
    lon.setncattr_string('long_name', cube.coord('longitude').long_name)
    lon.setncattr_string('units', cube.coord('longitude').units.origin)
    lat = ncout.createVariable(
        cube.coord('latitude').var_name, 'f8', 'polygon', zlib=True)
    lat.setncattr_string('standard_name', cube.coord('latitude').standard_name)
    lat.setncattr_string('long_name', cube.coord('latitude').long_name)
    lat.setncattr_string('units', cube.coord('latitude').units.origin)
    data = ncout.createVariable(
        cube.var_name, 'f4', ('time', 'polygon'), zlib=True)
    data.setncattr_string('standard_name', cube.standard_name)
    data.setncattr_string('long_name', cube.long_name)
    data.setncattr_string('units', cube.units.origin)
    for key, val in cube.metadata[-2].items():
        ncout.setncattr_string(key, val)
    times[:] = cube.coord('time').points
    lon[:] = plon
    lat[:] = plat
    polys[:] = polyid[:]
    data[:] = var[:]
    ncout.close()


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
