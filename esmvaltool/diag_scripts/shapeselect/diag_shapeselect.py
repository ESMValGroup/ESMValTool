"""Diagnostic to select grid points within a shapefile."""
from copy import deepcopy
import logging
import os
import sys
from netCDF4 import Dataset, num2date
import fiona
from shapely.geometry import shape, mapping, Point, Polygon, MultiPolygon
from shapely.geometry import MultiPoint
import shapefile
from shapely.ops import nearest_points
import numpy as np
import xlsxwriter
import csv

from matplotlib import pyplot as p 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import iris
from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))

def main(cfg):
    """Select grid points within shapefiles."""
    for filename, attributes in cfg['input_data'].items():
        logger.info("Processing variable %s from model %s",
                    attributes['standard_name'], attributes['model'])
        logger.debug("Loading %s", filename)
        cube = iris.load_cube(filename)
        ncts, nclon, nclat = shapeselect2(cfg, cube, filename)
        name = os.path.splitext(os.path.basename(filename))[0] + '_polygon'
        if cfg['write_xlsxcsv']:
            ncfile = Dataset(filename, 'r')
            otime = ncfile.variables['time']
            dtime = num2date(otime[:], otime.units, otime.calendar)
            wtime = []
            for tim in range(len(dtime)):
                wtime.append(str(dtime[tim]))
            pathx = os.path.join(cfg['work_dir'], name + '.xlsx')
            pathc = os.path.join(cfg['work_dir'], name + '.csv')
            workbook = xlsxwriter.Workbook(pathx)
            #worksheetMeta = workbook.add_worksheet()
            worksheet = workbook.add_worksheet()
            worksheet.write(0, 0, 'Date/lonlat')
            worksheet.write_row(0, 1, wtime)
            for row in range(ncts.shape[1]):
                # Better change to ID here
                #worksheet.write(row+1, 0, str(row))
                worksheet.write(row+1, 0, str(round(nclon[row], 3)) +
                                '_' + str(round(nclat[row], 3)))
                worksheet.write_row(row+1, 1, np.squeeze(ncts[:, row]))
            #### Convert to cls before closing
            #csvfile = open(pathc, 'w')
            #wr = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
            #for rownum in range(row):
            #    wr.writerow(worksheet.row_values(rownum))
            #csvfile.close()
            #workbook.close()

        # fixme: write nrows of header in first row
        # fixme: dump cfg
        # fixme: dump cube.metadata
        # fixme: write info about representative coordinates
        # fixme: write time coordinate
        # fixme: write lon and lat info in first two rows
                #timenp = np.zeros((ncts.shape[0]+2))
                #timenp[
            
        #
        #        file.write("id time \n")
        #        file.close()
        #        tvar = deepcopy(var)
        #        tvar = np.transpose(tvar)
        #        np.savetxt(path, tvar, delimiter=',')

        #if cfg['write_netcdf']:
        #    path = os.path.join(
        #        cfg['work_dir'],
        #        name + '.nc',
        #    )
        # fixme: add similar stuff as in csv to netcdf
        #    write_netcdf(path, polyid, var, cube, cfg)

def dates(otime):
    """ Converts iris time to human readable format """
    calendar = otime.units.calendar
    origin = otime.units.origin
    #otime.cell(0).point
    


def shapeselect2(cfg, cube, filename):
    """ New solution with fiona! """
    shppath = cfg['shppath']
    wgtmet = cfg['wgtmet']
    if ((cube.coord('latitude').ndim == 1 and
         cube.coord('longitude').ndim == 1)):
        coordpoints = [(x, y) for x in cube.coord('longitude').points
                  for y in cube.coord('latitude').points]
    else:
        logger.info("Support for 2-d coords not yet implemented!")
        sys.exit(1)
    points = MultiPoint(coordpoints)
    if cfg['evalplot']:
        # Set limits for map (This can definitely be improved!)
        shap = shapefile.Reader(shppath)
        llcrnrlon=shap.bbox[0]-1
        llcrnrlat=max((shap.bbox[1]-1,-90))
        urcrnrlon=shap.bbox[2]+1
        urcrnrlat=min((shap.bbox[3]+1,90))
        alons = []
        alats = []
        for lon, lat in coordpoints:
            if ( lat > llcrnrlat and lat < urcrnrlat
                 and lon > llcrnrlon and lon < urcrnrlon ):
                alons.append(lon)
                alats.append(lat)
        cols = ['go', 'bo', 'co']
        cnt=0
        for shapp in shap.shapeRecords():
            xm = [i[0] for i in shapp.shape.points[:]]
            ym = [i[1] for i in shapp.shape.points[:]]
            p.plot(xm, ym)
            p.plot(alons, alats, 'ro', markersize = 2)
            cnt += 1
        p.xlabel('Longitude')
        p.ylabel('Latitude')
    with fiona.open(shppath) as shp:
        gpx = []
        gpy = []
        cnt = -1
        ncts = np.zeros((cube.coord('time').shape[0], len(shp)))
        nclon = np.zeros((len(shp))) # Takes representative point
        nclat = np.zeros((len(shp)))
        for ishp, multipol in enumerate(shp):
            cnt += 1
            multi = shape(multipol['geometry'])
            if wgtmet == 'mean_inside':
                pth = 'o'
                gpx, gpy = mean_inside(gpx, gpy, points, multi, cube, cfg)
                if cfg['evalplot']:
                    p.plot(cube.coord('longitude').points[gpx],
                           cube.coord('latitude').points[gpy],
                           'bo',markersize=6)
                if len(gpx) == 0:
                    gpx, gpy = representative(gpx, gpy, points, multi, cube,
                                              cfg)
            elif wgtmet == 'representative':
                pth = 'x'
                gpx, gpy = representative(gpx, gpy, points, multi, cube, cfg)
                if cfg['evalplot']:
                    p.plot(cube.coord('longitude').points[gpx],
                           cube.coord('latitude').points[gpy],
                           'ro', markersize=3)
            print(len(gpx))
            if len(gpx) == 1:
                ncts[:,ishp] = np.reshape(cube.data[:, gpy, gpx],
                                          (cube.data.shape[0],))
                #for tt in range(cube.data.shape[0]):
                #    ncts[:,ishp] = tmpdat[tt] #cube.data[:, gpy, gpx][]
            else:
                print('lengths: ',len(gpx),len(gpy))
                #print(cube.data[:, gpy, gpx])
                #tmpdat = np.mean(cube.data[:, gpy, gpx],axis=1)
                #for tt in range(cube.data.shape[0]):
                ncts[:,ishp] = np.mean(cube.data[:, gpy, gpx],axis=1)
                    #tmpdat[tt] #np.mean(cube.data[:, gpy, gpx],axis=(1,2))
            print(ncts[:,ishp])
            print('***********************************')
            gx, gy = representative([], [], points, multi, cube, cfg)
            nclon[ishp] = cube.coord('longitude').points[gx]
            nclat[ishp] = cube.coord('latitude').points[gy]

    if cfg['evalplot']:
        name = os.path.splitext(os.path.basename(filename))[0]
        path = os.path.join(cfg['work_dir'],name + '.png')
        p.savefig(path)
        p.show()
    #print(ncts)
    return ncts, nclon, nclat

def mean_inside(gpx, gpy, points, multi, cube, cfg):
    for point in points:
        if point.within(multi):
            if point.x < 0 and np.min(cube.coord('longitude').points) >= 0:
                addx = 360.
            else:
                addx = 0.
            x, y = best_match(cube.coord('longitude').points,
                              cube.coord('latitude').points,
                              point.x + addx, point.y)
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
    x, y = best_match(cube.coord('longitude').points,
                      cube.coord('latitude').points,
                      npx, npy)
    gpx.append(x)
    gpy.append(y)
    return gpx, gpy

def best_match(i, j, px, py):
    """ Identifies the grid points in 2-d with minimum distance. """
    if i.shape != 2 or j.shape != 2:
        x = deepcopy(i)
        y = deepcopy(j)
        xx = np.zeros((len(x),len(y)))
        yy = np.zeros((len(x),len(y)))
        for i in range(0,len(y)):
            xx[:,i] = x
        for j in range(0,len(x)):
            yy[j,:] = y
    else:
        xx = deepcopy(i)
        yy = deepcopy(j)
    distance = ( (xx - px)**2 + (yy - py)**2 )**0.5
    ind = np.unravel_index(np.argmin(distance, axis=None), distance.shape)
    #print(np.min(distance), distance[ind[0],ind[1]])
    #print(xx[ind[0],ind[1]],yy[ind[0],ind[1]], px, py)
    return ind[0], ind[1]


def shapeselect(cfg, cube, filename):
    """
    Add some description here

    Two lines
    """
    shppath = cfg['shppath']
    shpidcol = cfg['shpidcol']
    wgtmet = cfg['wgtmet']
    if ((cube.coord('latitude').ndim == 1 and
         cube.coord('longitude').ndim == 1)):
        coord_points = [(x, y) for x in cube.coord('latitude').points
                        for y in cube.coord('longitude').points]
    # if lat and lon are matrices
    elif (cube.coord('latitude').ndim == 2 and
          cube.coord('longitude').ndim == 2):
        logger.info("Matrix coords not yet implemented with iris!")
        sys.exit(1)
    else:
        logger.info("Unexpected error: " +
                    "Inconsistency between grid lon and lat dimensions")
        sys.exit(1)
    points = MultiPoint(coord_points)
    # Import shapefile with catchments
    shp = shapefile.Reader(shppath)
    shapes = shp.shapes()
    fields = shp.fields[1:]
    records = shp.records()
    attr = [[records[i][j] for i in range(len(records))]
            for j in range(len(fields))]
    for shpid, xfld in enumerate(fields):
        if xfld[0] == shpidcol:
            index_poly_id = shpid
            break
    else:
        logger.info("%s not in shapefile!", shpid)
        sys.exit(1)
    poly_id = attr[index_poly_id]

    # This should be a loop over shapes instead
    # Then we are flexible to chose new method if one fails
    # if wgtmet == x:
    #    try:
    #        call method centroid_inside
    #    exception method fail:
    #        call method2 centroid
    #
    # --  Method: nearest_centroid ---
    if wgtmet == 'nearest_centroid':
        # get polygon centroids
        mulpol = MultiPolygon([shape(pol) for pol in shapes])
        cent = MultiPoint([pol.centroid for pol in mulpol])
        # find nearest point in netcdf grid for each catchment centroid
        selected_points = []
        for i in cent:
            #if i.coords[0] < -180 or i.coords[0] > 360:
            #    print('Shape file is likely not defined in lon-lat.')
            #    print('Conversion is not yet implemented. Aborting!')
            #    sys.exit(1)
            nearest = nearest_points(i, points)
            nearestgplon = np.where(cube.coord('longitude').points ==
                                    nearest[1].coords[0][1])
            nearestgplat = np.where(cube.coord('latitude').points ==
                                    nearest[1].coords[0][0])
            selected_points.append(list((nearestgplon[0], nearestgplat[0])))
            print(nearest[1].coords[0][1], nearest[1].coords[0][0],
                  nearest[0].coords[0][1], nearest[0].coords[0][0],
                  nearestgplon[0], nearestgplat[0],
                  cube.coord('longitude').points[nearestgplon[0]],
                  cube.coord('latitude').points[nearestgplat[0]]
                 )
            # Add a check if the point is inside polygon
            # when looping over shapes:
            # point = Point(lon, lat)
            # shapeX.contains(point)
            # Issue warning if outside
            # Implement forced inside (see https://stackoverflow.com/questions/33311616/find-coordinate-of-closest-point-on-polygon-shapely/33324058 )
        var = np.zeros((len(cube.coord('time').points), len(poly_id)))
        cnt = 0
        for point in selected_points:
            print(point[0],point[1])
            var[:, cnt] = np.squeeze(cube.data[:, point[0], point[1]])
            cnt += 1
    else:
        logger.info('ERROR: invalid weighting method %s', wgtmet)
        sys.exit(1)
    if cfg['evalplot']:
        shape_plot(selected_points, cube, filename, cfg)
    return poly_id, var




def write_netcdf(path, polyid, var, cube, cfg):
    """Write results to a netcdf file."""
    shppath = cfg['shppath']
    # wgtmet = cfg['wgtmet']
    ncout = Dataset(path, mode='w')
    ncout.createDimension('time', None)
    ncout.createDimension('polygon', len(polyid))
    times = ncout.createVariable('time', 'f8', ('time'), zlib=True)
    times.setncattr_string('standard_name', cube.coord('time').standard_name)
    times.setncattr_string('long_name', cube.coord('time').long_name)
    # if isinstance(cube.coord('time').units, str):
    #     times.setncattr_string('units',cube.coord('time').units)
    # else:
    #     tunit = cube.coord('time').units
    times.setncattr_string('calendar', cube.coord('time').units.calendar)
    times.setncattr_string('units', cube.coord('time').units.origin)
    polys = ncout.createVariable('polygon', 'f4', ('polygon'), zlib=True)
    polys.setncattr_string('standard_name', 'polygon')
    polys.setncattr_string('long_name', 'polygon')
    polys.setncattr_string('shapefile', shppath)
    data = ncout.createVariable(cube.var_name, 'f4', ('time', 'polygon'),
                                zlib=True)
    data.setncattr_string('standard_name', cube.standard_name)
    data.setncattr_string('long_name', cube.long_name)
    data.setncattr_string('units', cube.units.origin)
    for key, val in cube.metadata[-2].items():
        ncout.setncattr_string(key, val)
    times[:] = cube.coord('time').points
    polys[:] = polyid[:]
    data[:] = var[:]
    ncout.close()


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
