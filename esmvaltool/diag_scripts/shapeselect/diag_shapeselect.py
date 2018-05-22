"""Diagnostic to select grid points within a shapefile."""
import logging
import os
import iris
from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared.plot import quickplot

from netCDF4 import Dataset
import shapefile
from shapely.geometry import MultiPoint
from shapely.geometry import shape
from shapely.geometry.multipolygon import MultiPolygon
from shapely.ops import nearest_points

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Select grid points within shapefiles."""
    for filename, attributes in cfg['input_data'].items():
        logger.info("Processing variable %s from model %s",
                    attributes['standard_name'], attributes['model'],variables)
        #shapeselect(filename, variables, )
        print("hej!")

def shapeselect(ecvfile, varname, shppath, shpidcol, ncout,
               wgtmet='nearest_centroid'):
    """
    Add some description here

    Two lines
    """
    ecv = Dataset(ecvfile, mode='r')
    var_in = ecv.variables[varname]
    lat = ecv.variables['lat']
    lon = ecv.variables['lon']
    time = ecv.variables['time']

    # --- create Multipoint object representing the netcdf grid ---
    # if lat and lon are 1-dimensional: take all parwise combinations
    if len(lat.shape) == 1 and len(lon.shape) == 1:
        coord_points = [(x, y) for x in lat[:] for y in lon[:]]
    # if lat and lon are matrices
    if len(lat.shape) == 2:
        # check if lat.shape == lon.shape, issue error otherwise
        coord_points = [(lat[x, y], lon[x, y])
                        for x in range(lat.shape[0])
                        for y in range(lon.shape[1])]

    points = MultiPoint(coord_points)

    # Import shapefile with catchments
    shp = shapefile.Reader(shppath)
    shapes = shp.shapes()

    # extract field names from shapefile
    fields = shp.fields[1:]
    # extract attribute table from shapefile
    records = shp.records()
    attr = [[records[i][j] for i in range(len(records))]
            for j in range(len(fields))]

    # extract catchment id
    for shpid, xfld in enumerate(fields):
        if xfld[0] == shpidcol:
            index_catch_id = shpid
            break
    # ADD CHECK AND ERROR if shpidcol not in fields
    print(attr[index_catch_id])
    catch_id = attr[index_catch_id]

    # --  Method: nearest_centroid ---
    if wgtmet == 'nearest_centroid':
        # get catchments centroids
        mulpol = MultiPolygon([shape(pol) for pol in shapes])
        cent = MultiPoint([pol.centroid for pol in mulpol])

        # find nearest point in netcdf grid for each catchment centroid
        selected_points = []
        for i in cent:
            nearest = nearest_points(i, points)
            selected_points.append(list(nearest[1].coords)[0])

        # create table of results (time x one column for each catchment)
        var = []
        for point in selected_points:
            vari = var_in[:, point[0], point[1]]
            var.append(vari)

    else:
        print('ERROR: invalid weighting method')
        exit()

    # --- Write output netcdf ---
    # tries copy/creation netcdf
    outnetcdf = Dataset(ncout, "w", format="NETCDF4")

    # copy general attributes from original file and add prefix "orgfile"
    # outnetcdf.setncatts(ecv.__dict__)
    outnetcdf.setncatts({'orgfile_' + k: ecv.getncattr(k)
                         for k in ecv.ncattrs()})
    # add a "history" attribute with information on how the data was processed
    outnetcdf.setncatts({'history': 'data from the original file (orgfile:' +
                                    ecvfile + ') was processed by the catchment_selection script to produce summarized data for the cacthments provided in the shapefile ' +
                                    shppath + ' using the weigting method ' + wgtmet})

    # create dimensions
    time_out = outnetcdf.createDimension('time', len(time))
    bnds_out = outnetcdf.createDimension('bnds', ecv.dimensions['bnds'].size)
    catchment_id = outnetcdf.createDimension("catchment_id", len(mulpol))

    # create variables and set variable attributes
    for v_name in ('time', 'bnds', 'time_bnds'):
        varin = ecv.variables[v_name]
        outvar = outnetcdf.createVariable(v_name, varin.datatype,
                                          varin.dimensions)
        outvar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
        outvar[:] = varin[:]

    # or change back to 'i4' for integer naming
    catchment_id = outnetcdf.createVariable('catchment_id', 'i4',
                                            'catchment_id')
    var_out = outnetcdf.createVariable(varname, 'f4', ('time','catchment_id'))
    var_out.setncatts({'long_name':var_in.getncattr('long_name'),
                       'units':var_in.getncattr('units'),
                       'missing_value':var_in.getncattr('missing_value'),
                       '_FillValue':var_in.getncattr('_FillValue'),
                       'standard_name':var_in.getncattr('standard_name'),
                       'weighting_method':wgtmet})

    # fill with data
    catchment_id[:] = catch_id
    var_out[:] = var

    # Plot
    # from mpl_toolkits.basemap import Basemap
    # import matplotlib.pyplot as plt
    # map = Basemap(projection='mill',lon_0=180)
    # map.drawmapboundary(fill_color='aqua')
    # map.fillcontinents(color='#ddaa66',lake_color='aqua')
    # map.drawcoastlines()
    # map.readshapefile('./data/SUBID_WWH_1_0_1_WHISTv18_simpl90m_ForVisOnly',
    # 'SUBID')
    # plt.savefig('cdateRef.png', bbox_inches='tight')
    # plt.show()

if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)

