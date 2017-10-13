# -*- coding: utf-8 -*-

import numpy as np
import shapefile as shp
import os
import netCDF4
from geoval.core.netcdf import NetCDFHandler

from geoval.core.data import GeoData
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import collections

# TODO correct _set_cell_area
# TODO discuss: get_shape_statistics, get_regions

class GeoData(GeoData):

    def C_set_cell_area(self): #Overwritten due to error
        """
        set cell area size. If a cell area was already given (either by user or from file)
        nothing will happen. Otherwise it will be tried to calculate cell_area from
        coordinates using the CDO's
        The estimation of the cell area follows the following steps:
        1) try to estimate cellarea using cdo gridarea
        2) if this does not work:
            a) directory is write protected --> write to temporary directory
            b) unknown grid --> try to select another grid, using cdo selgrid
        3) it could be that the cdo's dont support the input filetype. In that case
           it is tried to read the lat/lon information using the netCDF4 library
           store these fields in a dummy nc3 file and then apply the cdo's
        4) if all this does not work, then cell_area is set to unity for all grid cells and a WARNING is raised
        """

        # TODO unittest implementation

        if hasattr(self, 'cell_area'):
            if self.cell_area is not None:
                return

        # in case that cell area shall explicitely NOT be caluclated
        if not self.cell_area: # TODO This line was wrong
            if self.ndim == 2:
                self.cell_area = np.ones(self.data.shape)
            elif self.ndim == 3:
                self.cell_area = np.ones(self.data[0, :, :].shape)
            else:
                raise ValueError('Invalid geometry!')
            return

        if (self.lat is None) or (self.lon is None):
            self._log_warning(
                "WARNING: cell area can not be calculated (missing coordinates)!")
            if self.ndim == 2:
                self.cell_area = np.ones(self.data.shape)
            elif self.ndim == 3:
                self.cell_area = np.ones(self.data[0, :, :].shape)
            else:
                raise ValueError('Invalid geometry!')
            return

        # calculate cell area from coordinates
        cell_file = self.filename[:-3] + '_cell_area.nc'

        if not os.path.exists(cell_file):  # calculate grid area using CDO's
            cdo = Cdo()
            try:
                cdo.gridarea(options='-f nc', output=cell_file,
                             input=self.filename)
            except:
                # occurs if you dont have write permissions
                print '   Seems that cell_area file can not be generated, try to generate in temporary directory'
                # generate some temporary filename
                cell_file = tempfile.mktemp(prefix='cell_area_', suffix='.nc')
                try:
                    cdo.gridarea(options='-f nc', output=cell_file,
                                 input=self.filename)
                    print '   Cell area file generated sucessfully in temporary file: ' + cell_file
                except:
                    # not sucessfull so far ... last try here by selecting an
                    # alternative grid (if available)
                    print(
                        "   Try to calculate gridarea using alternative grid")
                    # generate some temporary filename
                    cell_file = tempfile.mktemp(
                        prefix='cell_area_', suffix='.nc')
                    try:
                        cdo.gridarea(options='-f nc', output=cell_file,
                                     input='-selgrid,2 ' + self.filename)
                        print '   Cell area file generated sucessfully in temporary file: ' + cell_file
                    except:
                        try:
                            # store lat/lon coordinates in nc3 file and then
                            # apply the cdo's (last try)
                            coord_file = tempfile.mktemp(
                                prefix='coordinates_', suffix='.nc')
                            tmpdat = GeoData(None, None)
                            tmpdat.lat = self.lat
                            tmpdat.lon = self.lon
                            tmpxxx = np.zeros(self.lat.shape)
                            tmpdat.data = np.ma.array(
                                tmpxxx, mask=tmpxxx != tmpxxx)
                            tmpdat.save(coord_file, format='nc3')
                            del tmpxxx, tmpdat
                            cell_file = tempfile.mktemp(
                                prefix='cell_area_', suffix='.nc')
                            cdo.gridarea(options='-f nc',
                                         output=cell_file, input=coord_file)
                        except:
                            print('WARNING: no cell area could be generated!')

        # read cell_area file ---
        if os.path.exists(cell_file):
            # read cell area from file
            File = NetCDFHandler()
            File.open_file(cell_file, 'r')
            self.cell_area = File.get_variable('cell_area')
            File.close()

            # check geometries
            if self.data.ndim == 2:
                if self.cell_area.shape != self.data.shape:
                    raise ValueError(
                        'Invalid cell_area file: delete it manually and check again!')
            elif self.data.ndim == 1:
                if self.cell_area.shape != self.data.shape:
                    raise ValueError(
                        'Invalid cell_area file: delete it manually and check again!')
            elif self.data.ndim == 3:
                if self.cell_area.shape != self.data[0, :, :].shape:
                    raise ValueError(
                        'Invalid cell_area file: delete it manually and check again!')
        else:
            # no cell area calculation possible!!!
            # logger.warning('Can not estimate cell area! (setting all equal) ' + cell_file)

            self._log_warning(
                '*** WARNING: Can not estimate cell area! ' + cell_file)
            self._log_warning(' setting cell_area all to equal')

            if self.data.ndim == 2:
                self.cell_area = np.ones(self.data.shape)
            elif self.data.ndim == 3:
                self.cell_area = np.ones(self.data[0, :, :].shape)
            else:
                print 'actual geometry:  ', self.data.ndim, self.data.shape
                raise ValueError('Invalid geometry!')


    def get_shape_statistics(self,regions): #written before geoval was implemented
        """
        get statistical information for different polygons in shapefile
        Parameters
        ----------
        regions : masks for masked array
        """

        self.regionalized=dict()
        regname=regions.keys()

#        for s in np.arange(len(regions)):
#            loc_content=self.data.copy()
#            loc_content.mask=regions[regname[s]]
#
#            try:
#                a=np.nanmin(loc_content)
#            except:
#                a=np.nan
#            try:
#                b=np.nanmean(loc_content)
#            except:
#                b=np.nan
#            try:
#                c=np.nanmax(loc_content)
#            except:
#                c=np.nan
#            try:
#                d=np.nanstd(loc_content)
#            except:
#                d=np.nan
#            try:
#                e=loc_content.count()
#            except:
#                e=np.nan

        for s in np.arange(len(regions)):
            loc_content=self.copy()
            #loc_content.data.mask=regions[regname[s]]
            loc_content.data.mask=np.logical_or(loc_content.data.mask,regions[regname[s]])

# A_laue_ax+
            #print("==========================================================")
            #print(regname[s])
            #print(np.shape(loc_content.data.mask))
            #print(np.shape(regions[regname[s]]))
            #ncfile = NetCDFHandler()
            #ncfile.open_file('debug.nc', 'w')
            #dims = np.shape(loc_content.data.mask)
            #nx = dims[1]
            #ny = dims[0]
            #ncfile.create_dimension('x',size=nx)
            #ncfile.create_dimension('y',size=ny)
            #data = ncfile.create_variable('data','i4',('y','x'))
            #ncfile.assign_value('data', self.data.mask.astype(int))
            #ncfile.close()
            #print("==========================================================")
            #if (s == 1):
            #    exit()
# A_laue_ax-

            try:
                a=np.nanmin(loc_content.data)
            except:
                a=np.nan
            try:
                b=loc_content.fldmean(return_data=False)[0]
            except:
                b=np.nan
            try:
                c=np.nanmax(loc_content.data)
            except:
                c=np.nan
            try:
                lc2=loc_content.copy()
                lc2.data=lc2.data**2
                d=np.sqrt((lc2.fldmean(return_data=False)[0]-b**2))
            except:
                d=np.nan
            try:
                e=loc_content.data.count()
            except:
                e=np.nan

            self.regionalized[regname[s]]=[a,b,c,d,e]

    def get_regions(self,shape,column=0): #written before geoval was implemented
        """
        get setup for statistical information for different polygons in shapefile
        Parameters
        ----------
        shape : shp.Reader (shapefile.Reader)
            information on areas from a classic ESRI shapefile
        """
        assert isinstance(shape,shp.Reader)



        def point_in_poly(point,poly):

            testy=point[1]
            testx=point[0] if point[0]<180 else point[0]-360

            vertx=[item[0] for item in poly]
            verty=[item[1] for item in poly]

            nvert=len(poly)
            c=False

            j=nvert-1
            for i in range(0, nvert):
                if ( not((verty[i]>testy) == (verty[j]>testy)) and (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) ):
                    c = not c
                j=i

            return c


        regions=dict()
        regname=np.array(shape.records())[:,column]

        for s in range(len(shape.shapes())):
            breaks=shape.shapes()[s].parts
            breaks.append(0)
            loc_polys=[]
            loc_masks=[]
            for e in range(len(breaks)-1):

                loc_polys.append(shape.shapes()[s].points[breaks[e]:(breaks[e+1]-1)])

                if len(self.shape) == 3:
                    loc_mask=np.any(self.data,axis=0).mask.copy()
                elif len(self.shape) == 2:
                    loc_mask=self.data.mask.copy()
                else :
                    assert False, "wrong data dimensions"
                for i in np.arange(self.shape[0] if len(self.shape) == 2 else self.shape[1]):
                    for j in np.arange(self.shape[1] if len(self.shape) == 2 else self.shape[2]):
                        ll=[(self.lon[i,j]),self.lat[i,j]]

                        loc_mask[i,j]=point_in_poly(ll,loc_polys[-1])

                loc_masks.append(loc_mask)

#            if len(self.shape) == 3:
#                loc_mask=self.data.mask.copy()
#            elif len(self.shape) == 2:
#                loc_mask=self.data.mask.copy()
#            else :
#                assert False, "wrong data dimensions"

#            regions[regname[s]]=np.logical_or(loc_mask,np.logical_not(sum(loc_masks) % 2))
            regions[regname[s]]=np.logical_not(sum(loc_masks) % 2)

        return collections.OrderedDict(sorted(regions.items()))



#    def _save_netcdf(self, filename, varname=None, delete=False, compress=True, format='NETCDF4'):
#        """
#        saves the data object to a netCDF file
#
#        Parameters
#        ----------
#        filename : str
#            filename of output file
#        varname : str
#            name of output variable; this explicitely overwrites
#            self.varname, which is tried to be used as a first order
#        delete : bool
#            delete file if existing without asking
#        compress : bool
#            compress resulting file if supported by backend
#        format : str
#            output file format specifier as used by netCDF4 library
#        """
#
#        # check if output file already there
#        if os.path.exists(filename):
#            if delete:
#                os.remove(filename)
#            else:
#                raise ValueError('File already existing. Please delete \
#                                  manually or use DELETE option: \
#                                  %s' % filename)
#        # variable name
#        if varname is None:
#            if self.varname is None:
#                varname = 'var1'
#            else:
#                varname = self.varname
#
#        # create new file
#        File = NetCDFHandler()
#        File.open_file(filename, 'w', format=format)
#
#        # create dimensions
#        if hasattr(self, 'data'):
#            if self.data.ndim == 3:
#                if self.time is None:
#                    raise ValueError('No time variable existing! \
#                                      Can not write 3D data!')
#                nt, ny, nx = self.shape
#            elif self.data.ndim == 2:
#                ny, nx = self.shape
#            else:
#                raise ValueError('Current shape not supported here %s' %
#                                 self.shape)
#        else:
#            ny, nx = np.shape(self.lat)
#
#        # Create Dimensions
#        File.create_dimension('ny', size=ny)
#        File.create_dimension('nx', size=nx)
#
#        # Create variables
#        if hasattr(self, 'time'):
#            if self.time is not None:
#                File.create_dimension('time', size=len(self.time))
#                File.create_variable('time', 'd', ('time',))
#                File.F.variables['time'].units = self.time_str
#
#        if hasattr(self, 'data'):
#            if self.data.ndim == 3:
#                File.create_variable(
#                    varname, 'd', ('time', 'ny', 'nx'), zlib=compress)
#            elif self.data.ndim == 2:
#                File.create_variable(varname, 'd', ('ny', 'nx'), zlib=compress)
#
#        if self.ny is not None:
#            File.create_variable('ny', 'd', ('ny',))
#            File.F.variables['ny'].units = 'degrees_north'
#            File.F.variables['ny'].long_name = "y-coordinate"
#
#        if self.nx is not None:
#            File.create_variable('nx', 'd', ('nx',))
#            File.F.variables['nx'].units = 'degrees_east'
#            File.F.variables['nx'].long_name = "x_coordinate"
#
#        if self.lat is not None:
#            File.create_variable('lat', 'd', ('ny', 'nx'))
#            File.F.variables['lat'].units = 'degrees_north'
#            File.F.variables['lat'].axis = "Y"
#            File.F.variables['lat'].long_name = "latitude"
#
#        if self.lon is not None:
#            File.create_variable('lon', 'd', ('ny', 'nx'))
#            File.F.variables['lon'].units = 'degrees_east'
#            File.F.variables['lon'].axis = "X"
#            File.F.variables['lon'].long_name = "longitude"
#
#        if hasattr(self, 'data'):
#            if hasattr(self, 'cell_area'):
#                File.create_variable('cell_area', 'd', ('ny', 'nx'))
#
#        #/// write data
#        if hasattr(self, 'time'):
#            if self.time is not None:
#                File.assign_value('time', self.time)
#                File.F.variables['time'].calendar = self.calendar
#
#        if hasattr(self, 'data'):
#            File.assign_value(varname, self.data)
#
#        if self.lat is not None:
#            File.assign_value('lat', self.lat)
#        if self.lon is not None:
#            File.assign_value('lon', self.lon)
#        if self.ny is not None:
#            File.assign_value('ny', np.unique(self.lat))
#        if self.nx is not None:
#            File.assign_value('nx', np.unique(self.lon))
#
#
#        if hasattr(self, 'data'):
#            if hasattr(self, 'cell_area'):
#                if self.cell_area is not None:
#                    File.assign_value('cell_area', self.cell_area)
#
#        if hasattr(self, 'data'):
#            if hasattr(self, 'long_name'):
#                if self.long_name is not None:
#                    File.F.variables[varname].long_name = self.long_name
#            if hasattr(self, 'unit'):
#                if self.unit is not None:
#                    File.F.variables[varname].units = self.unit
#
#            File.F.variables[varname].scale_factor = 1.
#            File.F.variables[varname].add_offset = 0.
#            File.F.variables[varname].coordinates = "longitude latitude"
#
#        File.close()

    def _save_netcdf(self, filename, varname=None, delete=False, compress=True, format='NETCDF4'):
        """
        saves the data object to a netCDF file

        Parameters
        ----------
        filename : str
            filename of output file
        varname : str
            name of output variable; this explicitely overwrites
            self.varname, which is tried to be used as a first order
        delete : bool
            delete file if existing without asking
        compress : bool
            compress resulting file if supported by backend
        format : str
            output file format specifier as used by netCDF4 library
        """

        # check if output file already there
        if os.path.exists(filename):
            if delete:
                os.remove(filename)
            else:
                raise ValueError('File already existing. Please delete \
                                  manually or use DELETE option: \
                                  %s' % filename)
        # variable name
        if varname is None:
            if self.varname is None:
                varname = 'var1'
            else:
                varname = self.varname

        # create new file
        File = NetCDFHandler()
        File.open_file(filename, 'w', format=format)

        # create dimensions
        if hasattr(self, 'data'):
            if self.data.ndim == 3:
                if self.time is None:
                    raise ValueError('No time variable existing! \
                                      Can not write 3D data!')
                nt, ny, nx = self.shape
            elif self.data.ndim == 2:
                ny, nx = self.shape
            else:
                raise ValueError('Current shape not supported here %s' %
                                 self.shape)
        else:
            ny, nx = np.shape(self.lat)

        # Create Dimensions
        File.create_dimension('lat', size=ny)
        File.create_dimension('lon', size=nx)

        # Create variables
        if hasattr(self, 'time'):
            if self.time is not None:
                File.create_dimension('time', size=len(self.time))
                File.create_variable('time', 'd', ('time',))
                File.F.variables['time'].units = self.time_str

        if hasattr(self, 'data'):
            if self.data.ndim == 3:
                File.create_variable(
                    varname, 'd', ('time', 'lat', 'lon'), zlib=compress)
            elif self.data.ndim == 2:
                File.create_variable(varname, 'd', ('lat', 'lon'), zlib=compress)

        if self.lat is not None:
            File.create_variable('lat', 'd', ('lat',))
            File.F.variables['lat'].units = 'degrees_north'
            File.F.variables['lat'].long_name = "y-coordinate"

        if self.lon is not None:
            File.create_variable('lon', 'd', ('lon',))
            File.F.variables['lon'].units = 'degrees_east'
            File.F.variables['lon'].long_name = "x_coordinate"

        if self.lat is not None:
            File.create_variable('latitude', 'd', ('lat', 'lon'))
            File.F.variables['latitude'].units = 'degrees_north'
            File.F.variables['latitude'].axis = "Y"
            File.F.variables['latitude'].long_name = "latitude"

        if self.lon is not None:
            File.create_variable('longitude', 'd', ('lat', 'lon'))
            File.F.variables['longitude'].units = 'degrees_east'
            File.F.variables['longitude'].axis = "X"
            File.F.variables['longitude'].long_name = "longitude"

        if hasattr(self, 'data'):
            if hasattr(self, 'cell_area'):
                File.create_variable('cell_area', 'd', ('lat', 'lon'))

        #/// write data
        if hasattr(self, 'time'):
            if self.time is not None:
                File.assign_value('time', self.time)
                File.F.variables['time'].calendar = self.calendar

        if hasattr(self, 'data'):
            File.assign_value(varname, self.data)

        if self.lat is not None:
            File.assign_value('latitude', self.lat)
        if self.lon is not None:
            File.assign_value('longitude', self.lon)
        if self.lat is not None:
            File.assign_value('lat', np.unique(self.lat))
        if self.lon is not None:
            File.assign_value('lon', np.unique(self.lon))


        if hasattr(self, 'data'):
            if hasattr(self, 'cell_area'):
                if self.cell_area is not None:
                    File.assign_value('cell_area', self.cell_area)

        if hasattr(self, 'data'):
            if hasattr(self, 'long_name'):
                if self.long_name is not None:
                    File.F.variables[varname].long_name = self.long_name
            if hasattr(self, 'unit'):
                if self.unit is not None:
                    File.F.variables[varname].units = self.unit

            File.F.variables[varname].scale_factor = 1.
            File.F.variables[varname].add_offset = 0.
            File.F.variables[varname].coordinates = "longitude latitude"

        File.close()

    GeoData._set_cell_area=C_set_cell_area
    GeoData.get_regions=get_regions
    GeoData.get_shape_statistics=get_shape_statistics
    GeoData._save_netcdf=_save_netcdf






