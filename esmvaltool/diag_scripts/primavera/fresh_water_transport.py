# import glob
# import matplotlib
# from matplotlib import cm, colors, pyplot
# from matplotlib.animation import FuncAnimation
# # from mpl_toolkits.basemap import Basemap, maskoceans
# # from mpl_toolkits.basemap import shiftgrid
# import time
# import netCDF4
# import pandas as pd
# import numpy
# import numpy.ma as ma
# from numpy import logical_or, logical_and
# import matplotlib.pyplot as pyplot
# from scipy.interpolate import interp1d, interp2d
# import scipy.interpolate
# import sys
# import numbers
# import collections
# from matplotlib import colors as mcolors
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from geopy.distance import vincenty
# from geopy.distance import great_circle
# import xarray as xr
# import iris
# from iris.experimental.equalise_cubes import equalise_attributes
# import iris.pandas
# import cf_units
# import iris.coord_categorisation
# from iris.coords import DimCoord
# import iris.plot as iplot
# import iris.analysis.stats
# import iris.analysis.cartography
import os
import logging

import iris

from esmvalcore.preprocessor import extract_region

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata

logger = logging.getLogger(os.path.basename(__file__))

class FreshWaterTransport(object):
    def __init__(self, config):
        self.cfg = config
        self.filenames = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.output_name = None
        self.target_grid = self.cfg.get('target_grid')
        self.grid_cube = None
        self.sftlf = None

    def compute(self):
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        straits =  self.cfg.get(
            'straits',
            {
                'barents': {'latitude': (-174, -165), 'longitude': (63.9, 67.3)}
            }
        )
        for alias in data:
            var = group_metadata(data[alias], 'short_name')
            logger.info('Processing %s', alias)
            so = iris.load_cube(var['so'][0]['filename'])
            vo = iris.load_cube(var['vo'][0]['filename'])
            for strait, config in straits.items():
                self.plot_salinity_profile(so, vo, strait, config)

    def plot_salinity_profile(self, so, vo, strait, config, num=0, yearly=False):
        so = extract_region(so, config['latitude'][0], config['latitude'][1],
                            config['longitude'][0], config['longitude'][1])
        vo = extract_region(so, config['latitude'][0], config['latitude'][1],
                            config['longitude'][0], config['longitude'][1])
                            
        pass

        # sss=sal[0,:,:,:].data.squeeze()
        # uss=u[0,:,:,:].data.squeeze()
        # if (ind_latit==0):
        #     print 'lona: ',lona
        #     print 'lonmin:', lonmin
        #     print 'lonmax:', lonmax
        #     ind_lat, ind_lon=numpy.where((lona>=lonmin) & (lona<=lonmax) & (lata >= latmin) & (lata <= latmax))
        #     print 'ind_lat', ind_lat
        #     counter=collections.Counter(ind_lat)
        #     most_com=counter.most_common(1)
        #     ss=most_com[0]
        #     ind_latit=ss[0]

        # ind_lon=numpy.where((lona[ind_latit,:]>=lonmin) & (lona[ind_latit,:]<=lonmax))
        # #sss[0,ind_latit,:]=100
        # #uss[0,ind_latit,:]=100
        # #map_index(lona,lata,uss[0,:,:].squeeze())
        # il=ind_lon[0]
        # y=lata[ind_latit,ind_lon]
        # x=lona[ind_latit,ind_lon]
        # print 'lat', y
        # print 'lon', x
        # distxy=numpy.zeros(len(il))
        # d1=0
        # for dist in numpy.arange(0,len(il)):
        #     #print 'ind_latit', ind_latit
        #     p1=(lata[ind_latit,il[d1]], lona[ind_latit,il[d1]])
        #     p2=(lata[ind_latit,il[d1]+1], lona[ind_latit,il[d1]+1])
        #     #print 'points to follow', p1
        #     distxy[d1]=great_circle(p1, p2).m
        #     d1=d1+1
        # print distxy
        # print 'shape lat', lata.shape
        # sal=sal[:,:,ind_latit,ind_lon[0]].data.squeeze()
        # u=u[:,:,ind_latit,ind_lon[0]].data.squeeze()
        # print 'u',u[0,:]
        # u2=u
        # sal[sal<10]=0
        # sal[sal>45]=0
        # #pyplot.figure()
        # #pyplot.pcolormesh(sal[0],vmin=33, vmax=35)
        # #pyplot.title('salinity')
        # #pyplot.colorbar()
        # #pyplot.show()
        # #pyplot.figure()
        # #pyplot.pcolormesh(u[0],vmin=-0.5, vmax=0.5)
        # #pyplot.title('u')
        # #pyplot.colorbar()
        # #pyplot.show()

        # sal=sal*u
        # sizes=sal.shape
        # angle=numpy.arctan2(numpy.diff(y),numpy.diff(x))* 180 / numpy.pi
        # z2=numpy.arange(0,len(z)+2)
        # z2[0]=0
        # z2[1:-1]=z
        # #diffz=numpy.gradient(z2[0:-2])
        # diffz=numpy.gradient(z)
        # area=numpy.outer(distxy, diffz)
        # dz=numpy.ones(area.shape)
        # dx=numpy.ones(area.shape)
        # dx=distxy[0]*dx
        # dz=dz*diffz
        # area2=numpy.transpose(area)
        # #for ti in numpy.arange(0,sizes[0]):
        # #  for longs in numpy.arange(0,len(il)):
        # #print 'z levels', numpy.gradient(z2[0:-2])
        # #print 'len(z)', len(z)
        # stot=numpy.empty(sizes)
        # trans_t=numpy.empty(sizes[0])
        # utot=numpy.empty(sizes)
        # utrans_t=numpy.empty(sizes[0])
        # #print 'U2 SHAPEEEEE', u2.shape
        # for t in numpy.arange(0,sizes[0]):
        #     st=sal[t,...].squeeze()
        #     #st[st==-1000]=0
        #     #area2[st==-1000]=0
        #     mult=numpy.multiply(st,area2)
        #     stot[t,:,:]=mult
        #     trans=mult[mult!=0]
        #     trans_t[t]=numpy.sum(trans)

        #     ut=u2[t,...].squeeze()
        #     #print 'ut',ut
        #     #ut[ut==0]=numpy.nan
        #     umult=numpy.multiply(ut,area2)
        #     utot[t,:,:]=umult
        #     transu=umult[umult!=0]
        #     #print 'trans',trans
        #     utrans_t[t]=numpy.sum(transu)
        # ct=0
        # cm=0
        # cy=0
        # #print 'utrans_t',utrans_t
        # #yutran=numpy.empty([(sizes[0]/12)+1])
        # #ysaltran=numpy.empty([(sizes[0]/12)+1])
        # fwt=(utrans_t-(1/34.8)*trans_t)/1000
        # #print 'shape multiplication', stot.shape
        # hs = pd.Series(fwt,index=time)
        # #phony2period = lambda d: pd.Period(year=d.year, month=d.month, freq="M")
        # #hs.index = map(phony2period, time)

        # #print time
        # mult=ma.masked_invalid(mult)
        # #fig=pyplot.figure()
        # if yearly=="A":
        #     fig=pyplot.figure()
        #     hn=hs.resample("A")
        #     hs_freq =hs.resample("A").plot(title=name)
        #     hs_freq.set_xlabel("Time, years")
        #     hs_freq.set_ylabel("FWT transport")
        #     tit='freshwater_transport_'+name+'.pdf'
        #     fig.savefig(tit)

        # if yearly=="March":
        #     #hn=hs.groupby([hs.index.month==3]title=name).plot()
        #     hn=hs[hs.index.month==3].plot(title=name,xlabel='Time,years',ylabel='FWT')
        # if yearly=="November":
        #     #hn=hs.groupby([hs.index.month==3]).plot()
        #     hn=hs[hs.index.month==11].plot()
        # if yearly=="Monthly":
        #     #fig=pyplot.figure()
        #     hn=hs
        #     #hs_freq=hn.plot(title=name)
        #     #hn_freq = hn.resample("A").plot()
        # #hs.plot(style='g--')

        # #pyplot.title('u')
        # #pyplot.colorbar()

        # print 'mult', umult[0,:]
        # #pyplot.figure()
        # #pyplot.pcolormesh(area2)
        # #pyplot.colorbar()
        # #pyplot.title('dlon*dz')

        # #pyplot.figure()
        # #pyplot.plot(fwt)
        # #pyplot.title(name)
        # #pyplot.figure()
        # ##m = Basemap(llcrnrlon=-180.,llcrnrlat=-90.,urcrnrlon=180.,urcrnrlat=90, projection='cyl', resolution='l')
        # #m = Basemap(projection='nplaea',boundinglat=20,lon_0=0,resolution='l')
        # #x, y = m(lona,lata)
        # #fig_handle = m.pcolormesh(x,y,uss[0,:,:].squeeze()*sss[0,:,:].squeeze(),cmap='rainbow',rasterized=True,vmin=levels[0], vmax=levels[-1])
        # #m.drawcoastlines(linewidth=0.25)
        # #m.fillcontinents(color='darkkhaki')
        # #parallels = numpy.arange(0.,90.,10.)
        # #m.drawparallels(parallels,labels=[False,True,True,False])
        # #meridians = numpy.arange(-180.,180.,40.)
        # #m.drawmeridians(meridians,labels=[True,False,False,True])
        # #pyplot.colorbar()
        # newlon=lona[ind_latit,ind_lon].squeeze()

        # #print 'newlon', newlon
        # sal_col=sss[:,ind_latit,ind_lon].squeeze()
        # #pyplot.pcolormesh(newlon,-1*z,sal_col,cmap='YlOrRd',rasterized=True,vmin=31, vmax=36)
        # #pyplot.pcolormesh(newlon,-1*z,sal_col,cmap='bwr',rasterized=True,vmin=31.5, vmax=36.5)
        # #print 'shape z', z.shape
        # #negz=-1*z
        # negz=-1*z[z<prof]
        # #pyplot.figure()
        # #cMap = mcolors.ListedColormap(['saddlebrown','saddlebrown','saddlebrown'])
        # #pyplot.pcolormesh(newlon,-1*z[z<prof],sal_col[z<prof,:],cmap='rainbow',rasterized=True,vmin=levels[0], vmax=levels[-1])
        # #pyplot.colorbar()
        # ##pyplot.pcolormesh(newlon,-1*z[z<prof],land[z<prof],cmap=cMap,rasterized=True,vmin=-1000, vmax=-10)
        # ##CS = pyplot.contourf(newlon, -1*z[z<prof], sal_col[z<prof,:],levels = levels,cmap='rainbow',origin='lower',extend='neither',vmin=levels[0], vmax=levels[-1])
        # #pyplot.title(name)
        # #pyplot.axis([newlon.min(), newlon.max(), negz.min(), negz.max()])
        # return ind_latit,hn.values


def map_index(lona,lata,var):
    pyplot.figure()
    m = Basemap(projection='nplaea',boundinglat=20,lon_0=0,resolution='l')
    x, y = m(lona,lata)
    fig_handle = m.pcolormesh(x,y,var,cmap='rainbow',rasterized=True)
    m.drawcoastlines(linewidth=0.25)
    m.fillcontinents(color='darkkhaki')
    parallels = numpy.arange(0.,90.,10.)
    m.drawparallels(parallels,labels=[False,True,True,False])
    meridians = numpy.arange(-180.,180.,40.)
    m.drawmeridians(meridians,labels=[True,False,False,True])
    pyplot.colorbar()


def save_netcdf(transport,model, strait):

    namefile='mon.freshwater_transport_'+strait+'_'+model+'_1950-2014.nc'
    dataset= Dataset('trans_netcdf/'+namefile, 'w',format='NETCDF4_CLASSIC')
    time =dataset.createDimension('time',None)

    #latitudes = dataset.createVariable('latitude', numpy.float32,('lat',))

    times =dataset.createVariable('time', numpy.float64, ('time',))
    times.units= 'hours since 0001-01-01 00:00:00'
    times.calendar ='gregorian'
    dates=[]
    if 'ORAP' in model:
        dates=[]
        if 'A' in freq:

            for an in range(len(transport)):
                dates.append(datetime(1979, 6, 1) + an*timedelta(days=365))
        if 'Monthly' in freq:
            for an in range(len(transport)):
                dates.append(datetime(1979, 1, 15) + an*timedelta(days=30))
        times[:] = date2num(dates, units = times.units,calendar = times.calendar)
    else:
        for an in range(len(transport)):
                dates.append(datetime(1950, 1, 15) + an*timedelta(days=30))
        times[:] = date2num(dates, units = times.units,calendar = times.calendar)
        #times[:]=date2num(time1,units=times.units,calendar=times.calendar)
    fwt = dataset.createVariable('fwt', numpy.float32,('time'))
    fwt.units='mSv'
    #dataset.history = 'Created ' + time.ctime(time.time())
    fwt[0:len(transport)]=transport
    dataset.close()


def main():
    fwt_month_bering = []
    fwt_month_fram = []
    fwt_month_gin = []
    fwt_month_lab = []
    fwt_month_bar = []
    fwt_month_green = []
    fwt_month_ice = []
    fwt_month_baff = []
    for fn in filelist_so:

        ind_latit2,fwt_bering=plot_salinity_profile(lon,lat,z,cube_so,cube_vo,lonadj-174,lonadj-165,63.9,67.3,'BERINGSTRAIT',list(0.01*numpy.arange(-30.,35.,0.1)),4,0,80,"Monthly",time)
        #ind_latit1,fwt_fram=plot_salinity_profile(lon,lat,z,cube_so,cube_vo,lonadj-20,lonadj+14,77.,78.,'FRAMSTRAIT',list(0.1*numpy.arange(-44,44,0.2)),2,0,4000,"Monthly",time)
        #ind_latit6,fwt_gin=plot_salinity_profile(lon,lat,z,cube_so,cube_vo,lonadj-20,lonadj+20,+65,82,'GIN',list(0.1*numpy.arange(-44,44,0.2)),2,0,5000,"Monthly",time)
        ind_latit7,fwt_lab=plot_salinity_profile(lon,lat,z,cube_so,cube_vo,lonadj-70,lonadj-40,45,72,'LABRADOR',list(0.1*numpy.arange(-44,44,0.2)),2,0,6000,"Monthly",time)
        ind_latit1,fwt_bar=plot_salinity_profile(lon,lat,z,cube_so,cube_vo,lonadj+17.5,lonadj+52.5,79.7,80.3,'BARENTS',list(0.1*numpy.arange(-44,44,0.2)),6,0,450,"Monthly",time)
        ind_latit3,fwt_green=plot_salinity_profile(lon,lat,z,cube_so,cube_vo,lonadj-41,lonadj-21,64.6,66,'GREEN-ICE',list(0.1*numpy.arange(-44,44,0.2)),8,0,3000,'Monthly',time)
        #ind_latit3,fwt_ice=plot_salinity_profile(lon,lat,z,cube_so,cube_vo,lonadj-15.7,lonadj+12.4,64.6,65,'ICE-NORWAY',list(0.1*numpy.arange(-44,44,0.2)),10,ind_latit3,5000,'Monthly', time)
        ind_latit4,fwt_baff=plot_salinity_profile(lon,lat,z,cube_so,cube_vo,lonadj-78.4, lonadj-63, 78, 79,'NORTHBAFFIN',list(0.1*numpy.arange(-44,44,0.2)),12,0,1500, 'Monthly', time)

        fwt_month_bering.extend(fwt_bering)
        #fwt_month_fram.extend(fwt_fram)
        #fwt_month_gin.extend(fwt_gin)
        fwt_month_lab.extend(fwt_lab)
        fwt_month_bar.extend(fwt_bar)
        fwt_month_green.extend(fwt_green)
        #fwt_month_ice.extend(fwt_ice)
        fwt_month_baff.extend(fwt_baff)
    save_netcdf(fwt_month_bering,name, 'Bering_Strait')
    #save_netcdf(fwt_month_fram,name, 'Fram_Strait')
    #save_netcdf(fwt_month_gin,name, 'GIN')
    save_netcdf(fwt_month_lab,name, 'Labrador')
    save_netcdf(fwt_month_bar,name, 'Barents')
    save_netcdf(fwt_month_green,name, 'Green-Ice')
    #save_netcdf(fwt_month_ice,name, 'Ice-Norway')
    save_netcdf(fwt_month_baff,name, 'North_Baffin')






    #pyplot.figure()
    #pyplot.plot(fwt_month_fram)
    #pyplot.show()
if __name__ == "__main__":
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        FreshWaterTransport(config).compute()
