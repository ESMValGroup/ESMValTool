import os
import sys
import cf_units as unit
import iris
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from scipy.stats import pearsonr
from iris.coord_categorisation import add_month_number
from iris.coord_categorisation import add_day_of_year
from scipy.fftpack import rfft, irfft, fftfreq


def WinterExtract(cube):
    if not cube.coords('month_number'):
        iris.coord_categorisation.add_month_number(cube, 'time', name='month_number')
    return cube.extract(iris.Constraint(month_number=lambda cell: cell >= 11 or cell <= 4))

def SummerExtract(cube):
    if not cube.coords('month_number'):
        iris.coord_categorisation.add_month_number(cube, 'time', name='month_number')
    return cube.extract(iris.Constraint(month_number=lambda cell: cell >= 5 and cell <= 10))

def clmDayTLL(cube):
    if not cube.coords('day'):
        iris.coord_categorisation.add_day_of_year(cube, 'time', name='day')
    # Daily climatology
    return cube.aggregated_by('day', iris.analysis.MEAN)

def smthClmDayTLL(cube, nhar):
    # Now do a smooth climatology
    ntim, nlat, nlon = cube.shape
    ft = rfft(cube.data, axis=0)

    #ft[nhar] = 0.5*ft[nhar]
    #ft[nhar+1:] = 0.

    ft[nhar:-nhar] = 0.
    ft[nhar] = 0.5*ft[nhar]

    cube_sm = cube.copy()
    cube_sm.data = irfft(ft, axis=0)
    return cube_sm

def calcDayAnomTLL(cube, clim):
    if not cube.coords('day'):
        iris.coord_categorisation.add_day_of_year(cube, 'time', name='day')

    clim_days = clim.coord('day').points
    cube_days = cube.coord('day').points
    # Compute anomalies
    for n, day in enumerate(clim_days):
        inds = np.where(cube_days == day)
        cube.data[inds] -= clim.data[n]
    return cube

def regrid2obs(my_cubes):
    xdir = os.path.dirname(os.path.abspath(__file__))
    base_cube = iris.load_cube(os.path.join(xdir, 'obs_grid.nc'))[0]

    #print base_cube.regrid
    base_cube.coord('latitude').coord_system = None
    base_cube.coord('longitude').coord_system = None
    my_cubes.coord('latitude').coord_system = None
    my_cubes.coord('longitude').coord_system = None

    # For lat/lon regriding, make sure coordinates have bounds
    if my_cubes.coord('longitude').bounds is None:
        my_cubes.coord('longitude').guess_bounds()
    if my_cubes.coord('latitude').bounds is None:
        my_cubes.coord('latitude').guess_bounds()
    if base_cube.coord('longitude').bounds is None:
        base_cube.coord('longitude').guess_bounds()
    if base_cube.coord('latitude').bounds is None:
        base_cube.coord('latitude').guess_bounds()
    reg_cube = my_cubes.regrid(base_cube, iris.analysis.Linear())
    return reg_cube

def LanczosBandPassWeights(fc2, fc1, n):
    if n % 2 == 0:
        sys.exit("Number of weights should be an odd number!")
    else:
        n_2 = (n - 1) / 2
        if n_2 < (1.3 / (fc2 - fc1)):
            optWgt = (1.3 / (fc2 - fc1) + 1)
            sys.exit("Number of weights should be > " + str(optWgt))
        else:
            w = np.zeros([n])
            k = np.arange(-n_2, n_2 + 1, 1)
            k[n_2] = -999. # putting a dummy value
            sigma = np.sin(np.pi * k / n) / (np.pi * k / n)
            t01 = np.sin(2 * np.pi * fc2 * k) / (np.pi * k)
            t02 = np.sin(2 * np.pi * fc1 * k) / (np.pi * k)
            w = (t01 - t02) * sigma
            # When k == 0 then
            w[n_2] = 2 * (fc2 - fc1)
            return w

def RollingTimeWindowFilter(var, wgts):
    n = wgts.size
    n_2 = (n - 1) / 2
    ntime, nlat, nlon = var.shape
    vardata = var.data
    varfdata = vardata
    #varfdata[:,:,:] = -999.0
    times = np.arange(ntime)
    #print 'Filtering...'
    for la in range(nlat):
        #print 'Lat = ' + str(la)
        for lo in range(nlon):
            varfdata[n_2:ntime - n_2, la, lo] = np.convolve(vardata[:, la, lo], wgts, 'valid')
            varfdata[:, la, lo] = ma.masked_where(times <= n_2, varfdata[:, la, lo])
            varfdata[:, la, lo] = ma.masked_where(times >= ntime - n_2, varfdata[:, la, lo])
            #varfdata[:n_2+1,la,lo]  = -999.0
            #varfdata[ntime-n_2:,la,lo]  = -999.0
    var.data = varfdata
    return var

def Filter(var):
    # Filtering using predetermined filter coeffecients
    # FilterCoef.dat provided with the fortran code.
    filt_var = var.copy()
    mjo_dir = os.path.abspath(os.path.dirname(__file__))
    beta = np.loadtxt(os.path.join(mjo_dir, 'FilterCoef.dat'))
    lag = len(beta)
    ntime, nlat, nlon = var.shape
    print ntime, nlat, nlon
    sum = var.data[0]

    for it in range(ntime - 1):
        sum[:, :] = 0.
        for k in range(lag):
            k1 = ((ntime + it - k - 1) % (ntime))
            k2 = ((it + k) % (ntime))
            #print k, k1, k2
            sum += beta[k] * (var.data[k1] + var.data[k2])
        filt_var.data[it] = sum
    return filt_var

def anomRemoveClimAnnualCycle(var):
    # remove a climatological annual cycle from the data
    # to create anomaly.
    time = var.coord('time')
    year, month, day = getDates(time)
    mmdd = np.array(['%02d' % month[i] + '%02d' % day[i] for i in range(len(year))])

    anom = var.copy()
    # get unique values
    # taken from http://mattdickenson.com/2011/12/31/find-unique-values-in-list-python/
    mmdd_unique = np.array(sorted(list(set(mmdd))))

    # removing the climatological mean value for each day.
    for md in mmdd_unique:
        inds = np.where(mmdd == md)
        anom.data[inds] = var.data[inds] - var[inds].collapsed('time', iris.analysis.MEAN).data
    return anom

def splitTime(date):
    '''Split date string into yyyy, mm, dd integers'''
    d = str(date)
    year = int(d[0:4])
    month = int(d[5:7])
    day = int(d[8:10])
    return year, month, day

def getDates(time):
    '''splits the iris time coordinate into yyyy, mm, dd'''
    if time.units.calendar == 'gregorian':
        dates = unit.num2date(time.points, str(time.units), unit.CALENDAR_GREGORIAN)
    if time.units.calendar == 'standard':
        dates = unit.num2date(time.points, str(time.units), unit.CALENDAR_STANDARD)
    else:
        dates = unit.num2date(time.points, str(time.units), time.units.calendar)
    try:
        dates
    except NameError:
        print dates + " WASN'T defined after all!"
    else:
        year = np.zeros(len(dates), dtype=int)
        month = np.zeros(len(dates), dtype=int)
        day = np.zeros(len(dates), dtype=int)
        for i, date in enumerate(dates):
            year[i], month[i], day[i] = splitTime(date)
    return year, month, day

def region(coords):
    if len(coords) == 4 :
        loni, lati, lonf, latf = coords
        xy_slice = iris.Constraint(longitude=lambda cell:loni <= cell <= lonf,
                                 latitude=lambda cell: lati <= cell <= latf)
    else :
        print 'Wrong number of elements in the coordinate array.'
    return xy_slice
def pptocube(run, vardict, stashcode):
    import maverick.lib.loaddata as loaddata
    fname = vardict['fname']
    #print 'Dealing with PP stashsplit...'

    if not os.path.exists(fname):
        keys = dict() # This dictionary contains a series of pp header filters similiar to TIDLs ss() function including dt
        if vardict['varname'] == 'U850':
            keys['lblev'] = 850
        if vardict['varname'] == 'U200':
            keys['lblev'] = 200
        period = 'daily' # can be 'annual', 'seasonal', 'monthly'
        msi_stash_str = 'm01s{0[0]:02d}i{0[1]:03d}'.format(divmod(stashcode, 1000))
        cube = loaddata.load_run_ss(run, period, msi_stash_str, **keys)

        # Extract region
        coords = [0, -30, 360, 30]
        cube = cube.extract(region(coords))

        #print 'PP cube from Stashsplit'
        #print cube
        # Regridding data to 2.5 x 2.5 grid...
        print '*' * 50
        print 'Regridding data to 2.5 x 2.5 grid...'
        print '*' * 50
        cube = regrid2obs(cube)
        #print 'Regridded cube'
        #print cube
        iris.save(cube, fname)
        print '%s written. ' % fname
    else:
        print '%s exists. Reading Netcdf file...    ' % fname
        cube = iris.load_cube(fname)
        #print 'Netcdf cube'
        #print cube
    return cube
def AreaAverage(var, coords):
    # Generate area-weights array. As e1 and a1b are on the same grid we can
    # do this just once and re-use.
    # This method requires bounds on lat/lon coords, so first we must guess
    # these.
    ts = var.extract(region(coords))
    if ts.coord('latitude').bounds is None:
        ts.coord('latitude').guess_bounds()
    if ts.coord('longitude').bounds is None:
        ts.coord('longitude').guess_bounds()
    ts_grid_areas = iris.analysis.cartography.area_weights(ts)

    # Now perform an area-weighted collapse for each dataset:
    ts_area_mean = ts.collapsed(['latitude', 'longitude'],
                                                   iris.analysis.MEAN,
                                                   weights=ts_grid_areas)
    return ts_area_mean

def lagcorr(x, y, lag=None, verbose=False):
    '''Compute lead-lag correlations between 2 time series.

    <x>,<y>: 1-D time series.
    <lag>: lag option, could take different forms of <lag>:
          if 0 or None, compute ordinary correlation and p-value;
          if positive integer, compute lagged correlation with lag
          upto <lag>;
          if negative integer, compute lead correlation with lead
          upto <-lag>;
          if pass in an list or tuple or array of integers, compute
          lead/lag correlations at different leads/lags.

    Note: when talking about lead/lag, uses <y> as a reference.
    Therefore positive lag means <x> lags <y> by <lag>, computation is
    done by shifting <x> to the left hand side by <lag> with respect to
    <y>.
    Similarly negative lag means <x> leads <y> by <lag>, computation is
    done by shifting <x> to the right hand side by <lag> with respect to
    <y>.

    Return <result>: a (n*2) array, with 1st column the correlation
    coefficients, 2nd column correpsonding p values.

    Currently only works for 1-D arrays.
    '''

    if len(x) != len(y):
        print 'Input variables in lagcorr(x,y) are of different lengths.'
        sys.exit()

    #--------Unify types of <lag>-------------
    if np.isscalar(lag):
        if abs(lag) >= len(x):
            raise('Maximum lag equal or larger than array.')
        if lag < 0:
            lag = -np.arange(abs(lag) + 1)
        elif lag == 0:
            lag = [0, ]
        else:
            lag = np.arange(lag + 1)
    elif lag is None:
        lag = [0, ]
    else:
        lag = np.asarray(lag)

    #-------Loop over lags---------------------
    result = []
    if verbose:
        print '<lagcorr>: Computing lagged-correlations at lags:', lag

    for ii in lag:
        if ii < 0:
            result.append(pearsonr(y[:ii], x[-ii:]))
        elif ii == 0:
            result.append(pearsonr(x, y))
        elif ii > 0:
            result.append(pearsonr(y[ii:], x[:-ii]))

    result = np.asarray(result)

    return result

def LeadLagCorr(ts, var_tx):
    # lead-lag correlation
    ntime, nlon = var_tx.shape
    mxlag = 30
    lcorr = np.zeros([2 * mxlag + 1, nlon], dtype=float)
    #lcorr       = var_tx[0:2*mxlag+1, :]
    lags = np.arange(-mxlag, mxlag + 1)
    lons = var_tx.coord('longitude').points

    for lo in np.arange(nlon):
        dum = lagcorr(ts.data, var_tx[:, lo].data, lag=lags, verbose=False)
        lcorr[:, lo] = dum[:, 0]

    return lcorr, lags, lons
#def ceof(var):

def RMM_phases(pcs):
    ntime, neofs = pcs.shape

    pha = np.zeros(ntime, dtype=int)
    for it in np.arange(ntime):
        pc1 = pcs[it, 0]
        pc2 = pcs[it, 1]
        if pc1 != 0:
            tmp = np.arctan(pc2 / pc1)
            if pc2 * pc1 > 0:
                if pc1 < 0:
                    tmp = tmp + np.pi
            else:
                if pc2 > 0:
                    tmp = np.pi + tmp
                if pc2 < 0:
                    tmp = 2 * np.pi + tmp
        else:
            tmp = 10000
        if tmp != 10000:
            if tmp > np.pi and tmp < np.pi + np.pi / 4.:
                pha[it] = 1
            elif tmp > np.pi + np.pi / 4. and tmp < np.pi + 2 * np.pi / 4.:
                pha[it] = 2
            elif tmp > np.pi + 2 * np.pi / 4. and tmp < np.pi + 3 * np.pi / 4.:
                pha[it] = 3
            elif tmp > np.pi + 3 * np.pi / 4. and tmp < 2 * np.pi:
                pha[it] = 4
            elif tmp > 0. and tmp < np.pi / 4.:
                pha[it] = 5
            elif tmp > np.pi / 4. and tmp < 2 * np.pi / 4.:
                pha[it] = 6
            elif tmp > 2 * np.pi / 4. and tmp < 3 * np.pi / 4.:
                pha[it] = 7
            elif (tmp > 3 * np.pi / 4. and tmp < np.pi):
                pha[it] = 8
        else:
            pha[it] = -1
    return pha
