#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 16:19:40 2017

@author: bmueller
"""
import datetime
import math
from itertools import compress
import numpy as np
#from scipy import fft
from scipy import signal
from scipy import interpolate
from scipy.stats import t
import scipy.fftpack as fftpack
import scipy.optimize as optimize
from scipy.ndimage.filters import uniform_filter1d
import statsmodels.api as sm
from .models import LinearTrend, SineSeason3  # , SineSeasonk, SineSeason1
import matplotlib.pyplot as plt

#from bfast import BFAST # not working with rpy2

####
# Additional functions


def mysine(array, para1, para2, para3):
    """
    simple sine model
    """
    return para1 * np.sin(para2 * array + para3)


def multisine(array, para1, para2, para3):
    """
    multiple sine model
    """
    init = array*0.
    for i, _ in enumerate(para1):
        init += mysine(array, para1[i], para2[i], para3[i])
    return init


def wrapper_multisine(array, *args):
    """
    wrapper for multiple sine model
    """
    equal_len = int(1./3.*len(args))
    para1, para2, para3 = list(args[:equal_len]), \
        list(args[equal_len:2*equal_len]), \
        list(args[2*equal_len:3*equal_len])
    return multisine(array, para1, para2, para3)


# def mygauss(x, sigma):
#    """
#    simple gauss distribution model without defining mu
#    """
#    global mu
#    return mlab.normpdf(x, mu, sigma)


def nargmax(array, num=1):
    """
    get the n biggest values of x (np.array)
    """
    args = []
    array = array.astype(float)
    for _ in range(num):
        argmax = array.argmax()
        args.append(argmax)
        array[argmax] = -np.inf
    return args


####


class TempStab(object):
    """
    module to calculate measures of temporal stability
    """

    def __init__(self, dates, array, **kwargs):
        """
        Parameters
        ----------
        array: numpy nd.array
        """

        self.dates = dates[:]
        self.array = array.copy()
        self.prep = self.array.copy()
        self.season = []
        self.season_mod = self.__do_nothing__
        self.numdate = np.linspace(1, 100, num=len(dates))
        self.filler = None
        self.trend = None
        self.homogenized = None
        self.__season_removed__ = None
        self.__trend_removed__ = None
        self.__orig__ = self.array.copy()
        self.__prep_orig__ = self.array.copy()
        self.__numdate_orig__ = self.numdate.copy()
        self.__identified_gaps__ = []
        self.__min_time_step__ = None

        # constants:
        # numeric tolerance for shift
        self.frequency = 365.2425
        self.periods = np.array([1])

        # Initialize methods
        self.break_method = kwargs.get('breakpoint_method', None)
        self.__default_season__ = kwargs.get('default_season',
                                             [1., 0., 0., 0., 0., 0.])
        self.__homogenize__ = kwargs.get('homogenize', None)
#        self.__timescale__ = kwargs.get('timescale', 'months')
        self.__num_periods__ = kwargs.get('max_num_periods', 3)
        self.__periods_method__ = kwargs.get('periods_method', "autocorr")
        # TODO smoothing filter size is a hard question
        self.smoothing = kwargs.get('smoothing4periods', 3)
        self.temporal_res = kwargs.get('temporal_resolution', 1.)
        self.__detrend_bool__ = kwargs.get('detrend', False)
        self.__deseason_bool__ = kwargs.get('deseason', False)

        self.__run__ = kwargs.get('run', False)

        # Set methods, conversion and integritiy checks
#        self.__set_temp_fac__()
        self.__temp_fac__ = self.temporal_res
        self.__set_break_model__()
        self.__set_season_model__()
        self.__set_time__()
        self.__set_periods__()
        self.__scale_time__()
        self.__check__()

        # Run routine
        self.__preprocessing__()

        if self.__run__:
            self.analysis(homogenize=self.__homogenize__, **kwargs)

    def __do_nothing__(self, **kwargs):
        """
        does nothing!
        """
        pass

    def set_frequency(self, frequency):
        """
        set the annual frequency of numerical dates for finding gaps
        """
        self.frequency = frequency
        
    def __set_temp_fac__(self):
        """
        set the temporal resolution factor
        """
        
        options = {"daily" : self.__daily_res__,
                   "monthly" : self.__monthly_res__,
                   "yearly" : self.__yearly_res__,
#                   "nhrly" : self.__nhrly_res__() # not implemented properly
                   }
        options[self.temporal_res]()
        
    def __daily_res__(self):
        """
        set the temporal resolution factor for daily data
        """
        self.__temp_fac__ = 1.
        
    def __nhrly_res__(self, n):
        """
        set the temporal resolution factor for n-hourly data
        """
        self.__temp_fac__ = 1./24.*n
        
    def __monthly_res__(self):
        """
        set the temporal resolution factor for monthly data
        """
        self.__temp_fac__ = 365.2425/12.
        
    def __yearly_res__(self):
        """
        set the temporal resolution factor for yearly data
        """
        self.__temp_fac__ = 365.2425

#    def set_smoothing4periods(self, smoothing_window):
#        """
#        set smoothing window for simpler evaluation of periods (denoise)
#        """
#        self.smoothing = smoothing_window

    def __check__(self):
        """
        checks the integrity of the input
        """
        assert self.break_method is not None, \
            'Method for breakpoint estimation needs to be provided'
        assert len(self.dates) == len(self.array), \
            'Timeseries and data need to have the same dimension'
        assert isinstance(self.dates[0], datetime.datetime), \
            'Dates are not of type datetime'

    def __set_season_model__(self):
        # note that the seasonal model needs to be consistent with the one
        # used for the breakpoint estimation!!!
        # thus if BFAST is used, a similar seasonal model needs to be applied
        self.season_mod = SineSeason3
        # at the moment use the 3fold sine model like
        # in Verbesselt et al. (2010), Eq. 3/4

    def __set_time__(self):
        """
        convert datetime array into float array
        """

        def toordinaltime(this_datetime):
            """
            calculate ordinals plus time fraction
            """
            time = this_datetime.hour * 1. + \
                (this_datetime.minute * 1. +
                 (this_datetime.second * 1. +
                  this_datetime.microsecond/1000000.)/60.)/60.

            days = this_datetime.toordinal() + time/24.

            return days

        self.numdate = np.array([toordinaltime(d) for d in self.dates])


    def __set_periods__(self):
        """
        determine periodic frequencies of data sampling from data
        """
        
        print('Calculating the periods of the data ...')
        self.__detrend__()
        
        if self.__periods_method__ == "autocorr":
            self.__periods_autocorr__()
        elif self.__periods_method__ == "fft_optim":
            self.__periods_fft_optim__()
        else:
            assert False, "Method for determining periods undefined!"
            
        # correct for temporal resolution

        self.periods = [p*self.__temp_fac__ for p in self.periods]
            
        self.__trend_removed__ = None
        self.prep = self.array.copy()
        print('Calculating periods finished.')
            
        
    def __periods_autocorr__(self):
        """
        periods from autocorrelation
        """
        
        periods = []
        
        # normalized values
        yunbiased = self.prep-np.mean(self.prep)
        ynorm = np.sum(yunbiased**2)
        acor = np.correlate(yunbiased, yunbiased, "same")/ynorm
        
        x=np.arange(len(acor)/2)
    
        x=np.sort(np.concatenate((-x,np.array([0]),x)))

        while len(periods) < self.__num_periods__:
            
            maxs = signal.argrelmax(acor, order=3)[0]
            
            p=np.diff(x[maxs]).mean()
                
            if not math.isnan(p):
                periods.append(p)
            else:
                break
            
            acor = acor[maxs]
            
            x = x[maxs]
        
        self.periods = periods
        
        
    def __periods_fft_optim__(self):
        """
        periods from fft optimization
        """
        
        keep = [not math.isnan(pp) for pp in self.prep]
        loc_prep = np.array(list(compress(self.prep, keep)))
        loc_numdate = np.array(list(compress(self.numdate, keep)))

        # presmoothing needed
        loc_prep = uniform_filter1d(loc_prep, size=self.smoothing)
        self.prep[np.array(keep)] = loc_prep
        
        periods = []
      
        while len(periods) < self.__num_periods__:

            try:
                prephat = fftpack.rfft(loc_prep)
                idx = (prephat**2).argmax()

                freqs = fftpack.rfftfreq(prephat.size,
                                         d=np.min(np.abs(
                                                 np.diff(loc_numdate)
                                                 ))/(2*np.pi))
                frequency = freqs[idx]

                amplitude = loc_prep.max()
                guess = [amplitude, frequency, 0.]

                keep = [not math.isnan(pp) for pp in self.prep]
                loc_prep = np.array(list(compress(self.prep, keep)))
                loc_numdate = np.array(list(compress(self.numdate, keep)))

                (amplitude, frequency, phase), pcov = optimize.curve_fit(
                    mysine, loc_numdate, loc_prep, guess)

                period = 2*np.pi/frequency

                this_sine = mysine(self.numdate, amplitude, frequency, phase)
                self.prep -= this_sine

                periods.append(period)
                
                # reoccurences much longer than the time series
                # don't make sense
                keep = np.abs(periods) < len(self.prep)
                periods = list(compress(periods, keep))

            except RuntimeError:
                break
        
        self.periods = periods
        

    def __scale_time__(self):
        """
        scale to years, months, days
        """
        # TODO:
        # * change self.numdate
        # * change self.frequency
        # * change self.period
        # * maybe change self.dates
        pass

    def __preprocessing__(self, **kwargs):
        """
        perform preprocessing of the actual data
        The following options are available

        detrend : bool
            perform a linear detrending of the original data
        remove_season : bool
            remove the seasonality according to the model chosen
            if this option is used, then the overall linear trend
            is removed first and then the seasonality is removed thereafter
        """

        self.prep = self.__orig__.copy()

        # in case that season shall be removed do detrending first
        if self.__deseason_bool__:
            self.__detrend_bool__ = True

        #  remove linear trend res = x - (slope*t + offset)
        if self.__detrend_bool__:
            self.__detrend__()
            print("detrended")

        # remove seasonality
        if self.__deseason_bool__:
            count=0
            while count<1:
                try:
                    self.__deseason__()
                    count = 99
                except:
                    count += 1    
            
            print("deseasoned")

    def __deseason__(self):
        print('Deseasonalization of the data ...')
        keep = [not math.isnan(pp) for pp in self.prep]
        loc_prep = np.array(list(compress(self.prep, keep)))
        loc_numdate = np.array(list(compress(self.numdate, keep)))
#        sins = SineSeason1(loc_numdate, loc_prep, f=self.frequency)
#        sins.fit()  # estimate seasonal model parameters
#        self.__season_removed__ = sins.eval_func(self.numdate)
#        self.prep -= self.__season_removed__
        self.__season_removed__ = self.prep*0.
        self.__season_removed__[:] = 0.

        # presmoothing needed for better access on periods?
        # loc_prep = uniform_filter1d(loc_prep, size=self.smoothing)

        # setting best guess and bounds for seasons
        amplitudes = list(np.repeat((loc_prep.max()-loc_prep.min())/2.,
                                    len(self.periods)))
        freqs = [2*np.pi/p for p in self.periods]
        guess = amplitudes + freqs + list(np.repeat(1., len(self.periods)))

        ubound = list(np.repeat(np.inf, len(self.periods))) + \
            [f+10**(-16) for f in freqs] + \
            list(np.repeat(np.inf, len(self.periods)))
        lbound = list(np.repeat(0., len(self.periods))) + \
            [f-10**(-16)  for f in freqs] + \
            list(np.repeat(-np.inf, len(self.periods)))
            
        # fitting the curves to periods
        params, pcov = optimize.curve_fit(wrapper_multisine,
                                          loc_numdate,
                                          loc_prep,
                                          guess,
                                          bounds=(lbound, ubound))

        # updating periods
        self.periods = [2*np.pi/p for p in
                        list(params[len(self.periods):2*len(self.periods)])]

        self.__season_removed__ += wrapper_multisine(self.numdate,
                                                     *params)
        self.prep -= self.__season_removed__

        print('Deseasonalization finished.')

    def __detrend__(self):
        """
        substracting linear trend
        """
        print('Detrending of the data ...')
        keep = [not math.isnan(pp) for pp in self.prep]
        loc_prep = np.array(list(compress(self.prep, keep)))
        loc_numdate = np.array(list(compress(self.numdate, keep)))
        lint = LinearTrend(loc_numdate, loc_prep)
        lint.fit()
        self.__trend_removed__ = lint.eval_func(self.numdate)
        self.prep -= self.__trend_removed__
        print('Detrending finished.')

    def analysis(self, homogenize=None, **kwargs):
        """
        analyse for breakpoints
        """
        
        self.__prep_orig__ = self.prep.copy()
        if np.ma.isMaskedArray(self.prep):
            self.prep = self.prep.data
        self.__numdate_orig__ = self.numdate.copy()

        assert homogenize is not None, \
            'Homogenization argument need to be explicitely provided'

        # fill data gaps if existing
        self.__identify_gaps__()
        self.__fill_gaps__()

        # identify breakpoints in time based;
        # returns an array of indices with breakpoints
        self.breakpoints = np.sort(self.__calc_breakpoints__(self.array,
                                                             **kwargs))
        
        # if first element seems like breakpoint, delete it
        self.breakpoints=self.breakpoints[np.logical_not(self.breakpoints==1)]
        
        # if two breakpoints follow each other, delte first
        while np.any(np.diff(self.breakpoints) == 1):
            self.breakpoints = self.breakpoints[np.logical_not(
                    np.append((np.diff(self.breakpoints)==1)
                    ,False))]
            

        if len(self.breakpoints) > 0:
            # estimate linear trend parameters for
            # each section between breakpoints
            self.trend = self.__get_trend_parameters__(self.breakpoints,
                                                       self.prep)

            # in case that breakpoints occur, perform a normalizaion of the
            # timeseries which corresponds to a removal of the segmentwise
            # linear offsets of the linear trends, the slope is not corrected
            # for
            if homogenize:
                self.homogenized = self.__homogenization__()
#                fig = plt.figure()
#                ax = fig.add_subplot(111)
#                ax.plot(self.prep, label = "prep")
#                ax.plot(yn, label = "yn")
#                plt.legend()
#                plt.show()
#                assert False
            else:
                self.homogenized = self.prep.copy()
        else:
            self.trend = None
            self.homogenized = self.prep.copy()

        # perform final trend estimation
        keep = [not math.isnan(pp) for pp in self.__orig__]
        loc_orig = np.array(list(compress(self.__orig__, keep)))
        loc_numdate = np.array(list(compress(self.numdate, keep)))
        L_orig = LinearTrend(loc_numdate, loc_orig)  # TODO uncertatinties???
        L_orig.fit()  # should store also significance information if possible
        L_res = LinearTrend(self.numdate, self.homogenized)  # TODO uncertatinties???
        L_res.fit()  # should store also significance information if possible

        res = {}
        res.update({'homogenized_trend' : {'slope' : L_res.param[0], 'offset' : L_res.param[1]}})
        res.update({'original_trend' : {'slope' : L_orig.param[0], 'offset' : L_orig.param[1]}})
        res.update({'yn' : self.homogenized})
        res.update({'yorg' : self.array})
        res.update({'yraw' : self.__orig__})
        res.update({'season' : self.__season_removed__})
        res.update({'breakpoints' : self.breakpoints})
        res.update({'nbreak' : len(self.breakpoints)})
#        print(res)
#        import time
#        time.sleep(50)
        return res
    
    def __get_indices__(self, bp, i, n):
        """
        returns indices of boundaries so these can be
        directly used for indexing
        """
        i1 = bp[i]
        if i < len(bp)-1:
            i2 = bp[i+1]
        else:
            i2 = n
        return min([i1, i2]), max([i1, i2])

    def __get_trend_parameters__(self, bp, x, remove_seasonality=False):
        """
        calculate linear trend parameters for each segment between breakpoints

        TODO: significance of trends xxx, consideration of uncertainties of
        samples
        how to deal wih seasonality here ??? would need to be removed ???

        Parameters
        ----------
        bp : ndarray, list
            list with indices of breakpoints
        x : ndarray
            data array; this is explicitely provided as argument as it might be
            constructed from e.g. detrended or deseasonalized timeseries
        """
        trends = []

        if len(bp) == 0:
            return trends
        xx = self.numdate[0:bp[0]]
        yy = x[0:bp[0]]
        if remove_seasonality:
            assert False

        if bp[0]!=0:
        # estimate piecewise linear trend
            L = LinearTrend(xx, yy) ### TODO uncertainties???
            L.fit()
            L.__i1__ = 0
            L.__i2__ = bp[0]
            trends.append(L)
        else:
            L = LinearTrend(self.numdate[0:2], np.repeat(x[0],2)) ### TODO uncertainties???
            L.fit()
            L.__i1__ = 0
            L.__i2__ = bp[0]
            trends.append(L)
        n = len(x)
        for i in xrange(len(bp)):
            i1, i2 = self.__get_indices__(bp, i, n)
            if self.numdate[i1]==self.numdate[-1]:
                L = LinearTrend(self.numdate[(i1-1):i2], np.repeat(x[i1],2))   
                L.__i1__= i1
                L.__i2__ = i2
            else:
                L = LinearTrend(self.numdate[i1:i2], x[i1:i2])
                L.__i1__= i1
                L.__i2__ = i2
            L.fit()
            trends.append(L)
        return trends
    
    def __homogenization__(self):
        """
        perform a homogenization of the timeseries
        by removing the detected breakpoint impacts
        """
        # should use self.trend which has been estimated already before
#        res1 = self.__remove_resid_offset__(self.numdate, self.trend)
        # based on bias of residuals
        res2 = self.__remove_trend_offset__(self.numdate, self.trend)
        # based on linear trends of residuals
        
#        plt.plot(res1,res2)
#        plt.show()
        
        return res2
        
    
    def __remove_trend_offset__(self, x, trends):
        """
        remove offset of linear trends
        slope is not corrected for
        this routine is used to remove jumps
        caused by structural breakpoints in the timeseries
        """
        r = np.zeros_like(self.array) * np.nan

        O = 0.

        for i in xrange(len(trends)-1):
            T1=trends[i]  
            T2=trends[i+1]
            hlp = np.append(r[:T1.__i1__], self.array[T1.__i1__:T1.__i2__])
            F = T2.eval_func(x[T2.__i1__:T2.__i2__])
            if T1.__i1__==T1.__i2__:
                index=np.array(x[T1.__i1__])
                L_val = T1.eval_func(index)
            else:
                index=x[T1.__i1__:T1.__i2__]
                L_val = T1.eval_func(index)[-1]
            O = (L_val-F[0])
#            plt.plot(x[T1.__i1__:T2.__i2__],np.append(L,F))
#            plt.plot(x[T1.__i1__:T2.__i2__],np.append(L-O,F))
#            plt.show()
            hlp -= O
            r[:T1.__i2__] = hlp*1.
#            plt.plot(self.array + r*0, label = "orig")
#            plt.plot(r, label = "corr")
#            plt.legend()
#            plt.show()
        r[T2.__i1__:T2.__i2__] = self.array[T2.__i1__:T2.__i2__]
#        plt.plot(self.array + r*0, label = "orig")
#        plt.plot(r, label = "corr")
#        plt.legend()
#        plt.show()
        # correct for same mode (?)
#        print("mode")
#        print(mode(r).mode[0],mode(self.array).mode[0],mode(r).mode[0]-mode(self.array).mode[0])
#        print("mean")
#        print(np.mean(r),np.mean(self.array),np.mean(r)-np.mean(self.array))
#        print("median")
#        print(np.median(r),np.median(r),np.median(r)-np.median(r))
        return r
    
    
    def __remove_resid_offset__(self, x, trends):
        """
        remove offset of residual levels
        slope is not corrected for
        this routine is used to remove jumps
        caused by structural breakpoints in the timeseries
        """
        r = self.array.copy()
        
        O = 0.

        for i in xrange(len(trends)-1):
            T1=trends[i]  
            T2=trends[i+1]
            E2 = self.prep[T2.__i1__:T2.__i2__].mean()
            E1 = self.prep[T1.__i1__:T1.__i2__].mean()
            O += E1-E2
            r[T2.__i1__:T2.__i2__] = self.array[T2.__i1__:T2.__i2__]+O
            if i == len(trends)-2:
                r += self.array[-1] - r[-1]
        # correct for same mode (?)
#        print("mode")
#        print(mode(r).mode[0],mode(self.array).mode[0],mode(r).mode[0]-mode(self.array).mode[0])
#        print("mean")
#        print(np.mean(r),np.mean(self.array),np.mean(r)-np.mean(self.array))
#        print("median")
#        print(np.median(r),np.median(r),np.median(r)-np.median(r))
        return r


    def __fill_gaps__(self):
        """
        fill data gaps if existing
        """
        
        if sum(self.__identified_gaps__) > 0:
            keep = [not math.isnan(pp) for pp in self.prep]
            loc_prep = np.array(list(compress(self.prep, keep)))
            loc_numdate = np.array(list(compress(self.numdate, keep)))
            self.filler = self.numdate * 0. - 999.

            if self.__season_removed__ is not None:
                self.__linear__()
            elif self.__trend_removed__ is not None:
                self.__linear__()
                self.__gap_season__(loc_numdate, loc_prep)
            else:
                self.__linear__()
                self.__gap_season__(loc_numdate, loc_prep)
                self.__gap_trend__(loc_numdate, loc_prep)

            self.prep[self.__identified_gaps__] = \
                self.filler[self.__identified_gaps__]
        else:
            pass


    def __linear__(self):
        """
        produces linear filler without noise
        """
        self.filler = self.filler * 0. + np.mean(self.prep[np.logical_not(
            self.__identified_gaps__)])


    def __gap_season__(self, time, array):
        """
        produces linear filler with season (no noise)
        """
        season = self.season_mod(t=time,
                                 x=array,
                                 f=self.frequency)
        season.fit()
        self.filler = season.eval_func(self.numdate)


    def __gap_trend__(self, time, array):
        """
        produces linear filler with trend (no noise)
        """
        # calculate gaps environments (starts and stops)
        enlarged_gaps = signal.convolve(self.__identified_gaps__,
                                        np.array([1, 1, 1]))[1:-1] != 0
        starts_n_stops = np.logical_xor(enlarged_gaps,
                                        self.__identified_gaps__)
        sns_pos = self.numdate[starts_n_stops]

        # full linear model
        # (get that from detrend? externalize a single function?)
        lint = LinearTrend(time, array)
        lint.fit()

        for i in range(len(sns_pos)/2):
            if sns_pos[(i*2)] <= sns_pos[(i*2+1)]:
                this_gap = np.logical_and(self.numdate > sns_pos[(i*2)],
                                          self.numdate < sns_pos[(i*2+1)])
            else:
                this_gap = np.logical_and(self.numdate > sns_pos[(i*2+1)],
                                          self.numdate < sns_pos[(i*2)])
            self.filler[this_gap] = \
                lint.eval_func(self.numdate[this_gap])


    def __identify_gaps__(self):
        """
        identify data gaps in self.array
        """
        self.__set_up_new_ts__()
        self.__fill_with_na__()
        self.__calculate_na_indices__()
        self.__identified_gaps__ = self.__calculate_na_indices__()


    def __get_mean_timestep__(self):
        """
        calculates the minimum timestep
        """
        return (np.diff(self.numdate)).mean()
    
    def __get_min_timestep__(self):
        """
        calculates the minimum timestep
        """
        return (np.diff(self.numdate)).min()


    def __set_up_new_ts__(self):
        """
        set a new ts based on the minimum timestep if necessary
        """
        num = np.abs((self.numdate.max() - self.numdate.min()) /
                     self.__get_mean_timestep__()) + 1
        new_dates = np.linspace(self.numdate.min(),
                                self.numdate.max(),
                                num=np.ceil(num))

        if self.numdate[0] == min(self.numdate):
            self.numdate = new_dates
        else:
            self.numdate = new_dates[::-1]
            
        if np.all(np.isclose(self.numdate, self.__numdate_orig__)):
            self.numdate = self.__numdate_orig__
        else:
            num = np.abs((self.numdate.max() - self.numdate.min()) /
                         self.__get_min_timestep__()) + 1
            new_dates = np.linspace(self.numdate.min(),
                                    self.numdate.max(),
                                    num=np.ceil(num))
    
            if self.numdate[0] == min(self.numdate):
                self.numdate = new_dates
            else:
                self.numdate = new_dates[::-1]


    def __fill_with_na__(self):
        """
        fill any gap with nan
        """
        new_array = self.prep * np.NAN
        keep = np.isclose(self.numdate, self.__numdate_orig__)
        new_array[keep] = self.prep[keep]
        self.prep = new_array


    def __calculate_na_indices__(self):
        """
        get a boolean array as indices for nan values
        """
        gaps = np.isnan(self.prep)
        if not sum(gaps) == 0:
            print("gaps identified: " + str(sum(gaps)))
        return gaps


    def __calc_breakpoints__(self, x, **kwargs):
        """
        calculating breakpoints based on set function self.Break
        """
        return self.__chosen_break__(x,
                                     start=self.dates[0],
                                     frequency=self.periods,
                                     **kwargs)


    def __set_break_model__(self):
        """
        set break_method
        """

        if self.break_method is None:
            print("No breakpoint method assigned. Just gaps are filled.")
            self.__chosen_break__ = self.__break_none__
        elif self.break_method == 'olssum':
            self.__chosen_break__ = self.__break_olssum__
        elif self.break_method == 'CUMSUMADJ':
            self.__chosen_break__ = self.__break_CUMSUMADJ__
#        elif self.break_method == 'bfast':
#            self.__chosen_break__ = self.__break_bfast__
        elif self.break_method == 'dummy':
            self.__chosen_break__ = self.__break_dummy__
        elif self.break_method == 'wang':
            self.__chosen_break__ = self.__break_wang__
        else:
            assert False, 'ERROR: Unknown breakpoint method'


    def __break_none__(self, x, **kwargs):
        """
        no breakpoint analysis
        Returns
        res : ndarray
            array with breakpoint indices
        """
        return np.array([])


    def __break_wang__(self, x, **kwargs):
        assert False, \
            "Breakpoint method " + self.break_method + " not implemented yet!"


    def __break_dummy__(self, x, **kwargs):
        assert False, \
            "Breakpoint method " + self.break_method + " not implemented yet!"


#    def __break_bfast__(self, x, **kwargs):
#        """
#        calculate breakpoints in time using the BFAST method
#        Parameters
#        ----------
#        x : ndarray
#            data array
#
#        arguments in kwargs
#        start : datetime
#            start of timeseries
#        frequency : int
#            frequency of number of samples per year
#            e.g. 23 for 16-daily data = 365./16.
#
#        TODO: reasonable additional parameters!!!
#
#        Returns
#        -------
#        res : ndarray
#            array with indices of breakpoint occurence
#        """
##        assert False, \
##            "Breakpoint method " + self.break_method + " not working yet!"
#            
#        B = BFAST()
#        B.run(x, **kwargs)
#        # return detected breakpoints from last itteration
#        #print B.results[0]['output'][-1]['ciVt']  
#        # this contains the confidence of the breakpoint estimation
#        res = B.results[0]['output'][-1]['bpVt']
#        
#        try:  # capture that no breakpoints are detected
#            n = len(res)
#            del n
#        except:
#            if res == 0:
#                res = []
#            else:
#                assert False, 'CASE not covered yet'
#        return np.asarray(res).astype('int')


    def __break_olssum__(self, x, **kwargs):
        # threshold for change detected        
        # 1% 
        self.__thresh_change__ = 0.01
        
        # remove overall linear trend using OLS
        T = sm.add_constant(self.numdate)
        model = sm.OLS(x,T)
        results = model.fit()
        
        # initial result
        r=[]

        # estimate potential breakpoints from residual timeseries
        try:
            r = self.__get_breakpoints_spline__(self.numdate, results.resid)
        except:
            try:
                r = self.__get_breakpoints_spline__(np.flip(self.numdate,0),
                                                    np.flip(results.resid,0))
                r = [len(self.numdate)-ri for ri in r]
            except:
                print("There is something wrong with the datasets" +\
                      "for calculating breaking points." +\
                      " Flipping does not produce right order.")
        return(r)
        
        
    def __break_CUMSUMADJ__(self,x,**kwargs):
        # doi: 10.1016/j.jhydrol.2014.12.002
        
        alpha=0.05 #TODO make this an input
        
        m = len(x)/2
             
        csx = x.cumsum()
        
        CUMSUMADJ=[]
        
        for j in np.arange(len(x)):
            CUMSUMADJ.append(csx[j] - j/float(m) * csx[m])
            
        
        lin=np.polyfit(np.arange(len(CUMSUMADJ)),CUMSUMADJ,1)
        lin_fn=np.poly1d(lin)
        errors=lin_fn(np.arange(len(CUMSUMADJ)))-CUMSUMADJ
        std_error=np.sqrt(1./(len(x)-2.)*sum(errors**2))
        quantile = t.ppf(1.-alpha, len(x)-2)
        deltaC=std_error*quantile
        pot_bps=abs(errors)>deltaC
        pot_bps2=self.__single_peaks_CUMSUMADJ__(errors,pot_bps)
        r=list(compress(xrange(len(pot_bps2)), pot_bps2))
        return(r)

    def __single_peaks_CUMSUMADJ__(self,x,choice)  :
        #if there are 2 breakpoints following each other, it's the bigger one!
        new_choice=np.repeat(False,len(x))
        i=0
        while i<len(x):
            v=[]
            if choice[i]:
                j=i
                while (choice[min([j,len(x)-1])] and j<len(x)):
                    v.append(x[j])
                    j+=1
                new_choice[i+np.argmax(np.abs(v))]=True
                i=j
            else:
                i+=1
        return(new_choice)
                    
        
    def __get_breakpoints_spline__(self, x, y):
        # estimate breakpoints using splines
        #http://stackoverflow.com/questions/29382903/
        #       how-to-apply-piecewise-linear-fit-in-python
        yn = interpolate.splrep(x, y, k=1, s=0)
        d1 = np.abs(interpolate.splev(x, yn, der=1))  
        # absolute of first derivative
        
        # now find the results that are not similar to others
        # first sort results
        sidx = d1.argsort()
        idx = np.arange(len(d1))[sidx]
        
        ds = d1[d1.argsort()]
        
        r = []
        for i in xrange(len(ds)-1,-1,-1):
            if ds[i] > 1.E-6:
                # check if ratio between subsequent slope values changes 
                # by more than self.__thresh_change__ --> does not work if 
                # jumps with same magnitude occur!
                if np.abs(1.-(ds[i]/ds[i-1])) > self.__thresh_change__:   
                    r.append(idx[i]+1)
                else:
                    break
                
        return(r)
        
        
    def reanalysis(self):
        """
        redo all calulations after gapfilling and homogenization
        """
        # reinitialize
        self.array = None
        
        if self.homogenized is not None:
            homogenize = True
            self.array = self.homogenized.copy()
        else:
            if self.__season_removed__ is not None:
                self.array = self.prep + \
                    self.__trend_removed__ + self.__season_removed__
            elif self.__trend_removed__ is not None:
                self.array = self.prep + self.__trend_removed__
            else:
                self.array = self.prep.copy()
                
        self.prep = self.array.copy()
        self.__orig__ = self.array.copy()
        self.__prep_orig__ = self.array.copy()
                
        # resetting preprocessing
        self.prep = self.array.copy()
        self.periods = np.array([1])
        self.__trend_removed__ = None
        self.__season_removed__ = None
        
        # preprocessing
        self.__set_periods__()
        self.__preprocessing__()
        
        # resetting analysis
        self.homogenized = None
        
        # analysis
        self.analysis(homogenize=homogenize)
            
        
        
