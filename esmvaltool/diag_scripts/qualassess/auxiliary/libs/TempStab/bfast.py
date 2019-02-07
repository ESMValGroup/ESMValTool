"""
BFAST wrapper

example session

library(bfast)
f = bfast(harvest, max.iter=1000)

should give

> f$jump
$x
[1] 2004.652 2004.652

$y
[1] 0.8282047 0.5740273



using this class it would look like the following


"""

# todo how to handle gaps ???

import numpy as np

from rpy2.robjects.packages import importr
from rpy2 import robjects

import numpy as np
import multiprocessing

from multiprocessing import Process

from datetime import datetime

def worker_wrapper(arg):
    args, kwargs = arg
    return _run_bfast(*args, **kwargs)

def _run_bfast(x, maxiter=1000, season='harmonic', frequency=1, start=None, maxbreaks=100):
    assert season in ['harmonic','dummy','none']
    assert start is not None
    try:
        rbfast = importr('bfast')
    except:
        print('ERROR: Package BFAST is not installed yet!')
        return None
    rts = robjects.r['ts']

    # convert to time series
    TS = rts(robjects.FloatVector(x), frequency=365.25, start=start)  # todo set frequencies etc

    # perform structural change analysis
    print('Running BFAST with maxiter=%s and maxbreak=%s' % (maxiter, maxbreaks))
    fit = rbfast.bfast(TS, season=season, max_iter=maxiter, breaks=maxbreaks)  # todo adapt parameters
    return _extract_results(fit)

def _map_output(o):
    Tt = np.array(o[0])
    St = np.array(o[1])
    Nt = np.array(o[2])
    Vt = np.array(o[3])
    # the breakpoint information is extraced in R by a call to the breakpoint function
    # here we only extract the breakdates at the moment
    # much more detail like e.g. significance of breakdate estimation is existing
    # but still need to check how to best extract TODO
    bpVt = np.array(o[4][0])   # dates of breakpoints
    # breakpoints within confidence limits of 2.5 ... 97.5%
    ciVt = np.array(o[6][0])
    Wt = np.array(o[7])
    # bpwkt  TODO bpWt is still missing!!!
    return {'Tt' : Tt, 'St' : St, 'Nt':Nt, 'Vt':Vt, 'Wt' :Wt, 'bpVt' : bpVt, 'ciVt' : ciVt}

def _extract_results(fit):

    res = {}
    Yt = fit.rx(1)
    res.update({'Yt' : np.array(Yt)})

        # the output contains results for each itteration
    output = fit.rx(2)[0]
    routput = []
    for o in output:
        routput.append(_map_output(o))
    res.update({'output' : routput})

    nobp = fit.rx(3)
    res.update({'nobp' : {'Vt' : np.array(nobp[0][0]).astype('bool'), 'Wt' : np.array(nobp[0][1]).astype('bool')} })

    Magnitude = np.array(fit.rx(4)[0])
    res.update({'Magnitude' : Magnitude})

    Mags = np.array(fit.rx(5)[0])
    res.update({'Mags' : Mags})

    Time = np.array(fit.rx(6)[0])
    res.update({'Time' : Time})

    jump = {}
    jump.update({'x' : np.array(fit.rx(7)[0][0])})
    jump.update({'y' : np.array(fit.rx(7)[0][1])})
    res.update({'jump' : jump})

    return res


class BFAST(object):
    def __init__(self, **kwargs):
        self._check()

        # import R objects and functions
        # do this only during initialization to enhance speed

        self.ncpu = kwargs.get('ncpu', multiprocessing.cpu_count())
        #print 'NCPU: ', self.ncpu

    def _check(self):
        pass

    def run(self, x, start=None, frequency=None, maxiter=1000, maxbreaks=100):
        """
        perform actual analysis

        Parameters
        ----------
        start : datetime
            specifies the start time, currently the year and month information is used
        """
        #~ assert len(t) == len(x)
        #~ assert t.ndim == 1
        if x.ndim == 1:
            parallel = False
        elif x.ndim == 2:
            parallel = True
        else:
            assert False

        assert start is not None
        assert frequency is not None

        if parallel:
            print('Doing BFAST analysis in parallel with %s CPUs' % self.ncpu)
            if True:
                # some wrapping is needed to be able to treat the kwargs
                # http://stackoverflow.com/questions/34031681/passing-kwargs-with-multiprocessing-pool-map
                jobs = []
                pool=multiprocessing.Pool(processes=self.ncpu)
                start_date = self._set_date(start)
                keywords = {'start': start_date, 'frequency': frequency, 'maxiter' : maxiter, 'maxbreaks' : maxbreaks}

                for i in xrange(len(x)):
                    jobs.append([x[i]])
                arg = [(j, keywords) for j in jobs]
                self.results = pool.map(worker_wrapper, arg)
            else:
                self.results = []
                for i in xrange(len(x)):
                    print '**', i
                    res = self._run_bfast(x[i], start=self._set_date(start), frequency=frequency, maxiter=maxiter, maxbreaks=maxbreaks)
                    self.results.append(res)
        else:  # no parallel processing
            self.results = [_run_bfast(x, start=self._set_date(start), frequency=frequency, maxiter=maxiter, maxbreaks=maxbreaks)]


    def _set_date(self, d):
        # d should be datetime object
        assert type(d) is datetime
        return robjects.FloatVector(np.array([d.year,d.month]))

    def harvest(self):
        """ load harvest data as an example """
        return np.loadtxt('harvest.txt')


def _xxx():


    #~ import matplotlib.pyplot as plt


    #t = np.arange(12*31) #generate a 30 year timeseries
    #x = 0.2*np.sin(2.*np.pi*1./12.*t)

    #x[100:] += 0.05



    #~ plt.close('all')


    import time


    # run example using the 'harvest' data from BFAST package
    B = BFAST()
    x = B.harvest()
    #~ t = np.arange(len(x))

    # stack data to mimic multiple timeseries
    X = np.vstack([x,x,x,x,x,x,x,x])


    if False:
        f=plt.figure()
        ax=f.add_subplot(111)
        ax.plot(x)
        plt.show()

    # test with individual samples
    print 'Starting ...'
    t = time.time()
    for i in xrange(len(X)):
        print i
        B.run(X[i], start=datetime(2000,4,1), frequency=23)
    elapsed = time.time() - t
    print 'Time for individual processing: ', elapsed

    print 'Starting again ...'
    t = time.time()
    B.run(X, start=datetime(2000,4,1), frequency=23)
    elapsed = time.time() - t
    print 'Time for loop within BFAST: ', elapsed


