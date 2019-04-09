
# todo consider correldated noise!!!
# take into account the uncertainties when doing the fitting


from scipy.optimize import curve_fit
import numpy as np


class Models(object):
    def __init__(self, t, x, **kwargs):
        self.t = t
        self.x = x
        self.e = kwargs.get('e', None)
        assert len(self.t) == len(self.x), \
            'ERROR: unequal size of timeseries and data'
        if self.e is not None:
            assert len(self.x) == len(self.e), \
                'ERROR: data and uncertainties need to have same size'

    def fit(self):
        self.param, self.paramcov = curve_fit(self.func, self.t, self.x)

    def _eval_subset(self, t):
        """
        evaluate function for a subset in time
        requires that these have been set previously
        """
        return self.eval_func(t[self._i1:self._i2])


class LinearTrend(Models):
    def __init__(self, t, x, e=None):
        super(LinearTrend, self).__init__(t, x, e=e)
        self.numParam = 2

    def func(self, t, slope, intercept):
        # print 't: ', t
        # print 'slope: ', slope
        # print 'intercept: ', intercept
        # print 'res: ', type(t*slope + intercept)
        # print type(t)
        return t*slope + intercept

    def eval_func(self, t):
        return self.func(t, self.param[0], self.param[1])


class SineSeason(Models):
    def __init__(self, t, x, **kwargs):
        super(SineSeason, self).__init__(t, x, **kwargs)
        self.f = kwargs.get('f', None)
        assert self.f is not None, 'Frequency not provided!'

    def _sine(self, a, p0, k):
        # sine curve description like in Eq. 4
        p = 2.*np.pi*self.t*k/self.f
        return a*np.sin(p+p0)


class SineSeason1(SineSeason):
    def __init__(self, t, x, **kwargs):
        super(SineSeason1, self).__init__(t, x, **kwargs)
        self.numParam = 2

    def func(self, t, a1, p0_1):
        self.t = t
        return self._sine(a1, p0_1, 1.)

    def eval_func(self, t):
        return self.func(t, self.param[0], self.param[1])


class SineSeasonk(SineSeason):
    def __init__(self, t, x, **kwargs):
        super(SineSeasonk, self).__init__(t, x, **kwargs)
        self.numParam = 3

    def func(self, t, a1, p0_1, k):
        self.t = t
        return self._sine(a1, p0_1, k)

    def eval_func(self, t):
        return self.func(t, self.param[0], self.param[1], self.param[2])


class SineSeason3(SineSeason):
    def __init__(self, t, x, **kwargs):
        """
        References
        ----------
        [1] Verbesselt et al. (2010): doi: 10.1016/j.rse.2010.08.003
        """
        super(SineSeason3, self).__init__(t, x, **kwargs)
        self.numParam = 6

    def func(self, t, a1, a2, a3, p0_1, p0_2, p0_3):
        # Eq 4 in [1]
        self.t = t
        return self._sine(a1, p0_1, 1.) + \
            self._sine(a2, p0_2, 2.) + \
            self._sine(a3, p0_3, 3.)

    def eval_func(self, t):
        return self.func(t, self.param[0], self.param[1], self.param[2],
                         self.param[3], self.param[4], self.param[5])
