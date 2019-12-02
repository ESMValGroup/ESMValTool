import numpy as np
from scipy.stats import pearsonr
from statsmodels.tsa.stattools import acf
from scipy.stats import norm
from scipy.stats.mstats import linregress
from scipy.signal import detrend

def standardize1d(x):
    return (x-np.nanmean(x))/np.nanstd(x)

def pearsonr1d(a,b):
    commonmask = (np.isfinite(a) & np.isfinite(b))
    a,b = a[commonmask],b[commonmask]
    if a.size < 2:
        return np.nan, np.nan
    else:
        return pearsonr(a,b)

def rmsd1d(a,b):
    result = np.sqrt(np.nanmean((a-b)**2))
    return result

def absdiffaxismean1d(a,b):
    """Absolute difference of along-axis mean values"""
    return np.nanmean(a)-np.nanmean(b)

def reldiffaxismean1d(a,b):
    """Relative difference of along-axis mean values"""
    return ((np.nanmean(a)-np.nanmean(b))/np.nanmean(b))*100.

def weatherhead1d(inputdata,trend_magnitude=None):
    """This framework follows  Weatherhead et al. 1998 [1] for estimating the amount of years needed to detect a trend of certain magnitude with a probability of 90%.
    Data has to be provided with yearly frequency (e.g. yearly means, means for DJF, etc.)
    Parameters
    ----------
    trend_name : str
        The name of the trend method used for detrending the data
    trend_magnitude : float
        The magnitude of the trend in [data units/decade] in which one would be interested
    Returns
    -------
        A dictionary containing:
            trend_magnitude : same as above
            std_res : the standard deviation of the residuals
            acf_res : the 1-lag autocorrelation
            n_star : the estimated number of years
    References
    -----------
    [1] Weatherhead, Betsy & C. Reinsel, Gregory & C. Tiao, George & Meng, Xiao-Li & Choi, Dongseok & Cheang, Wai-Kwong & Keller, Teddie & DeLuisi, John & Wuebbles, Donald & Kerr, J & J. Miller, Alvin & Oltmans, Samuel. (1998). Factors affecting the detection of trends: Statistical considerations and applications to environmental data. Journal of Geophysical Research. 1031. 17149-17162. 10.1029/98JD00995.
    """
    # Return nan if any nan in array
    if np.any(~np.isfinite(inputdata)):
        return np.nan, np.nan, np.nan
    
    # Calculate according to Weatherhead et al.
    std_res = np.std(inputdata)
    # This is a workaround, related to the following Dask issue: https://github.com/dask/dask/pull/3742
    if len(inputdata)==1:
        acf_res = inputdata
    else:
        acf_res = acf(inputdata, nlags=1)[1]
    n_star = ((3.3 * std_res / np.abs(trend_magnitude)) * (
        (1 + acf_res) / (1 - acf_res))**.5)**(2. / 3.)
    return std_res, acf_res, n_star

def mannkendall1d(x, alpha=0.05):
    """
    This function is derived from code originally posted by Sat Kumar Tomer
    (satkumartomer@gmail.com)
    See also: http://vsp.pnnl.gov/help/Vsample/Design_Trend_Mann_Kendall.htm
    The purpose of the Mann-Kendall (MK) test (Mann 1945, Kendall 1975, Gilbert
    1987) is to statistically assess if there is a monotonic upward or downward
    trend of the variable of interest over time. A monotonic upward (downward)
    trend means that the variable consistently increases (decreases) through
    time, but the trend may or may not be linear. The MK test can be used in
    place of a parametric linear regression analysis, which can be used to test
    if the slope of the estimated linear regression line is different from
    zero. The regression analysis requires that the residuals from the fitted
    regression line be normally distributed; an assumption not required by the
    MK test, that is, the MK test is a non-parametric (distribution-free) test.
    Hirsch, Slack and Smith (1982, page 107) indicate that the MK test is best
    viewed as an exploratory analysis and is most appropriately used to
    identify stations where changes are significant or of large magnitude and
    to quantify these findings.
    #############################################################################
    MIT License
    Copyright (c) 2017 Michael Schramm
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
    #############################################################################    
    Input:
        x:   a vector of data
        alpha: significance level (0.05 default)
    Output:
        trend: tells the trend (increasing, decreasing or no trend)
        h: True (if trend is present) or False (if trend is absence)
        p: p value of the significance test
        z: normalized test statistics
    Examples
    --------
      >>> x = np.random.rand(100)
      >>> trend,h,p,z = mk_test(x,0.05)
    """
    x = x[~np.isnan(x)]
    
    n = len(x)
    
    if n <= 1:
        return np.nan

    # calculate S
    s = 0
    for k in range(n - 1):
        for j in range(k + 1, n):
            s += np.sign(x[j] - x[k])

    # calculate the unique data
    unique_x = np.unique(x)
    g = len(unique_x)

    # calculate the var(s)
    if n == g:  # there is no tie
        var_s = (n * (n - 1) * (2 * n + 5)) / 18
    else:  # there are some ties in data
        tp = np.zeros(unique_x.shape)
        for i in range(len(unique_x)):
            tp[i] = sum(x == unique_x[i])
        var_s = (n * (n - 1) * (2 * n + 5) - np.sum(tp * (tp - 1) *
                                                    (2 * tp + 5))) / 18

    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:  # s == 0:
        z = 0

    # calculate the p_value
    p = 2 * (1 - norm.cdf(abs(z)))  # two tail test
    h = abs(z) > norm.ppf(1 - alpha / 2)

    if (z < 0) and h:
        trend = -1
    elif (z > 0) and h:
        trend = 1
    else:
        trend = 0
    return trend #, h, p, z

def validfrac1d(x):
    return np.isfinite(x).sum()/x.size


def lineartrend1d(y,x=None, alpha=0.05):
    y = np.array(y).flatten()
    if x is None:
        x = np.arange(len(y), dtype=float)
    else:
        x = np.array(x, dtype=float).flatten()
    
    # Do own masking of missing values
    isfinite_mask = np.isfinite(y)
    y = y[isfinite_mask]
    x = x[isfinite_mask]
    # Catching the case of less than two valid points
    if y.size > 1:
        linoutput = linregress(x, y)
        return linoutput.slope, linoutput.pvalue
    else:
        return np.nan, np.nan

def theilslopes1d(y,x=None):
    '''
    Adapted from scipy.stats.theilslopes, leaving out calculation of confidence intervals and allowing for nan values.
    '''
    y = np.array(y).flatten()
    if x is None:
        x = np.arange(len(y), dtype=float)
    else:
        x = np.array(x, dtype=float).flatten()
    if len(x) != len(y):
        raise ValueError("Incompatible lengths ! (%s<>%s)" % (len(y), len(x)))
    # Do own masking of missing values
    isfinite_mask = np.isfinite(y)
    y = y[isfinite_mask]
    x = x[isfinite_mask]

    # Catching the case of less than two valid points
    if y.size > 1:
        # Compute sorted slopes only when deltax > 0
        deltax = x[:, np.newaxis] - x
        deltay = y[:, np.newaxis] - y
        slopes = deltay[deltax > 0] / deltax[deltax > 0]
        medslope = np.nanmedian(slopes)
        return medslope
    else:
        return np.nan
