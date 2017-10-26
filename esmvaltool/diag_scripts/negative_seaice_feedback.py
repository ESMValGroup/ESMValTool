#!/usr/bin/python
# Author:    Francois Massonnet
# Date:      October 2017
# Purpose:   Set of functions to evaluate the negative ice growth-
#            ice thickness feedback
# Reference: Massonnet et al., in prep.
# Other remarks: The codes presented here is a copy of
#                TECLIM's GitHub code developed by F. Massonnet:
#                http://www.climate.be:3000/TECLIM/ClimateData.git
#                branch develop-fmasson
#
# Structure of the script
#
# 1. A function to evaluate sea-ice volume from spatial fields of
#    sea ice thickness ("compute_volume")
# 2. A function to detrend time series ("detrend")
# 3. A function to evaluate the negative sea ice feedback and its
#    statistics ("negative_seaice_feedback")
# 4. An example of use

# Module imports
from netCDF4 import Dataset
import numpy as np
import sys
import scipy.stats
import esmvaltool.interface_scripts.preprocess
import iris
import logging
logger = logging.getLogger(__name__)


# 1. Function to evaluate sea-ice volume
def compute_volume(avgthickness, cellarea, mask=1):
    """ Input: - avgthickness: sea ice or snow volume per unit cell area, in meters
                               numpy array. if 3-D, time is assumed to be 1st
               - cellarea: array of grid cell areas (sq. meters)
               - mask (1 on ocean, 0 on continent)

        Output: Sea ice or snow volume in the region defined by the mask
    """

    import sys
    import numpy as np

    if np.max(mask) != 1.0 or np.min(mask) < 0.0:
        sys.exit("(compute_volume): mask not between 0 and 1")

    if np.max(avgthickness) > 20.0:
        logger.warning("(compute_volume): large sea ice thickness")
        logger.warning("(compute_volume): np.max(avgthickness) = " + str(np.max(avgthickness)))

    if len(avgthickness.shape) == 3:
        nt, ny, nx = avgthickness.shape
        vol = np.asarray([np.sum(avgthickness[jt, :, :] * cellarea * mask) / 1e12 for jt in range(nt)])
    elif len(avgthickness.shape) == 2:
        vol = np.sum(avgthickness * cellarea * mask) / 1e12
    else:
        sys.exit("(compute_volume): avgthickness has not 2 nor 3 dimensions")

    return vol


# 2. Function to detrend time series
def detrend(data, order=1, period=None):
    import numpy as np
    import sys

    """ Input: data: 1-D numpy array of size N, assumed to be sampled at
                     evenly spaced times
               order: order of the polynomial for detrending
               period: possible existing periodicity of the signal, coming e.g.
                       from external forcing,  expressed in units time steps. 
                       That is, data[i] and data[i + period] correspond
                       to two realizations of the process at times where the
                       forcing might be similar. Common examples include
                       the seasonal cycle forcing, or the diurnal forcing.
   
                       If "period" is not None, the detrending is performed
                       separately for each time step (e.g., all 1st of January,
                       all 2nd of January, ..., all 31st of December in case
                       of annual cycle.
  
                       If "period" is None, the detrending is performed on
                       the given time series
                                            
        Output: the signal detrended using a least-square polynomial
                regression of order "order"
  
    """

    if len(data.shape) != 1:
        sys.exit("(detrend): non-conform input data")

    # Remove possible nans from the data. All the regression (ie polyfit)
    # parameters will be estimated based on the no-nan data but the
    # residuals will be computed from the original data in order
    # to keep the same size and to restitute NaNs where they appeared

    data_nonan = data[~np.isnan(data)]

    N = len(data)
    N_nonan = len(data_nonan)

    # If the signal has no periodicity, we just make a linear regression
    if period is None:
        time_nonan = np.arange(N_nonan)
        time = np.arange(N)
        p = np.polyfit(time_nonan, data_nonan, order)
        residuals = data - np.sum([p[i] * time ** (order - i) \
                                   for i in range(order + 1)], axis=0)

    # If the signal contains a periodical component, we do the regression
    # time step per time step
    else:
        residuals = np.empty([N])

        # For each time step of the period, detrend
        for jP in np.arange(period):
            raw = data[np.arange(jP, N, period)]
            raw_nonan = raw[~np.isnan(raw)]
            time = np.arange(len(raw))
            time_nonan = np.arange(len(raw_nonan))
            p = np.polyfit(time_nonan, raw_nonan, order)
            residuals[np.arange(jP, N, period)] = \
                raw - np.sum([p[i] * time ** (order - i) \
                              for i in range(order + 1)], axis=0)

            # Note that another common option is to first remove a seasonal cycle
            # and then detrend the anomalies. However this assumes that a cycle can
            # be estimated, which in presence of a trend is tricky because the
            # trend component interferes with the mean. I have tried that and it
            # gives ugly step-wise anomalies. Detrending day per day seems the
            # the most natural way to do, at least as long as we assume that the
            # raw signal at some time is the result of a seasonal cycle depending
            # on the position of the time step in the period, plus a common trend,
            # plus some noise.

    return residuals


# 3. Function to estimate the negative feedback
def negative_seaice_feedback(volume, period, order=1):
    """ Function to estimate the negative ice-thickness ice growth feedback
        and its significance.

        INPUTS
          (1) volume = 1-D numpy array containing time series of sea ice volume
          (2) period = period of the signal (period expressed
                       in time steps: 12 for monthly, 365 for daily...)
          (3) order  = order of the polynomial detrending (integer >=0)

        OUTPUTS
          (1) Feedback parameter expressed as the regression
              between dV on V_min (see (3))
          (2) Correlation between those two and the p-value under the null
              hypothesis of no correlation between dV and V_min
          (3) [V_min, dV]: detrended time series of annual minimum of sea ice
              volume, detrended series of wintertime volume production
    """

    if len(volume.shape) != 1:
        sys.exit("(negative_seaice_feedback) volume is not 1-D")

    nt = len(volume)
    if nt // period != 1.0 * nt / period:
        sys.exit("(negative_seaice_feedback) length of volume series is \
              not multiple of period")

    # 1. Locate the minima for each year
    imin = [t + np.nanargmin(volume[t:t + period]) \
            for t in np.arange(0, nt, period)]

    # 2. Locate the maxima for each year
    imax = [t + np.nanargmax(volume[t:t + period]) \
            for t in np.arange(0, nt, period)]

    # 3. Detrend series. A one-year shift is introduced to make sure we
    #    compute volume production *after* the summer minimum
    Vmin = detrend(volume[imin[:-1]], order=order)
    dV = detrend(volume[imax[1:]] - volume[imin[:-1]], order=order)

    # 4. Compute diagnostics
    # If all Vmins are zero or all dVs are zero, return Nan (pathological case)
    if np.max(Vmin) == 0.0 or np.max(dV == 0.0):
        nf = np.nan
        r = np.nan
        pval = np.nan
        sd = np.nan
    else:
        r = np.corrcoef(Vmin, dV)[0, 1]
        N = len(Vmin)
        tstat = r / np.sqrt((1 - r ** 2) / (N - 2))  # The t-statistic.
        # Under the null hypothesis
        # of no correlation,
        # tstat follows a student's
        # law with  N - 2 dof.
        pval = 1.0 - scipy.stats.t.cdf(np.abs(tstat), N - 2)

        if pval > 0.05:
            logger.warning("(negative_seaice_feedback) Check the scatterplot of dV versus V_min, it is most likely "
                           "suspicious, and the feedback  factor likely meaningless: p-value: " + str(pval))

        try:
            fit, cov = np.polyfit(Vmin, dV, 1, cov=True)
            nf = fit[0]  # Fit parameter
            sd = np.sqrt(cov[0, 0])  # Standard deviation on it
        except ValueError:
            logger.error("(negative_seaice_feedback) PROBLEM, series badly conditioned: "
                         "Input volume: {0} Vmin: {1} dv: {2}".format(volume, Vmin, dV))
            raise
    # 5. Return
    return [nf, [r, pval, sd], [Vmin, dV]]


def main(project_info):
    var = project_info['RUNTIME']['currDiag'].variables[0]

    # Open data
    # Root directory to the data

    for model_info in project_info['ALLMODELS']:
        sit_path = esmvaltool.interface_scripts.preprocess.get_cf_fullpath(project_info, model_info, var)

        # Load cell area
        filearea = "/esnas/autosubmit/con_files/mesh_mask_nemo.Ec2.3_O1L42.nc"
        f = Dataset(filearea, mode="r")
        e1t = f.variables["e1t"][:]
        e2t = f.variables["e2t"][:]
        cellarea = e1t[0, :, :] * e2t[0, :, :]
        f.close()

        sit = iris.load_cube(sit_path)
        # Compute integrated volume over domain and append to existing time series
        volume = compute_volume(sit.data, cellarea, mask=1.0 * (sit.coord('latitude').points > 80.0))

        nf, stats, _ = negative_seaice_feedback(volume, period=12, order=2)

        logger.info("Negative feedback: ".ljust(20) + str(nf))
        logger.info("P-Value: ".ljust(20) + str(stats[1]))
