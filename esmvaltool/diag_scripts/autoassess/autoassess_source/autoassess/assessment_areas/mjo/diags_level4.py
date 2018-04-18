import matplotlib.pyplot as plt
import numpy as np
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import os

import mjo_utils as mu

def diagnos_level4(runid, out_dir):
    print 'Starting Level 4 diagnostics...'
    #
    metrics = dict()
    undef = np.nan



    # TODO move to level 3, requires Obs data to be processed
    # Compute pattern correlation between OBS lead-lag composites
    varnames = ['OLR', 'U850', 'U200']
    for varname in varnames:
        xdir = os.path.dirname(os.path.abspath(__file__))
        obs_lcor_name = os.path.join(xdir, "OBS_" + varname + "_LeadLagCorr.nc")
        model_lcor_name = os.path.join(out_dir, runid + "_" + varname + "_LeadLagCorr.nc")
        if os.path.isfile(obs_lcor_name):
            obs_cube = iris.load_cube(obs_lcor_name)
            model_cube = iris.load_cube(model_lcor_name)
            nlag, nlon = obs_cube.shape
            rcorr = np.corrcoef(np.reshape(obs_cube.data, nlag*nlon), np.reshape(model_cube.data, nlag*nlon))
            metrics['Lead-lag pattern corr. ' + varname] = round(rcorr[0,1], 2)
        else:
            print '%s does not exist.' % obs_lcor_name
            metrics['Lead-lag pattern corr. ' + varname] = undef

    for x in sorted(metrics.keys()):
        print x + ' = ' + str(metrics[x])

    print '*' * 50
    print 'Level 4 diagnostics completed.'
    print '*' * 50
    return metrics

