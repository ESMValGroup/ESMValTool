import numpy as np
import os
    
def load_betas_intercept(stat, inpath):
    if not type(stat) == list:
        stat = [stat]
    filename_int = 'monthly_anom_intercept_uv_SLP_gradSLP_'
    filename_bs = 'monthly_anom_betas_uv_SLP_gradSLP_'
    filename_end = '_075x075_19790101-20000101_1000.npy'
    if len(stat) <= 68:
        betas = {}
        intercept = {}
        for ifname in stat:
            betas[str(ifname)] = np.load(inpath + '/' + filename_bs +
                                         'allStats' +  #...
                                         filename_end).item()[str(ifname)]
            intercept[str(ifname)] = np.load(inpath + '/' + filename_int +
                                             'allStats' +  #...
                                             filename_end).item()[str(ifname)]
    else:
        betas = np.load(inpath + '/' + filename_bs + 'allStats' +
                        filename_end).item()
        intercept = np.load(inpath + '/' + filename_int + 'allStats' +
                            filename_end).item()

    return betas, intercept
