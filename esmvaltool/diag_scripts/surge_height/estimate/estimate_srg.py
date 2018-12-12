import os

import numpy as np
import pandas as pd

#from .load import load_betas_intercept as llbi
from load.load_monmean_srgclim import load_monmean_srgclim


def estimate_srg(X, dates, stat, betas, intercept, data_dir):
    if not type(stat) == list:
        stat = [stat]
    #
    #global srg_est_full, srg_est
    srg_est = {}
    srg_est_full = {}
    monanom_srg = {}
    for s in stat:
        srg_est[s] = np.zeros(len(dates)).tolist()
        for t in range(len(dates)):
            srg_est[s][t] = sum(X[s][t] * betas[s][:len(X[s][t])]) + float(
                intercept[s])
        # Add seasonal cycle back to surge
        srg_dir = os.path.join(data_dir + '/srgclim/')
        monanom_srg[s] = load_monmean_srgclim(s, srg_dir)
        srg_est_full[s] = []
        for t in range(len(dates)):
            srg_est_t = srg_est[s][t] + monanom_srg[s][dates[t].month - 1]
            srg_est_full[s].append(srg_est_t)

    return srg_est_full
