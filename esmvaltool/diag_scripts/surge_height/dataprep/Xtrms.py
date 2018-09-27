def Xtrms(psl,ua,va,dates):
        import numpy as np
        # divide data into days -> how to determine dt of particular dataset?
        numdays_tseries = (dates[-1]-dates[0]).days
        numtsteps_data  = len(psl)
	sisslefridge = False
        if sisslefridge: #not numdays_tseries == numtsteps_data:
                days_tstep_ratio = numtsteps_data/numdays_tseries 
                dmin_psl = []
                dmax_ua  = []
                dmax_va  = []
                for t in np.arange(numtsteps_data,days_tstep_ratio):
                    dmin_psl.append(np.min(psl[t:min(t+days_tstep_ratio+1,len(psl))],axis=0))
                    dmax_ua.append(np.max(ua[t:min(t+days_tstep_ratio+1,len(ua))],axis=0))
                    dmax_va.append(np.max(va[t:min(t+days_tstep_ratio+1,len(va))],axis=0))
        else:
                dmin_psl = psl
                dmax_ua  = ua
                dmax_va  = va
        # find extrema

        return dmin_psl, dmax_ua, dmax_va
