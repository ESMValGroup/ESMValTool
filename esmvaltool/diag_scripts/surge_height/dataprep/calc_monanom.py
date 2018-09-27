def calc_monanom(psl, ua, va):
	import numpy as np
        # ------------------------------
        # Calculate monthly climatology
        # ------------------------------
	#climatology = ds.groupby('time.month').mean('time')
	# PROBLEM: This only caluclates the climatology for the input file not the full model run!!!
	climatology_va  = va.groupby('time.month').mean('time')
	climatology_ua  = ua.groupby('time.month').mean('time')
	climatology_psl = psl.groupby('time.month').mean('time')

        # ----------------------------
	# Determine monthly anomalies
        # ----------------------------
	#anomalies = ds.groupby('time.month') - climatology
	apsl_array = xr.apply_ufunc(
			lambda x, m: x - m,
			psl.groupby('time.day'),climatology_psl)
	aua_array  = xr.apply_ufunc(
			lambda x, m: x - m,
			ua.groupby('time.day'),climatology_ua)
	ava_array = xr.apply_ufunc(
			lambda x, m: x - m,
			va.groupby('time.day'),climatology_va)

        # ----------------------------
	# Convert into list object
        # ----------------------------
	apsl = apsl_array.psl.values
	aua  = aua_array.ua.values
	ava  = ava_array.va.values


        return apsl, aua, ava
