def build_predictX(dates,pseudo_pcs_SLP, pseudo_pcs_gradlatSLP, pseudo_pcs_gradlonSLP, pseudo_pcs_u, pseudo_pcs_v):
	m1day_wgt = 0.5
	m2day_wgt = 0.25
	xnpcs     = 1000 
	#
	X = []
	for t in range(len(dates)):
		X.append([])
	#
	# I. same day
	for t in range(len(dates)):
		X[t].extend([x for x in pseudo_pcs_SLP[t].tolist()[:xnpcs]])
		X[t].extend([x for x in pseudo_pcs_gradlatSLP[t].tolist()[:xnpcs]])
		X[t].extend([x for x in pseudo_pcs_gradlonSLP[t].tolist()[:xnpcs]])
		X[t].extend([x for x in pseudo_pcs_u[t].tolist()[:xnpcs]])
		X[t].extend([x for x in pseudo_pcs_v[t].tolist()[:xnpcs]])
	#
	# II. previous day
	X[0].extend([m1day_wgt*x for x in pseudo_pcs_SLP[0].tolist()[:xnpcs]])
	X[0].extend([m1day_wgt*x for x in pseudo_pcs_gradlatSLP[0].tolist()[:xnpcs]])
	X[0].extend([m1day_wgt*x for x in pseudo_pcs_gradlonSLP[0].tolist()[:xnpcs]])
	X[0].extend([m1day_wgt*x for x in pseudo_pcs_u[0].tolist()[:xnpcs]])
	X[0].extend([m1day_wgt*x for x in pseudo_pcs_v[0].tolist()[:xnpcs]])
	#
	for t in range(0,len(dates)-1):
		X[t+1].extend([m1day_wgt*x for x in pseudo_pcs_SLP[t].tolist()[:xnpcs]])
		X[t+1].extend([m1day_wgt*x for x in pseudo_pcs_gradlatSLP[t].tolist()[:xnpcs]])
		X[t+1].extend([m1day_wgt*x for x in pseudo_pcs_gradlonSLP[t].tolist()[:xnpcs]])
		X[t+1].extend([m1day_wgt*x for x in pseudo_pcs_u[t].tolist()[:xnpcs]])
		X[t+1].extend([m1day_wgt*x for x in pseudo_pcs_v[t].tolist()[:xnpcs]])
	#
#	# II. 2 days previously
#	X[0].extend([m2day_wgt*x for x in pseudo_pcs_SLP[0].tolist()[:xnpcs]])
#	X[0].extend([m2day_wgt*x for x in pseudo_pcs_gradlatSLP[0].tolist()[:xnpcs]])
#	X[0].extend([m2day_wgt*x for x in pseudo_pcs_gradlonSLP[0].tolist()[:xnpcs]])
#	X[0].extend([m2day_wgt*x for x in pseudo_pcs_u[0].tolist()[:xnpcs]])
#	X[0].extend([m2day_wgt*x for x in pseudo_pcs_v[0].tolist()[:xnpcs]])
#	#
#	X[1].extend([m2day_wgt*x for x in pseudo_pcs_SLP[1].tolist()[:xnpcs]])
#	X[1].extend([m2day_wgt*x for x in pseudo_pcs_gradlatSLP[1].tolist()[:xnpcs]])
#	X[1].extend([m2day_wgt*x for x in pseudo_pcs_gradlonSLP[1].tolist()[:xnpcs]])
#	X[1].extend([m2day_wgt*x for x in pseudo_pcs_u[1].tolist()[:xnpcs]])
#	X[1].extend([m2day_wgt*x for x in pseudo_pcs_v[1].tolist()[:xnpcs]])
#	#
#	for t in range(0,len(dates)-2):
#		X[t+2].extend([m2day_wgt*x for x in pseudo_pcs_SLP[t].tolist()[:xnpcs]])
#		X[t+2].extend([m2day_wgt*x for x in pseudo_pcs_gradlatSLP[t].tolist()[:xnpcs]])
#		X[t+2].extend([m2day_wgt*x for x in pseudo_pcs_gradlonSLP[t].tolist()[:xnpcs]])
#		X[t+2].extend([m2day_wgt*x for x in pseudo_pcs_u[t].tolist()][:xnpcs])
#		X[t+2].extend([m2day_wgt*x for x in pseudo_pcs_v[t].tolist()][:xnpcs])

	return X
