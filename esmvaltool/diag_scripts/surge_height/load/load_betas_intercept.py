def load_betas_intercept(stat):
	import numpy as np
	global betas, intercept
	if not type(stat) == list:
		stat = [stat]
	#?import data (to remove path to files)?
	PATHin = '/usr/people/ridder/Documents/0_models/scripts/python/0_MAGIC/ESMValTool_scripts/v3/data/'
	filename_int = 'monthly_anom_intercept_uv_SLP_gradSLP_'
	filename_bs  = 'monthly_anom_betas_uv_SLP_gradSLP_'
	filename_end = '_075x075_19790101-20000101_1000.npy'
	if len(stat) <= 68:
		betas = {}
		intercept = {}
		for ifname in stat:
			betas[str(ifname)] =  np.load(PATHin + filename_bs + 'allStats' + 	#...
							filename_end).item()[str(ifname)]
			intercept[str(ifname)] = np.load(PATHin + filename_int + 'allStats' + 	#...
							filename_end).item()[str(ifname)]
	else:
		betas     = np.load(PATHin + filename_bs + 'allStats' + filename_end).item()
		intercept = np.load(PATHin + filename_int + 'allStats' + filename_end).item()
