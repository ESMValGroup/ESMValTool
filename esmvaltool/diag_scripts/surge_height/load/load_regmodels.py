def load_regmodels(stat):
	import pickle
	EOFinPATH  = '/usr/people/ridder/nobackup_2/python_data/python_databases/0_MAGIC/solver/'
	#EOFinPATH = '/usr/people/ridder/Documents/0_models/scripts/python/0_MAGIC/ESMValTool_scripts/v3/data/'
	fname     = 'monthly_anom_solver_uv_SLP_gradSLP_'
	fname_end = '_075x075_19790101-20000101_1000.pkl'
	global model
	#
	if not type(stat) == list:
		stat = [stat]
	model = {}
	for ifname in stat:
		with open(EOFinPATH + fname + ifname + fname_end,'rb') as input:
			model[i] = pickle.load(input)

