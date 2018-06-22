def load_EOFs():
	import pickle
	#EOFinPATH  = '/usr/people/ridder/nobackup_2/python_data/python_databases/0_MAGIC/EOFs/final/'
	#EOFinPATH = '/usr/people/ridder/Documents/0_models/scripts/python/0_MAGIC/ESMValTool_scripts/v3/data/'
	EOFinPATH = '/usr/people/ridder/Documents/0_models/ESMValTool/diag_scripts/lib/python/surge_estimator/data/'
	global SLPsolver, gradlonsolver, gradlatsolver, usolver, vsolver
	#
	with open(EOFinPATH + 'monthly_anom_EOF_SLPsolver_075x075_19790101-20000101.pkl','rb') as input:
		SLPsolver = pickle.load(input)
	#
	with open(EOFinPATH + 'monthly_anom_EOF_gradlatsolver_075x075_19790101-20000101.pkl','rb') as input:
		gradlatsolver = pickle.load(input)
	#
	with open(EOFinPATH + 'monthly_anom_EOF_gradlonsolver_075x075_19790101-20000101.pkl','rb') as input:
		gradlonsolver = pickle.load(input)
	#
	with open(EOFinPATH + 'monthly_anom_EOF_usolver_075x075_19790101-20000101.pkl','rb') as input:
		usolver =  pickle.load(input)
	#
	with open(EOFinPATH + 'monthly_anom_EOF_vsolver_075x075_19790101-20000101.pkl','rb') as input:
		vsolver =  pickle.load(input)

