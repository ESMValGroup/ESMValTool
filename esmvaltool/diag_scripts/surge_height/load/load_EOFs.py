def load_EOFs():
	import pickle
	EOFinPATH = './esmvaltool/diag_scripts/surge_height/data/'
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

