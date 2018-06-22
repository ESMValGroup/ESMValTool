def grad_psl(msl):
	import numpy as np
	global gradlatSLP, gradlonSLP
	gradlatSLP = np.gradient(msl,axis = 1)
	gradlonSLP = np.gradient(msl,axis = 2)
