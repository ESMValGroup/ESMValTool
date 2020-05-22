#################################################################################
# A set of functions for the detection and attribution of climate signals
# Megan Kirchmeier-Young
#
# Many functions adapted from other code as noted.
#
# References:
#	Allen, M. R., and S. F. B. Tett, 1999: Checking for model consistency in optimal fingerprinting. Climate Dyn., 15, 419-434, doi:10.1007/s003820050291.
#	Allen, M. R., and P. A. Stott, 2003: Estimating signal amplitudes in optimal fingerprinting, part I: Theory. Climate Dyn., 21, 477-491, doi:10.1007/s00382-003-0313-9.
#	Ribes, A., S. Planton, and L. Terray, 2013: Application of regularised optimal fingerprinting to attribution. Part I: Method, properties and idealised analysis. Climate Dyn., 41, 2817-2836, doi:10.1007/s00382-013-1735-7.
#################################################################################

def regC(xL, center_flag=1):

	# given the control data xL formed into a nxp matrix, calculates the
	#     regularized covariance matrix cl_hat
	# see Ribes et al. (2009) following Ledoit and Wolf (2004)
	#
	# gives the option to center the input matrix first (default); would set
	#     flag to 0 if input matrix has already been multiplied by the projection
	#     matrix to guarantee full rank of the cov matrix
	
	import numpy as np
	
	if np.size(np.shape(xL))>2:
		raise ValueError('input data has too many dimensions')
	n,p = np.shape(xL)

	#center the input data	
	if center_flag == 1:
		xLc = xL - np.tile(np.mean(xL,axis=0),(n,1))
	else:
		xLc = xL	
	
	#calculate covariance matrix estimate
	#c_hat = np.dot(xLc.T,xLc)/n
	c_hat = np.dot(xLc,xLc.T)/p


		
	v_hat = np.trace(c_hat)/n
	
	y1 = c_hat - v_hat*np.eye(n)
	del_hat2 = np.trace(np.dot(y1.T,y1))/n
	
	y2a = 0
	for i in np.arange(p):
		w = xLc[:,i]
		w = np.expand_dims(w,axis=1)
		y2b = np.dot(w,w.T) - c_hat
		y2c = np.trace(np.dot(y2b.T, y2b))/n	
		y2a += y2c

	y2 = y2a/(p**2)
	
	beta_hat2 = min(del_hat2,y2)

	alpha_hat2 = del_hat2 - beta_hat2


	gamma_hat = alpha_hat2/del_hat2
	
	rho_hat = beta_hat2*v_hat/del_hat2
	
	cl_hat = gamma_hat*c_hat + rho_hat*np.eye(n)

	return cl_hat

#################################################################################
def reduce_dim(x,y,noise1,noise2,ns=1,ind='None',return_flag=0,order='F'):
	
	#reduces the dimension of the input variables so Cov matrix will be full-rank
	#used when the data are centered (because one entry is linear combination 
	#     of others)
	#
	#ns is the length of the spatial dimension (default is 1)
	#assumes order s=1,t=1; s=1,t=2;...; s=2,t=1; s=2,t=2;...
	#   if grouped by time instead, change order to 'C' (only if no missing data)
	#
	#if any spatial dimensions are missing time steps, supply ind with list of
	#number of time steps for each spatial dimensions
	#
	#from Ribes code and ECOF code
	
	import numpy as np
	from scipy import linalg	
	
	n1 = np.size(y)	#assumes y is a single column of data (or row, I guess)
	if np.mod(float(n1),ns)!= 0:
		#raise ValueError('please check your dimensions')
		if ind is 'None':
			raise ValueError('check dimensions or supply indices')

	if ind is 'None':	
		n = int(n1/ns)
		m = np.eye(n) - 1./n
		
		u,s,v1 = linalg.svd(m)
			
		p = u[:,:-1].T

		if ns==1:
			xr = np.dot(p,x)
			yr = np.dot(p,y)
			cn1 = np.dot(p,noise1)
			cn2 = np.dot(p,noise2)
		else:

			yr = np.dot(p,np.reshape(y,(-1,ns),order=order)).flatten(order=order)[:,None]
			xr = np.zeros(((n-1)*ns,len(x[0,:])))
			for j in xrange(len(x[0,:])):
				xr[:,j] = np.dot(p,np.reshape(x[:,j],(-1,ns),order=order)).flatten(order=order)
			cn1 = np.zeros(((n-1)*ns,len(noise1[0,:])))
			for j in xrange(len(noise1[0,:])):
				cn1[:,j] = np.dot(p,np.reshape(noise1[:,j],(-1,ns),order=order)).flatten(order=order)	
			cn2 = np.zeros(((n-1)*ns,len(noise2[0,:])))
			for j in xrange(len(noise2[0,:])):
				cn2[:,j] = np.dot(p,np.reshape(noise2[:,j],(-1,ns),order=order)).flatten(order=order)					

	else:
		
		if len(ind) != ns:
			raise ValueError('please provide nt for each spatial entry')
		if order is 'C':
			raise ValueError('can only handle missing data when grouped by spatial dimension')
		if return_flag!=1:
			raise ValueError('getting the p matrix probably will not be useful here')
		
		yr = np.zeros((sum(ind)-ns))
		xr = np.zeros((sum(ind)-ns,len(x[0,:])))
		cn1 = np.zeros((sum(ind)-ns,len(noise1[0,:])))
		cn2 = np.zeros((sum(ind)-ns,len(noise2[0,:])))
		
		q = 0;w=0
		for i in xrange(ns):
			n = ind[i]
			m = np.eye(n) - 1./n
			u,s,v1 = linalg.svd(m)
			p = u[:,:-1].T

			yr[q:q+n-1,:] = np.dot(p,y[w:w+n])[:,None]
			if len(x[0,:])==1:
				xr[q:q+n-1,:] = np.dot(p,x[w:w+n])[:,None]
			else:
				for j in xrange(len(x[0,:])):
					xr[q:q+n-1,j] = np.dot(p,x[w:w+n,j])
			for j in xrange(len(noise1[0,:])):
				cn1[q:q+n-1,j] = np.dot(p,noise1[w:w+n,j])
			for j in xrange(len(noise2[0,:])):
				cn2[q:q+n-1,j] = np.dot(p,noise2[w:w+n,j])
			w=w+n;q=q+n-1	

	
	if return_flag==0:
		return (xr,yr,cn1,cn2)
	else:
		return p
 
#################################################################################
def rct_mc(C, x0, ne, n1, n2, nmc=10000, flag_plotden=0):
	
	# following Ribes et al. (2013) and accompanying code, uses Monte Carlo
	# 	approach for the residual consistency test
	#
	# input C is covariance matrix, x0 is response matrix, nmc is the number
	# 	of Monte Carlo simulations, n1 and n2 are the dimensions of the
	#    climate noise matrices (from the control runs, probably), and ne is
	#    the number of ensemble members that went into each ensemble mean
	#
	# returns set of nmc estimates of the test statistic for the residual 
	# 	consistency test
	
	import numpy as np
	from scipy import linalg
	
	n,m = np.shape(x0)
	
	if np.size(np.shape(x0))>2:
		raise ValueError('input x0 has too many dimensions')
	if np.size(np.shape(C))>2:
		raise ValueError('input C has too many dimensions')
	if np.shape(C)[0]!=n or np.shape(C)[1]!=n:
		raise ValueError('dimensions of C must be consistent with dimensions of x0')

	#if ne is a single value, assume all signals have same ensemble size 
	#   and expand ne to size m
	if np.size(ne)==1:
		ne = np.repeat(float(ne),m)
	else:
		ne = np.array(ne,dtype='float')	
		
	beta0 = np.ones((m,1))

	#c12 = np.real(linalg.inv(linalg.sqrtm(C)))
	c12 = np.real(linalg.sqrtm(C))
	
	r2s = np.zeros((nmc,1))
	
	for i in np.arange(nmc):
		#simulate obs and response
		y =  np.dot(x0,beta0) + np.dot(c12,np.random.randn(n,1))
		x = x0 + np.dot(c12,np.random.randn(n,m))/(np.ones((n,1))*np.sqrt(ne))
		x1 = x * (np.ones((n,1))*np.sqrt(ne))
		
		#simulate climate noise and calculate the reg cov matrix
		cn1 = np.dot(c12,np.random.randn(n,n1))
		cn2 = np.dot(c12,np.random.randn(n,n2))
		c1reg = regC(cn1)
		c1reg12 = np.real(linalg.inv(linalg.sqrtm(c1reg)))

		#prewhiten with the reg cov matrix
		xc = np.dot(c1reg12,x1)
		yc = np.dot(c1reg12,y)
		cn2c = np.dot(c1reg12,cn2)
		
		#calculate test statistic for residual consistency test, using TLS
		z = np.hstack((xc,yc))
		u,s,v1 = linalg.svd(z)
		#v = v1.T
		
		lam = s**2
		um = np.matrix(u[:,m]).T
		r2s[i] = lam[-1]/(np.dot(np.dot(np.dot(um.T,cn2c),cn2c.T),um)/n2)
		
	if flag_plotden==1:
		import matplotlib.pyplot as plt
		plt.figure		
		plt.hist(r2s, bins=30)
		plt.show()
		
	return r2s
		
#################################################################################

def tls(x,y,cn1,cn2='None',ne=1,rof_flag=1,CI_flag=1,alpha=0.10,flag_2S=0,flag_3S=0,flag_4S=0,error_flag=0,RCT_flag=True):
	
	# returns the coefficients from a total least squares regression
	# expects x and y to have column data (and should be centered)
	#
	# x is a nxm matrix of response data where n is time and m is the signal
	# y is a nx1 array of the corresponding observations
	# ne is the number of members that went into the calculation of the ensemble
	#     mean; should be of length m; if all m signals have same size ensemble
	#     a single input is allowed
	# cn11 and cn22 are sections of the data (probably from a control run) 
	#     used for calculating the covariance matrix; cn=climate noise
	#     cn11 is used for calculating the regularized covariance matrix and
	#     in the the residual consistency test Monte Carlo simulations
	#	 if cn2 input, will be used for the residual consistency test
	#     --does not need to be regularized
	# rof_flag is a flag indicating whether to use the ROF approach (1)
	# CI_flag is a flag indicating whether the confidence intervals of beta 
	#     should be calculated; default is yes (1)
	# alpha is the significance level used to construct the CI on beta
	# flag_2S is a flag indicating whether to adjust the output for a combination
	#     of forcing signals (e.g., if 1 with inputs ALL+NAT, switch to ANT+NAT)
	# error_flag determines if errors on x and y will be calculated/output
	#
	# References:
	#      Ribes et al. (2013), Allen and Stott (2003), Van Huffel and 
	#      Vandewalle (1991)
	#      And code from Ribes (in scilab) and Feng (in R)
	
	
	import numpy as np
	from scipy import linalg, stats
	import warnings

	if np.shape(x)[0] != np.shape(y)[0]:
		raise ValueError('x and y must have same first dimension')


	n,m = np.shape(x)	#time dimension, number of signals
	
	if np.shape(cn1)[0]!=n:
		raise ValueError('dimension mismatch: x and cn1')
	p1 = np.shape(cn1)[1]
	if cn2 is not 'None':
		if np.shape(cn2)[0]!=n:
			raise ValueError('dimension mismatch: x and cn2')
		p2 = np.shape(cn2)[1]
	
	#if ne is a single value, assume all signals have same ensemble size 
	#   and expand ne to size m
	if np.size(ne)==1:
		ne = np.repeat(float(ne),m)
	else:
		ne = np.array(ne,dtype='float')	

	#scale variance or something
	if m==1:
		x1 = x*np.sqrt(ne)
	else:
		x1 = np.dot(x,np.diag(np.sqrt(ne)))

	creg = regC(cn1,center_flag=0)
	
	#if scaling by creg
	if rof_flag==1:
		cl12 = np.real(linalg.inv(linalg.sqrtm(creg)))
		x1 = np.dot(cl12,x1)
		y = np.dot(cl12,y)
	
	z = np.hstack((x1,y))		#combine x and y into one matrix
	u,s,v1 = linalg.svd(z)
	v = v1.T
	
	#method for computation of TLS coefs from Van Huffel and Vandewalle (1991)
	beta = np.dot(linalg.inv(np.dot(x1.T,x1)-(s[m]**2)*np.identity(m)),np.dot(x1.T,y))
	#scale results by sqrt(ensemble size) because other people do this	
	beta = np.array(beta)*np.expand_dims(np.sqrt(ne),axis=1)
	
	###########################
	#residual consistency test#
	###########################
	if RCT_flag:
		lam = s**2
		if rof_flag==0:
			w2 = cn2
		else:
			w2 = np.dot(cl12,cn2)
			
		um = np.expand_dims(u[:,m], axis=1)	  #want (m+1)th column of u
		r2 = lam[-1]/(np.dot(np.dot(np.dot(um.T,w2),w2.T),um)/p2)	#from AS03
		print (r2)

		#cna = np.hstack((cn1,cn2))
		#chat = np.dot(cna,cna.T)/p1
		#r2s = rct_mc(chat,x,ne=ne,n1=p1,n2=p2)
		r2s = rct_mc(creg,x,ne=ne,n1=p1,n2=p2)
		nr = np.size(r2s)
		h = 1.06*np.std(r2s,ddof=1)*nr**(-1./5)	#following Ribes2013 (see function gke)
		pvi = stats.norm.sf(r2*np.ones((nr,1)),loc=r2s,scale=h*np.ones((nr,1)))
		rc_pvalue = np.sum(pvi)/nr
		rc_mcCI = np.percentile(r2s,[5,95])	#non-parametric CI from MC results
		rc_mcCI = np.append(rc_mcCI,r2)	#also save test stat for comparison
	
	###############
	#calculate CIs#
	###############
	if CI_flag==1:
		if rof_flag==0: #NPG - copied here, so code will run even if RCT_flag not set.
			w2 = cn2
		else:
			w2 = np.dot(cl12,cn2)
		betaCI = np.empty((m,2))

		lam2_hat = np.empty((m+1,1))
		lam = s**2
		for j in np.arange(m+1):
			um = np.matrix(u[:,j]).T
			lam2_hat[j] = lam[j]/(np.dot(np.dot(um.T,np.dot(w2,w2.T)),um)/p2)
			
		#use m-sphere or whatever, as in AS03 and Ribes code
		npt = 10000 #Original version. Replaced with 2000 to speed up code for testing.
#		npt = 2000 

		if m==1:
			pts = np.array([1,-1])
			pts = np.expand_dims(pts,axis=1)
		else:
			pts1 = np.random.randn(npt,m)
			pts = pts1 / (np.expand_dims(np.sqrt(np.sum(pts1**2,axis=1)),axis=1)*np.ones((1,m)))
			
		dlam2 = lam2_hat - np.min(lam2_hat)
		critf = (stats.f.ppf(.9,1,p2))	#because Ribes code
		a = np.sqrt(critf)*pts
		if np.sum(dlam2[:-1])!=0:
			b_m1 = a / np.dot(np.ones((np.size(pts[:,0]),1)),np.sqrt(dlam2[:-1]).T)
			b_m2 = np.lib.scimath.sqrt(1-np.sum(b_m1**2,axis=1)) #this sqrt function returns complex values if input is (-)
			
			b_m2 = np.expand_dims(b_m2,axis=1)	#put the needed dimension back in; no idea why this pos removed it in the first place
			cont = 1 	#continue flag
		else:
			betaCI[:] = np.nan
			cont = 0
		if cont==1:
			if np.any(b_m2<=0):		#must be strictly positive
				betaCI[:] = np.nan
			else:
				vpts = np.dot(np.hstack((b_m1, b_m2)),v.T)
				vpts2 = np.dot(vpts,np.diag(np.sqrt(np.append(ne,1))))
				if flag_2S==1:
					vpts2[:,1] = vpts2[:,0] + vpts2[:,1]	#if using P to switch to ANT+NAT, this gives the appropriate CI for this NAT (from Reza)
				if flag_3S==1:
					vpts2[:,1] = vpts2[:,0] + vpts2[:,1]	#updated to work for ALL+NAT+ANT1 to ANT2+NAT+ANT1)
					vpts2[:,2] = vpts2[:,0] + vpts2[:,2]
				if flag_4S==1:
					vpts2[:,1] = vpts2[:,0] + vpts2[:,1]	#updated to work for ALL+NAT+ANT1+ANT2 to ANT3+NAT+ANT1+ANT2)
					vpts2[:,2] = vpts2[:,0] + vpts2[:,2]
					vpts2[:,3] = vpts2[:,0] + vpts2[:,3]
				for i in np.arange(m):
					vc_2d_pts = vpts2[:,i]+1j*vpts2[:,m]	#1j gives complex i because python
					vs = v1[:,-1]					
					vc_2d_ref = vs[i]*np.sqrt(ne[i])+1j*vs[m]
					vprod_2d = vc_2d_pts/vc_2d_ref
					arg = np.sort(np.imag(np.log(vprod_2d)))
					delta_arg_min = arg[0]
					delta_arg_max = arg[-1]
					delta_max_1 = np.amax(arg[1:]-arg[:-1])
					delta_max = np.amax([delta_max_1, arg[0]-arg[-1]+2*np.pi])
					if delta_max < np.pi:
						betaCI[i,:] = np.nan
					else:
						if  np.argmax([delta_max_1, arg[0]-arg[-1]+2*np.pi])!=1:
							warnings.warn('Apparently the second entry here should be larger...')
						arg_ref = np.imag(np.log(vc_2d_ref))
						arg_min = delta_arg_min + arg_ref
						arg_max = delta_arg_max + arg_ref
						betaCI[i,0] = -1/np.tan(arg_min)
						betaCI[i,1] = -1/np.tan(arg_max)

	if error_flag==1:
		um = np.expand_dims(u[:,m], axis=1)
		vm = np.expand_dims(v[:,m], axis=1)
		err = s[m]*np.dot(um,vm.T)

						
	if flag_2S==1:	#switch from ALL+NAT to ANT+NAT
		p = np.array([[1,0],[1,1]])
		beta = np.dot(p,beta)
	
	if flag_3S==1:
		p = np.array([[1,0,0],[1,1,0],[1,0,1]])
		beta = np.dot(p,beta)
		
	if flag_4S==1:
		p = np.array([[1,0,0,0],[1,1,0,0],[1,0,1,0],[1,0,0,1]])
		beta = np.dot(p,beta)
	
	#function outputs (depends on inputs)
	out = {}
	out['beta'] = beta
	
	if RCT_flag==1:
		out['rc_pvalue'] = rc_pvalue
		out['rc_mcCI'] = rc_mcCI
	if CI_flag==1:
		out['betaCI'] = betaCI
	if error_flag==1:
		out['err'] = err

	return out
	
#################################################################################
def adjust_beta(xa,y,cn1,xn='None',weights='None'):
	
	#returns the time series in columns of xa (ALL) and xn (NAT) that have 
	#    been adjusted using the scaling factors
	#the inputs xa, xn, and y should be for the same years
	#assumes xa and xn have same number of ensemble members
	
	import numpy as np
	from scipy import linalg
	
	#add checks for appropriate size of inputs	
	
	n,ne = np.shape(xa)	
	
	#calculate ensemble mean
	if weights is 'None':
		ema = np.mean(xa,axis=1)
		ne2 = ne
		if xn is not 'None':		
			emn = np.mean(xn,axis=1)
	else:
		ema = np.squeeze(np.sum(xa*np.tile(weights[None,:], (np.shape(xa)[0],1)),axis=1))
		ne2 = 1/np.sum(weights**2)
		if xn is not 'None':	
			emn = np.squeeze(np.sum(xn*np.tile(weights[None,:], (np.shape(xn)[0],1)),axis=1))
	
	#calculate anomalies
	mt = np.mean(ema)
	ema1 = ema - mt
	xa1 = xa - mt
	y1 = y - np.mean(y)
	if xn is not 'None':
		xn1 = xn - mt
		emn1 = emn - np.mean(emn)	
	
	#scale variance
	ema2 = ema1*np.sqrt(ne2)
	if xn is not 'None':
		emn2 = emn1*np.sqrt(ne2)
		
	#reduce dimensions

		
	#pre-whiten with regularized cov. matrix
	creg = regC(cn1,center_flag=0)
	cl12 = np.real(linalg.inv(linalg.sqrtm(creg)))
	xac = np.dot(cl12,xa1)
	emac = np.dot(cl12,ema2)
	if xn is not 'None':	
		xnc = np.dot(cl12,xn1)
		emnc = np.dot(cl12,emn2)

	#fit scaling factors
	outdict1 = tls(ema1[:,None],y1,cn1,ne=ne2,CI_flag=0,error_flag=1)
	beta1 = outdict1['beta']; err1 = outdict1['err']
	if xn is not 'None':	
		outdict2 = tls(np.vstack((ema1,emn1)).T,y1,cn1,ne=ne2,CI_flag=0,flag_2S=1,error_flag=1)
		beta2 = outdict2['beta']; err2 = outdict1['err']

	#adjust with beta
	emat = np.tile(emac[:,None],(1,ne))
	era = err1[:,0]
	xab = xac - emat + beta1*(emat - era[:,None])
	if xn is not 'None':
		emnt = np.tile(emnc[:,None],(1,ne))
		ern = err2[:,1]
		xnb = xnc - emnt + beta2[1]*(emnt - ern[:,None])
	
	#un-pre-whiten
	q = linalg.inv(cl12)
	xab1 = np.dot(q,xab)
	if xn is not 'None':
		xnb1 = np.dot(q,xnb)
		
	#un-reduce dims
	
	
	if xn is not 'None':
		return xab1,xnb1
	else:
		return xab1
	
#################################################################################
def att_trend(emA,y,cn1,cn2,ne,yrs,emN='None',start_year='None',end_year='None',o_flag=0):
	
	#returns the trends attributable to ALL (and NAT if supplied) forcings over
	#  the period provided (or the whole time series, which is default)
	#uses the CI on beta to provide a CI on the trend
	#should input the ensemble means
	#if o_flag is set to 1, then will also output rescaled ensemble means
	
	import numpy as np
	from scipy import linalg
	import statfunc_mk as sf
	
	#add checks for appropriate size of inputs	
	
	n = np.size(emA)
	
	if np.shape(emA) != np.shape(y):
		raise ValueError('emA and y must have same shape')
	if emN is not 'None':
		if np.shape(emA) != np.shape(emN):
			raise ValueError('emA and emN must have same shape')
	if np.size(yrs) != n:
		raise ValueError('yrs must match size of emA')
	
	if np.shape(cn1)[0]!=n:
		raise ValueError('dimension mismatch: x and cn1')	
	
	ne = float(ne)
	
	if np.size(np.shape(y))==1:
		y = y[:,None]
	if np.size(np.shape(emA))==1:
		emA = emA[:,None]
	if emN is not 'None':
		if np.size(np.shape(emN))==1:
			emN = emN[:,None]
	
	#make sure time series are centered
	mt = np.mean(emA)
	ema1 = emA - mt
	y1 = y - np.mean(y)
	if emN is not 'None':
		emn1 = emN - np.mean(emN)	
	
	#scale variance
	ema2 = ema1*np.sqrt(ne)
	if emN is not 'None':
		emn2 = emn1*np.sqrt(ne)
		
	#pre-whiten with regularized cov. matrix
	creg = regC(cn1,center_flag=0)
	cl12 = np.real(linalg.inv(linalg.sqrtm(creg)))
	emac = np.dot(cl12,ema2)
	if emN is not 'None':	
		emnc = np.dot(cl12,emn2)

	#fit scaling factors
	outdict1 = tls(ema1,y1,cn1,cn2,ne=ne,CI_flag=1,error_flag=1)
	beta1 = outdict1['beta']; beta1CI=outdict1['beta1CI']; err1=outdict1['err']
	if emN is not 'None':	
		outdict2 = tls(np.hstack((ema1,emn1)),y1,cn1,cn2,ne=ne,CI_flag=1,flag_2S=1,error_flag=1)
		beta2 = outdict2['beta']; beta2CI=outdict2['beta1CI']; err2=outdict2['err']
		
	#adjust with beta
	xab = beta1*(emac - err1[:,0][:,None])
	xabl = beta1CI[0,0]*(emac - err1[:,0][:,None])
	xabu = beta1CI[0,1]*(emac - err1[:,0][:,None])
	if emN is not 'None':
		xnb = beta2[1]*(emnc - err2[:,1][:,None])
		xnbl = beta2CI[1,0]*(emac - err2[:,1][:,None])
		xnbu = beta2CI[1,1]*(emac - err2[:,1][:,None])
	
	#un-pre-whiten
	q = linalg.inv(cl12)
	xab1 = np.dot(q,xab)
	xabl1 = np.dot(q,xabl)
	xabu1 = np.dot(q,xabu)
	if emN is not 'None':
		xnb1 = np.dot(q,xnb)
		xnbu1 = np.dot(q,xnbu)
		xnbl1 = np.dot(q,xnbl)
		
	#un-reduce dims
	
	#un-scale var
	xab1 = xab1 / np.sqrt(ne)
	xabu1 = xabu1 / np.sqrt(ne)
	xabl1 = xabl1 / np.sqrt(ne)
	if emN is not 'None':
		xnb1 = xnb1 / np.sqrt(ne)
		xnbu1 = xnbu1 / np.sqrt(ne)
		xnbl1 = xnbl1 / np.sqrt(ne)		
	
	#calculate trend
	if start_year is 'None':
		start_year = np.min(yrs)
	if end_year is 'None':
		end_year = np.max(yrs)
	
	ind = np.zeros((n,))
	ind[(yrs>=start_year)&(yrs<=end_year)] = 1
	ind = ind.astype('bool')
	ny = np.sum(ind)
	
	tA = []
	tA.append(sf.regr(yrs[ind],xabl1[ind,:],add_int=1)[1]*ny)
	tA.append(sf.regr(yrs[ind],xab1[ind,:],add_int=1)[1]*ny)
	tA.append(sf.regr(yrs[ind],xabu1[ind,:],add_int=1)[1]*ny)
	tA = np.squeeze(tA)	
	
	if emN is not 'None':
		tN = []
		tN.append(sf.regr(yrs[ind],xnbl1[ind,:],add_int=1)[1]*ny)
		tN.append(sf.regr(yrs[ind],xnb1[ind,:],add_int=1)[1]*ny)
		tN.append(sf.regr(yrs[ind],xnbu1[ind,:],add_int=1)[1]*ny)
		tN = np.squeeze(tN)
	
	if emN is not 'None':
		if o_flag==1:
			emA1 = [xabl1,xab1,xabu1]
			emN1 = [xnbl1,xnb1,xnbu1]
			return tA,tN,emA1,emN1
		else:
			return tA,tN
	else:
		if o_flag==1:
			emA1 = [xabl1,xab1,xabu1]
			return tA,emA1
		else:
			return tA	
	
	
	
	
