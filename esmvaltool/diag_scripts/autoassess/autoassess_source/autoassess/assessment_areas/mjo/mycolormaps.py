import matplotlib as mpl
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

from scipy import interpolate
from pylab import *
from numpy import outer

def _spreadColorbar(C):
    x = np.linspace(0,256,len(C))
    xnew = np.arange(256)
    fr = interpolate.interp1d(x, C[:,0])
    fg = interpolate.interp1d(x, C[:,1])
    fb = interpolate.interp1d(x, C[:,2])
    C = np.array([fr(xnew).astype(int), fg(xnew).astype(int), fb(xnew).astype(int)]).T
    return C
    

def getcolors(name,  colorReverse = False):
    '''
    try:  
        os.environ["NCARG_ROOT"]
    except KeyError: 
        print "Please set the environment variable NCARG_ROOT"
        sys.exit(1)
    '''
    #NCLcolormapdir = os.environ['NCARG_ROOT']+'/lib/ncarg/colormaps/'
    #NCLcolormapdir = '../colormaps/'
    colorbarfilename = os.path.join(os.path.dirname(os.path.abspath(__file__)), name+'.rgb')
    
    #NCLcolormapdir+name+'.rgb'
    
    if os.path.isfile(colorbarfilename):
        C = np.loadtxt(colorbarfilename, skiprows=1)
        if len(C)<256:
            C = _spreadColorbar(C)
                        
        if not name=='test_default':
            cm = mpl.colors.ListedColormap(C/256.)
            if colorReverse :
                cm = mpl.colors.ListedColormap(C[::-1, :]/256.)
        else:
            cm = mpl.colors.ListedColormap(C)
            if colorReverse :
                cm = mpl.colors.ListedColormap(C[::-1, :])
        
        return cm
    else: 
        print colorbarfilename+ " does not exist."
        sys.exit(1)

def main():
    c = getcolors('ncl_default',  colorReverse= True)
    a=outer(arange(0,1,0.01),ones(10))
    imshow(a,aspect='auto',cmap=c,origin="lower")
    plt.show()
if __name__ == '__main__':
	main()
