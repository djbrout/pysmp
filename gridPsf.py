import numpy as np
import os
import scipy.signal
import scipy.ndimage as nd
import sys
from sys import getsizeof
import time

#@profile
def gridOne(psffile=None, pixstep=None, numpixsteps=None,
         x=None, y=None, gridsize=None, verbose=True):

    bigpsfarr = np.zeros((gridsize,gridsize,numpixsteps*2,numpixsteps*2,100),dtype='float32')
    for i in range(gridsize):
        for j in range(gridsize):
            for k in range(numpixsteps*2):
                for l in range(numpixsteps*2):
                    bigpsfarr[i,j,k,l] = .00001
    print bigpsfarr.shape
    print bigpsfarr[20,20,10,9]
    print 'size of one psf grid:',getsizeof(bigpsfarr)/1000./1000.,getsizeof(0.00001)

    return bigpsfarr


def interpolate():
    x = np.arange(20)
    y = np.arange(20)
    from scipy.interpolate import UnivariateSpline
    spl = UnivariateSpline(x, y, k=1)
    psf = np.zeros(30*30*100)
    vfunc = np.vectorize(spl)

    t1 = time.time()
    #for i in range(30*30*100):
    #psf =  spl(10.101)
    psf = vfunc(np.ones(30*30*100)/10.)
    t2 = time.time()
    print 'one mcmc step time:',t2-t1

def gridAll(psffiles=[], pixstep=None, numpixsteps=None,
         xs=[], ys=[], gridsize=None, verbose=True):

    return

if __name__ == "__main__":
    psffile = ''
    pixstep = .01
    numpixsteps = 10
    x = 200
    y = 200
    gridsize=30
    #gridOne(psffile, pixstep, numpixsteps, x, y, gridsize)
    interpolate()