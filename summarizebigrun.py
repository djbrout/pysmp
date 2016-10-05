
import numpy as np
import exceptions
import os
import sys
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import pyfits as pf
import scipy.signal
from copy import copy
import time
from astropy.io import fits
from scipy.interpolate import UnivariateSpline
import sigma_clip
import meanclip
import dilltools as dt


resultsdir = '/pnfs/des/scratch/pysmp/smp_02/'
isfermigrid = True
cacheddata = False

def go(resultsdir,isfermigrid=False):

    if isfermigrid:
        useifdh = True
    else:
        useifdh = False
    tmpwriter = dt.tmpwriter(useifdh=useifdh)

    if not cacheddata:
        data = grabdata(tmpwriter,resultsdir)
    else:
        data = np.load(os.path.join(resultsdir,'Summary','sumdata.npz'))

    print data.shape

def grabdata(tmpwriter,resultsdir):

    files = os.listdir(os.path.join(resultsdir, 'lightcurves'))
    smpfiles = []
    for f in files:
        if '.smp' in f:
            smpfiles.append(os.path.join(resultsdir, 'lightcurves', f))

    print "Found " + str(len(smpfiles)) + " .smp files"

    if not os.path.exists(os.path.join(resultsdir,'Summary')):
        os.makedirs(os.path.join(resultsdir,'Summary'))
    outfile = os.path.join(resultsdir,'Summary','sumdata.npz')

    for f in smpfiles:
        #print open(f,'r').readlines()[0]
        data = dt.readcol(f)
        print data.shape
        sys.exit()

if __name__ == "__main__":
    go(resultsdir)