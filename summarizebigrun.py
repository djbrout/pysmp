
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

def go(resultsdir):
    files = os.listdir(os.path.join(resultsdir,'lightcurves'))
    smpfiles = []
    for f in files:
        if '.smp' in f:
            smpfiles.append(os.path.join(resultsdir,'lightcurves',f))
            print smpfiles[-1]

if __name__ == "__main__":
    go(resultsdir)