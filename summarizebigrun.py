
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


resultsdir = '/pnfs/des/scratch/pysmp/smp_02_simnosnnoskyerr'
isfermigrid = True
cacheddata = True

def go(resultsdir,isfermigrid=False):

    if isfermigrid:
        useifdh = True
    else:
        useifdh = False
    tmpwriter = dt.tmpwriter(useifdh=useifdh)

    if not cacheddata:
        data = grabdata(tmpwriter,resultsdir)
    else:
        #data = np.load(os.path.join(resultsdir,'Summary','sumdata.npz'))
        data = np.load('tmp.npz')
    print data.keys()
    print len(data['Flux'])

    plotpercentageresid(data['Flux'],data['FakeMag'],data['FitZPT'],data['FakeZPT'])
    plotsigmaresid(data['Flux'],data['Fluxerr'],data['FakeMag'], data['FitZPT'], data['FakeZPT'])


def grabdata(tmpwriter,resultsdir):

    files = os.listdir(os.path.join(resultsdir, 'lightcurves'))
    smpfiles = []
    for f in files:
        if '.smp' in f:
            smpfiles.append(os.path.join(resultsdir, 'lightcurves', f))

    print "Found " + str(len(smpfiles)) + " .smp files"

    if not os.path.exists(os.path.join(resultsdir,'Summary')):
        os.makedirs(os.path.join(resultsdir,'Summary'))
    #outfile = os.path.join(resultsdir,'Summary','sumdata.npz')
    outfile = 'tmp.npz'
    bigdata = {'Flux':[],'Fluxerr':[],'FakeMag':[],'FitZPT':[],'FakeZPT':[]}

    for f in smpfiles:
        data = dt.readcol(f)
        try:
            print len(data['FLUX']),len(data['FLUXERR']),len(data['FAKEMAG']),len(data['ZPT']),(data['FAKEZPT'])
            bigdata['Flux'].extend(data['FLUX'])
            bigdata['Fluxerr'].extend(data['FLUXERR'])
            bigdata['FakeMag'].extend(data['FAKEMAG'])
            bigdata['FitZPT'].extend(data['ZPT'])
            bigdata['FakeZPT'].extend(data['FAKEZPT'])
            print f,'read in'
        except:
            print 'Columns missing in file '+f
    print 'saving to cachfile'
    np.savez(outfile,**bigdata)
    print 'saved'
    #tmpwriter.savez(outfile,*bigdata)
    return bigdata


def plotpercentageresid(flux,fakemag,fitzpt,fakezpt):
    flux = np.asarray(flux)
    fakemag = np.asarray(fakemag)
    print fakemag.shape
    print flux.shape
    #print fakemag[0].shape
    #sys.exit()
    fitzpt = np.asarray(fitzpt)
    fakezpt = np.asarray(fakezpt)

    fakeflux = 10**(.4*(31. - fakemag))
    fakeflux *= 10**(-1*.4*(fitzpt - fakezpt))

    ww = fakemag < 99.

    fig = plt.figure(figsize=(15, 10))
    plt.scatter(fakemag[ww],(flux[ww]-fakeflux[ww])/fakeflux[ww],alpha=.5)
    ax, ay, aystd = bindata(fakemag[ww],(flux[ww]-fakeflux[ww])/fakeflux[ww],
                            np.arange(min(fakemag[ww]),max(fakemag[ww]), .5))
    plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')

    plt.plot([20,27],[0,0])
    plt.xlim(20,25)
    plt.ylim(-.2,.2)
    plt.xlabel('Fake Mag')
    plt.ylabel('Percentage Flux Difference')
    plt.savefig('percentagefluxdiff.png')
    print 'saved png'

def plotsigmaresid(flux,fluxerr,fakemag,fitzpt,fakezpt):
    flux = np.asarray(flux)
    fluxerr = np.asarray(fluxerr)
    fakemag = np.asarray(fakemag)
    fitzpt = np.asarray(fitzpt)
    fakezpt = np.asarray(fakezpt)
    fakeflux = 10 ** (.4 * (31. - fakemag))
    fakeflux *= 10 ** (-1 * .4 * (fitzpt - fakezpt))

    ww = fakemag < 99.
    print flux[0:10]
    print fakeflux[0:10]
    print flux.shape
    print fakeflux.shape
    d = (flux - fakeflux) / fluxerr

    chisq = (flux - fakeflux) ** 2 / fluxerr ** 2
    chisq = np.nanmean(chisq[abs(d) < 3])

    plt.clf()

    dc = d[abs(d) < 3]
    rms = np.sqrt(np.nanmean(np.square(dc)))

    plt.hist(d, bins=np.arange(-10, 10, .25), normed=True,
             label='RMS: ' + str(round(rms, 3)) + '\nChiSq (3sig cut) ' + str(round(chisq, 3)) + '\nMedian ' + str(
                 round(np.median(d), 3)) + ' +- ' + str(round(np.std(d), 3)))
    import matplotlib.mlab as mlab
    import math
    mean = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.arange(-5, 5, .1)
    plt.plot(x, mlab.normpdf(x, mean, sigma), color='black', label='Gaussian Normal')

    plt.xlim(-5, 5)
    plt.ylim(0,.5)
    plt.xlabel('STDEV')
    plt.ylabel('Normalized Count')
    plt.legend()
    plt.savefig('stdresid.png')

    plt.clf()
    plt.scatter(fakemag,d,alpha=.3)
    ax, ay, aystd = bindata(fakemag, d, np.arange(min(fakemag), max(fakemag), .5))
    plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')

    plt.plot([20, 27], [0, 0])
    plt.xlim(20, 25)
    plt.ylim(-4., 4.)
    plt.xlabel('Fake Mag')
    plt.ylabel('STD')
    plt.savefig('stdvsmag.png')



    print 'saved stdresid.png'



def bindata(x, y, bins, returnn=False):
    medians = np.zeros(len(bins) - 1)
    mads = np.zeros(len(bins) - 1)
    nums = np.zeros(len(bins) - 1)

    for i in np.arange(len(bins) - 1):
        bs = bins[i]
        bf = bins[i + 1]
        ww = [(x > bs) & (x < bf)]
        yhere = y[ww]
        yhere = yhere[np.isfinite(yhere)]
        ss = [abs(yhere) < 3. * np.std(yhere)]
        try:
            nums[i] = len(yhere[ss])
            medians[i] = np.median(yhere[ss])
            mads[i] = 1.48 * np.median(abs(yhere[ss] - medians[i])) * 1 / np.sqrt(len(yhere[ss]))
        except IndexError:
            print 'excepted'
            nums[i] = 0.
            medians[i] = np.nan
            mads[i] = np.nan
    xvals = (bins[1:] + bins[:-1]) / 2.
    if returnn:
        return xvals, medians, mads, nums
    return xvals, medians, mads


if __name__ == "__main__":
    go(resultsdir)