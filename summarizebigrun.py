
import numpy as np
import exceptions
import os
import sys
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import NullFormatter

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
#
# fakedir='/pnfs/des/scratch/pysmp/DESY1_imgList_fake/'
# resultsdir = '/pnfs/des/scratch/pysmp/smp_02_simnosnnoskyerr'
# isfermigrid = True
# cacheddata = True
# cd = 'tmp_snse.npz'
def go(fakedir,resultsdir,cacheddata,cd,isfermigrid=False):

    if isfermigrid:
        useifdh = True
    else:
        useifdh = False
    tmpwriter = dt.tmpwriter(useifdh=useifdh)

    if not cacheddata:
        data = grabdata(tmpwriter,resultsdir)
    else:
        #data = np.load(os.path.join(resultsdir,'Summary','sumdata.npz'))
        data = np.load(cd)
    print data.keys()
    print len(data['Flux'])

    plotpercentageresid(data['Flux'],data['FakeMag'],data['FitZPT'],data['FakeZPT'])
    plotsigmaresid(data['Flux'],data['Fluxerr'],data['FakeMag'], data['FitZPT'], data['FakeZPT'],data['HostMag'])


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
    outfile = cd
    bigdata = {'Flux':[],'Fluxerr':[],'FakeMag':[],'FitZPT':[],'FakeZPT':[],'HostMag':[]}

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

        fakef = f.split('/')[-1][:17]
        filt = f.split('/')[-1][18]
        fakefile = os.path.join(fakedir,fakef+'.dat')
        ff = open(fakefile,'r').readlines()
        print 'fileter',filt
        hostmag = -999
        #raw_input()
        for l in ff:
            key = l.split(':')[0]
            if key == 'HOSTGAL_SB_FLUXCAL':
                if filt == 'r':
                    hostmag = 27.5 - 2.5 * np.log10(float(l.split()[2]))
        print 'hostmag',hostmag
        bigdata['HostMag'].extend(data['FLUX']*0 + hostmag)
        #raw_input()
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

    ww = fakemag < 24.5

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

def plotsigmaresid(flux,fluxerr,fakemag,fitzpt,fakezpt,hostmag):
    flux = np.asarray(flux)
    fakemag = np.asarray(fakemag)
    fitzpt = np.asarray(fitzpt)
    fakezpt = np.asarray(fakezpt)
    fakeflux = 10 ** (.4 * (31. - fakemag))
    fakeflux *= 10 ** (-1 * .4 * (fitzpt - fakezpt))
    fluxerr = np.asarray(fluxerr)+(abs(flux)/4.)**.5
    hostmag = np.array(hostmag)


    d = (flux - fakeflux) / fluxerr
    ww = fakemag < 28.

    #fakemag[fakemag==99] = 29.5
    # flux = flux[ww]
    # fakemag = fakemag[ww]
    # fitzpt = fitzpt[ww]
    # fakezpt = fakezpt[ww]
    # fakeflux= fakeflux[ww]
    # fluxerr=fluxerr[ww]
    # d = d[ww]
    # hostmag = hostmag[ww]

    #print flux[0:10]
    #print fakeflux[0:10]
    #print flux.shape
    #print fakeflux.shape
    d = (flux - fakeflux) / fluxerr

    chisq = (flux - fakeflux) ** 2 / fluxerr ** 2
    chisq = np.nanmean(chisq[abs(d) < 3])

    plt.clf()

    dc = d[abs(d) < 3]

    dc99 = d[fakemag == 99]
    rms99 = np.sqrt(np.nanmean(np.square(dc99[abs(dc99) < 3.])))
    dcr = d[fakemag != 99]
    rmsr = np.sqrt(np.nanmean(np.square(dcr[abs(dcr) < 3.])))

    #f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

    # fig = plt.figure(figsize=(16, 12))
    # gs = gridspec.GridSpec(1, 3, width_ratios=[4, 1])
    # ax1 = plt.subplot(gs[0])
    # ax2 = plt.subplot(gs[1])

    nullfmt = NullFormatter()  # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom+height/2., width, height/2.]
    rect_scatterflux = [left, bottom, width, height/2.]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom+height/2., 0.2, height/2.]
    rect_histyflux = [left_h, bottom, 0.2, height/2.]

    # start with a rectangular Figure
    plt.figure(1, figsize=(8, 8))

    ax1 = plt.axes(rect_scatter)
    ax3 = plt.axes(rect_histx)
    ax2 = plt.axes(rect_histy)
    ax4 = plt.axes(rect_scatterflux)
    ax5 = plt.axes(rect_histyflux)

    # no labels
    ax2.yaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax5.yaxis.set_major_formatter(nullfmt)


    ax2.hist(d, bins=np.arange(-10, 10, .25), normed=True,label='RMS Fakemag = 99: ' + str(round(rms99, 3))+
                                                                '\nRMS Fakemag < 99: '+ str(round(rmsr, 3))
             ,orientation='horizontal')
             #label='RMS: ' + str(round(rms, 3)) + '\nChiSq (3sig cut) ' + str(round(chisq, 3)) + '\nMedian ' + str(
             #   round(np.median(d), 3)) + ' +- ' + str(round(np.std(d), 3)),

    import matplotlib.mlab as mlab
    import math
    mean = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.arange(-5, 5, .1)
    ax2.plot(mlab.normpdf(x, mean, sigma),x, color='black', label='Gaussian Normal')

    ax2.set_ylim(-4, 4)
    ax2.set_xlim(0,.5)
    #.xlabel('STDEV')
    #plt.ylabel('Normalized Count')
    ax2.legend(fontsize='small')
    #plt.savefig('stdresid.png')

    #plt.clf()
    ax1.scatter(fakemag,d,alpha=.3,color='blue')
    ax, ay, aystd = bindata(fakemag, d, np.arange(min(fakemag), max(fakemag), .5))
    ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

    ax1.plot([20, 27], [0, 0])
    ax1.set_xlim(20, 27)
    ax1.set_ylim(-3., 3.)
    ax1.set_xlabel('Fake Mag')
    ax1.set_ylabel('STD')

    ax, ayrms= dt.binrms(fakemag, d, np.arange(min(fakemag), max(fakemag), .1),.5)
    ax3.plot(ax, ayrms, color='blue',label='RMS',linewidth=3)


    ax3.plot(ax,ax*0+1.,linestyle='--',color='black')
    ax3.set_ylim(.7,1.5)
    ax3.legend(fontsize='small')

    fresid = (flux - fakeflux) / fakeflux
    ax5.hist(fresid, bins=np.arange(-.15,.15,.01),color='blue', orientation='horizontal')

    ax4.scatter(fakemag,fresid,alpha=.3,color='blue')
    ax, ay, aystd = dt.bindata(fakemag,fresid,
                            np.arange(min(fakemag), max(fakemag), .1),window=.5)
    ax4.plot(ax, ay, linewidth=3, color='green')
    ax4.plot(ax, ay+aystd, linewidth=2, color='green',linestyle='--')
    ax4.plot(ax, ay-aystd, linewidth=2, color='green',linestyle='--')
    ax4.set_xlim(ax1.get_xlim())
    ax4.set_ylim(-.15,.15)
    ax4.set_xlabel('Fake Mag')
    ax5.set_xlabel('Counts')

    ax3.set_xlim(ax1.get_xlim())
    ax2.set_ylim(ax1.get_ylim())
    ax5.set_ylim(ax4.get_ylim())
    ax2.xaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    #plt.tight_layout()
    plt.subplots_adjust(wspace=0.001,hspace=0.001)
    plt.savefig('std.png')


    #--------------------------------------------------------------------------------------
    plt.clf()
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    nullfmt = NullFormatter()  # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    plt.figure(1, figsize=(8, 8))

    ax1 = plt.axes(rect_scatter)
    ax3 = plt.axes(rect_histx)
    ax2 = plt.axes(rect_histy)

    # no labels
    ax2.yaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)

    ax2.hist(d, bins=np.arange(-10, 10, .25), normed=True, label='RMS Fakemag = 99: ' + str(round(rms99, 3))+
                                                                '\nRMS Fakemag < 99: '+ str(round(rmsr, 3)),
             # label='RMS: ' + str(round(rms, 3)) + '\nChiSq (3sig cut) ' + str(round(chisq, 3)) + '\nMedian ' + str(
             #   round(np.median(d), 3)) + ' +- ' + str(round(np.std(d), 3)),
             orientation='horizontal')
    import matplotlib.mlab as mlab
    import math
    mean = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.arange(-5, 5, .1)
    ax2.plot(mlab.normpdf(x, mean, sigma), x, color='black', label='Gaussian Normal')

    #ax2.set_ylim(-3, 3)
    ax2.set_xlim(0, .5)
    # .xlabel('STDEV')
    # plt.ylabel('Normalized Count')
    ax2.legend(fontsize='small')
    # plt.savefig('stdresid.png')

    # plt.clf()
    ax1.scatter(hostmag, d, alpha=.3)
    ax, ay, aystd = bindata(hostmag, d, np.arange(min(hostmag), max(hostmag), .5))
    ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

    ax1.plot([21.5, 27], [0, 0])
    ax1.set_xlim(21.5, 27)
    ax1.set_ylim(-3., 3.)
    ax1.set_xlabel('27.5 - 2.5 * log( Host SB Fluxcal )')
    ax1.set_ylabel('STD')

    ax, ayrms = dt.binrms(hostmag, d, np.arange(min(hostmag), max(hostmag), .1), .5)
    ax3.plot(ax, ayrms, color='blue', label='ALL SNe', linewidth=3)
    ax3.plot(ax, ax * 0 + 1., linestyle='--',color='black')

    ww = fakemag == 99
    ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag), max(hostmag), .1), .5)
    ax3.plot(ax, ayrms, color='red', label='FakeMag = 99', linewidth=3)

    ww = fakemag < 22.
    ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag), max(hostmag), .1), .5)
    ax3.plot(ax, ayrms, color='green', label='FakeMag < 22', linewidth=3)

    ax3.set_ylim(.7, 1.5)
    ax3.set_ylabel('RMS')
    ax3.legend(fontsize='small')

    ax3.set_xlim(ax1.get_xlim())
    ax2.set_ylim(ax1.get_ylim())

    # plt.tight_layout()
    plt.subplots_adjust(wspace=0.001, hspace=0.001)
    plt.savefig('hostmagstd.png')




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
    fakedir = '/pnfs/des/scratch/pysmp/DESY1_imgList_fake/'
    resultsdir = '/pnfs/des/scratch/pysmp/smp_02_simnosnnoskyerr'
    isfermigrid = False
    cacheddata = True
    cd = 'tmp_snse.npz'

    import sys, getopt

    try:
        args = sys.argv[1:]

        opt, arg = getopt.getopt(
            args, "fd:rd:cd:cdf",
            longopts=["fakedir=", "resultsdir=", "cacheddata", "cashedfile="])

    except getopt.GetoptError as err:
        print "No command line arguments"


    for o,a in opt:
        if o in ["-h","--help"]:
            print __doc__
            sys.exit(0)
        elif o in ["-fd","--fakedir"]:
            fakedir = a
        elif o in ["-rd","--resultsdir"]:
            resultsdir = a
        elif o in ["-cd","--cacheddata"]:
            cacheddata = True
        elif o in ["-cdf","--cachedfile"]:
            cd = a

    go(fakedir,resultsdir,cacheddata,cd)