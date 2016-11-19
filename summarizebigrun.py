
import numpy as np
import exceptions
import os
import sys
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import NullFormatter
import iterstat
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
from copy import copy
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
        #grabstardata("/pnfs/des/persistent/smp/v5/","/pnfs/des/persistent/smp/v5/stardatav5.npz")
        #sys.exit()
        data = grabdata(tmpwriter,resultsdir,cd)
        sys.exit()
    else:
        #data = np.load(os.path.join(resultsdir,'Summary','sumdata.npz'))
        data = np.load(cd)
        dostars = False
        if dostars:
            stardata = np.load('/pnfs/des/persistent/smp/v2/stardata2.npz')
            plotstarrms(stardata['starflux'], np.sqrt(stardata['starfluxerr'] ** 2), stardata['starzpt'],
                        stardata['catmag'], stardata['chisq'], stardata['rmsaddin'], stardata['sky'], stardata['skyerr'],
                        stardata['poisson'],
                        title='rmsaddin_')
    print data.keys()
    print len(data['Flux'])
    print np.unique(data['field'])
    #raw_input()
    plotpercentageresid(data['Flux'],data['FakeMag'],data['FitZPT'],data['FakeZPT'],resultsdir)
    plotsigmaresid(data['Flux'],data['Fluxerr'],data['FakeMag'], data['FitZPT'], data['FakeZPT'],data['HostMag'],
                   data['Chisq'],data['rmsaddin'],data['field'],resultsdir)
    #starmag = stardata['starzpt'] - 2.5*np.log10(stardata['starflux'])
    #starmagerr = - 2.5*np.log10(stardata['starflux']) + 2.5*
    #err = 10**(.4*(data['starzpt']-2.5*np.log10()))

def lookup_rms_addin(smpfile,obsid):
    pass

def grabstardata(imagedir,outfile):
    bigdata = {'starflux': [], 'starfluxerr': [], 'starzpt': [], 'diffimzpt':[], 'catmag': [], 'chisq': [], 'rmsaddin': [],
               'sky':[], 'skyerr': [],'psf':[],'poisson':[]}
    zptfiles = []
    cntr = 0
    for dirName, subdirList, fileList in os.walk(imagedir):
        if cntr > 200.: continue
        #print('Found directory: %s' % dirName)
        for fname in fileList:

            if 'globalstar.npz' in fname:
                #print('\t%s' % fname)
                print os.path.join(imagedir,dirName,fname)
                if not 'SN-S2' in fname:
                    if not 'SN-S1' in fname: continue
                try:
                    zptdata = np.load(os.path.join(imagedir,dirName,fname))
                except:
                    print 'could not load'
                    continue
                #zptdata = np.load('/pnfs/des/persistent/smp/v2/20131119_SN-S2/r_21/SNp1_256166_SN-S2_tile20_r_21+fakeSN_rband_dillonzptinfo_globalstar.npz')
                print zptdata.keys()
                if not fname in zptfiles:
                    try:
                        if max(zptdata['cat_mag'])>21.1:
                            continue
                        #if True:
                        bigdata['skyerr'].extend(zptdata['skyerr'])
                        bigdata['sky'].extend(zptdata['sky'])
                        bigdata['starflux'].extend(zptdata['flux_starh'])
                        bigdata['starzpt'].extend(zptdata['flux_starh']*0. + zptdata['fit_zpt'])
                        bigdata['catmag'].extend(zptdata['cat_mag'])
                        bigdata['chisq'].extend(zptdata['chisqu'])
                        bigdata['diffimzpt'].extend(zptdata['fakezpt'])
                        bigdata['starfluxerr'].extend(zptdata['flux_star_std'])
                        psfs = zptdata['psfs']
                        for i in range(len(psfs)):
                            bigdata['psf'].append(psfs[i,:,:])
                            bigdata['poisson'].append(np.sqrt(np.sum(psfs[i,:,:].ravel()**2*zptdata['flux_starh'][i])))
                            #print zptdata['flux_starnormm'][i],zptdata['flux_star_std'][i],bigdata['poisson'][-1]
                            #raw_input()
                        cm = zptdata['cat_mag']
                        fs = zptdata['flux_starh']
                        zp = zptdata['fit_zpt']
                        ww = (cm < 18.) & (cm > 16.)

                        #plt.scatter(cm[ww],float(zp) - cm[ww] - 2.5*np.log10(fs[ww]))
                        # plt.scatter(cm[ww],- 2.5*np.log10(fs[ww]))
                        # plt.savefig('testzpt.png')
                        md, std = iterstat.iterstat(float(zp) - cm[ww] - 2.5*np.log10(fs[ww]),
                                                     startMedian=True, sigmaclip=3, iter=10)
                        print 'worked now std',std
                        bigdata['rmsaddin'].extend(zptdata['flux_starh']*0. + std)
                        #print 'read in ',fname
                        zptfiles.append(fname)
                        cntr += 1
                        print 'CNTR',cntr
                    except:
                        print 'FAILED', fname
                        pass
    np.savez(outfile, **bigdata)

def grabdata(tmpwriter,resultsdir,cd):

    files = os.listdir(os.path.join(resultsdir, 'lightcurves'))
    smpfiles = []
    for f in files:
        if '.smp' in f:
            smpfiles.append(os.path.join(resultsdir, 'lightcurves', f))

    print "Found " + str(len(smpfiles)) + " .smp files"

    if not os.path.exists(os.path.join(resultsdir,'Summary')):
        os.makedirs(os.path.join(resultsdir,'Summary'))
    os.system('rm '+cd+' -f')
    #outfile = os.path.join(resultsdir,'Summary','sumdata.npz')
    outfile = cd
    bigdata = {'Flux':[],'Fluxerr':[],'FakeMag':[],'FitZPT':[],'FakeZPT':[],'HostMag':[],'Chisq':[],
               'starflux':[],'starfluxerr':[],'starzpt':[],'catmag':[],'rmsaddin':[],'field':[]}
    zptfiles = []
    #deep = 0
    tot = len(smpfiles)
    cntr = 0
    for f in smpfiles[:]:
        cntr += 1
        print cntr, 'of',tot
        deep = 0
        data = dt.readcol(f)
        '''
        sn = f.split('/')[-1][0:17]+'.dat'
        snd = open('/pnfs/des/scratch/pysmp/DESY1_imgList_fake/'+sn,'r').read()
        if not '-S1' in snd:
            if not '-S2' in snd:
                deep = 1
                print 'deep'
                continue
        '''
                #raw_input()
        #print data.keys()
        #raw_input()
        try:
            #print data['ID_OBS']
            #raw_input()
            print len(data['FLUX']),len(data['FLUXERR']),len(data['FAKEMAG']),len(data['ZPT']),(data['FAKEZPT'])
            data2 = dt.readcol('./working/lightcurves/' + f.split('/')[-1])
            rms = np.mean(data2['RMSADDIN'][data2['RMSADDIN'] > 0.0])
            bigdata['rmsaddin'].extend(data['CHI2'] * 0. + rms)

            bigdata['Flux'].extend(data['FLUX'])
            bigdata['Fluxerr'].extend(data['FLUXERR'])
            bigdata['FakeMag'].extend(data['FAKEMAG'])
            bigdata['FitZPT'].extend(data['ZPT'])
            bigdata['FakeZPT'].extend(data['FAKEZPT'])
            bigdata['Chisq'].extend(data['CHI2'])

            for m, faz, fiz in zip(data['MJD'],data['FAKEZPT'], data['ZPT']):
                if abs(faz - fiz) > 1:
                    print f
                    print m
                    raw_input('STOPPED')

                    # try:
            #bigdata['rmsaddin'].extend(data['RMSADDIN'])
            #     #print data['RMSADDIN']
            #     #print np.mean(data['RMSADDIN'])
            #     #raw_input()
            # except:
            #     data2 = dt.readcol('./working/lightcurves/'+f.split('/')[-1])
            #     rms = np.mean(data2['RMSADDIN'][data2['RMSADDIN'] > 0.0])
            #     print rms
            #     raw_input()
            #     bigdata['rmsaddin'].extend(data['CHI2']*0. + rms)
            bigdata['field'].extend(data['CHI2']*0 + np.float(deep))
            print f,'read in'
        except:
            print 'Columns missing in file '+f

        # for sf in data['ZPTFILE']:
        #     zptdata = np.load(sf)
        #     if not sf in zptfiles:
        #         try:
        #             bigdata['starfluxerr'].extend(zptdata['flux_star_std'])
        #             bigdata['starflux'].extend(zptdata['flux_star'])
        #             bigdata['starzpt'].extend(zptdata['fit_zpt'])
        #             bigdata['catmag'].extend(zptdata['cat_mag'])
        #             print 'read in ',sf
        #             zptfiles.append(sf)
        #         except:
        #             print 'Missing flux_star_std'
        '''
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
                    if float(l.split()[2]) <= 0.:
                        hgf = 1.
                    else:
                        hgf = float(l.split()[2])
                    hostmag = 27.5 - 2.5 * np.log10(hgf)
        print 'hostmag',hostmag
        bigdata['HostMag'].extend(data['FLUX']*0 + hostmag)
        '''
        #raw_input()
    print 'saving to cachfile'
    #np.savez(outfile,**bigdata)
    print 'saved'
    #tmpwriter.savez(outfile,*bigdata)
    return bigdata


def plotpercentageresid(flux,fakemag,fitzpt,fakezpt,outdir):
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
    plt.savefig(outdir+'/percentagefluxdiff.png')
    print 'saved png'

def plotsigmaresid(flux,fluxerr,fakemag,fitzpt,fakezpt,hostmag,chisqarr,rmsaddin,deep,outdir):
    flux = np.asarray(flux)
    fakemag = np.asarray(fakemag)
    fluxerr = np.asarray(fluxerr)
    fitzpt = np.asarray(fitzpt)
    fakezpt = np.asarray(fakezpt)
    rmsaddin = np.asarray(rmsaddin)
    fakeflux = 10 ** (.4 * (31. - fakemag))
    fakeflux *= 10 ** (-1 * .4 * (fitzpt - fakezpt))
    chisqarr = np.asarray(chisqarr)
    #print np.sqrt(10**(.4*(fitzpt - hostmag))/3.)
    #am = np.argmax(np.sqrt(10**(.4*(fitzpt - hostmag))/3.))
    #for a,f,fe in zip(np.sqrt(10**(.4*(fitzpt - hostmag))/3.),flux,fluxerr):
    #    print a,f,fe
    #print ,flux[am],fluxerr[am]
    #raw_input()

    #fluxerr = np.sqrt(np.asarray(fluxerr)**2+(abs(flux)/3.) + 10**(.4*(fitzpt - hostmag))/3.)
    #fitmag = 31-2.5*np.log10(flux)


    #MY CALCULATION TO GET FLUX RMS ADDIN
    #.1 = -2.5*np.log10(flux) + 2.5*np.log10(flux+rmsfluxerr)
    #rmsaddin + 2.5*log(flux) = 2.5*np.log10(flux+rmsfluxerr)
    #(rmsaddin + 2.5*log(flux))/2.5 = np.log10(flux+rmsfluxerr)
    #10**(.4*(rmsaddin + 2.5*log(flux))) = flux + rmsfluxerr

    fluxrmsaddin = 10 ** (.4 * (rmsaddin + 2.5 * np.log10(flux))) - flux
    #fluxerr = (fluxerr**2 + fluxrmsaddin**2)**.5

    frms = 10**(.4*(2.5*np.log10(abs(flux))+rmsaddin)) - abs(flux)

    #print frms[0:100]
    #print flux[0:100]*rmsaddin[0:100]
    #raw_input()

    fitmag = 31-2.5*np.log10(flux)
    fitmagerr = ((-2.5*np.log10(flux)+2.5*np.log10(flux+fluxerr))**2+(rmsaddin)**2)**.5
    fakemag = fakemag + fakezpt - fitzpt

    hostmag = np.array(hostmag)

    fifx = 10**(.4*(31-fitmag))
    fafx = 10**(.4*(31-fakemag))
    fime = 10**(.4*(31-fitmag+fitmagerr)) - 10**(.4*(31-fitmag))

    #d = (fifx-fafx)/fime
    #d = (fitmag - fakemag)/(fitmagerr*1.08)


    ff = copy(fakeflux)
    ff[ff < 1.] = 1.
    plt.hist(fluxerr/ff,bins=20,normed=True)
    plt.savefig('sinosn.png')
    plt.clf()
    np.savez('simnosn.npz',flux=flux,fakeflux=ff,fluxerr=np.sqrt(fluxerr**2 + abs(flux)/3.8))

    d = (flux - fakeflux) / ((fluxerr**2 - abs(flux)/3.8 + frms**2)**.5*1.08)

    ww = (flux != 0.) #& (deep == 0)

    #fakemag[fakemag==99] = 29.5
    flux = flux[ww]
    fakemag = fakemag[ww]
    fitzpt = fitzpt[ww]
    fakezpt = fakezpt[ww]
    fakeflux= fakeflux[ww]
    fluxerr=fluxerr[ww]
    d = d[ww]
    hostmag = hostmag[ww]
    chisqarr = chisqarr[ww]

    #print flux[0:10]
    #print fakeflux[0:10]
    #print flux.shape
    #print fakeflux.shape
    #d = (flux - fakeflux) / fluxerr

    chisq = (flux - fakeflux) ** 2 / fluxerr ** 2
    chisq = np.nanmean(chisq[abs(d) < 3])

    plt.clf()

    dc = d[abs(d) < 3]

    dc99 = d[fakemag > 90]
    rms99 = np.sqrt(np.nanmean(np.square(dc99[abs(dc99) < 3.])))
    dcr = d[fakemag < 28]
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
    plt.figure(1, figsize=(16, 12))

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
    fakemag[fakemag == 99] = 28.5

    ax1.scatter(fakemag,d,alpha=.3,color='blue')
    ax, ay, aystd = dt.bindata(fakemag, d, np.arange(19.5, max(fakemag), .1),window=.5)
    ax1.plot([19, 28.7], [0, 0],color='grey')
    ax1.plot(ax, ay, linewidth=3, color='orange', label='SMP')
    ax1.plot(ax, ay+aystd, linewidth=2, color='orange',linestyle='--', label='SMP')
    ax1.plot(ax, ay-aystd, linewidth=2, color='orange',linestyle='--', label='SMP')

    #ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

    ax1.set_xlim(19, 28.7)
    ax1.set_ylim(-3., 3.)
    ax1.set_xlabel('Fake Mag')
    ax1.set_ylabel('STD')

    #ax, ayrms= dt.binrms(fakemag, d, np.arange(19.5, max(fakemag), .1),.5)
    #ax3.plot(ax, ayrms, color='blue',label='RMS',linewidth=3)


    ax3.plot([0,100],[1.,1.],linestyle='--',color='black')
    ax3.set_ylim(.7,1.5)
    ax3.legend(fontsize='small')

    fresid = np.zeros(flux.shape)
    for i,f,ff in zip(range(len(flux)),flux,fakeflux):
        if f == 0.:
            fresid[i] = np.nan
        else:
            fresid[i] = (f - ff) / max([abs(ff),1.])
    #fresid[abs(fakeflux) < 1.] = flux[abs(fakeflux) < 1.] - fakeflux[abs(fakeflux) < 1.]

    ax5.hist(fresid, bins=np.arange(-.155,.15,.01),color='blue', orientation='horizontal')


    ax4.scatter(fakemag,fresid,alpha=.3,color='blue')
    axa, aya, aystd = dt.bindata(fakemag,fresid,
                            np.arange(19.5, 26., .1),window=2.)
    ax4.plot([19, 28.7], [0, 0],color='grey')

    ax, ayrms = dt.binrms(fakemag, d, np.arange(19.5, max(fakemag), .1), 1.5)
    ax3.plot(ax, ayrms, color='blue', label='ALL SNe', linewidth=3)
    ax3.plot(ax, ax * 0 + 1., linestyle='--', color='black')

    # ww = hostmag > 25.
    # ax, ayrms = dt.binrms(fakemag[ww], d[ww], np.arange(19.5, max(fakemag), .1), .5)
    # ax3.plot(ax, ayrms, color='red', label='HostMag > 25.', linewidth=3)
    #
    # ww = hostmag < 23.
    # ax, ayrms = dt.binrms(fakemag[ww], d[ww], np.arange(19.5, max(fakemag), .1), .5)
    # ax3.plot(ax, ayrms, color='green', label='HostMag < 23', linewidth=3)
    # ax3.legend(fontsize='small')

    ax4.plot(axa, aya, linewidth=3, color='orange')
    ax4.plot(axa, aya+aystd, linewidth=2, color='orange',linestyle='--')
    ax4.plot(axa, aya-aystd, linewidth=2, color='orange',linestyle='--')
    ax4.set_xlim(ax1.get_xlim())
    ax4.set_ylim(-.1,.1)
    ax4.set_xlabel('Fake Mag')
    ax5.set_xlabel('Counts')
    ax3.set_ylabel('RMS')
    ax4.set_ylabel('(fitflux - fakeflux)/fakeflux')

    ax3.set_xlim(ax1.get_xlim())
    ax2.set_ylim(ax1.get_ylim())
    ax5.set_ylim(ax4.get_ylim())
    ax2.xaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax1.xaxis.set_major_formatter(nullfmt)

    plt.subplots_adjust(wspace=0.001,hspace=0.001)
    plt.savefig(outdir+'/std.png')


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

    rect_scatter = [left, bottom+height/2., width, height/2.]
    rect_scatterflux = [left, bottom, width, height/2.]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom+height/2., 0.2, height/2.]
    rect_histyflux = [left_h, bottom, 0.2, height/2.]

    # start with a rectangular Figure
    plt.figure(1, figsize=(32, 24))

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
    ww = fakemag > 28.
    ax1.scatter(hostmag,d,alpha=.3,color='blue')
    ax, ay, aystd = dt.bindata(hostmag[ww], d[ww], np.arange(min(hostmag), 27.5, .1),window=1.5)
    ax1.plot([min(hostmag), max(hostmag)], [0, 0],color='grey')
    ax1.plot(ax, ay, linewidth=3, color='orange', label='SMP')
    ax1.plot(ax, ay+aystd, linewidth=2, color='orange',linestyle='--', label='SMP')
    ax1.plot(ax, ay-aystd, linewidth=2, color='orange',linestyle='--', label='SMP')

    #ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

    ax1.set_xlim(21,28)
    ax1.set_ylim(-3., 3.)
    ax1.set_xlabel('Host Mag')
    ax1.set_ylabel('STD')

    ax, ayrms= dt.binrms(hostmag, d, np.arange(min(hostmag), 27.5, .1),.5)
    #ax3.plot(ax, ayrms, color='blue',label='RMS',linewidth=3)


    ax3.plot([0,100],[1.,1.],linestyle='--',color='black')
    ax3.set_ylim(.7,1.5)

    fresid = np.zeros(flux.shape)
    for i,f,ff in zip(range(len(flux)),flux,fakeflux):
        if f == 0.:
            fresid[i] = np.nan
        else:
            fresid[i] = (f - ff) / max([abs(ff),1.])
    #fresid[abs(fakeflux) < 1.] = flux[abs(fakeflux) < 1.] - fakeflux[abs(fakeflux) < 1.]

    ax5.hist(fresid, bins=np.arange(-.155,.15,.01),color='blue', orientation='horizontal')

    ax4.scatter(hostmag,fresid,alpha=.3,color='blue')
    ax, ay, aystd = dt.bindata(hostmag,fresid,
                            np.arange(min(hostmag), 27.5, .1),window=1.)
    ax4.plot([min(hostmag), max(hostmag)], [0, 0],color='grey')

    ax4.plot(ax, ay, linewidth=3, color='orange')
    ax4.plot(ax, ay+aystd, linewidth=2, color='orange',linestyle='--')
    ax4.plot(ax, ay-aystd, linewidth=2, color='orange',linestyle='--')
    ax4.set_xlim(ax1.get_xlim())
    ax4.set_ylim(-.1,.1)
    ax4.set_xlabel('Host Mag')
    ax5.set_xlabel('Counts')
    ax3.set_ylabel('RMS')
    ax4.set_ylabel('(fitflux - fakeflux)/fakeflux')

    ax3.set_xlim(ax1.get_xlim())
    ax3.set_ylim(.8,1.4)
    ax2.set_ylim(ax1.get_ylim())
    ax5.set_ylim(ax4.get_ylim())
    ax2.xaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax1.xaxis.set_major_formatter(nullfmt)
    plt.subplots_adjust(wspace=0.001,hspace=0.001)

    ww = fakemag < 100.
    ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag), max(hostmag), .1), 1.5)
    ax3.plot(ax, ayrms, color='blue', label='ALL SNe w/ Light', linewidth=3)
    ax3.plot(ax, ax * 0 + 1., linestyle='--',color='black')

    # ww = fakemag > 28.
    # ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag), max(hostmag), .1), 1.5)
    # ax3.plot(ax, ayrms, color='red', label='FakeMag = 99', linewidth=3)
    #
    # ww = fakemag < 22.
    # ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag), max(hostmag), .1), 1.5)
    # ax3.plot(ax, ayrms, color='green', label='FakeMag < 22', linewidth=3)
    # ax3.legend(fontsize='small')
    #
    # ww = (fakemag > 22.) & (fakemag < 28)
    # ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag), max(hostmag), .1), 1.5)
    # ax3.plot(ax, ayrms, color='purple', label='FakeMag < 22', linewidth=3)
    # ax3.legend(fontsize='small')

    plt.savefig(outdir+'/hostmagstd.png')







    # ----------------------------------------------------------------------------------------------------
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

    rect_scatter = [left, bottom + height / 2., width, height / 2.]
    rect_scatterflux = [left, bottom, width, height / 2.]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom + height / 2., 0.2, height / 2.]
    rect_histyflux = [left_h, bottom, 0.2, height / 2.]

    # start with a rectangular Figure
    plt.figure(1, figsize=(16, 12))

    ax1 = plt.axes(rect_scatter)
    ax3 = plt.axes(rect_histx)
    ax2 = plt.axes(rect_histy)
    ax4 = plt.axes(rect_scatterflux)
    ax5 = plt.axes(rect_histyflux)

    # no labels
    ax2.yaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax5.yaxis.set_major_formatter(nullfmt)

    ax2.hist(d, bins=np.arange(-10, 10, .25), normed=True, label='RMS Fakemag = 99: ' + str(round(rms99, 3)) +
                                                                 '\nRMS Fakemag < 99: ' + str(round(rmsr, 3))
             , orientation='horizontal')
    # label='RMS: ' + str(round(rms, 3)) + '\nChiSq (3sig cut) ' + str(round(chisq, 3)) + '\nMedian ' + str(
    #   round(np.median(d), 3)) + ' +- ' + str(round(np.std(d), 3)),

    import matplotlib.mlab as mlab
    import math
    mean = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.arange(-5, 5, .1)
    ax2.plot(mlab.normpdf(x, mean, sigma), x, color='black', label='Gaussian Normal')

    ax2.set_ylim(-4, 4)
    ax2.set_xlim(0, .5)
    # .xlabel('STDEV')
    # plt.ylabel('Normalized Count')
    ax2.legend(fontsize='small')
    # plt.savefig('stdresid.png')

    # plt.clf()
    ax1.scatter(chisqarr, d, alpha=.3, color='blue')
    ax, ay, aystd = dt.bindata(chisqarr, d, np.arange(0.8, 1.2, .001), window=.01)
    ax1.plot([min(chisqarr), max(chisqarr)], [0, 0], color='grey')
    ax1.plot(ax, ay, linewidth=3, color='orange', label='SMP')
    ax1.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP')
    ax1.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP')

    # ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

    ax1.set_xlim(0.8, 1.2)
    ax1.set_ylim(-3., 3.)
    ax1.set_xlabel('Chi Sq')
    ax1.set_ylabel('STD')

    ax, ayrms = dt.binrms(chisqarr, d, np.arange(0.8, 1.2, .001), .01)
    # ax3.plot(ax, ayrms, color='blue',label='RMS',linewidth=3)


    #ax3.plot([0, 100], [1., 1.], linestyle='--', color='black')
    ax3.set_ylim(.7, 1.8)

    fresid = np.zeros(flux.shape)
    for i, f, ff in zip(range(len(flux)), flux, fakeflux):
        if f == 0.:
            fresid[i] = np.nan
        else:
            fresid[i] = (f - ff) / max([abs(ff), 1.])
    # fresid[abs(fakeflux) < 1.] = flux[abs(fakeflux) < 1.] - fakeflux[abs(fakeflux) < 1.]

    ax5.hist(fresid, bins=np.arange(-.155, .15, .01), color='blue', orientation='horizontal')

    ax4.scatter(chisqarr, fresid, alpha=.3, color='blue')
    ax, ay, aystd = dt.bindata(chisqarr, fresid,
                               np.arange(0.8, 1.2, .001), window=.01)
    ax4.plot([min(chisqarr), max(chisqarr)], [0, 0], color='grey')

    ax4.plot(ax, ay, linewidth=3, color='orange')
    ax4.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--')
    ax4.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--')
    ax4.set_xlim(ax1.get_xlim())
    ax4.set_ylim(-.1, .1)
    ax4.set_xlabel('Chi Sq')
    ax5.set_xlabel('Counts')
    ax3.set_ylabel('RMS')
    ax4.set_ylabel('(fitflux - fakeflux)/fakeflux')

    ax3.set_xlim(ax1.get_xlim())
    ax2.set_ylim(ax1.get_ylim())
    ax5.set_ylim(ax4.get_ylim())
    ax2.xaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax1.xaxis.set_major_formatter(nullfmt)
    plt.subplots_adjust(wspace=0.001, hspace=0.001)

    ax, ayrms = dt.binrms(chisqarr, d, np.arange(0.8, 1.2, .005), .02)
    ax3.plot(ax, ayrms, color='blue', label='ALL SNe', linewidth=3)
    ax3.plot(ax, ax * 0 + 1., linestyle='--', color='black')

    ww = hostmag > 25
    ax, ayrms = dt.binrms(chisqarr[ww], d[ww], np.arange(0.8, 1.2, .005), .02)
    ax3.plot(ax, ayrms, color='red', label='Hostmag > 25', linewidth=3)

    ww = hostmag < 23.
    ax, ayrms = dt.binrms(chisqarr[ww], d[ww], np.arange(0.8, 1.2, .005), .02)
    ax3.plot(ax, ayrms, color='green', label='HostMag < 23', linewidth=3)
    ax3.legend(fontsize='small',loc=2)

    plt.savefig(outdir+'/chisqstd.png')


    print 'saved stdresid.png'




def plotstarrms(flux,fluxerr,zpt,catmag,chisq,rmsaddin,sky,skyerr,poisson,title=''):
    catflux = 10 ** (.4 * (zpt - catmag))
    ff = (flux - catflux) / catflux
    st = np.std(ff)
    ww = (catmag < 21.) & (rmsaddin < 1.) & (abs(ff) < 5*st)

    flux = flux[ww]
    fluxerr = fluxerr[ww]
    zpt = zpt[ww]
    catmag = catmag[ww]
    skyerr= skyerr[ww]
    sky = sky[ww]
    rmsaddin = rmsaddin[ww]
    poisson = poisson[ww]
    #print -2.5*np.log10(skyerr)+zpt
    #raw_input('skyerr in mags')
    #starmagerr = np.sqrt((-2.5*np.log10(sky)+2.5*np.log10(skyerr))**2+rmsaddin**2)
    starmagerr = rmsaddin
    fluxerro = copy(fluxerr)
    #fluxerr = np.sqrt(fluxerr**2)
    catflux = 10**(.4*(zpt-catmag))
    chisq = chisq[ww]
    #chisq = chisq[ww]*1/np.sqrt((abs(flux) / 3.))
    # plt.clf()
    # plt.scatter(catmag,(flux-catflux)/catflux)
    # plt.ylim(-.5,.5)
    # plt.savefig('starresid.png')
    # plt.clf()
    # plt.scatter(catmag,(flux-catflux)/fluxerr)
    # plt.ylim(-5,5)
    # plt.savefig('starstd.png')

    #fluxerr = np.sqrt(fluxerr**2 + )

    starmag = -2.5*np.log10(flux) + zpt

    #print 'fluxerr vs rmsadding' ,np.median((-2.5*np.log10(flux) + 2.5*np.log10(flux+fluxerr))), np.median(rmsaddin)
    #raw_input()
    starmagerr2 = ((-2.5*np.log10(flux) + 2.5*np.log10(flux+fluxerr))**2 + rmsaddin**2 + (-2.5*np.log10(flux) + 2.5*np.log10(flux+poisson))**2 )**.5
    #starmagerr3 = ((-2.5*np.log10(sky) + 2.5*np.log10(sky+skyerr))**2 + rmsaddin[ww]**2)**.5
    skymagerr = -2.5*np.log10(sky) + 2.5*np.log10(sky+skyerr)

    print starmag[0:10]
    print catmag[0:10]
    dmz = (starmag - catmag) / starmagerr
    dmam = (starmag - catmag) / starmagerr2

    #dmas = (starmag - catmag) / starmagerr3
    dsss = (starmag - catmag) / skymagerr


    #raw_input('printing mags')
    d = (flux - catflux) / fluxerr
    ds = (flux - catflux) / skyerr
    dp = (flux-catflux) / poisson

    #chisq = (flux - catflux) ** 2 / catflux

    plt.clf()
    plt.scatter(catmag,chisq)
    plt.ylim(0,15)
    plt.xlim(17,21)
    plt.axhline(1)
    plt.xlabel('Cat Mag')
    plt.ylabel('Chi Sq')
    plt.savefig(title+'chivscat.png')
    #chisq = np.nanmean(chisq[abs(d) < 3])

    plt.clf()

    dc = dmam[abs(dmam) < 3]

    rms = np.sqrt(np.nanmean(np.square(dc[abs(dc) < 3.])))

    plt.clf()

    nullfmt = NullFormatter()  # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom + height / 2., width, height / 2.]
    rect_scatterflux = [left, bottom, width, height / 2.]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom + height / 2., 0.2, height / 2.]
    rect_histyflux = [left_h, bottom, 0.2, height / 2.]

    # start with a rectangular Figure
    plt.figure(1, figsize=(16, 12))

    ax1 = plt.axes(rect_scatter)
    ax3 = plt.axes(rect_histx)
    ax2 = plt.axes(rect_histy)
    ax4 = plt.axes(rect_scatterflux)
    ax5 = plt.axes(rect_histyflux)

    # no labels
    ax2.yaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax5.yaxis.set_major_formatter(nullfmt)

    ax2.hist(dmam, bins=np.arange(-10.3, 10, .3), normed=True, label='RMS: ' + str(round(rms, 3))
             , orientation='horizontal')
    # label='RMS: ' + str(round(rms, 3)) + '\nChiSq (3sig cut) ' + str(round(chisq, 3)) + '\nMedian ' + str(
    #   round(np.median(d), 3)) + ' +- ' + str(round(np.std(d), 3)),

    import matplotlib.mlab as mlab
    import math
    mean = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.arange(-5, 5, .1)
    ax2.plot(mlab.normpdf(x, mean, sigma), x, color='black', label='Gaussian Normal')

    ax2.set_ylim(-4, 4)
    ax2.set_xlim(0, .5)
    # .xlabel('STDEV')
    # plt.ylabel('Normalized Count')
    ax2.legend(fontsize='small')
    # plt.savefig('stdresid.png')

    # plt.clf()

    ax1.scatter(catmag, dmam, alpha=.02, color='blue')
    ax, ay, aystd = dt.bindata(catmag, dmam, np.arange(min(catmag), max(catmag), .1), window=.5)
    ax1.plot([min(catmag), max(catmag)], [0, 0], color='grey')
    ax1.plot(ax, ay, linewidth=3, color='orange', label='SMP')
    ax1.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP')
    ax1.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP')

    # ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

    ax1.set_xlim(16., max(catmag))
    ax1.set_ylim(-3., 3.)
    ax1.set_xlabel('Cat Mag')
    ax1.set_ylabel('STD')

    # ax, ayrms= dt.binrms(fakemag, d, np.arange(19.5, max(fakemag), .1),.5)
    # ax3.plot(ax, ayrms, color='blue',label='RMS',linewidth=3)


    ax3.plot([0, 100], [1., 1.], linestyle='--', color='black')
    ax3.set_ylim(0., 3.)
    ax3.legend(fontsize='small')

    fresid = np.zeros(flux.shape)
    for i, f, ff in zip(range(len(flux)), flux, catflux):
        if f == 0.:
            fresid[i] = np.nan
        else:
            fresid[i] = (f - ff) / max([abs(ff), 1.])
    # fresid[abs(fakeflux) < 1.] = flux[abs(fakeflux) < 1.] - fakeflux[abs(fakeflux) < 1.]

    ax5.hist(fresid, bins=np.arange(-.155, .15, .001), color='blue', orientation='horizontal')

    ax4.scatter(catmag, fresid, alpha=.02, color='blue')
    ax, ay, aystd = dt.bindata(catmag, fresid,
                               np.arange(16., max(catmag), .1), window=1.)
    ax4.plot([19, 28.7], [0, 0], color='grey')

    ax, ayrms = dt.binrms(catmag, d, np.arange(16., max(catmag), .1), .5)
    ax3.plot(ax, ayrms, color='blue', label='Chisq Min Err', linewidth=3,alpha=.7)
    # ax, ayrms = dt.binrms(catmag, dsss, np.arange(16., max(catmag), .1), .5)
    # ax3.plot(ax, ayrms, color='green', label='Skyerr', linewidth=3,alpha=.4)
    ax, ayrms = dt.binrms(catmag, dmz, np.arange(16., max(catmag), .1), .5)
    # print ayrms
    # raw_input('zpt scatter err')
    ax3.plot(ax, ayrms, color='red', label='ZPT Scatter Err', linewidth=3,alpha=.7)

    ax, ayrms = dt.binrms(catmag, dp, np.arange(16., max(catmag), .1), .5)
    ax3.plot(ax, ayrms, color='green', label='Poisson Err', linewidth=3,alpha=.7)

    # ax, ayrms = dt.binrms(catmag, dmas, np.arange(16., max(catmag), .1), .5)
    # ax3.plot(ax, ayrms, color='orange', label='ZPT Scatter Err and Sky Err', linewidth=3,alpha=.4)
    ax, ayrms = dt.binrms(catmag, dmam, np.arange(16., max(catmag), .1), .5)
    ax3.plot(ax, ayrms, color='orange', label='ZPT Scat + Shot + ChisMin', linewidth=3,alpha=.95)
    ax3.plot(ax, ax * 0 + 1., linestyle='--', color='black')
    ax3.legend(fontsize='x-small',loc = 'center right', bbox_to_anchor = (1.25, 0.8))
    # ww = hostmag > 25.
    # ax, ayrms = dt.binrms(catmag[ww], d[ww], np.arange(19.5, max(fakemag), .1), .5)
    # ax3.plot(ax, ayrms, color='red', label='HostMag > 25.', linewidth=3)
    #
    # ww = hostmag < 23.
    # ax, ayrms = dt.binrms(fakemag[ww], d[ww], np.arange(19.5, max(fakemag), .1), .5)
    # ax3.plot(ax, ayrms, color='green', label='HostMag < 23', linewidth=3)
    # ax3.legend(fontsize='small')

    ax4.plot(ax, ay, linewidth=3, color='orange')
    ax4.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--')
    ax4.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--')
    ax4.set_xlim(ax1.get_xlim())
    ax4.set_ylim(-.025, .025)
    ax4.set_xlabel('Cat Mag')
    ax5.set_xlabel('Counts')
    ax3.set_ylabel('RMS')
    ax4.set_ylabel('(fitflux - catflux)/catflux')

    ax3.set_xlim(ax1.get_xlim())
    ax3.set_ylim(.5,2.)
    ax2.set_ylim(ax1.get_ylim())
    ax5.set_ylim(ax4.get_ylim())
    ax2.xaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax1.xaxis.set_major_formatter(nullfmt)

    plt.subplots_adjust(wspace=0.001, hspace=0.001)
    plt.savefig(title+'starstd.png')



    #------------------------------------------------------------------------------------------------





    plt.clf()

    dc = d[abs(d) < 3]

    rms = np.sqrt(np.nanmean(np.square(dc[abs(dc) < 3.])))

    plt.clf()

    nullfmt = NullFormatter()  # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom + height / 2., width, height / 2.]
    rect_scatterflux = [left, bottom, width, height / 2.]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom + height / 2., 0.2, height / 2.]
    rect_histyflux = [left_h, bottom, 0.2, height / 2.]

    # start with a rectangular Figure
    plt.figure(1, figsize=(16, 12))

    ax1 = plt.axes(rect_scatter)
    ax3 = plt.axes(rect_histx)
    ax2 = plt.axes(rect_histy)
    ax4 = plt.axes(rect_scatterflux)
    ax5 = plt.axes(rect_histyflux)

    # no labels
    ax2.yaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax5.yaxis.set_major_formatter(nullfmt)

    ax2.hist(d, bins=np.arange(-10, 10, .25), normed=True, label='RMS: ' + str(round(rms, 3))
             , orientation='horizontal')
    # label='RMS: ' + str(round(rms, 3)) + '\nChiSq (3sig cut) ' + str(round(chisq, 3)) + '\nMedian ' + str(
    #   round(np.median(d), 3)) + ' +- ' + str(round(np.std(d), 3)),

    import matplotlib.mlab as mlab
    import math
    mean = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.arange(-5, 5, .1)
    ax2.plot(mlab.normpdf(x, mean, sigma), x, color='black', label='Gaussian Normal')

    ax2.set_ylim(-4, 4)
    ax2.set_xlim(0, .5)
    # .xlabel('STDEV')
    # plt.ylabel('Normalized Count')
    ax2.legend(fontsize='small')
    # plt.savefig('stdresid.png')

    # plt.clf()

    ax1.scatter(chisq, d, alpha=.3, color='blue')
    ax, ay, aystd = dt.bindata(chisq, d, np.arange(.1, 5, .01), window=.1)
    ax1.plot([.1, 5.], [0, 0], color='grey')
    ax1.plot(ax, ay, linewidth=3, color='orange', label='SMP')
    ax1.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP')
    ax1.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP')

    # ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

    ax1.set_xlim(.1, 5.)
    ax1.set_ylim(-3., 3.)
    ax1.set_xlabel('Chisq')
    ax1.set_ylabel('STD')

    # ax, ayrms= dt.binrms(fakemag, d, np.arange(19.5, max(fakemag), .1),.5)
    # ax3.plot(ax, ayrms, color='blue',label='RMS',linewidth=3)


    ax3.plot([0, 100], [1., 1.], linestyle='--', color='black')
    ax3.set_ylim(.7, 1.7)
    ax3.legend(fontsize='small')

    fresid = np.zeros(flux.shape)
    for i, f, ff in zip(range(len(flux)), flux, catflux):
        if f == 0.:
            fresid[i] = np.nan
        else:
            fresid[i] = (f - ff) / max([abs(ff), 1.])
    # fresid[abs(fakeflux) < 1.] = flux[abs(fakeflux) < 1.] - fakeflux[abs(fakeflux) < 1.]

    ax5.hist(fresid, bins=np.arange(-.155, .15, .01), color='blue', orientation='horizontal')

    ax4.scatter(chisq, fresid, alpha=.3, color='blue')
    ax, ay, aystd = dt.bindata(chisq, fresid,
                               np.arange(.1, 5, .01), window=.1)
    ax4.plot([.1, 5], [0, 0], color='grey')

    ax, ayrms = dt.binrms(chisq, d, np.arange(.1, 5., .01), .1)
    ax3.plot(ax, ayrms, color='blue', label='ALL Standard Stars in r Band', linewidth=3)
    ax3.plot(ax, ax * 0 + 1., linestyle='--', color='black')

    # ww = hostmag > 25.
    # ax, ayrms = dt.binrms(catmag[ww], d[ww], np.arange(19.5, max(fakemag), .1), .5)
    # ax3.plot(ax, ayrms, color='red', label='HostMag > 25.', linewidth=3)
    #
    # ww = hostmag < 23.
    # ax, ayrms = dt.binrms(fakemag[ww], d[ww], np.arange(19.5, max(fakemag), .1), .5)
    # ax3.plot(ax, ayrms, color='green', label='HostMag < 23', linewidth=3)
    # ax3.legend(fontsize='small')

    ax4.plot(ax, ay, linewidth=3, color='orange')
    ax4.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--')
    ax4.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--')
    ax4.set_xlim(ax1.get_xlim())
    ax4.set_ylim(-.025, .025)
    ax4.set_xlabel('Chisq')
    ax5.set_xlabel('Counts')
    ax3.set_ylabel('RMS')
    ax4.set_ylabel('(fitflux - catflux)/catflux')

    ax3.set_xlim(ax1.get_xlim())
    ax2.set_ylim(ax1.get_ylim())
    ax5.set_ylim(ax4.get_ylim())
    ax2.xaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax1.xaxis.set_major_formatter(nullfmt)

    plt.subplots_adjust(wspace=0.001, hspace=0.001)
    plt.savefig(title+'starchi.png')





    print 'saved starstd.png'

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
    resultsdir = '/pnfs/des/scratch/pysmp/smp_04_modelerrors'
    resultsdir = '/pnfs/des/scratch/pysmp/smp_02_simnosnnoskyerr'
    resultsdir= './working/'
    resultsdir= '/export/scratch0/ps1sn1/data/v10.0/GPC1v3/eventsv1/smpworkspace/PS_TEST1/'
    #resultsdir = './workingsimnosn'
    isfermigrid = False
    cacheddata = False
    cd = '/pnfs/des/scratch/pysmp/smp_04_modelerrors/np_data/summary_results.npz'
    #cd = '/pnfs/des/scratch/pysmp/smp_02_simnosnnoskyerr/np_data/summary_results.npz'
    import sys, getopt

    try:
        args = sys.argv[1:]

        opt, arg = getopt.getopt(
            args, "fd:rd:cd:cdf",
            longopts=["fakedir=", "resultsdir=", "cacheddata", "cashedfile="])

    except getopt.GetoptError as err:
        print "No command line argument s"


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