import os
#import matplotlib as m

#m.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(threshold=np.nan)
import sys
import dilltools as dt

default_checkstar_file = '/Volumes/ExtraSpace/pysmp_downloads/des_fake_00229567_r_standardstarfits.txt'


# def checkstars(checkstarfile=default_checkstar_file):
#     cols = dt.readcol(checkstarfile,delim='\t')
#     print cols.keys()
#     plt.scatter(-2.5*np.log10(cols['Fit Flux']),(cols['Fit Flux']-cols['Galsim Fit Flux'])/cols['Fit Flux'],
#                 alpha=.1,color='black')
#     xvals, medians, mads = dt.bindata(-2.5*np.log10(cols['Fit Flux']),(cols['Fit Flux']-cols['Galsim Fit Flux'])/cols['Fit Flux'],
#                                       np.arange(-14,-9,.25))
#     plt.errorbar(xvals,medians,mads,color='blue',fmt='o')
#     plt.axhline(0,color='blue')
#     plt.xlabel('-2.5*np.log10(fitflux)')
#     plt.ylabel('myPSF Flux - GalsimPSF Flux / myPSF Flux')
#     plt.ylim(-.001,.001)
#     plt.title('Star Fits, 1CCD, All Epochs, All Stars')
#     plt.savefig('/Volumes/ExtraSpace/pysmp_downloads/starfit_resids.png')
#     #print cols.keys()
#
#
#     plt.clf()
#     plt.hist(cols['Fit Flux Chisq'], label='my PSF Model median='+str(round(np.median(cols['Fit Flux Chisq']),4)),
#              bins=np.arange(.525,2,.05),alpha=.75,normed=True)
#     plt.hist(cols['Galsim Fit Flux Chisq'], label='Galsim PSF Model median='+str(round(np.median(cols['Galsim Fit Flux Chisq']),4)),
#              bins=np.arange(.525,2,.05),alpha=.75,normed=True)
#     plt.xlabel('Reduced Chisq')
#     plt.legend(prop={'size':10})
#     plt.savefig('/Volumes/ExtraSpace/pysmp_downloads/starfit_chisqhist.png')
#
#     plt.clf()
#     plt.hist(cols['Fit Flux DMS'], label='my PSF Model', bins=np.arange(-10500,5000,1000),alpha=.75,normed=True)
#     plt.hist(cols['Galsim Fit Flux DMS'], label='Galsim PSF Model', bins=np.arange(-10500,5000,1000),alpha=.75,normed=True)
#     plt.xlabel('Data - Sim')
#     plt.legend()
#     plt.savefig('/Volumes/ExtraSpace/pysmp_downloads/starfit_dmshist.png')
#
#     plt.clf()
#     plt.scatter(cols['Galsim Fit Flux Chisq'], (cols['Fit Flux'] - cols['Galsim Fit Flux']) / cols['Fit Flux'],
#                 alpha=.1, color='black')
#
#     xvals, medians, mads = dt.bindata(cols['Galsim Fit Flux Chisq'],
#                                       (cols['Fit Flux'] - cols['Galsim Fit Flux']) / cols['Fit Flux'],
#                                       np.arange(.5, 3, .1))
#     plt.errorbar(xvals, medians, mads, color='blue', fmt='o')
#     plt.axhline(0, color='blue')
#     plt.xlim(.5,3)
#     plt.ylim(-.001,.001)
#     plt.xlabel('Galsim Chisq')
#     plt.ylabel('myPSF Flux - GalsimPSF Flux / myPSF Flux')
#     plt.title('Star Fits, 1CCD, All Epochs, All Stars')
#     plt.savefig('/Volumes/ExtraSpace/pysmp_downloads/starfit_vs_galsimchisq.png')
#
#     plt.clf()
#     x = np.sqrt((np.absolute(cols['xstar']-np.round(cols['xstar'])) - .5)**2+
#                 (np.absolute(cols['ystar']-np.round(cols['ystar'])) -.5)**2)
#     plt.scatter(x, (cols['Fit Flux'] - cols['Galsim Fit Flux']) / cols['Fit Flux'],
#                 alpha=.1, color='black')
#
#     xvals, medians, mads = dt.bindata(x, (cols['Fit Flux'] - cols['Galsim Fit Flux']) / cols['Fit Flux'],
#                                       np.arange(0, .7, .1))
#     plt.errorbar(xvals, medians, mads, color='blue', fmt='o')
#     plt.axhline(0, color='blue')
#     plt.xlim(0, .7)
#     plt.ylim(-.001, .001)
#     plt.xlabel('distance from center of pixel')
#     plt.ylabel('myPSF Flux - GalsimPSF Flux / myPSF Flux')
#     plt.title('Star Fits, 1CCD, All Epochs, All Stars')
#     plt.savefig('/Volumes/ExtraSpace/pysmp_downloads/starfit_vs_dtc.png')
#
#     return cols


def checkstars(smpfile):
    data = dt.readcol(smpfile)
    #print data.keys()
    #print data['ZPTFILE']
    zptfiles = data['ZPTFILE'][:-1]
    mjd = []
    fitmag = []
    catmag = []
    fitzpt = []
    ra = []
    dec = []
    c = 0
    for z in zptfiles:
        c += 1
        #print z
        #print z
        #print zd.keys()
        try:
            print z
            zd = np.load(z)

            mjd.extend(zd['cat_mag']*0.+ zd['mjd'])
            fitmag.extend(zd['fit_mag'])
            #print zd['mpfit_mag']
            #print zd['cat_mag']
            #raw_input(  )
            catmag.extend(zd['cat_mag'])
            fitzpt.extend(zd['cat_mag']*0. + zd['fit_zpt'])
            #ra.extend()
            #dec.extend()
            #raw_input()
        except:
            print 'missing column'
    print len(fitmag),len(catmag),len(fitzpt),len(mjd)
    fitmag = np.array(fitmag)
    catmag = np.array(catmag)
    fitzpt = np.array(fitzpt)
    mjd = np.array(mjd)
    ww = catmag > -99.
    resid = fitmag - catmag + fitzpt
    #resid = resid[ww]
    #mjd = mjd[ww]

    # md, std, num = dt.iterstat(resid,startMedian=True, sigmaclip=1.5, iter=10)
    # plt.hist(resid,bins=np.arange(-.2025,.2,.005),label='Median:'+str(round(md,5))+'\nSTD: '+str(round(std,3)))
    # plt.xlim(-.2,.2)
    # plt.xlabel('Magnitude Residual From Zpt Fit')
    # plt.ylabel('Counts')
    # plt.title('CAT MAG > 20.')
    # plt.legend()
    # plt.savefig('zpttestgt20.png')
    # print 'saved zpttest.png'

    ww = abs(resid) < .5

    catmag = catmag[ww]
    resid = resid[ww]

    std = np.std(resid)
    ww = abs(resid) < 2*std

    catmag = catmag[ww]
    resid = resid[ww]

    plt.scatter(catmag,resid,alpha=.1)
    plt.plot([min(catmag),max(catmag)],[0,0],color='black')
    ax, ay, aystd = dt.bindata(catmag, resid,
                            np.arange(min(catmag), max(catmag), .1))
    plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')
    plt.xlim(min(catmag),max(catmag))
    plt.ylim(min(ay)-.01,max(ay)+.01)
    plt.xlabel('Catalog Mag')
    plt.ylabel('Fit Mag - (Cat Mag + Fit Zpt)')
    plt.savefig('zptresid.png')
    #print np.unique(mjd)
    #print np.unique(mjd[resid>.06])



if __name__ == '__main__':
    a = checkstars('/pnfs/des/scratch/pysmp/smp_02/lightcurves/des_fake_00211042_r.smp')
