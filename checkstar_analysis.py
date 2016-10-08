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
    print data.keys()
    print data['ZPTFILE']
    zptfiles = data['ZPTFILE'][:-1]
    #mjd = []
    fitmag = []
    catmag = []
    fitzpt = []
    ra = []
    dec = []
    for z in zptfiles:
        zd = np.load(z)
        print z
        #print zd.keys()
        #mjd.append(zd['mjd'])
        fitmag.extend(zd['mpfit_mag'])
        catmag.extend(zd['cat_mag'])
        fitzpt.extend(zd['cat_mag']*0. + zd['mpfit_zpt'])
        #ra.extend()
        #dec.extend()
        #raw_input()
    print len(fitmag),len(catmag),len(fitzpt)
    fitmag = np.array(fitmag)
    catmag = np.array(catmag)
    fitzpt = np.array(fitzpt)
    resid = fitmag - catmag + fitzpt
    plt.hist(resid,bins=np.arange(-.1,.1,.005),normed=1)
    plt.savefig('zpttest.png')
    print 'saved zpttest.png'

if __name__ == '__main__':
    a = checkstars('/pnfs/des/scratch/pysmp/smp_02/lightcurves/des_fake_00212904_r.smp')
