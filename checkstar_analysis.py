import os
#import matplotlib as m

#m.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(threshold=np.nan)
import sys
import dilltools as dt

default_checkstar_file = '/Volumes/ExtraSpace/pysmp_downloads/des_fake_00229567_r_standardstarfits.txt'


def checkstars(checkstarfile=default_checkstar_file):
    cols = dt.readcol(checkstarfile,delim='\t')
    print cols.keys()
    plt.scatter(-2.5*np.log10(cols['Fit Flux']),(cols['Fit Flux']-cols['Galsim Fit Flux'])/cols['Fit Flux'],
                alpha=.1,color='black')
    xvals, medians, mads = dt.bindata(-2.5*np.log10(cols['Fit Flux']),(cols['Fit Flux']-cols['Galsim Fit Flux'])/cols['Fit Flux'],
                                      np.arange(-14,-9,.25))
    plt.errorbar(xvals,medians,mads,color='blue',fmt='o')
    plt.axhline(0,color='blue')
    plt.xlabel('-2.5*np.log10(fitflux)')
    plt.ylabel('myPSF Flux - GalsimPSF Flux / myPSF Flux')
    plt.title('Star Fits, 1CCD, All Epochs, All Stars')
    plt.savefig('/Volumes/ExtraSpace/pysmp_downloads/starfit_resids.png')
    #print cols.keys()


    plt.clf()
    plt.hist(cols['Fit Flux Chisq'], label='my PSF Model', bins=np.arange(.5,2,.05),alpha=.75)
    plt.hist(cols['Galsim Fit Flux Chisq'], label='Galsim PSF Model', bins=np.arange(.5,2,.05),alpha=.75)
    plt.xlabel('Reduced Chisq')
    plt.legend()
    plt.savefig('/Volumes/ExtraSpace/pysmp_downloads/starfit_chisqhist.png')

    plt.clf()
    plt.hist(cols['Fit Flux DMS'], label='my PSF Model', bins=np.arange(-10000,5000,1000),alpha=.75)
    plt.hist(cols['Galsim Fit Flux DMS'], label='Galsim PSF Model', bins=np.arange(-10000,5000,1000),alpha=.75)
    plt.xlabel('Data - Sim')
    plt.legend()
    plt.savefig('/Volumes/ExtraSpace/pysmp_downloads/starfit_dmshist.png')

    plt.clf()
    plt.scatter(cols['Galsim Fit Flux Chisq'], (cols['Fit Flux'] - cols['Galsim Fit Flux']) / cols['Fit Flux'],
                alpha=.1, color='black')


    xvals, medians, mads = dt.bindata(cols['Galsim Fit Flux Chisq'],
                                      (cols['Fit Flux'] - cols['Galsim Fit Flux']) / cols['Fit Flux'],
                                      np.arange(.5, 3, .1))
    plt.errorbar(xvals, medians, mads, color='blue', fmt='o')
    plt.axhline(0, color='blue')
    plt.xlim(.5,3)
    plt.xlabel('Galsim Chisq')
    plt.ylabel('myPSF Flux - GalsimPSF Flux / myPSF Flux')
    plt.title('Star Fits, 1CCD, All Epochs, All Stars')
    plt.savefig('/Volumes/ExtraSpace/pysmp_downloads/starfit_vs_galsimchisq.png')

    return cols




if __name__ == '__main__':
    a = checkstars()
