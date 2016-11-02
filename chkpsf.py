# !/usr/bin/env python
# D. Jones - 6/5/16

from txtobj import txtobj
import numpy as np
import pyfits
import rdpsf
import dao_value
import rebin

import matplotlib as m

m.use('Agg')
import matplotlib.pyplot as plt

import pyfits
import os
from astropy import wcs

import scipy.optimize as opti


def resid(param, psf, im, sigma, fitrad, sky, psfmag):
    model = psf * param / 10 ** (-0.4 * (psfmag - 25)) + sky
    residsig = (im - model) / sigma
    return np.array(residsig.ravel())

def simstamp(param, psf, im, sigma, fitrad, sky, psfmag):
    model = psf * param / 10 ** (-0.4 * (psfmag - 25)) + sky
    return model
def fit(
        fileroot='/export/scratch0/ps1sn1/data/v10.0/GPC1v3/eventsv1/workspace/PSc560121/g/PSc560121.md01s043.g.ut090831e.1917665_14.sw',
        xpos=None, ypos=None, radius=8, pdf_pages=None, ra=None, dec=None, title='', returnstamps = False, maskfile=None):
    # xpos = xpos +1
    # ypos = ypos +1
    # from matplotlib.backends.backend_pdf import PdfPages
    # pdf_pages = PdfPages('daophot_resid.pdf')
    p = pyfits.open('%s.fcmp' % fileroot)
    p.verify("fix")

    if os.path.exists('test.fcmp'):
        os.remove('test.fcmp')
    p.writeto('test.fcmp', output_verify='fix')
    # fcmp = p[0].header
    # print p[1]

    fcmp = txtobj('test.fcmp', cmpheader=True)
    # print fcmp.__dict__['class']
    # print fcmp['class']
    # raw_input()
    im = pyfits.getdata('%s.fits' % fileroot)
    mask = pyfits.getdata(maskfile)

    w = wcs.WCS('%s.fits' % fileroot)
    #results2 = w.wcs_world2pix(np.array([[ra, dec]]), 0)
    # xpos,ypos =results2[0][0], results2[0][1]

    imhdr = pyfits.getheader('%s.fits' % fileroot)
    fullpsf, hpsf = rdpsf.rdpsf('%s.dao.psf.fits' % fileroot)
    impsf = pyfits.getdata('%s.dao.psf.fits' % fileroot)

    psfsize = np.shape(impsf)[0]

    fcmp.Xpos = fcmp.Xpos[1:].astype(float)
    fcmp.Ypos = fcmp.Ypos[1:].astype(float)
    fcmp.__dict__['class'] = fcmp.__dict__['class'][1:].astype(float)

    fcmp.flux = fcmp.flux[1:].astype(float)
    fcmp.dflux = fcmp.dflux[1:].astype(float)
    # for x,y,flux,fluxerr in zip(fcmp.Xpos,fcmp.Ypos,
    #                            fcmp.flux,fcmp.dflux):

    x = xpos
    y = ypos
    # print fcmp.Xpos-xpos
    # print fcmp.Ypos-ypos
    # raw_input()
    #print fcmp.__dict__['class']
    #print fcmp.Xpos
    #print xpos
    #raw_input()
    ww = (abs(fcmp.Xpos - xpos) < 1.) & (abs(fcmp.Ypos - ypos) < 1.)
    thisclass = fcmp.__dict__['class'][ww]
    #print 'THIS CLASS IS', thisclass
    #print 'all classes', fcmp.__dict__['class']
    # flux = fcmp.flux
    # fluxerr = fcmp.dflux

    ny, nx = np.shape(im)
    psfy, psfx = np.shape(impsf)
    ixlo, iylo = int(x - radius), int(y - radius)
    if ixlo < 0: ixlo = 0
    if iylo < 0: iylo = 0
    ixhi = int(x + radius) + 1
    iyhi = int(y + radius) + 1
    if ixhi > (nx - 1): ixhi = nx - 1
    if iyhi > (ny - 1): iyhi = ny - 1
    ixx = ixhi - ixlo + 1
    iyy = iyhi - iylo + 1
    dx = np.arange(ixx) + ixlo - x
    dy = np.arange(iyy) + iylo - y
    psf1d = impsf.reshape(np.shape(impsf)[0] ** 2.)
    gauss = [hpsf['GAUSS1'], hpsf['GAUSS2'], hpsf['GAUSS3'],
             hpsf['GAUSS4'], hpsf['GAUSS5']]
    dx = dx.reshape(1, len(dx))
    dy = dy.reshape(len(dy), 1)
    dx = rebin.rebin(dx, [np.shape(dx)[1], np.shape(dx)[1]])
    dy = rebin.rebin(dy, [len(dy), len(dy)])
    model = dao_value.dao_value(dx, dy, gauss,
                                impsf,  # psf1d=psf1d,
                                deriv=False)  # ,ps1d=False)
    subim = im[iylo - 1:iyhi, ixlo - 1:ixhi]
    submask = mask[iylo - 1:iyhi, ixlo - 1:ixhi]
    submask[submask != 0] = 9
    submask[submask == 0 ] = 1
    submask[submask == 9 ] = 0
    # scaledpsf = model+impsf[psfy/2+1-radius:psfy/2+1+radius+1,
    #                        psfx/2+1-radius:psfx/2+1+radius+1]
    # print model.shape
    # print flux.shape
    # print hpsf['PSFMAG']
    # print imhdr['SKYADU']
    chisqvec = []
    fluxvec = []
    substamp = model.shape[0]
    fitrad = np.zeros([substamp, substamp])
    radius = 4
    for x in np.arange(substamp):
        for y in np.arange(substamp):
            if np.sqrt((substamp / 2. - x) ** 2 + (substamp / 2. - y) ** 2) < radius:
                fitrad[int(x), int(y)] = 1.
    '''
    for flux in range(1,500000,200):
        scaledpsf = model*flux/10**(-0.4*(hpsf['PSFMAG']-25)) + imhdr['SKYADU']
        chisq = np.sum(fitrad*(subim-scaledpsf)**2/imhdr['SKYSIG']**2)
        chisqvec.append(chisq)
        fluxvec.append(flux)

    chisqvec = np.array(chisqvec)
    fluxvec = np.array(fluxvec)
    flux = fluxvec[np.argmin(chisqvec)]
    scaledpsf = model*flux/10**(-0.4*(hpsf['PSFMAG']-25)) + imhdr['SKYADU']

    #resid(param,psf,im,sigma,fitrad,sky,psfmag)

    #print model, subim, imhdr['SKYSIG']
    '''
    fluxls, cov = opti.leastsq(resid, 100000,
                               args=(model, subim, imhdr['SKYSIG'], fitrad, imhdr['SKYADU'], hpsf['PSFMAG']),full_output=False)

    # print 'flux fit comparo',flux,fluxls,
    scaledpsf = model*fluxls/10**(-0.4*(hpsf['PSFMAG']-25)) + imhdr['SKYADU']

    fluxerr = 100.
    chisq = 1.
    dms = 1.
    good = False
    if len(thisclass) == 1:
        if thisclass[0] == 1:
            good = True

    if not pdf_pages is None:
        fig = plt.figure()
        plt.clf()
        axim = plt.subplot(131)
        axpsf = plt.subplot(132)
        axdiff = plt.subplot(133)
        for ax,title in zip([axim,axpsf,axdiff],['image','model','difference']):
            ax.set_title(title)
            axim.imshow(subim,
            cmap='gray',interpolation='nearest')
            axpsf.imshow(model,cmap='gray',interpolation='nearest')
            axdiff.imshow(subim-scaledpsf,cmap='gray',interpolation='nearest')
            #plt.colorbar()
            axim = plt.subplot(131)
            axpsf = plt.subplot(132)
            axdiff = plt.subplot(133)
        #for ax,title in zip([axim,axpsf,axdiff],['image','model','difference']):
            if good:
                ax.set_title(title + 'GOOD')
            else:
                ax.set_title(title + 'BADD')
            #ax.set_title(title)
        axim.imshow(subim,cmap='gray',interpolation='nearest')
        axpsf.imshow(scaledpsf,cmap='gray',interpolation='nearest')
        ax = axdiff.imshow(subim-scaledpsf,cmap='gray',interpolation='nearest')
        cbar = fig.colorbar(ax)
        #plt.imshow((subim-scaledpsf)/imhdr['SKYSIG'],cmap='gray',interpolation='nearest')
        #plt.colorbar()
        if good:
            plt.title(title + 'GOOD' )
        else:
            plt.title(title + 'BADD' )
        pdf_pages.savefig(fig)

        #pdf_pages.close()
        #plt.savefig('')




    if returnstamps:
        rpsf = model
        good = True#we know this wont be in the starcat file so set to good is true
        print 'fluxls', fluxls
        print 'maxpsf', np.max(rpsf)
        return fluxls,fluxerr,chisq,dms,good,subim, rpsf, imhdr['SKYSIG'], fitrad, imhdr['SKYADU'], hpsf['PSFMAG'], submask
    print fluxls
    print np.max(model)
    #raw_input('fluxls')
    sstamp = simstamp(fluxls,model, subim, imhdr['SKYSIG'], fitrad, imhdr['SKYADU'], hpsf['PSFMAG'])

    return fluxls, fluxerr, chisq, dms, good, subim, sstamp
