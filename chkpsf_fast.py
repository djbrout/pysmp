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
    residsig = (im - model) / sigma * fitrad
    return np.array(residsig.ravel())

def simstamp(param, psf, im, sigma, fitrad, sky, psfmag):
    model = psf * param / 10 ** (-0.4 * (psfmag - 25)) + sky
    return model
def fit( fileroot,im,mask, impsf,fullpsf,imhdr,hpsf,fitrad,
        xpos=None, ypos=None, radius=4, pdf_pages=None, ra=None, dec=None, title='', returnstamps = False, maskfile=None):

    good = False

    fluxerr = 100.
    chisq = 1.
    dms = 1.

    x = xpos
    y = ypos

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
    try:
        model = dao_value.dao_value(dx, dy, gauss,
                                    impsf,  # psf1d=psf1d,
                                    deriv=False)  # ,ps1d=False)
    except:
        return 1, 1, 0, 0, False, 0, 0, 0


    # subim = im[iylo - 1:iyhi, ixlo - 1:ixhi]
    # #print 'modelshape', model.shape, 'imshape', subim.shape
    # #raw_input()
    # submask = mask[iylo - 1:iyhi, ixlo - 1:ixhi]
    # submask[submask != 0] = 9
    # submask[submask == 0 ] = 1
    # submask[submask == 9 ] = 0
    # # scaledpsf = model+impsf[psfy/2+1-radius:psfy/2+1+radius+1,
    #                        psfx/2+1-radius:psfx/2+1+radius+1]
    # print model.shape
    # print flux.shape
    # print hpsf['PSFMAG']
    # print imhdr['SKYADU']
    chisqvec = []
    fluxvec = []
    #substamp = model.shape[0]

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
    # fluxls, cov = opti.leastsq(resid, 100000,
    #                            args=(model, subim, imhdr['SKYSIG'], fitrad, imhdr['SKYADU'], hpsf['PSFMAG']),full_output=False)
    #print cov.shape
    #print fluxls, cov
    #raw_input('covshape')
    # print 'flux fit comparo',flux,fluxls,
    # scaledpsf = model*fluxls/10**(-0.4*(hpsf['PSFMAG']-25)) + imhdr['SKYADU']





    return model
