# !/usr/bin/env python
# D. Jones - 6/5/16

import numpy as np
import dao_value
import rebin

def fit(im, impsf,hpsf,xpos=None, ypos=None, radius=32):

    x = xpos
    y = ypos

    ny, nx = np.shape(im)
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
    gauss = [hpsf['GAUSS1'], hpsf['GAUSS2'], hpsf['GAUSS3'],
             hpsf['GAUSS4'], hpsf['GAUSS5']]
    dx = dx.reshape(1, len(dx))
    dy = dy.reshape(len(dy), 1)
    dx = rebin.rebin(dx, [np.shape(dx)[1], np.shape(dx)[1]])
    dy = rebin.rebin(dy, [len(dy), len(dy)])
    try:
        model = dao_value.dao_value(dx, dy, gauss, impsf, deriv=False)
    except:
        return 1, 1, 0, 0, False, 0, 0, 0

    subim = im[iylo - 1:iyhi, ixlo - 1:ixhi]

    model = model/ 10 ** (-0.4 * (hpsf['PSFMAG'] - 25))

    return model,subim
