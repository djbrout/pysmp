import os
from copy import copy
thereisgalsim = True
try:
    import galsim
    import galsim.des
except ImportError:
    thereisgalsim = False

def build_psfex(psffile, x, y, imfile, fitrad, dogalsim=False):
    '''
    Inputs from dump_psfex output file:

    PSF: Xcoord, Ycoord, dx, dy, psfval

    Returns psfout such that psfout[Ycoord+5, Xcoord+5] = psfval

    e.g. if a line reads "PSF: 17 18 1.266 0.341 1.649823e-02"

    print psfout[23, 22] will return .01649823


    PSF_CENTER: 17 17        # PSF coords
    PSF_MAX:    17 18        # PSF coords
    IMAGE_CENTER: 554 3914   # IMAGE coords (PSF_CENTER here)
    IMAGE_CORNER: 1 1      # pixel index at corner

    Center of PSF in image coordinates is 553 3913
    This corresponds to psfout[22,22]


    '''
    psf = os.popen("dump_psfex -inFile_psf %s -xpix %s -ypix %s -gridSize %s" % (psffile, x, y,
                                                                                 35)).readlines()
    xin = copy(x)
    yin = copy(y)
    IMAGE_CENTERX = -9
    IMAGE_CENTERY = -9
    ix, iy, psfval = [], [], []
    for line in psf:
        line = line.replace('\n', '')
        if line.startswith('PSF:'):
            linelist = line.split()
            ix += [int(linelist[1])];
            iy += [int(linelist[2])];
            psfval += [float(linelist[5])]
        elif line.startswith("IMAGE_CENTER"):
            linelist = line.split()
            IMAGE_CENTERX = float(linelist[1]);
            IMAGE_CENTERY = float(linelist[2])

    # IMAGE_CENTERX -= IMAGE_CORNERX; IMAGE_CENTERY -= IMAGE_CORNERY
    ix, iy, psfval = np.array(ix), np.array(iy), np.array(psfval)
    psfout = np.zeros((2 * fitrad + 1, 2 * fitrad + 1))
    for x, y, p in zip(ix, iy, psfval):
        if x >= (35 - 2 * fitrad - 1) / 2 and y >= (35 - 2 * fitrad - 1) / 2 and x < (2 * fitrad + 1) and y < (2 * fitrad + 1):
            psfout[y - (35 - 2 * fitrad - 1) / 2, x - (35 - 2 * fitrad - 1) / 2] = p

    if dogalsim and thereisgalsim and IMAGE_CENTERX != -9.:
        #worldpsf = wcs.toWorld(gs_psfex,image_pos=psfcenter)  # convert coords
        wcs = galsim.FitsWCS(imfile)  # read in wcs
        psfstamp = galsim.ImageF(array=copy(psfout * 0.))  # make emtpy stamp
        des_psfex = galsim.des.DES_PSFEx(psffile)  # read in psf file
        psfcenter = galsim.PositionD(xin, yin)  # set galsim position
        gs_psfex = des_psfex.getPSF(psfcenter)  # get psf at position
        off = galsim.PositionD(-IMAGE_CENTERX + xin,
                               -IMAGE_CENTERY + yin)  # set subpixel offset
        gs_psfex.drawImage(image=psfstamp, offset=off)  # draw psf to stamp at offset
        psfout = psfstamp.array

    return (psfout), (IMAGE_CENTERX, IMAGE_CENTERY)
