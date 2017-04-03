import os
from copy import copy
import numpy as np

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
    return
def build(psffile, x, y, stampsize):


    print "dump_psfex -inFile_psf %s -xpix %s -ypix %s -gridSize %s" % (psffile, x, y,
                                                                        stampsize)

    psf = os.popen("dump_psfex -inFile_psf %s -xpix %s -ypix %s -gridSize %s" % (psffile, x, y,
                                                                                 stampsize)).readlines()

    ix, iy, psfval = [], [], []
    for line in psf:
        #print line
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

    ix, iy, psfval = np.array(ix), np.array(iy), np.array(psfval)
    psfout = np.zeros((stampsize, stampsize))
    for x, y, p in zip(ix, iy, psfval):
        psfout[y, x] = p

    return (psfout), (IMAGE_CENTERX, IMAGE_CENTERY)

