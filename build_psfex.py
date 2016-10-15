import os
import numpy as np
import sys

def build(psffile, x, y, stampsize):

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


def buildall(psffile, x, y, stampsize):
    pstring = ''
    for p,xi,yi in zip(psffile,x,y):
        pstring += "dump_psfex -inFile_psf %s -xpix %s -ypix %s -gridSize %s; " % (p, xi, yi, stampsize)

    psf = os.popen(pstring).readlines()
    print psf
    sys.exit()
    ix, iy, psfval = [], [], []
    for line in psf:
        # print line
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