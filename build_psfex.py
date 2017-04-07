import os
import numpy as np
import sys
try:
    import psfex
except:
    pass

def build(psffile, x, y, stampsize,psfexworked=True):
    #print 'psffile'
    #raw_input()
    try:
        a = psfex.PSFEx(psffile)
        im = a.get_rec(y, x)[3:-4, 3:-4]
        im /= np.sum(im.ravel())
    except:

        print "dump_psfex -inFile_psf %s -xpix %s -ypix %s -gridSize %s" % (psffile, x, y,
                                                                            stampsize)

        psf = os.popen("dump_psfex -inFile_psf %s -xpix %s -ypix %s -gridSize %s" % (psffile, x, y,
                                                                                     stampsize)).readlines()

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
        im = np.zeros((stampsize, stampsize))
        for x, y, p in zip(ix, iy, psfval):
            im[y, x] = p


    return im, (round(x), round(y))



def buildall(psffile, x, y, stampsize):
    #print psffile
    #raw_input()
    pstring = ''
    badvec = []
    for i,p,xi,yi in zip(np.arange(len(x)),psffile,x,y):
        #if i > 5: continue
        #print p
        if 'SN-' in p:
            pstring += "dump_psfex -inFile_psf "+p+" -xpix "+str(xi)+" -ypix "+str(yi)+" -gridSize "+str(stampsize)+"; "
        else:
            badvec.append(i)
    #print pstring
    psf = os.popen(pstring).readlines()
    #print psf
    #sys.exit()
    badvec = np.array(badvec)
    psfouts = []
    imagecenterouts = []
    keepgoing = True
    isnewpsf = True
    l = -1
    vec = 0
    for line in psf:
        if not keepgoing: continue
        l+=1
        if isnewpsf:
            isnewpsf = False
            ix, iy, psfval = [], [], []

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
            isnewpsf = True
            ix, iy, psfval = np.array(ix), np.array(iy), np.array(psfval)
            psfout = np.zeros((stampsize, stampsize))
            for x, y, p in zip(ix, iy, psfval):
                psfout[y, x] = p
            psfouts.append(psfout)
            imagecenterouts.append([IMAGE_CENTERX, IMAGE_CENTERY])
            vec += 1
            if vec in badvec:
                psfouts.append(psfout*0.)
                imagecenterouts.append([0,0])

    return psfouts, imagecenterouts