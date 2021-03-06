'''

Dillon Brout
dbrout@physics.upenn.edu


Python function to grab
sextractor sky and skyerr
values from any image

USAGE:
im = '/global/cscratch1/sd/dbrout/v3/20130902_SN-S2/r_21/SNp1_230168_SN-S2_tile20_r_21.fits'
background, rms = runsextractor.getsky_and_skyerr(im)

'''


import sewpy
import logging
import pyfits as pf
import dilltools as dt
import os
import numpy as np
import sys
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, NullLocator
from matplotlib.colors import LinearSegmentedColormap, colorConverter
from matplotlib.ticker import ScalarFormatter
from matplotlib import rcParams, rc
rcParams['font.family'] = 'serif'  # no serif on palomides?
rcParams['savefig.dpi'] = 300  # higher res outputs
rcParams['legend.numpoints'] = 1
rcParams['legend.markerscale'] = 0.7
rcParams['legend.handlelength'] = 0.5
rcParams['legend.handletextpad'] = 0.5
rcParams['legend.borderpad'] = 0.5
rcParams['legend.borderaxespad'] = 0.2
rcParams['legend.columnspacing'] = 1.0


def azimuthalAverage(image, center=None):
    """
    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is
             None, which then uses the center of the image (including
             fracitonal pixels).

    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if not center:
        center = np.array([(x.max() - x.min()) / 2.0, (x.max() - x.min()) / 2.0])

    r = np.hypot(x - center[0], y - center[1])

    # Get sorted radii
    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    # Get the integer part of the radii (bin size = 1)
    r_int = r_sorted.astype(int)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    rind = np.where(deltar)[0]  # location of changed radius
    nr = rind[1:] - rind[:-1]  # number of radius bin

    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    radial_prof = tbin / nr

    return radial_prof


def run(imagefilename,weightfilename,survey='DES',index='',bigreturn=False):
    print 'inside getsky and skyerr'
    if survey == 'DES':
        sexpath = "sex"
    if survey == 'PS1':
        sexpath = "/export/scratch0/ps1sn1/pipe/v10.0gpc1/photpipe/Cfiles/bin/linux/sex"

    newfilename = '/global/cscratch1/sd/dbrout/sewpy_logs/'+imagefilename.split('/')[-1]

    # im = pf.getdata(imagefilename)
    # dt.save_fits_image(im, newfilename,go=True)
    logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.DEBUG)
    sew = sewpy.SEW(
            workdir='/global/cscratch1/sd/dbrout/sewpy_logs/'
            , sexpath=sexpath
            , loglevel="CRITICAL"
            , params = ["X_IMAGE", "Y_IMAGE","ISOAREA_IMAGE"]
            , config={"WEIGHT_TYPE":"NONE,MAP_WEIGHT","WEIGHT_IMAGE":weightfilename
                      ,"back_size":"256"
                      ,"catalog":"test.cat"
                      ,"DETECT MINAREA":"4"
                      ,"THRESH_TYPE":"RELATIVE"
                      ,"DETECT_THRESH":"1.1"
                      }

        )
    out = sew(imagefilename)
    path = out['logfilepath']
    log = open(path, 'r')
    background = -9
    rms = -9
    print 'running S-Extractor'
    for line in log.readlines():
        print line
    print '-'*100
    print out["table"]

    im = pf.getdata(imagefilename)
    wgt = pf.getdata(weightfilename)
    wgt[wgt<1e-5] = np.nan
    wgt[wgt>0] = 1.
    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    ax.imshow(np.log10(im*wgt), cmap="Greys")
    from matplotlib.patches import Ellipse
    ells = [Ellipse(xy=(x,y),
                    width=xa ,height=ya,
                    angle=ang)
            for x,y,xa,ya,ang in zip(out["table"]['XWIN_IMAGE'],out["table"]['YWIN_IMAGE'],
                                 out["table"]['AWIN_IMAGE']*np.log10(out["table"]['FLUX_AUTO']+10)*4.5+2,
                                 out["table"]['BWIN_IMAGE']*np.log10(out["table"]['FLUX_AUTO']+10)*4.5+2,
                                 out["table"]['THETAWIN_IMAGE'])]

    for e in ells:
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
        #e.set_alpha(np.random.rand())
        e.set_facecolor('none')
        e.set_edgecolor('red')
        e.set_linewidth(.2)

    #ax.set_ylim(1000,4000)
    #asdf

    plt.savefig('testext.png',dpi=1000)
    #os.popen('upload testext.png')

    import skimage.draw
    for x, y, xa, ya, ang in zip(out["table"]['XWIN_IMAGE'], out["table"]['YWIN_IMAGE'],
                                         out["table"]['AWIN_IMAGE'] * np.log10(out["table"]['FLUX_AUTO']+10) * 3. + 2,
                                         out["table"]['BWIN_IMAGE'] * np.log10(out["table"]['FLUX_AUTO']+10) * 3. + 2,
                                         out["table"]['THETAWIN_IMAGE']):

        rr, cc = skimage.draw.ellipse(y, x, ya,xa,shape=im.shape, rotation=2*3.14-np.deg2rad(ang))
        im[rr, cc] = np.nan
    plt.clf()
    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    ax.imshow(np.log10(im * wgt), cmap="Greys")
    plt.savefig('testextmask.png',dpi=1000)
    #os.popen('upload testextmask.png')



    # from scipy import fftpack
    # F1 = fftpack.fft(im[np.isfinite(im)].ravel().astype(float))
    # F2 = fftpack.fftshift(F1)
    # #psd2D = np.abs(F2) ** 2
    # #psd1D = azimuthalAverage(psd2D)
    # #print psd2D
    # #print psd1D
    # print F1
    # print F2
    # plt.clf()
    # plt.semilogy(F2)
    # plt.xlabel('Spatial Frequency')
    # plt.ylabel('Power')
    # #plt.ylim(-.2e8,.2e8)
    # plt.savefig('plots/noisepower.png', dpi=100)
    # print os.popen('source ~/.bash_profile.ext; upload plots/noisepower.png').read()
    # sys.exit()
    nx,ny = im.shape[0],im.shape[1]

    groupings = [1.,2.,3.,4.,6.,8.,16.,32.,64.]
    resultsdict = {}
    for g in groupings:
        print 'calculating grouping',g
        resultsdict[g] = []
        for x in np.arange(0,nx-64,g):
            if x%100 == 0: print x
        #for x in np.arange(0, 512, g):
            #for y in np.arange(0, 512, g):
            for y in np.arange(0,ny-64,g):
                resultsdict[g].append(np.mean(im[int(x):int(x+g),int(y):int(y+g)]))

        resultsdict[g] = np.array(resultsdict[g])
        print ''
    plt.clf()
    for g in groupings:
        hist, bin_edges = np.histogram(resultsdict[g][np.isfinite(resultsdict[g])], bins=np.arange(-705.,700,10), normed=True)
        bin_centers = (bin_edges[1:] + bin_edges[:-1])/2.
        #print bin_centers
        #raw_input()
        plt.plot(bin_centers,hist,label='Group %d'%g,linewidth=3.,alpha=.75)
        #plt.hist(resultsdict[g][np.isfinite(resultsdict[g])], bins=np.arange(-505,500,10),
        #         type='step', label='Group %d'%g)
    plt.legend()
    plt.xlim(-300,300)
    plt.xlabel('Grouped Pixel Values')
    plt.savefig('plots/correlatednoise.png',dpi=100)
    print os.popen('source ~/.bash_profile.ext; upload plots/correlatednoise.png').read()

    plt.clf()
    stds = []
    for g in groupings:
        stds.append(np.std(resultsdict[g][np.isfinite(resultsdict[g])]))
    stds = np.array(stds)
    groupings = np.array(groupings)
    #plt.plot(np.log(groupings),stds,label='DES IMAGE',linewidth=3)
    #plt.plot(np.log(groupings),stds[0]/np.sqrt(groupings**2),label='1/sqrt(N^2)',linewidth=3)
    plt.plot(np.log(groupings),stds/(stds[0]/np.sqrt(groupings**2)),label='Ratio',linewidth=3)
    plt.grid(True)
    plt.legend()
    plt.xticks(np.log(groupings),groupings)
    plt.xlabel('Pixel Groupings')
    #lt.xlim(-1,7)
    plt.ylabel('Ratio of DES STD / ROOT(Nside^2)')
    plt.savefig('plots/correlatednoisestds.png', dpi=100)
    print os.popen('source ~/.bash_profile.ext; upload plots/correlatednoisestds.png').read()

    print stds
    print stds[0]/np.sqrt(groupings**2)
    return

im = '/global/cscratch1/sd/masao/diffim/output/FPH_V8/20151008_SN-C3/z_05/SNY3_483208_SN-C3_tile81_z_05.fits'
weight = '/global/cscratch1/sd/masao/diffim/output/FPH_V8/20151008_SN-C3/z_05/SNY3_483208_SN-C3_tile81_z_05.weight.fits'
run(im,weight)
#print 'bbb', b, 'rms', r