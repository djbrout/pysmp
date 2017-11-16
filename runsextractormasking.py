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
                      ,
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
    wgt[wgt<1e-6] = 0.
    wgt[wgt>0] = 1.
    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    import numpy as np
    ax.imshow(np.log10(im*wgt), cmap="Greys")
    from matplotlib.patches import Ellipse
    ells = [Ellipse(xy=(x,y),
                    width=xa ,height=ya,
                    angle=ang)
            for x,y,xa,ya,ang in zip(out["table"]['XWIN_IMAGE'],out["table"]['YWIN_IMAGE'],
                                 out["table"]['AWIN_IMAGE']*np.log10(out["table"]['FLUX_AUTO'])*1.5,
                                 out["table"]['BWIN_IMAGE']*np.log10(out["table"]['FLUX_AUTO'])*1.5,
                                 out["table"]['THETAWIN_IMAGE'])]

    for e in ells:
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
        #e.set_alpha(np.random.rand())
        e.set_facecolor('none')
        e.set_edgecolor('red')
        e.set_linewidth(1.)
    #ax.set_ylim(1000,4000)

    plt.savefig('testext.png',dpi=1000)

    return

im = '/global/cscratch1/sd/masao/diffim/output/FPH_V8/20151008_SN-C3/z_05/SNY3_483208_SN-C3_tile81_z_05.fits'
weight = '/global/cscratch1/sd/masao/diffim/output/FPH_V8/20151008_SN-C3/z_05/SNY3_483208_SN-C3_tile81_z_05.weight.fits'
run(im,weight)
#print 'bbb', b, 'rms', r