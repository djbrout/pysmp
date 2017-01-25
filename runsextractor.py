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

def getsky_and_skyerr(imagefilename,im,xlow,xhi,ylow,yhi,survey='DES'):
    if survey == 'DES':
        sexpath = "sex"
        fermigrid  = True
    if survey == 'PS1':
        sexpath = "/export/scratch0/ps1sn1/pipe/v10.0gpc1/photpipe/Cfiles/bin/linux/sex"
        fermigrid = False

    #im = pf.getdata(imagefilename)
    #hdr = pf.getheader(imagefilename)
    im = im[ylow:yhi,xlow:xhi]
    if not os.path.exists('sewpy_logs/'):
        os.makedirs('sewpy_logs/')
    newfilename = 'sewpy_logs/trimmed_'+imagefilename.split('/')[-1]
    #dt.savefits(im, newfilename,fermigrid=fermigrid)
    dt.save_fits_image(im, newfilename,go=True)

    logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.DEBUG)
    sew = sewpy.SEW(
            workdir='sewpy_logs/'
            , sexpath=sexpath
            , loglevel="CRITICAL"
            , config={"checkimage_type":"BACKGROUND,BACKGROUND_RMS","checkimage_name":imagefilename+'.background, '+imagefilename+'.background_rms'
                      ,"back_size":"256"}

        )
    out = sew(imagefilename)
    path = out['logfilepath']
    log = open(path, 'r')
    background = -9
    rms = -9
    for line in log.readlines():
        if 'Background:' in line.split(' '):
            background = line.split('Background: ')[1].split(' ')[0]
            rms = line.split('RMS: ')[1].split(' ')[0]

    os.remove(newfilename)
    return float(background), float(rms)

#im = '/global/cscratch1/sd/dbrout/v3/20130902_SN-S2/r_21/SNp1_230168_SN-S2_tile20_r_21.fits'
#background, rms = getsky_and_skyerr(im)
#print 'bbb', b, 'rms', r