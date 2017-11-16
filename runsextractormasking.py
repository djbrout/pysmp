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
            , params = ["X_IMAGE", "Y_IMAGE", "FLUX_APER(3)","ISOCOR", "FLAGS"]
            , config={"WEIGHT_TYPE":"NONE,MAP_WEIGHT","WEIGHT_IMAGE":weightfilename
                      # "checkimage_type":"BACKGROUND,BACKGROUND_RMS",
                      # "checkimage_name":'/global/cscratch1/sd/dbrout/sewpy_logs/'+index+'_'+imagefilename.split('/')[-1]+
                      #       '.background, '+
                      #       '/global/cscratch1/sd/dbrout/sewpy_logs/'+index+'_'+imagefilename.split('/')[-1]+
                      #       '.background_rms'
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

    return

im = '/global/cscratch1/sd/masao/diffim/output/FPH_V8/20151008_SN-C3/z_05/SNY3_483208_SN-C3_tile81_z_05.fits'
weight = '/global/cscratch1/sd/masao/diffim/output/FPH_V8/20151008_SN-C3/z_05/SNY3_483208_SN-C3_tile81_z_05.weight.fits'
run(im,weight)
#print 'bbb', b, 'rms', r