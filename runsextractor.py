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

def getsky_and_skyerr(imagefilename,weightfilename,im,xlow,xhi,ylow,yhi,survey='DES',index='',bigreturn=False):
    print 'inside getsky and skyerr'
    if survey == 'DES':
        sexpath = "sex"
        #sexpath = "/global/common/cori/contrib/hpcosmo/hpcports_gnu-4.0/sextractor-2.19.5_b1614e12-4.0/bin/sex"
        fermigrid  = True
    if survey == 'PS1':
        sexpath = "/export/scratch0/ps1sn1/pipe/v10.0gpc1/photpipe/Cfiles/bin/linux/sex"
        fermigrid = False

    #im = pf.getdata(imagefilename)
    #hdr = pf.getheader(imagefilename)
    #im = im[ylow:yhi,xlow:xhi]
    if not os.path.exists('/global/cscratch1/sd/dbrout/sewpy_logs/'):
        os.makedirs('/global/cscratch1/sd/dbrout/sewpy_logs/')
    newfilename = '/global/cscratch1/sd/dbrout/sewpy_logs/'+index+'trimmed_'+imagefilename.split('/')[-1]
    print newfilename
    #dt.savefits(im, newfilename,fermigrid=fermigrid)
    dt.save_fits_image(im, newfilename,go=True)
    print 'saved trimmed file'
    logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.DEBUG)
    sew = sewpy.SEW(
            workdir='/global/cscratch1/sd/dbrout/sewpy_logs/'
            , sexpath=sexpath
            , loglevel="CRITICAL"
            , config={"WEIGHT_TYPE":"NONE,MAP_WEIGHT","WEIGHT_IMAGE":weightfilename,"checkimage_type":"BACKGROUND,BACKGROUND_RMS","checkimage_name":'/global/cscratch1/sd/dbrout/sewpy_logs/'+index+'_'+imagefilename.split('/')[-1]+
                                                                                      '.background, '+
                                                                                      '/global/cscratch1/sd/dbrout/sewpy_logs/'+index+'_'+imagefilename.split('/')[-1]+
                                                                                      '.background_rms'
                      ,"back_size":"256"}

        )
    out = sew(newfilename)
    path = out['logfilepath']
    log = open(path, 'r')
    background = -9
    rms = -9
    print 'getting background'
    for line in log.readlines():
        if 'Background:' in line.split(' '):
            background = line.split('Background: ')[1].split(' ')[0]
            rms = line.split('RMS: ')[1].split(' ')[0]

    if bigreturn:
        bkgrnd = pf.getdata(
            '/global/cscratch1/sd/dbrout/sewpy_logs/' + index + '_' + imagefilename.split('/')[-1] + '.background')
        bkgrndrms = pf.getdata(
            '/global/cscratch1/sd/dbrout/sewpy_logs/' + index + '_' + imagefilename.split('/')[-1] + '.background_rms')

    print 'removing files'
    try:
        os.remove(newfilename)
        os.remove(newfilename.split('.fits')[0]+'.cat.txt')
        os.remove(newfilename.split('.fits')[0]+'.log.txt')

        os.remove('/global/cscratch1/sd/dbrout/sewpy_logs/'+index+'_'+imagefilename.split('/')[-1]+'.background')
        os.remove('/global/cscratch1/sd/dbrout/sewpy_logs/'+index+'_'+imagefilename.split('/')[-1]+'.background_rms')
    except:
        pass

    print 'done'*100

    if bigreturn:
        return float(background), float(rms), bkgrnd, bkgrndrms
    #raw_input('sextractor stop')
    return float(background), float(rms)

#im = '/global/cscratch1/sd/dbrout/v3/20130902_SN-S2/r_21/SNp1_230168_SN-S2_tile20_r_21.fits'
#background, rms = getsky_and_skyerr(im)
#print 'bbb', b, 'rms', r