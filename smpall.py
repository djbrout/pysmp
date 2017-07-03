#!/usr/bin/env python

'''

Scene Modeling Photometric Pipeline

Dillon Brout
6/16/2016
dbrout@physics.upenn.edu

'''


import numpy as np
import exceptions
import os
import sys
sys.path.append("./mpfit/")
import mpfitexpr
sys.path.append("/global/homes/d/dbrout/GalSim-1.3.0")
sys.path.append("/global/homes/d/dbrout/GalSim-1.3.0/lib")
from iminuit import Minuit
from scipy.odr import *

#import scipy.ndimage
import matplotlib as m
import mcmc as mcmc3
#import mcmc3galsimpixshift as mcmc3galsimpixshift
m.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import pyfits as pf
import scipy.signal
from copy import copy
try:
    import galsim
    import galsim.des
except:
    print 'Could not import galsim, galsim features not functional...'
import time
from astropy.io import fits
from scipy.interpolate import UnivariateSpline
import sigma_clip
import meanclip
import cntrd,aper,getpsf,rdpsf
import runsextractor
import pkfit_norecent_noise_smp
import dilltools as dt
import chkpsf

inflateskyerrworked = False
try:
    import inflateskyerr
    inflateskyerrworked = True
except:
    pass
from lmfit import Minimizer, Parameters

import scipy.optimize as opti
psfexworked = False
try:
    import psfex
    psfexworked = True
except:
    pass
snkeywordlist = {'SURVEY':'string','SNID':'string','FILTERS':'string',
                 'PIXSIZE':'float','NXPIX':'float','NYPIX':'float',
                 'ZPFLUX':'float','RA':'string',
                 'DECL':'string','PEAKMJD':'float','WEIGHT_BADPIXEL':'string',
                 'HOSTGAL_SB_FLUXCAL':[],
                 'STARCAT':'string', 'PSF_UNIT':'string', 'NOBS':'float'}
snkeywordlist_fake = {'SURVEY':'string','SNID':'string','FILTERS':'string',
                 'PIXSIZE':'float','NXPIX':'float','NYPIX':'float',
                 'ZPFLUX':'float','RA':'string', 'FAKE_RA':'string','FAKE_DEC':'string',
                 'DECL':'string','PEAKMJD':'float','WEIGHT_BADPIXEL':'string',
                 'HOSTGAL_SB_FLUXCAL':[],
                 'STARCAT':'string', 'PSF_UNIT':'string', 'NOBS':'float'}
snvarnameslist_fake = {'ID_OBS': 'string','ID_COADD': 'string','MJD':'float','BAND':'string',
                  'IMAGE_NAME_SEARCH':'string','IMAGE_NAME_WEIGHT':'string',
                  'FILE_NAME_PSF':'string','FAKE_TRUEMAG':'float','ZP':'float',
                  'FLUX':'float','FLUXERR':'float','PHOTFLAG':'string','SKYSIG':'float'}
snvarnameslist = {'ID_OBS': 'string','ID_COADD': 'string','MJD':'float','BAND':'string',
                  'IMAGE_NAME_SEARCH':'string','IMAGE_NAME_WEIGHT':'string',
                  'FILE_NAME_PSF':'string','ZP':'float',
                  'FLUX':'float','FLUXERR':'float','PHOTFLAG':'string','SKYSIG':'float'}

snvarnameslist_ps = {'ID_OBS': 'string', 'ID_COADD': 'string', 'MJD': 'float', 'BAND': 'string',
                  'IMAGE_NAME_SEARCH': 'string', 'IMAGE_NAME_MASK': 'string', 'IMAGE_NAME_NOISE': 'string',
                  'FILE_NAME_PSF': 'string', 'FAKE_TRUEMAG': 'float', 'ZP': 'float',
                  'FLUX': 'float', 'FLUXERR': 'float', 'PHOTFLAG': 'string','SKYSIG':'float'}

paramkeywordlist = {'STAMPSIZE':'float ','RADIUS1':'float',
                    'RADIUS2':'float','SUBSTAMP':'float',
                    'MAX_MASKNUM':'float','RDNOISE_NAME':'string',
                    'GAIN_NAME':'string','FWHM_MAX':'float',
                    'PSF_MAX':'float','WEIGHT_TYPE':'string',
                    'MASK_TYPE':'string','MJDPLUS':'float','MJDMINUS':'float',
                    'BUILD_PSF':'string','CNTRD_FWHM':'float','FITRAD':'float',
                    'FORCERADEC':'string','FRA':'float','FDEC':'float',
                    'FIND_ZPT':'string','PIXELATION_FACTOR':'float','SKYERR_RADIUS':'float',
                    'NEARBY_STARS_PIXEL_RAD':'float','GALAXY_MODEL_STEPS':'float','SN_PLUS_GALMODEL_STEPS':'float',
                    'SN_PLUS_GALMODEL_STEPS_GALSIM':'float',
                    'SN_SHIFT_STD':'float','HDR_PLATESCALE_NAME':'string','HDR_AIRMASS_NAME':'string',
                    'HDR_PSF_FWHM':'string','MINZPTSTARS':'float','FLUX_STD_DIV':'float','GALMODEL_DIV':'float',
                    'CUTINFOBANDS':[],'ZPTMIN':[],'SKYMIN':[],'SKYMAX':[],'SKYERRMAX':[]
                    }

def save_fits_image(image,filename):
    hdu = pf.PrimaryHDU(image)
    if os.path.exists(filename):
        os.remove(filename)
    hdu.writeto(filename)
    return

def gaussian(x, a, b, c):
    val = a * np.exp(-(x - b)**2 / (2*c**2))
    return val

def bindata(x,y,bins):
        medians = np.zeros(len(bins)-1)
        mads = np.zeros(len(bins)-1)
        xvals = (bins[:-1]+bins[1:])/2.

        for i in np.arange(len(bins)-1):
                bs = bins[i]
                bf = bins[i+1]
                ww = [(x>bs)&(x<bf)]
                try:
                        medians[i] = np.median(y[ww])
                        mads[i] = 1.48*np.median(abs(y[ww]-medians[i]))*1/np.sqrt(len(y[ww]))
                except IndexError:
                        medians[i] = np.nan
                        mads[i] = np.nan
        return xvals,medians,mads
def parabola(x,a,b,c):
    y = a*x**2+b*x+c
    return y

class get_snfile:
    def __init__(self,snfile, rootdir, useweights):
        varnames = ''
        fin = open(snfile,'r')
        survey = ''
        for line in fin:
            line = line.replace('\n','')
            if not line.startswith('#') and line.replace(' ',''):
                if not line.replace(' ','').startswith('OBS:') and \
                        not line.replace(' ','').startswith('VARNAMES:'):
                    key,val = line.split('#')[0].split(':')
                    key = key.replace(' ','')
                    #print key
                    #raw_input()
                    if key.upper() == 'SURVEY':
                        survey = val
                    if key.upper() == 'HOSTGAL_SB_FLUXCAL':
                        val = val.split()
                        self.__dict__[key.lower()] = val
                    elif key.upper() != 'WEIGHT_BADPIXEL' and (key.upper() != 'STARCAT'):
                        val = val.split()[0]
                        val = val.replace(' ','')
                        self.__dict__[key.lower()] = val
                    elif key.lower() == 'starcat' and 'des' in snfile:
                        #print 'here1'
                        catfilter = val.split()[0]                        
                        if filt.lower() == catfilter.lower():
                            #print val
                            self.__dict__["starcat"] = {catfilter.lower(): os.path.join(rootdir,val.split()[1])}
                        elif filt.lower() == 'all':
                            if "starcat" in self.__dict__:
                                self.__dict__["starcat"][val.split()[0]] = os.path.join(rootdir,val.split()[1])
                            else:
                                self.__dict__["starcat"] = {}
                                self.__dict__["starcat"][val.split()[0]] = os.path.join(rootdir,val.split()[1])
                    elif key.lower() == 'starcat' and 'PS1' in survey:
                        print 'ps1'
                        #raw_input()
                        self.__dict__["starcat"] = {filt.lower(): val.strip()}
                    elif key.lower() == 'starcat':
                        print 'here33'
                        catfilter = val.split()[0]
                        if filt.lower() == catfilter.lower():
                            # print val
                            self.__dict__["starcat"] = {catfilter.lower(): os.path.join(rootdir, val.split()[1])}
                    else:
                        #print 'here2',key.lower()
                        try:
                            self.__dict__[key.lower()] = np.array(val.split()).astype('float')
                        except:
                            raise exceptions.RuntimeError("Error : WEIGHT_BADPIXEL cannot be parsed!")
                
                elif line.split(' ')[0] == 'VARNAMES:':
                    varnames = filter(None,line.split('VARNAMES:')[-1].split(' '))
                    for v in varnames:
                        self.__dict__[v.lower()] = np.array([])
                elif line.replace(' ','').startswith('OBS:'):
                    vals = filter(None,line.split(' '))[1:]
                    if not varnames:
                        raise exceptions.RuntimeError("Error : Variable names are not defined!!")
                    elif len(varnames) != len(vals):
                        raise exceptions.RuntimeError("Error : Number of variables provided is different than the number defined!!!")
                    for var,val in zip(varnames,vals):
                        self.__dict__[var.lower()] = np.append(self.__dict__[var.lower()],val)

        catalog_exists = True
        for p in snkeywordlist.keys():
            if not self.__dict__.has_key(p.lower()):
                if p.lower() != 'starcat':
                    print "Error : keyword %s doesn't exist in supernova file!!!"%p
                    #raise exceptions.RuntimeError("Error : keyword %s doesn't exist in supernova file!!!"%p)
                else:
                    catalog_exists = False
            if snkeywordlist[p] == 'float':
                try:
                    self.__dict__[p.lower()] = float(self.__dict__[p.lower()])
                except:
                    raise exceptions.RuntimeError('Error : keyword %s should be set to a number!'%p)

        if useweights:
            for p in snvarnameslist.keys():
                if not self.__dict__.has_key(p.lower()):
                    if p.lower() != 'starcat':
                        print "Error : field %s doesn't exist in supernova file!!!"%p
                        #raise exceptions.RuntimeError("Error : field %s doesn't exist in supernova file!!!"%p)
                    elif catalog_exists == False:
                        raise exceptions.RuntimeError("Error : field %s doesn't exist in supernova file!!!"%p)
                if snvarnameslist[p] == 'float':
                    try:
                        self.__dict__[p.lower()] = self.__dict__[p.lower()].astype('float')
                    except:
                        raise exceptions.RuntimeError('Error : keyword %s should be set to a number!'%p)
        else:
            for p in snvarnameslist_ps.keys():
                if not self.__dict__.has_key(p.lower()):
                    if p.lower() != 'starcat':
                        raise exceptions.RuntimeError("Error : field %s doesn't exist in supernova file!!!" % p)
                    elif catalog_exists == False:
                        raise exceptions.RuntimeError("Error : field %s doesn't exist in supernova file!!!" % p)
                if snvarnameslist_ps[p] == 'float':
                    try:
                        self.__dict__[p.lower()] = self.__dict__[p.lower()].astype('float')
                    except:
                        raise exceptions.RuntimeError('Error : keyword %s should be set to a number!' % p)

class get_params:
    def __init__(self,paramfile):
        print os.listdir('.')
        print 'opening paramfile',paramfile
        fin = open(paramfile,'r')
        for line in fin:
            print line
            line = line.replace('\n','')
            if not line.startswith('#') and line.replace(' ',''):
                try:
                    key,val = line.split('#')[0].split(':')
                except:
                    raise exceptions.RuntimeError('Invalid format!  Should be key: value')
                key = key.replace(' ','')
                val = val.replace(' ','')
                if len(val.split(',')) > 1:
                    self.__dict__[key.lower()] = val.split(',')
                else:
                    self.__dict__[key.lower()] = val
        for p in paramkeywordlist.keys():

            if not self.__dict__.has_key(p.lower()):
                raise exceptions.RuntimeError("Error : keyword %s doesn't exist in parameter file!!!"%p)
            if paramkeywordlist[p] == 'float':
                try:
                    self.__dict__[p.lower()] = float(self.__dict__[p.lower()])
                except:
                    raise exceptions.RuntimeError('Error : keyword %s should be set to a number!'%p)

class smp:
    def __init__(self,snparams,params,rootdir,psf_model):
        self.snparams = snparams
        self.params = params
        self.rootdir = rootdir
        self.psf_model = psf_model

    #@profile
    def main(self,nodiff=False,nozpt=False,rootdir='',outdir='' ,
             nomask=False,outfile='',debug=False,
             verbose=False, clear_zpt=False,clear_checkstars=True,mergeno=0,
             mpfit_or_mcmc='mpfit',usefake=False,
             snfile='/test.dat',gal_model=None,stardumppsf=True,
             dogalfit=True,dosnfit=True,dogalsimfit=True, dogalsimpixfit=True,dosnradecfit=True,
             usediffimzpt=False,useidlsky=False,fixgalzero=True,floatallepochs=False,dailyoff=False,
             doglobalstar=True,exactpos=True,bigstarcatalog=None,
             stardeltasfolder=None, zptfoldername=None, galaxyfoldername=None,dobigstarcat=False,useweights=True,
             dosextractor=True,fermigrid=False,zptoutpath='./zpts/',fermigriddir=None,worker=False,
             savezptstamps=False,fermilog=False,isdonedir=None,oldformat=False,continu=False,continudir=None
             ):

        print 'snfile',snfile

        sstime = time.time()

        if fermigrid & worker:
            #if not os.path.exists(os.path.join(outdir,SNfoldername)):
            print 'ifdh mkdir ',outdir
            os.system('ifdh mkdir '+outdir)
        else:
            if not os.path.exists(outdir):
                os.mkdir(outdir)


        #print filt
        #print self.snparams.photflag
        #raw_input()
        #if 'x' in self.snparams.photflag[0]:
        #    self.snparams.photflag = ~(self.snparams.photflag == '0x00')
        #else:
        #    self.snparams.photflag = (self.snparams.photflag != '1')

        #print self.snparams.photflag
        print 'Starting Scene Modeling Photometry'
        if fermigrid and worker:
            useifdh = True
        else:
            useifdh = False
        self.tmpwriter = dt.tmpwriter(tmp_subscript=snfile.split('/')[-1].split('.')[0]+'_'+filt,useifdh=useifdh)
        print 'done with tmpwriter line 267'
        tstart = time.time()
        from txtobj import txtobj
        from astropy import wcs
        import astropy.io.fits as pyfits
        self.outfile = outfile
        print 'line 275'
        oldoutdir = copy(outdir)
        self.outdir = oldoutdir
        #print oldoutdir
        #sys.exit()
        if fermigrid & worker: #NEED TO ZIP AND COPY ALL DATA BACK TO OLDOUTFULE AFTER SMP IS DONE
            oldoutfile = copy(outfile)
            outfile = ''
            outdir = ''
        mainoutdir = copy(outdir)

        # cspath = os.path.join(outdir,foldername+'/SNe/starfits/')
        # if fermigrid and worker:
        #     if not os.path.exists(cspath)
        #     os.popen('ifdh mkdir '+cspath)
        #
        # if not os.path.exists(cspath):
        #     os.makedirs(cspath)
        self.checkstarfile = os.path.join(outdir,'SNe/starfits/'+snfile.split('/')[-1].split('.')[0]
                                          +'_'+filt+'_standardstarfits.txt')
        #print self.checkstarfile
        #print 'checkstarfile'
        #raw_input()
        self.snfn = snfile.split('/')[-1].split('.')[0]
        if nozpt:
            self.zpt_fits = './zpts/zpt_plots.txt'
            self.big_zpt = './zpts/big_zpt'
            #self.checkstarfile = os.path.join(outfile,foldername+'/SNe/'+snfile.split('/')[-1].split('.')[0] + '/'+filt+'/standardstarfits.txt')
            #print 'checkstarfile',self.checkstarfile
            #self.checkstarfile = self.params.checkstarfile
            if not os.path.exists('/'.join(self.checkstarfile.split('/')[:-1])):
                os.makedirs('/'.join(self.checkstarfile.split('/')[:-1]))
            if not os.path.exists('./zpts'):
                os.makedirs('./zpts/')
            #a = open(self.zpt_fits,'w')
            #a.write('ZPT FILE LOCATIONS\n')
            #a.close()
            self.tmpwriter.writefile('ZPT FILE LOCATIONS\n',self.zpt_fits)
            #print self.zpt_fits
            #raw_input('first instance of tmpwriter')
            if clear_zpt:
                #big = open(self.big_zpt+'.txt','w')
                #big.write('Exposure Num\tRA\tDEC\tCat Zpt\tMPFIT Zpt\tMPFIT Zpt Err\tMCMC Zpt\tMCMC Zpt Err\tMCMC Model Errors Zpt\tMCMC Model Errors Zpt Err\tCat Mag\tMP Fit Mag\tMCMC Fit Mag\tMCMC Model Errors Fit Mag\tMCMC Analytical Simple\tMCMC Analytical Weighted\n')
                #big.close()
                self.tmpwriter.writefile('Exposure Num\tRA\tDEC\tCat Zpt\tMPFIT Zpt\tMPFIT Zpt Err\tMCMC Zpt\tMCMC Zpt Err\tMCMC Model Errors Zpt\tMCMC Model Errors Zpt Err\tCat Mag\tMP Fit Mag\tMCMC Fit Mag\tMCMC Model Errors Fit Mag\tMCMC Analytical Simple\tMCMC Analytical Weighted\n',self.big_zpt)

            #b = open(self.checkstarfile, 'w')
            #b.write('Exposure Num\tMJD\tRA\tDEC\txstar\tystar\tCat Zpt\tMPFIT Zpt\tMPFIT Zpt Err\tFit Flux\tFit '
            #        'Flux Err\tFit Flux Chisq\tFit Flux DMS\tGalsim Fit Flux\tGalsim Fit Flux Err\t'
            #        'Galsim Fit Flux Chisq\tGalsim Fit Flux DMS\tCat Mag\n')
            #b.close()
            self.tmpwriter.writefile('Exposure Num\tMJD\tRA\tDEC\txstar\tystar\tCat Zpt\tMPFIT Zpt\tMPFIT Zpt Err\tFit Flux\tFit '
                    'Flux Err\tFit Flux Chisq\tFit Flux DMS\tGalsim Fit Flux\tGalsim Fit Flux Err\t'
                    'Galsim Fit Flux Chisq\tGalsim Fit Flux DMS\tCat Mag\n',self.checkstarfile)
            #if clear_checkstars:
            #    big = open(self.checkstarfile,'w')
            #    big.write('Exposure Num\tMJD\tRA\tDEC\tCat Zpt\tMPFIT Zpt\tMPFIT Zpt Err\tFit Flux\tFit Flux Err\tCat Mag\n')
            #    big.close()

        self.verbose = verbose
        params,snparams = self.params,self.snparams
        #print "FAKE TRUE MAGS"
        #print snparams.fake_truemag
        #print self.psf_model
        #raw_input()
        snparams.psf_model = self.psf_model
        snparams.snfile = snfile
        self.usefake = usefake
        self.gal_model = gal_model
        self.stardumppsf = stardumppsf
        self.dogalfit = dogalfit
        self.dosnfit = dosnfit
        self.dogalsimfit = dogalsimfit
        self.dogalsimpixfit = dogalsimpixfit
        self.usediffimzpt = usediffimzpt
        self.fixgalzero = fixgalzero
        self.floatallepochs = floatallepochs
        self.dosnradecfit = dosnradecfit
        self.rickfakestarfile = ''
        self.dosextractor = dosextractor
        self.worker=worker
        #self.lcfilepath=lcfilepath
        self.savezptstamps = savezptstamps
        self.fermilog = fermilog
        self.isdonedir = isdonedir
        self.oldformat = oldformat
        self.continu = continu
        self.continudir = continudir

        self.useweights = useweights
        if not self.useweights:
            self.snparams.image_name_weight = zip(self.snparams.image_name_noise,self.snparams.image_name_mask)


        self.exactpos = exactpos

        self.ras = []
        self.decs = []
        self.deltaras = []
        self.fakestarfluxes = []
        self.fakestarfluxerrs = []
        self.fakestarzpts = []
        self.deltadecs = []
        self.deltamjds = []
        self.airmasses = []
        self.x_stars = []
        self.y_stars = []
        badflags = []
        pixstart = None
        self.fermigrid = fermigrid
        self.zptoutpath = zptoutpath

        if snparams.psf_model == 'psfex' and not snparams.__dict__.has_key('psf'):
            raise exceptions.RuntimeError('Error : PSF must be provided in supernova file!!!')
        if filt != 'all':
            snparams.nvalid = 0
          
            for b in snparams.band:
                if b in filt:
                    snparams.nvalid +=1
        else:
            snparams.nvalid = snparams.nobs
        self.zptstamps = os.path.join(oldoutdir,'zptstamps')
        self.lcfilepath = os.path.join(oldoutdir,'lightcurves')
        if self.fermigrid and self.worker:
            os.popen('ifdh mkdir '+self.lcfilepath).read()
        else:
            os.popen('mkdir '+self.lcfilepath)

        if not os.path.exists(self.zptstamps):
            if fermigrid & worker:
                if self.zptstamps.split('/')[1] != 'pnfs':
                    raise ValueError(
                        '--zptoutpath must be located at /pnfs/des/persistent/desdm/ for fermigrid running')
                os.system('ifdh mkdir ' + self.zptstamps)
            else:
                os.makedirs(self.zptstamps)

        #print snparams.nvalid,params.substamp
        params.substamp = int(params.substamp)
        smp_bkg = np.zeros([snparams.nvalid,params.substamp,params.substamp])

        smp_im = np.zeros([snparams.nvalid,params.substamp,params.substamp])
        smp_noise = np.zeros([snparams.nvalid,params.substamp,params.substamp])
        smp_psf = np.zeros([snparams.nvalid,params.substamp,params.substamp])
        smp_mask = np.zeros([snparams.nvalid, params.substamp, params.substamp])
        smp_starind = np.zeros([snparams.nvalid, 500])
        smp_starra = np.zeros([snparams.nvalid, 500])
        smp_stardec = np.zeros([snparams.nvalid, 500])
        smp_dict = {'scale':np.zeros(snparams.nvalid),
                    'scale_err':np.zeros(snparams.nvalid),
                    'mcmc_scale':np.zeros(snparams.nvalid),
                    'mcmc_scale_err':np.zeros(snparams.nvalid),
                    'image_scalefactor':np.zeros(snparams.nvalid),
                    'snx':np.zeros(snparams.nvalid),
                    'sny':np.zeros(snparams.nvalid),
                    'scalefactor':np.zeros(snparams.nvalid),
                    'fwhm_arcsec':np.zeros(snparams.nvalid),
                    'sky':np.zeros(snparams.nvalid),
                    'skyerr':np.zeros(snparams.nvalid),
                    'flag':np.ones(snparams.nvalid),
                    'descriptiveflag': np.ones(snparams.nvalid),
                    'fitflag':np.ones(snparams.nvalid),
                    'psf':np.zeros(snparams.nvalid),
                    'psfcenter':[],
                    'psf_fwhm':np.zeros(snparams.nvalid),
                    'fakepsf':np.zeros(snparams.nvalid),
                    'zpt':np.zeros(snparams.nvalid),
                    'zpterr':np.zeros(snparams.nvalid),   
                    'skysig':np.zeros(snparams.nvalid),
                    'sexrms':np.zeros(snparams.nvalid),
                    'sexsky':np.zeros(snparams.nvalid),
                    'total_skyerr':np.zeros(snparams.nvalid),
                    'mjd':np.zeros(snparams.nvalid),
                    'mjd_flag':np.zeros(snparams.nvalid),
                    'gain': np.zeros(snparams.nvalid),
                    'mjdoff':[],
                    'rmsaddin': np.zeros(snparams.nvalid),
                    'mjdslopeinteroff':[],
                    'cat_mag':np.zeros(snparams.nvalid),
                    'image_filename':np.chararray(snparams.nvalid,itemsize=200),
                    'psf_filename':np.chararray(snparams.nvalid,itemsize=200),
                    'weight_filename':np.chararray(snparams.nvalid,itemsize=200),
                    'zpt_file':np.chararray(snparams.nvalid,itemsize=200),
                    'imwcs':[],
                    'mask':[],
                    'hostgal_mag':np.zeros(snparams.nvalid),
                    'hostgal_sbmag':np.zeros(snparams.nvalid),
                    'fakemag':np.zeros(snparams.nvalid),
                    'fakezpt':np.zeros(snparams.nvalid),
                    'diffim_flux':np.zeros(snparams.nvalid),
                    'diffim_fluxerr':np.zeros(snparams.nvalid),
                    'id_obs':np.zeros(snparams.nvalid),
                    'id_coadd':np.zeros(snparams.nvalid),
                    'snra':np.zeros(snparams.nvalid),
                    'sndec':np.zeros(snparams.nvalid),
                    'xoff': np.zeros(snparams.nvalid),
                    'yoff': np.zeros(snparams.nvalid),
                    'raoff': np.zeros(snparams.nvalid),
                    'decoff': np.zeros(snparams.nvalid),
                    'notbrightflag':np.ones(snparams.nvalid),
                    'fullims':[],
                    'impsfs':[],
                    'hpsfs':[],
                    'fileroots':np.chararray(snparams.nvalid,itemsize=200),

                    'psf_centerx':np.zeros(snparams.nvalid),
                    'psf_centery':np.zeros(snparams.nvalid),
                    'expnum':np.zeros(snparams.nvalid),

                    }
        smp_scale = np.zeros(snparams.nvalid)
        smp_sky = np.zeros(snparams.nvalid)
        smp_flag = np.zeros(snparams.nvalid)
        for i in np.arange(snparams.nvalid):
            smp_dict['image_filename'][i] = 'na'


        """
        band = 'r'
        if os.path.exists(snparams.starcat[band]):
            starcat = txtobj(snparams.starcat[band],useloadtxt=True, des=True)                                                                                                                                                       
            if not starcat.__dict__.has_key('mag_%s'%band):                                                                                                                                                                          
                try:                                                                                                                                                                                                                 
                    print starcat.__dict__                                                                                                                                                                                           
                    starcat.mag = starcat.__dict__[band]                                                                                                                                                                             
                    starcat.dmag = starcat.__dict__['d%s'%band]                                                                                                                                                                      
                except:                                                                                                                                                                                                             
                    raise exceptions.RuntimeError('Error : catalog file %s has no mag column!!'%snparams.starcat[band])                                                                                                              
        """

        i = 0

        index = 0
        orig_nozpt = copy(nozpt)

        starids = []
        starras = []
        #starcatras = []
        stardecs = []
        cntrs = 0
        cols = None

        filename = snparams.snfile.split('/')[-1].split('.')[0] +'_'+ filt

        if not nozpt:
            npoutdir = os.path.join(oldoutdir, 'np_data/' + filt + '/')
            datf = os.path.join(npoutdir, filename + '_withSn.npz')
            print datf
            ls = os.popen('ifdh ls ' + datf).read()
            #print 'ls',ls
            print len(ls) > 0
            #raw_input('hhhh')
            if len(ls) > 0:
                os.system('ifdh cp --force=xrootd '+datf+' .')
                dat = np.load(filename + '_withSn.npz')
                print dat['modelvec_nphistory'].shape[0]
                if dat['modelvec_nphistory'].shape[0] == 1000:
                    print 'already fitted. Exiting gracefully'

        staroutdir = os.path.join(oldoutdir,'staroffsets',filt)
        if self.fermigrid and self.worker:
            print os.popen('ifdh mkdir '+os.path.join(oldoutdir,'staroffsets')).read()
            print os.popen('ifdh mkdir '+staroutdir).read()
        else:
            if not os.path.exists(staroutdir):
                try:
                    os.mkdir(os.path.join(oldoutdir,'staroffsets'))
                except:
                    print 'path already exists', os.path.join(oldoutdir,'staroffsets')
                try:
                    os.mkdir(staroutdir)
                except:
                    print 'path already exists', staroutdir
        star_offset_file = os.path.join(staroutdir,filename+'band_starGlobalOffsets.npz')
        print 'loading',star_offset_file

        if not nozpt:
            if fermigrid and worker:
                #lsa = os.popen('ifdh ls '+staroutdir).read()
                #print lsa
                ls = os.popen('ifdh ls '+star_offset_file).read()
                print ls
                if len(ls) > 0:
                    print star_offset_file
                    #sys.exit()
                    os.system('ifdh cp --force=xrootd '+star_offset_file+' .')
                    star_offset_file = filename+'band_starGlobalOffsets.npz'
                else:
                    print 'Could not find star offset file. Calculating...'
                    #raw_input('hereweare')
                    nozpt = True
            else:
                try:
                    print 'loading',star_offset_file
                    staroffsets = np.load(star_offset_file)
                    print 'Found and loaded star offset file...'
                except:
                    if not doglobalstar:
                        pass
                    else:
                        print 'Could not find star offset file. Calculating...'
                        nozpt = True

        
        #############################################################################################################################
        ################################################# GET STAR GLOBAL OFFSETS ###################################################

        print 'getting star global offsets'
        print nozpt
        #print self.usefake
        didglobalstar = False
        #sys.exit()

        if self.fermilog:
            self.fermilogfile = os.path.join(oldoutdir,'fermilogs', snfile.split('/')[-1].split('.')[0] + '.smplog')
            os.system('ifdh mkdir '+os.path.join(oldoutdir,'fermilogs'))

            print self.fermilogfile
            if os.path.exists(self.fermilogfile):
                print os.popen('ifdh rm '+self.fermilogfile).read()
            self.tmpwriter.writefile('Starting SMP\n',self.fermilogfile)
            print 'Follow .smplog at ', self.fermilogfile
            self.tmpwriter.appendfile('running globalstar\n',self.fermilogfile)
            print 'DOING FERMILOG'*100
            #sys.exit()
        else:
            self.fermilogfile = None
        #print self.snparams.__dict__
        if self.snparams.survey == 'PS1':
            print 'reading in starcaftile'
            self.starcatfile = self.snparams.starcat[self.snparams.band[0]]
            print self.snparams.band[0]
            print self.starcatfile
            #raw_input()
            #print self.starcatfile['g']
            #starcat = dt.readcol)
            self.starcatmagfile = self.starcatfile.replace('/abswcs','')
            self.starcatmagfile = self.starcatmagfile.replace('full.abswcs.','')
            #'/export/scratch0/ps1sn1/data/v10.0/GPC1v3/eventsv1/absphotcat/0x1615/PSc000010.0x1615.dougpv3.full.abswcs.phot1'
            #'/export/scratch0/ps1sn1/data/v10.0/GPC1v3/eventsv1/absphotcat/0x1615/PSc000010.0x1615.dougpv3.phot1'
            print self.starcatmagfile
            starcat = txtobj(self.starcatfile, useloadtxt=True)
            starcatmag = txtobj(self.starcatmagfile, useloadtxt=True)

            rastarr = []
            decstarr = []
            smag = []
            #print starcatmag.__dict__
            #raw_input()
            for iii, ram,decm,m in zip(range(len(starcatmag.__dict__['ra'])),starcatmag.__dict__['ra'],
                               starcatmag.__dict__['dec'],starcatmag.__dict__[filt.lower()]):
                #print starcat.__dict__['ra'][0:10]
                #print xm
                #print np.sqrt((starcat.__dict__['Xpos']-xm)**2 + (starcat.__dict__['Ypos']-ym)**2)
                wwww = np.sqrt((starcat.__dict__['ra']-ram)**2 + (starcat.__dict__['dec']-decm)**2) < .001
                try:
                    print starcat.__dict__['ra'][wwww]
                    if len(starcat.__dict__['ra'][wwww]) >0 :
                        rastarr.append(float(starcat.__dict__['ra'][wwww][0]))
                        decstarr.append(float(starcat.__dict__['dec'][wwww][0]))
                        smag.append(m)
                    #raw_input()
                except:
                    print 'could not find associated star'
                #raw_input()
            rastarr = np.array(rastarr)
            decstarr = np.array(decstarr)
            smag = np.array( smag )
            # print starcatmag.__dict__
            # print starcatmag.__dict__['i'] - starcat.__dict__['i']
            # print 'done reading in starcatfile'
            # raw_input()
            #print starcat.__dict__
            # print self.starcat.__dict__
            # print self.starcat.__dict__['RA']
            # print self.starcat.__dict__['DEC']
            # print self.starcat.__dict__['MAG_PSF_MEAN_%s'%filt.upper()]
            #print starcat.__dict__
            #raw_input()
            #starcat.bigra = np.array(starcat.__dict__['ra'][1:], dtype='float')
            #starcat.bigdec = np.array(starcat.__dict__['dec'][1:], dtype='float')
            print rastarr
            print rastarr.shape
            starcat.bigra = rastarr[1:]
            starcat.bigdec = decstarr[1:]
            starcat.bigmag = smag[1:]
            # starcat.x = xstarr
            # starcat.y = ystarr
            # starcat.bigmag = np.array(starcat.__dict__[filt.lower()][1:], dtype='float')
            starcat.mag = starcat.bigmag
            starcat.bigid = np.arange(len(starcat.bigra))
            starcat.objid = np.arange(len(starcat.bigra))
            starcat.ra = starcat.bigra
            starcat.dec = starcat.bigdec
            newra = starcat.bigra
            newdec = starcat.bigdec
        elif self.snparams.survey == 'DES':
            print 'reading in starcaftile'
            wehavestarcat = False
            self.starcatfile = 'catalogs/des/RykoffY3A1Catalog_AB_Beta.tab'
            if nozpt:
                starcat = txtobj(self.starcatfile, useloadtxt=True)
                print 'done reading in starcatfile'
                wehavestarcat = True
                #print self.starcat.__dict__
                #print self.starcat.__dict__['RA']
                #print self.starcat.__dict__['DEC']
                #print self.starcat.__dict__['MAG_PSF_MEAN_%s'%filt.upper()]

                starcat.bigra = np.array(starcat.__dict__['RA'][1:],dtype='float')
                starcat.bigdec = np.array(starcat.__dict__['DEC'][1:],dtype='float')
                starcat.bigmag = np.array(starcat.__dict__['MAG_PSF_MEAN_%s'%filt.upper()][1:],dtype='float')
                starcat.bigid = np.array(starcat.__dict__['MATCH_OBJECT_ID'][1:],dtype='float')
                starcat.id = starcat.bigid
                starcat.ra = starcat.bigra
                starcat.dec = starcat.bigdec



        else:
            raise Exception('Unknown survey. Not supported. Please contact Dillon Brout at dbrout@phsyics.upenn.edu')
        #print self.starcat.ra
        #sys.exit()
        badindices = []
        for imfile,noisefile,psffile,band, j in \
                zip(snparams.image_name_search,snparams.image_name_weight,snparams.file_name_psf,snparams.band, range(len(snparams.band))):
            #print doglobalstar, snparams.mjd[j], nozpt
            if not doglobalstar:
                continue
            if snparams.mjd[j] == 0:
                continue
            #if round(snparams.mjd[j]) != 56559:
            #    continue
            #if float(snparams.mjd[j]) < 57045:
            #    continue
            if not nozpt:
                continue
            if not band == filt:
                continue

            if nozpt:
                if self.snparams.survey == 'DES':
                    if not wehavestarcat:
                        starcat = txtobj(self.starcatfile, useloadtxt=True)
                        print 'done reading in starcatfile'
                        wehavestarcat = True


            skysig=np.nan
            if cntrs > 1:
               continue
            #if snparams.mjd[j] != 56636.:
            #    if snparams.mjd[j] < 57000.:
            #        continue
            didglobalstar = True
            #nozpt = copy(orig_nozpt)
            if self.usefake:
                imfile = ''.join(imfile.split('.')[:-1]) + '+fakeSN.fits'

            if not self.oldformat:
                if self.usefake:
                    if not self.snparams.survey == 'PS1':
                        imfile = imfile.replace('p1', 'Y1')
                        noisefile = noisefile.replace('p1', 'Y1')
                        psffile = psffile.replace('p1', 'Y1')
                    #noisefile = ''.join(noisefile.split('.')[:-1])+'+fakeSN.fits'
                else:
                    if not self.snparams.survey == 'PS1':
                        imfile = imfile.replace('p1','Y1')
                        noisefile = noisefile.replace('p1','Y1')
                        psffile = psffile.replace('p1', 'Y1')

            imfile = os.path.join(rootdir,imfile)
            longimfile = copy(imfile)
            longpsffile = copy(psffile)


            #print longimfile
            #raw_input()
            self.impath = '/'.join(imfile.split('/')[:-1])
            try:
                noisefile = os.path.join(rootdir,noisefile)
            except:
                noisefile = [os.path.join(rootdir,noisefile[0]),os.path.join(rootdir,noisefile[1])]

            if self.snparams.survey == 'DES':
                if self.usefake:
                    longnoisefile = copy(imfile.replace('+fakeSN.fits','.weight.fits'))
                else:
                    longnoisefile = copy(imfile.replace('.fits','.weight.fits'))
            else:
                longnoisefile = copy(noisefile)

            if self.fermilog:
                self.tmpwriter.appendfile('running globalstar on '+longimfile+'\n', self.fermilogfile)

            psffile = os.path.join(rootdir,psffile)
            #print imfile
            #raw_input()
            print 'hereeeee'
            if self.fermigrid & self.worker:
                #print imfile
                #os.system('IFDH_CP_MAXRETRIES=1; ifdh cp ' + imfile + ' .')
                #print 'line 529 copied image files to here'
                #sys.exit()
                #print 'ifdh cp '+imfile+' .'

                ifdhls = os.popen('ifdh lss '+imfile).read()
                #print ifdhls
                #print 'line 576'
                #raw_input()
                if (len(ifdhls) > 1):
                    if (int(ifdhls.split()[-1]) > 0) :
                        print 'Copying over',imfile
                        os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp --force=xrootd '+imfile+' .').read()
                        #imfilel = copy(imfilel)
                        imfile = imfile.split('/')[-1]
                        print 'imfile',imfile
                        # if self.usefake:
                        #     #if '.gz' in imfile:
                        #     print 'ifdh','IFDH_CP_MAXRETRIES=1; ifdh cp ' + imfilel.split('.fits.gz')[0]+ '+fakeSN.fits.gz' + ' .'
                        #     os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp ' + imfilel.split('.fits.gz')[0]+ '+fakeSN.fits.gz' + ' .').read()
                        #     #imfile = imfilel.split('/')[-1]
                        #     #else:
                        #     os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp ' + imfilel.split('.fits')[
                        #         0] + '+fakeSN.fits' + ' .').read()
                        #     imfile = imfilel.split('/')[-1]
                        #print 'IFDH_CP_MAXRETRIES=1; ifdh cp '+noisefile+' .'
                        os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp --force=xrootd '+noisefile+' .').read()
                        noisefile = noisefile.split('/')[-1]
                        weightsfile = noisefile
                        #print 'ifdh cp ' + psffile + ' .'
                        os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp --force=xrootd ' + psffile + ' .').read()
                        psffile = psffile.split('/')[-1]
                        #print 'copied all files'
                        #print os.popen('ifdh ls .').read()
                        #sys.exit()
                        print 'here2'
                        #raw_input()
                else:
                    #print 'file not found',imfile
                    #continue

                    ifdhls = os.popen('ifdh ls ' + imfile+'.fz').read()
                    #print ifdhls
                    #print 'line 610'
                    #raw_input()
                    if len(ifdhls) > 0:
                        a = os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp --force=xrootd ' + imfile + '.fz .').read()
                        #print a
                        a = os.popen('funpack '+imfile.split('/')[-1]+'.fz').read()
                        #print a
                        os.popen('rm '+imfile.split('/')[-1]+'.fz').read()
                        # imfilel = copy(imfilel)
                        imfile = imfile.split('/')[-1]
                        #print 'imfile', imfile
                        # if self.usefake:
                        #     #if '.gz' in imfile:
                        #     print 'ifdh','IFDH_CP_MAXRETRIES=1; ifdh cp ' + imfilel.split('.fits.gz')[0]+ '+fakeSN.fits.gz' + ' .'
                        #     os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp ' + imfilel.split('.fits.gz')[0]+ '+fakeSN.fits.gz' + ' .').read()
                        #     #imfile = imfilel.split('/')[-1]
                        #     #else:
                        #     os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp ' + imfilel.split('.fits')[
                        #         0] + '+fakeSN.fits' + ' .').read()
                        #     imfile = imfilel.split('/')[-1]
                        # print 'IFDH_CP_MAXRETRIES=1; ifdh cp '+noisefile+' .'
                        os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp --force=xrootd ' + noisefile + '.fz .').read()
                        os.popen('funpack '+noisefile.split('/')[-1]+'.fz')
                        os.popen('rm '+noisefile.split('/')[-1]+'.fz')
                        noisefile = noisefile.split('/')[-1]
                        weightsfile = noisefile
                        # print 'ifdh cp ' + psffile + ' .'
                        os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp --force=xrootd ' + psffile + ' .').read()
                        #os.popen('funpack '+psffile)
                        psffile = psffile.split('/')[-1]
                        # print 'copied all files'
                        # print os.popen('ifdh ls .').read()
                        # sys.exit()
                        #print 'here2'
                        #raw_input()
                    else:
                        print 'file not found', imfile
                        continue
                    #raw_input('cpied fz filehere')
                print 'grabbed sn files'
                #sys.exit()
            try:
                self.field = imfile.split('/')[-1].split('-')[1].split('_')[0]
                self.ccdnum = imfile.split('/')[-1].split('_')[-1][:2]
                self.expnum = imfile.split('/')[-1].split('_')[1]

            except:
                self.ccdnum = np.nan
                self.field = np.nan
                self.expnum = np.nan

            # self.laskerstarcat = 'CCD'+self.ccdnum+'_'+filt+'band.starcat'
            # print self.laskerstarcat
            #sys.exit()

            if filt != 'all' and band not in filt:
                #print('filter %s not in filter list %s for image file %s'%(band,filt,imfile))
                #print 'filter %s,%s not in filter list for image file %s'%(band,filt,imfile)
                continue

            #imfileloc = copy(imfile)

            #if fermigrid and worker:
            #    self.rootdir = '.'

            #print imfile
            #print imfileloc


            if not worker:
                if self.useweights:
                    weightsfile = os.path.join(self.rootdir,noisefile)
                else:
                    noisefile, maskfile = os.path.join(self.rootdir,noisefile[0]),os.path.join(self.rootdir,noisefile[1])

            if not self.oldformat:
                if self.usefake:
                    if not self.snparams.survey == 'PS1':
                        print imfile
                        changename = True
                        if changename:
                            if int(imfile.replace(self.rootdir,'')[:8]) < 20140601:
                                imfile = imfile.replace('p1', 'Y1')
                                noisefile = noisefile.replace('p1', 'Y1')
                                psffile = psffile.replace('p1', 'Y1')
                            elif int(imfile.replace(self.rootdir,'')[:8]) < 20150601:
                                imfile = imfile.replace('p1', 'Y2')
                                noisefile = noisefile.replace('p1', 'Y2')
                                psffile = psffile.replace('p1', 'Y2')
                            elif int(imfile.replace(self.rootdir,'')[:8]) < 20160601:
                                imfile = imfile.replace('p1', 'Y3')
                                noisefile = noisefile.replace('p1', 'Y3')
                                psffile = psffile.replace('p1', 'Y3')
                            elif int(imfile.replace(self.rootdir,'')[:8]) < 20170601:
                                imfile = imfile.replace('p1', 'Y4')
                                noisefile = noisefile.replace('p1', 'Y4')
                                psffile = psffile.replace('p1', 'Y4')
                else:
                    if not self.snparams.survey == 'PS1':
                        print imfile

                        #raw_input('imh')
                        if int(imfile.replace(self.rootdir,'')[:8]) < 20140601:
                            imfile = imfile.replace('p1', 'Y1')
                            noisefile = noisefile.replace('p1', 'Y1')
                            weightsfile = weightsfile.replace('p1', 'Y1')
                            psffile = psffile.replace('p1', 'Y1')
                        elif int(imfile.replace(self.rootdir,'')[:8]) < 20150601:
                            imfile = imfile.replace('p1', 'Y2')
                            imfile = imfile.replace('Y1', 'Y2')
                            noisefile = noisefile.replace('p1', 'Y2')
                            noisefile = noisefile.replace('Y1', 'Y2')
                            psffile = psffile.replace('p1', 'Y2')
                            psffile = psffile.replace('Y1', 'Y2')
                            weightsfile = weightsfile.replace('p1', 'Y2')
                            weightsfile = weightsfile.replace('Y1', 'Y2')

                        elif int(imfile.replace(self.rootdir,'')[:8]) < 20160601:
                            imfile = imfile.replace('p1', 'Y3')
                            imfile = imfile.replace('Y1', 'Y3')
                            noisefile = noisefile.replace('p1', 'Y3')
                            noisefile = noisefile.replace('Y1', 'Y3')
                            psffile = psffile.replace('p1', 'Y3')
                            psffile = psffile.replace('Y1', 'Y3')
                            weightsfile = weightsfile.replace('p1', 'Y3')
                            weightsfile = weightsfile.replace('Y1', 'Y3')

                        elif int(imfile.replace(self.rootdir,'')[:8]) < 20170601:
                            imfile = imfile.replace('p1', 'Y4')
                            imfile = imfile.replace('Y1', 'Y4')
                            noisefile = noisefile.replace('p1', 'Y4')
                            noisefile = noisefile.replace('Y1', 'Y4')
                            psffile = psffile.replace('p1', 'Y4')
                            psffile = psffile.replace('Y1', 'Y4')
                            weightsfile = weightsfile.replace('p1', 'Y4')
                            weightsfile = weightsfile.replace('Y1', 'Y4')

            # print imfile
            # print imfile.replace(self.rootdir, '')[:8]
            # raw_input('stopped')
            if not os.path.exists(imfile):
                #print os.popen('ls -ltr *.fz').read()
                #print 'funpack %s.fz' % imfile
                #os.system('funpack %s.fz' % imfile)
                #sys.exit()
                if not os.path.exists(imfile+'.fz'):
                    print('Error : file %s does not exist'%imfile)
                    continue
                    print('Error : file %s does not exist'%imfile)
                    raise exceptions.RuntimeError('Error : file %s does not exist'%imfile)
                else:
                    print os.popen('ls -ltr *.fz').read()
                    print 'funpack %s.fz'%imfile
                    os.system('funpack %s.fz'%imfile)
                    #sys.exit()

            if self.useweights:
                if not os.path.exists(weightsfile):
                    os.system('gunzip %s.gz'%weightsfile)
                    if not os.path.exists(weightsfile):
                        os.system('funpack %s.fz'%weightsfile)
                        if not os.path.exists(weightsfile):
                            raise exceptions.RuntimeError('Error : file %s does not exist'%weightsfile)
            else:

                if not os.path.exists(noisefile):
                    os.system('gunzip %s.gz' % noisefile)
                    if not os.path.exists(noisefile):
                        os.system('funpack %s.fz' % noisefile)
                        if not os.path.exists(noisefile):
                            raise exceptions.RuntimeError('Error : file %s does not exist' % noisefile)
                if not os.path.exists(maskfile):
                    os.system('gunzip %s.gz' % maskfile)
                    if not os.path.exists(maskfile):
                        os.system('funpack %s.fz' % maskfile)
                        if not os.path.exists(maskfile):
                            raise exceptions.RuntimeError('Error : file %s does not exist' % maskfile)

            if not os.path.exists(psffile):
                if os.path.exists(psffile+'.gz'):
                    os.system('gunzip %s.gz' % psffile)
                elif not os.path.exists(psffile+'.fz'):
                    raise exceptions.RuntimeError('Error : file %s does not exist'%psffile)
                else:
                    os.system('/global/u1/d/dbrout/cfitsio/funpack %s.fz'%psffile)

            # print maskfile
            # mask = pyfits.getdata(maskfile)
            # print mask.shape
            # raw_input()
            if not nomask:
                if useweights:

                    maskfile = os.path.join(self.rootdir,snparams.image_name_search[j])
                    mask = pyfits.getdata(maskfile)

            print 'here3'
            #raw_input()
            # if self.usefake:
            #     fakeim = ''.join(imfile.split('.')[:-1]) + '+fakeSN.fits'
            #     print 'IFDH_CP_MAXRETRIES=1; ifdh cp ' + ''.join(longimfile.split('.')[:-1]) + '+fakeSN.fits' + ' .'
            #     print 'ifdh', 'IFDH_CP_MAXRETRIES=1; ifdh cp ' + longimfile.split('.fits.gz')[
            #         0] + '+fakeSN.fits.gz' + ' .'
            #     sys.exit()
            #     os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp ' + ''.join(longimfile.split('.')[:-1]) + '+fakeSN.fits' + ' .')\
            #         .read()
            #     if not os.path.exists(fakeim):
            #         print 'ifdh', 'IFDH_CP_MAXRETRIES=1; ifdh cp ' + longimfile.split('.fits.gz')[
            #             0] + '+fakeSN.fits.gz' + ' .'
            #         os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp ' + longimfile.split('.fits.gz')[
            #             0] + '+fakeSN.fits.gz' + ' .').read()
            #         os.system('funpack %s.fz' % fakeim)
            #         os.system('gunzip %s.gz' % fakeim)
            #         print longimfile
            #         print longimfile.split('.fits.gz')[0] + '+fakeSN.fits.gz' + ' .'
            #         sys.exit()
            #     imfile = fakeim
            try:
                im = pyfits.getdata(imfile)
                hdr = pyfits.getheader(imfile)
                print 'running globalstars on ',imfile
                #print 'image shape',im.shape
            except:
                print 'Image is EMPTY, skipping star...'
                continue
            #raw_input()
            print 'got image data!'
            #sys.exit()
            #fakeim_hdr = pyfits.getheader(fakeim)
            #snparams.cat_zpts[imfile] = fakeim_hdr['HIERARCH DOFAKE_ZP']

            snparams.platescale = hdr[self.params.hdr_platescale_name]
            snparams.airmass = hdr[self.params.hdr_airmass_name]

            if self.useweights:
                weights = pyfits.getdata(weightsfile)
            else:
                noise = pyfits.getdata(noisefile)
                mask = pyfits.getdata(maskfile)


            psf = pyfits.getdata(psffile)

            if params.weight_type.lower() == 'ivar':
                print 'ivar'
                #raw_input()
                noise = np.sqrt(1/noise)
            elif params.weight_type.lower() != 'noise':
                raise exceptions.RuntimeError('Error : WEIGHT_TYPE value %s is not a valid option'%params.WEIGHT_TYPE)

            if nomask:
                if self.useweights:
                    mask = np.zeros(np.shape(weights))
                    maskcols = np.where((weights < 1e-8) |
                                        (np.isfinite(weights) == False))
                    mask[maskcols] = 100.0

            wcsworked = True
            try:
                w =wcs.WCS(imfile)
            except:
                wcsworked = False
                hl = pf.open(imfile)
                import starlink.Ast as Ast
                import starlink.Atl as Atl
                fitschan = Ast.FitsChan( Atl.PyFITSAdapter(hl[ 0 ]) )
                encoding = fitschan.Encoding
                wcsinfo = fitschan.read()
                w = wcs.WCS(imfile+'old')
            
            if wcsworked:
                results =  w.wcs_pix2world(np.array([[0,0]]), 0)
                ra1, dec1 = results[0][0], results[0][1]
                results2 =  w.wcs_pix2world(np.array([[snparams.nxpix-1,
                                   snparams.nypix-1]]), 0)
                ra2, dec2 =results2[0][0], results2[0][1]
            else:
                xpix = [0]
                ypix = [0]
                radtodeg = 360/(2*3.14159)
                results =  wcsinfo.tran([xpix,ypix])
                ra1, dec1 = results[0]*radtodeg, results[1]*radtodeg
                results2 =  wcsinfo.tran([[snparams.nxpix-1], [snparams.nypix-1]])
                ra2, dec2 =results2[0]*radtodeg, results2[1]*radtodeg

            ra_high = np.max([ra1,ra2])
            ra_low = np.min([ra1,ra2])
            dec_high = np.max([dec1,dec2])
            dec_low = np.min([dec1,dec2])
            #print 'starcat'*50
            #print snparams.starcat
            #sys.exit()


            #
            # if type(snparams.starcat) == np.array:
            #     if os.path.exists(snparams.starcat[j]):
            #         starcat = txtobj(snparams.starcat[j],useloadtxt=True)
            #         if not starcat.__dict__.has_key('mag'):
            #             try:
            #                 starcat.mag = starcat.__dict__[band]
            #                 starcat.dmag = starcat.__dict__['d%s'%band]
            #             except:
            #                 raise exceptions.RuntimeError('Error : catalog file %s has no mag column!!'%snparams.starcat[j])
            #     else:
            #         raise exceptions.RuntimeError('Error : catalog file %s does not exist!!'%snparams.starcat[j])
            # elif type(snparams.starcat) == dict and 'des' in snfile:
            #     starcatfile = None
            #     starcatloc = '/'.join(longimfile.split('/')[0:-1])+'/'
            #     #print starcatloc
            #     #print imfile
            #     #print longimfile
            #     #sys.exit()
            #     if fermigrid and worker:
            #         starcatloc = '/'.join(longimfile.split('/')[0:-1])+'/'
            #         #starcatloc = '/'.join(longimfile.split('/')[0:-2]) + '/g'+longimfile.split('/')[-2][1:]
            #         print starcatloc
            #         ifdhls = os.popen('ifdh ls ' + starcatloc + '/').read()
            #         #print ifdhls
            #         #print 'ls on imfileloc'
            #         ifdhls = os.popen('ifdh ls ' + starcatloc + '/STARCAT*.LIST').read()
            #         #print ifdhls
            #         #print 'ls on imfileloc/STARCAT*.LIST'
            #         #print 'len starcat',len(ifdhls)
            #         #sys.exit()
            #         if len(ifdhls) > 0:
            #             os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp ' + ifdhls.strip() + ' .').read()
            #             a = os.popen('ls STARCAT*').read()
            #             #print '743',a
            #             starcatfile = ifdhls.strip().split('/')[-1]
            #             #print '745',starcatfile
            #             #sys.exit()
            #             starcatloc = ''
            #             ifdhls = os.popen('ifdh ls  ./STARCAT*.LIST').read()
            #             #print ifdhls
            #             #print 'sssssssssss'
            #             starcatloc = ''
            #             #raw_input()
            #         else:
            #             continue
            #     else:
            #         for fl in os.listdir(starcatloc):
            #             if 'STARCAT' in fl:
            #                 starcatfile = fl
            #     #raw_input()
            #     print starcatloc+starcatfile
            #     if os.path.exists(starcatloc+starcatfile):
            #         starcat = txtobj(starcatloc+starcatfile,useloadtxt=True, des=True)
            #         if not starcat.__dict__.has_key('mag_%s'%band):
            #             try:
            #                 print starcat.__dict__
            #                 starcat.mag = starcat.__dict__[band]
            #                 starcat.dmag = starcat.__dict__['d%s'%band]
            #             except:
            #                 raise exceptions.RuntimeError('Error : catalog file %s has no mag column!!'%snparams.starcat[band])
            #     else:
            #         raise exceptions.RuntimeError('Error : catalog file %s does not exist!!'%snparams.starcat[band])
            # else:
            #     print snparams.starcat[filt]
            #     if os.path.exists(snparams.starcat[filt]):
            #         starcat = txtobj(snparams.starcat[filt],useloadtxt=True)
            #         if not starcat.__dict__.has_key('mag'):
            #             try:
            #                 starcat.mag = starcat.__dict__[band]
            #                 starcat.dmag = starcat.__dict__['d%s'%band]
            #             except:
            #                 print snparams.starcat[filt]
            #                 raise exceptions.RuntimeError('Error : catalog file %s has no mag column!!'%snparams.starcat[filt])
            #     else:
            #         raise exceptions.RuntimeError('Error : catalog file %s does not exist!!'%snparams.starcat[filt])
            #

            #if not self.snparams.survey == 'PS1':

            if nozpt:

                #self.rdnoise = hdr[params.rdnoise_name]
                #self.gain = hdr[params.gain_name]
                # print params.gain_name
                # print 'self.gain',self.gain
                # raw_input()

                # print starcat.ra
                # print 'ralow',ra_low
                # print 'rahi',ra_high
                # print 'min(starcat.ra)',min(starcat.ra)
                # print 'max(starcat.ra)',max(starcat.ra)
                # sys.exit()
                # print starcat.__dict___
                if self.snparams.survey == 'DES':
                    cols = np.where((starcat.bigra > ra_low) &
                                    (starcat.bigra < ra_high) &
                                    (starcat.bigdec > dec_low) &
                                    (starcat.bigdec < dec_high))[0]

                    starcat.ra = starcat.bigra[cols]
                    starcat.dec = starcat.bigdec[cols]
                    starcat.mag = starcat.bigmag[cols]
                    starcat.objid = starcat.bigid[cols]


                    #self.tmpwriter.savez(self.zptoutpath+'/'+str(self.field)+'_CCD'+str(self.ccdnum)+'_standards.npz',
                    #         ra=starcat.ra,dec=starcat.dec,mag=starcat.mag,id=starcat.objid)

                    if not len(cols):
                        # print "Error : No stars in image!!"
                        # badindices.append(j)
                        # continue
                        raise exceptions.RuntimeError("Error : No stars in image!!")

                if wcsworked:
                    coords = zip(*w.wcs_world2pix(np.array(zip(starcat.ra, starcat.dec)), 0))
                    x_star, y_star = [], []
                    for xval, yval in zip(*coords):
                        x_star += [xval]
                        y_star += [yval]
                    print starcat.ra
                    print x_star
                else:
                    coords = wcsinfo.tran([starcat.ra / radtodeg, starcat.dec / radtodeg], False)

                x_star, y_star = [], []
                for xval, yval in zip(*coords):
                    x_star += [xval]
                    y_star += [yval]

                x_star1, y_star1 = np.array(x_star), np.array(y_star)
                x_star, y_star = cntrd.cntrd(im, x_star1, y_star1, params.cntrd_fwhm)
                newra, newdec = zip(*w.wcs_pix2world(np.array(zip(x_star, y_star)), 0))
                # try:
                #     starcat.objid += 0.
                # except:
                #     starcat.objid = np.arange(len(starcat.mag))

                # for rrr in starcat.objid:
                #    starids.append(rrr)
                # for rrr,zzz in zip(newra,newdec):
                #    starras.append(rrr)
                #    stardecs.append(zzz)
                # for rrr in starcat.ra:
                #    starcatras.append(rrr)
                cntrs += 1
                # print len(starids)
                # print len(starras)
                # print starras

                # raw_input()

        if nozpt:
            starids = np.array(starcat.objid)
            starras = np.array(newra)
            stardecs = np.array(newdec)
            starmags = np.array(starcat.mag)

            print starras
            print star_offset_file
            self.tmpwriter.savez(star_offset_file, starras=np.array(newra), stardecs=np.array(newdec),
                                 starids=np.array(starcat.objid), starmags=np.array(starcat.mag))
        else:
            staroffsets = np.load(star_offset_file)
            starras = np.array(staroffsets['starras'])
            stardecs = np.array(staroffsets['stardecs'])
            starids = np.array(staroffsets['starids'])
            starmags = np.array(staroffsets['starmags'])
            # sys.exit()
        starglobalras = []
        starglobaldecs = []
        starglobalids = []
        starglobalmags = []
        print starids.shape
        print np.unique(starids)
        #raw_input()
        for ide in np.unique(starids):
            #print starids
            ww = (starids == ide)
            starglobalids.append(ide)
            #print starras[ww]
            try:
                #print 'aaa', starras
                #raw_input()
                print ide,np.median(starras[ww]),np.median(stardecs[ww]),np.median(starmags[ww])
                starglobalras.append(np.median(starras[ww]))
                starglobaldecs.append(np.median(stardecs[ww]))
                starglobalmags.append(np.median(starmags[ww]))
            except:
                starglobalras.append(np.nan)
                starglobaldecs.append(np.nan)
                starglobalmags.append(np.nan)

        starglobalids = np.array(starglobalids)
        starglobalras = np.array(starglobalras)
        starglobaldecs = np.array(starglobaldecs)
        starglobalmags = np.array(starglobalmags)

        # if self.dobigstarcat:
        #     # scampra,scampdec,mag_starwrong = self.getProperCatRaDec(starglobalras,starglobaldecs)
        #     # offsetra = np.array(starglobalras) - np.array(scampra)
        #     # offsetdec = np.array(starglobaldecs) - np.array(scampdec)
        #     #jra,jdec,jmag = self.getJamesCat()
        #     print 'we are here inside bigstarcat'
        # else:
        offsetra = np.array(starglobalras) * 0.
        offsetdec = np.array(starglobalras) * 0.
        print 'Done with globalstars'
        if self.fermilog:
            self.tmpwriter.appendfile('\nDone with globalstars\n ', self.fermilogfile)

        #############################################################################################################################
        #############################################################################################################################
        #print offsetra
        #print offsetdec
        print 'Done with centroiding!!'
        #sys.exit()
        #orig_nozpt = copy(nozpt)
        print nozpt
        #raw_input()
        startedstarcat = False
        cccc = 0
        for imfile,noisefile,psffile,band, j in \
                zip(snparams.image_name_search,snparams.image_name_weight,snparams.file_name_psf,snparams.band, range(len(snparams.band))):
            nozpt = copy(orig_nozpt)

            # if not '20151029' in imfile:
            #     continue

            #raw_input('stopped')
            # if round(float(snparams.mjd[j])) != 56935:
            #    continue
            try:
                self.field = imfile.split('/')[-1].split('-')[1].split('_')[0]
                self.ccdnum = imfile.split('/')[-1].split('_')[-1][:2]
                self.expnum = imfile.split('/')[-1].split('_')[1]

            except:
                self.ccdnum = np.nan
                self.field = np.nan
                self.expnum = np.nan

            # print imfile
            # raw_input('immmmmm')
            # if 'Y4' in imfile:
            #     print 'yfour inside now skipping'
            #     continue


            if self.fermigrid:
                self.laskerstarcat = os.path.join(self.outdir, 'stardata', 'lasker',
                                                  'CCD' + self.ccdnum + '_' + filt + 'band.starcat')
                print self.laskerstarcat
                if not startedstarcat:
                    try:
                        os.system('ifdh  mkdir '+os.path.join(self.outdir,'stardata'))
                    except:
                        pass
                    try:
                        os.system('ifdh mkdir '+os.path.join(self.outdir,'stardata','lasker'))
                    except:
                        pass
                    fout = open('tmp.txt', 'w')
                    print >> fout, '# ID\t RA\t DEC\t CATMAG\t FILTER\t EXPNUM\n'
                    fout.close()
                    startedstarcat = True

                    ls = os.popen('ifdh ls ' + self.laskerstarcat).read()
                    if len(ls) > 0:
                        self.laskerstarcat = 'lsc.txt'
                        #os.popen('ifdh rm ' + self.laskerstarcat)
                    else:
                        os.popen('ifdh cp --force=xrootd tmp.txt ' + self.laskerstarcat)
            #sys.exit()



            #smp_dict['mjd'][j] = float(snparams.mjd[j])
            #smp_dict['mjd'][i] = float(snparams.mjd[j])
            print imfile
            epochtime  = time.time()

            if j in badindices:
                continue
            if round(snparams.mjd[j]) == 0:
                #raw_input('mjdddd')
                continue

            #if round(snparams.mjd[j]) != 56636.:
            #    if snparams.mjd[j] < 57000.:
            #        continue
            #if round(snparams.mjd[j],2) != 56030.33:
            #    continue
            #raw_input('passed')

            if filt != 'all' and band not in filt:
                # print('filter %s not in filter list %s for image file %s'%(band,filt,imfile))
                # print 'filter %s,%s not in filter list for image file %s'%(band,filt,imfile)
                continue
            print imfile

            if cccc > 6000:
                #cccc += 1
                continue


            #raw_input()
            skysig=np.nan
            badflag = 0
            descriptiveflag = 0
            nozpt = copy(orig_nozpt)

            # longimfile = copy(imfile )
            # if self.fermigrid & self.worker:
            #     imfile = imfile.split('/')[-1]
            #     noisefile = noisefile.split('/')[-1]
            #     weightsfile = noisefile
            #     psffile = psffile.split('/')[-1]
            #
            #
            # try:
            #     self.ccdnum = imfile.split('/')[1].split('_')[1]
            #     self.field = imfile.split('/')[0].split('-')[1]
            # except:
            #     self.ccdnum = np.nan
            #     self.field = np.nan

            ########NEWWWWWWWWWWWWWWWW#######################
            if self.usefake:
                imfile = ''.join(imfile.split('.')[:-1]) + '+fakeSN.fits'

            if not self.oldformat:
                if self.usefake:
                    if not self.snparams.survey == 'PS1':
                        print imfile
                        changename = True
                        if changename:
                            if int(imfile[:8]) < 20140601:
                                imfile = imfile.replace('p1', 'Y1')
                                noisefile = noisefile.replace('p1', 'Y1')
                                psffile = psffile.replace('p1', 'Y1')
                            elif int(imfile[:8]) < 20150601:
                                imfile = imfile.replace('p1', 'Y2')
                                noisefile = noisefile.replace('p1', 'Y2')
                                psffile = psffile.replace('p1', 'Y2')
                            elif int(imfile[:8]) < 20160601:
                                imfile = imfile.replace('p1', 'Y3')
                                noisefile = noisefile.replace('p1', 'Y3')
                                psffile = psffile.replace('p1', 'Y3')
                            elif int(imfile[:8]) < 20170601:
                                imfile = imfile.replace('p1', 'Y4')
                                noisefile = noisefile.replace('p1', 'Y4')
                                psffile = psffile.replace('p1', 'Y4')
                else:
                    if not self.snparams.survey == 'PS1':
                        print imfile

                        #raw_input('imh')
                        if int(imfile[:8]) < 20140601:
                            imfile = imfile.replace('p1', 'Y1')
                            noisefile = noisefile.replace('p1', 'Y1')
                            psffile = psffile.replace('p1', 'Y1')
                        elif int(imfile[:8]) < 20150601:
                            imfile = imfile.replace('p1', 'Y2')
                            noisefile = noisefile.replace('p1', 'Y2')
                            psffile = psffile.replace('p1', 'Y2')
                        elif int(imfile[:8]) < 20160601:
                            imfile = imfile.replace('p1', 'Y3')
                            noisefile = noisefile.replace('p1', 'Y3')
                            psffile = psffile.replace('p1', 'Y3')
                        elif int(imfile[:8]) < 20170601:
                            imfile = imfile.replace('p1', 'Y4')
                            noisefile = noisefile.replace('p1', 'Y4')
                            psffile = psffile.replace('p1', 'Y4')
                # noisefile = ''.join(noisefile.split('.')[:-1])+'+fakeSN.fits'


            if 'Y4' in imfile:
                print 'yfour inside now skipping'
                continue



            imfile = os.path.join(rootdir, imfile)
            longimfile = copy(imfile)
            self.impath = '/'.join(imfile.split('/')[:-1])
            try:
                noisefile = os.path.join(rootdir, noisefile)
            except:
                noisefile = [os.path.join(rootdir, noisefile[0]), os.path.join(rootdir, noisefile[1])]

            psffile = os.path.join(rootdir, psffile)
            longpsffile = copy(psffile)
            if self.snparams.survey == 'DES':
                if self.usefake:
                    longnoisefile = copy(imfile.replace('+fakeSN.fits','.weight.fits'))
                else:
                    longnoisefile = copy(imfile.replace('.fits','.weight.fits'))
            else:
                longnoisefile = copy(noisefile)

            # print imfile
            # raw_input()
            print 'running zeropoints on ' + longimfile
            if self.fermilog:
                self.tmpwriter.appendfile('running zeropoints on ' + longimfile+'\n', self.fermilogfile)

            #print 'hereeeee'
            if self.fermigrid & self.worker:
                # print imfile
                # os.system('IFDH_CP_MAXRETRIES=1; ifdh cp ' + imfile + ' .')
                # print 'line 529 copied image files to here'
                # sys.exit()
                # print 'ifdh cp '+imfile+' .'

                ifdhls = os.popen('ifdh lss ' + imfile).read()
                #print ifdhls
                #print 'line 576'
                # raw_input()
                print ifdhls
                print ifdhls
                if (len(ifdhls) > 0):
                    if int(ifdhls.split()[-1]) == 0:
                        print 'EMTPY FILE!'*10
                        continue
                    ifdhlss = os.popen('ifdh lss ' + imfile.split('/')[-1]).read()
                    if not (len(ifdhlss)>0):
                        os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp --force=xrootd ' + imfile + ' .').read()
                    # imfilel = copy(imfilel)
                    imfile = imfile.split('/')[-1]
                    print 'imfile', imfile
                    # if self.usefake:
                    #     #if '.gz' in imfile:
                    #     print 'ifdh','IFDH_CP_MAXRETRIES=1; ifdh cp ' + imfilel.split('.fits.gz')[0]+ '+fakeSN.fits.gz' + ' .'
                    #     os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp ' + imfilel.split('.fits.gz')[0]+ '+fakeSN.fits.gz' + ' .').read()
                    #     #imfile = imfilel.split('/')[-1]
                    #     #else:
                    #     os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp ' + imfilel.split('.fits')[
                    #         0] + '+fakeSN.fits' + ' .').read()
                    #     imfile = imfilel.split('/')[-1]
                    # print 'IFDH_CP_MAXRETRIES=1; ifdh cp '+noisefile+' .'
                    print 'reading in noisefile',noisefile
                    ifdhlss = os.popen('ifdh lss ' + noisefile.split('/')[-1]).read()
                    if not (len(ifdhlss)>0):
                        os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp --force=xrootd ' + noisefile + ' .').read()
                    noisefile = noisefile.split('/')[-1]
                    weightsfile = noisefile
                    # print 'ifdh cp ' + psffile + ' .'
                    print 'reading in psffile'
                    ifdhlss = os.popen('ifdh lss ' + psffile.split('/')[-1]).read()
                    if not (len(ifdhlss)>0):
                        os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp --force=xrootd ' + psffile + ' .').read()
                    psffile = psffile.split('/')[-1]
                    # print 'copied all files'
                    # print os.popen('ifdh ls .').read()
                    # sys.exit()
                    #print 'here2'
                    # raw_input()
                else:
                    print 'file not found', imfile
                    # continue

                    ifdhls = os.popen('ifdh lss ' + imfile + '.fz').read()
                    #print ifdhls
                    #print 'line 610'
                    # raw_input()
                    if (len(ifdhls) > 0) :
                        if int(ifdhls.split()[-1] == 0) :
                            print 'EMPTY FILE!'*10
                            continue
                        a = os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp --force=xrootd ' + imfile + '.fz .').read()
                        print a
                        a = os.popen('funpack ' + imfile.split('/')[-1] + '.fz').read()
                        print a
                        # imfilel = copy(imfilel)
                        imfile = imfile.split('/')[-1]
                        print 'imfile', imfile
                        # if self.usefake:
                        #     #if '.gz' in imfile:
                        #     print 'ifdh','IFDH_CP_MAXRETRIES=1; ifdh cp ' + imfilel.split('.fits.gz')[0]+ '+fakeSN.fits.gz' + ' .'
                        #     os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp ' + imfilel.split('.fits.gz')[0]+ '+fakeSN.fits.gz' + ' .').read()
                        #     #imfile = imfilel.split('/')[-1]
                        #     #else:
                        #     os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp ' + imfilel.split('.fits')[
                        #         0] + '+fakeSN.fits' + ' .').read()
                        #     imfile = imfilel.split('/')[-1]
                        # print 'IFDH_CP_MAXRETRIES=1; ifdh cp '+noisefile+' .'
                        noisefile = longimfile.split('.fits')[0]+'.weight.fits'

                        print 'reading in noisefile',noisefile,'asdfasdf'
                        lfz = os.popen('ifdh lss ' + noisefile + '.fz').read()
                        if (len(lfz) > 0):
                            print 'there is fz'
                            os.popen('ifdh cp --force=xrootd ' + noisefile + '.fz .').read()
                            print 'funpacking noisefile'
                            os.popen('funpack ' + noisefile.split('/')[-1] + '.fz')
                        else:
                            print 'there is no fz'
                            os.popen('ifdh cp --force=xrootd ' + noisefile + ' .').read()
                        noisefile = noisefile.split('/')[-1]
                        weightsfile = noisefile
                        # print 'ifdh cp ' + psffile + ' .'
                        print 'reading in psf file'
                        os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp --force=xrootd ' + psffile + ' .').read()
                        # os.popen('funpack '+psffile)
                        psffile = psffile.split('/')[-1]
                        print 'done grabbing images'
                        # print 'copied all files'
                        # print os.popen('ifdh ls .').read()
                        # sys.exit()
                        #print 'here2'
                        # raw_input()
                    else:
                        print 'file not found', imfile
                        continue
                        # print 'grabbed sn files'
                        # sys.exit()
            try:
                self.ccdnum = imfile.split('/')[-1].split('_')[5].split('+')[0]
                self.field = imfile.split('/')[-1].split('_')[2].split('-')[1]
            except:
                self.ccdnum = np.nan
                self.field = np.nan

            try:
                print 'reading in'
                im = pyfits.getdata(imfile)
                hdr = pyfits.getheader(imfile)
                print 'done reading in'
                #print 'image shape', im.shape
            except:
                print 'Image is EMPTY, skipping star...'

                continue

            snparams.platescale = hdr[self.params.hdr_platescale_name]
            snparams.airmass = hdr[self.params.hdr_airmass_name]

            if filt != 'all' and band not in filt:
                # print('filter %s not in filter list %s for image file %s'%(band,filt,imfile))
                # print 'filter %s,%s not in filter list for image file %s'%(band,filt,imfile)
                continue

            # imfileloc = copy(imfile)

            # if fermigrid and worker:
            #    self.rootdir = '.'

            # print imfile
            # print imfileloc

            if not worker:
                if self.useweights:
                    weightsfile = os.path.join(self.rootdir, noisefile)
                else:
                    noisefile, maskfile = os.path.join(self.rootdir, noisefile[0]), os.path.join(self.rootdir,
                                                                                                 noisefile[1])

            if not os.path.exists(imfile):
                # print os.popen('ls -ltr *.fz').read()
                # print 'funpack %s.fz' % imfile
                # os.system('funpack %s.fz' % imfile)
                # sys.exit()
                if not os.path.exists(imfile + '.fz'):
                    print('Error : file %s does not exist' % imfile)
                    continue
                    print('Error : file %s does not exist' % imfile)
                    raise exceptions.RuntimeError('Error : file %s does not exist' % imfile)
                else:
                    print os.popen('ls -ltr *.fz').read()
                    print 'funpack %s.fz' % imfile
                    os.system('funpack %s.fz' % imfile)
                    # sys.exit()

            if self.useweights:
                if not os.path.exists(weightsfile):
                    os.system('gunzip %s.gz' % weightsfile)
                    if not os.path.exists(weightsfile):
                        os.system('funpack %s.fz' % weightsfile)
                        if not os.path.exists(weightsfile):
                            raise exceptions.RuntimeError('Error : file %s does not exist' % weightsfile)
            else:

                if not os.path.exists(noisefile):
                    os.system('gunzip %s.gz' % noisefile)
                    if not os.path.exists(noisefile):
                        os.system('funpack %s.fz' % noisefile)
                        if not os.path.exists(noisefile):
                            if not os.path.exists(noisefile[:-3]):
                                raise exceptions.RuntimeError('Error : file %s does not exist' % noisefile)
                            else:
                                noisefile = noisefile[:-3]
                if not os.path.exists(maskfile):
                    os.system('gunzip %s.gz' % maskfile)
                    if not os.path.exists(maskfile):
                        os.system('funpack %s.fz' % maskfile)
                        if not os.path.exists(maskfile):
                            if not os.path.exists(maskfile[:-3]):
                                raise exceptions.RuntimeError('Error : file %s does not exist' % maskfile)
                            else:
                                maskfile = maskfile[:-3]
            if not os.path.exists(psffile):
                if os.path.exists(psffile + '.gz'):
                    os.system('gunzip %s.gz' % psffile)
                elif not os.path.exists(psffile + '.fz'):
                    raise exceptions.RuntimeError('Error : file %s does not exist' % psffile)
                else:
                    os.system('/global/u1/d/dbrout/cfitsio/funpack %s.fz' % psffile)

            # if not nomask:
            #     if useweights:
            #         maskfile = os.path.join(self.rootdir, snparams.image_name_search[j])
            #         mask = pyfits.getdata(maskfile)



            ########NEWWWWWWWWWWWWWWWWENDDDDDDDDDD############

            try:

                self.rickfakestarfile = 'data/fixmagCoords_SN-'+self.field+'.dat'
            except:
                self.rickfakestarfile = ''
            #print filt
            #raw_input('filttt')
            if filt != 'all' and band not in filt:
                if verbose: print('filter %s not in filter list for image file %s'%(band,filt,imfile))
                print 'filter %s,%s not in filter list for image file %s'%(band,filt,imfile)
                continue
            cccc += 1


            if self.useweights:
                weights = pyfits.getdata(weightsfile)
                maskfile = weightsfile
                mask = weights
                mask[mask<0.000000001] = 0
                mask[mask>=0.000000001]= 1
            else:
                noise = pyfits.getdata(noisefile)
                mask = pyfits.getdata(maskfile)
                mask[mask>0] = -1000
                mask[mask==0] = 1
                mask[mask==-1000] = 0
                weights = 1./noise**2 * mask
                #weights = pyfits.getdata(weightsfile)


            psf = pyfits.getdata(psffile)

            if params.weight_type.lower() == 'ivar':
                print 'ivar'
                #raw_input()
                noise = np.sqrt(1/noise)
            elif params.weight_type.lower() != 'noise':
                raise exceptions.RuntimeError('Error : WEIGHT_TYPE value %s is not a valid option'%params.WEIGHT_TYPE)
            # if nomask:
            #     if self.useweights:
            #         mask = np.zeros(np.shape(weights))
            #         maskcols = np.where((weights < 0) |
            #                             (np.isfinite(weights) == False))
            #         mask[maskcols] = 100.0

            wcsworked = True
            try:
                w =wcs.WCS(imfile)
            except:
                wcsworked = False
                hl = pf.open(imfile)
                import starlink.Ast as Ast
                import starlink.Atl as Atl
                fitschan = Ast.FitsChan( Atl.PyFITSAdapter(hl[ 0 ]) )
                encoding = fitschan.Encoding
                wcsinfo = fitschan.read()
                w = wcs.WCS(imfile+'old')

            
            if wcsworked:
                results =  w.wcs_pix2world(np.array([[0,0]]), 0)
                ra1, dec1 = results[0][0], results[0][1]
                results2 =  w.wcs_pix2world(np.array([[snparams.nxpix-1,
                                   snparams.nypix-1]]), 0)
                ra2, dec2 =results2[0][0], results2[0][1]
            else:
                xpix = [0]
                ypix = [0]
                radtodeg = 360/(2*3.14159)
                results =  wcsinfo.tran([xpix,ypix])

                ra1, dec1 = results[0]*radtodeg, results[1]*radtodeg
                results2 =  wcsinfo.tran([[snparams.nxpix-1], [snparams.nypix-1]])
                ra2, dec2 =results2[0]*radtodeg, results2[1]*radtodeg

            ra_high = np.max([ra1,ra2])
            ra_low = np.min([ra1,ra2])
            dec_high = np.max([dec1,dec2])
            dec_low = np.min([dec1,dec2])

            if not self.usefake:
                #try:
                snparams.RA = float(snparams.ra)
                snparams.DECL = float(snparams.decl)
                #except:
                    #snparams.RA = astCoords.hms2decimal(snparams.ra, ':')
                    #snparams.DECL = astCoords.dms2decimal(snparams.decl, ':')

            else:
                #try:
                if self.exactpos:
                    snparams.RA = float(snparams.fake_ra)
                    snparams.DECL = float(snparams.fake_dec)
                else:
                    snparams.RA = float(snparams.ra)
                    snparams.DECL = float(snparams.decl)
                # except:
                #     if self.exactpos:
                #         snparams.RA = astCoords.hms2decimal(snparams.fake_ra,':')
                #         snparams.DECL = astCoords.dms2decimal(snparams.fake_dec,':')
                #     else:
                #         snparams.RA = astCoords.hms2decimal(snparams.ra,':')
                #         snparams.DECL = astCoords.dms2decimal(snparams.decl,':')


            if params.forceradec.lower() == 'true':
                #print 'forcingggg'
                #raw_input()
                #params.fra = 161.1183929
                #params.fdec = 57.8055115
                print float(params.fra),float(params.fdec)
                xsn,ysn = zip(*w.wcs_world2pix(np.array([[float(params.fra),float(params.fdec)]]), 0))
            else:
                xsn,ysn = zip(*w.wcs_world2pix(np.array([[snparams.RA,snparams.DECL]]), 0))

            xsn = xsn[0]
            ysn = ysn[0]
            if self.snparams.survey == 'PS1':
                xsn += 1
                ysn += 1
            #if self.snparams.survey == 'DES':
            #    xsn += 1
            #    ysn += 1
            #print xsn,ysn
            #raw_input('stopppppppppp')

            #if self.snparams.psf_model.lower() == 'psfex':
            #    xsn = xsn
            #    ysn = ysn
            #print 'usefake',self.usefake
            #print 'snra sndec',snparams.RA,snparams.DECL
            #print 'xsn ysnnnnnnnnnnnnnnnnn', xsn, ysn
            #print 'snparams.npix ypix',snparams.nxpix,snparams.nypix
            testoffset = False
            offsetx = .6
            offsety = .3
            if testoffset:
                xsn += offsetx
                ysn += offsety



            i += 1
            smp_dict['image_filename'][i] = imfile
            smp_dict['zpt_file'][i] = 'na'
            smp_dict['psf_filename'][i] = psffile
            smp_dict['weight_filename'][i] = weightsfile

            if xsn < 0 or ysn < 0 or xsn > snparams.nxpix-1 or ysn > snparams.nypix-1:
                print "Error : SN Coordinates %s,%s are not within image"%(snparams.ra,snparams.decl)
                badflag = 1
                descriptiveflag = 1 #
                #raw_input()


            #print 'about to do zeropoints'
            #sys.exit()
            if snparams.psf_model.lower() == 'daophot':
                #self.psf = rdpsf.rdpsf(psffile)[0]/10.**(0.4*(25.-magzpt))
                self.psf = rdpsf.rdpsf(psffile)[0]
                self.psf = self.psf/np.sum(self.psf)
                self.psfcenter = None
                if params.build_psf == 'yes':
                    #self.rdnoise = hdr[params.rdnoise_name]
                    #self.gain = hdr[params.gain_name]  # 1
                    #cols = (starglobalras > ra_low) & (starglobalras < ra_high) & (starglobaldecs > dec_low) & (
                    #starglobaldecs < dec_high)

                    #if not len(cols):
                    #    raise exceptions.RuntimeError("Error : No stars in image!!")
                    #print starcat.__dict__
                    #raw_input()
                    mag_star = starcat.mag
                    #mag_star = starcat.mag[cols]
                    coords = zip(*w.wcs_world2pix(np.array(zip(starglobalras, starglobaldecs)), 0))
                    #coords = zip(*w.wcs_world2pix(np.array(zip(starcat.ra[cols], starcat.dec[cols])), 0))
                    x_star, y_star = [], []

                    for xval, yval in zip(*coords):
                        print xval,yval
                        x_star += [xval]
                        y_star += [yval]

                    x_star, y_star = np.array(x_star), np.array(y_star)
                    cols = (x_star > 0) & (y_star > 0)
                    x_star = x_star[cols]
                    y_star = y_star[cols]
                    mag_star = mag_star[cols]

                    print x_star

                    try:
                        mag,magerr,flux,fluxerr,sky,skyerr,badflagx,outstr = \
                            aper.aper(im,x_star,y_star,apr = params.fitrad,verbose=False)
                    except:
                        print sys.exc_info()
                        if not os.path.exists(self.isdonedir):
                            os.mkdirs(self.isdonedir)
                        os.system(
                            'touch ' + os.path.join(self.isdonedir, self.snparams.snfile.split('/')[-1].split('.')[0] + '.done'))
                        sys.exit()

                    # self.rdnoise = hdr[params.rdnoise_name]
                    # self.gain = hdr[params.gain_name]
                    # #if not os.path.exists(psffile) or params.clobber_psf == 'yes':
                    # gauss,psf,magzpt = getpsf.getpsf(im,x_star,y_star,mag,sky,
                    #                                  self.rdnoise,self.gain,
                    #                                  range(len(x_star)),params.fitrad,params.fitrad-1.,
                    #                                  psffile)
                    # hpsf = pyfits.getheader(psffile)
                    # self.psf = psf
                    # self.psfcenter = None
                    #else:
                    #    print('PSF file exists.  Not clobbering...')
                    #    hpsf = pyfits.getheader(psffile)
                    #    magzpt = hpsf['PSFMAG']

                elif nozpt:
                    print 'getting stars'
                    #self.rdnoise = hdr[params.rdnoise_name]
                    #self.gain = hdr[params.gain_name] #1
                    cols = (starglobalras > ra_low) & (starglobalras < ra_high) & (starglobaldecs > dec_low) & (starglobaldecs < dec_high)

                    #print starglobalras[cols]


                    if not len(cols):
                        print 'Error: No stars in image!'
                        continue
                        raise exceptions.RuntimeError("Error : No stars in image!!")
                    
                    #mag_star = starcat.mag[cols]
                    mag_star = starglobalmags[cols]
                    coords = zip(*w.wcs_world2pix(np.array(zip(starglobalras[cols], starglobaldecs[cols])), 0))
                    x_star,y_star = [],[]

                    for xval,yval in zip(*coords):
                        x_star += [xval]
                        y_star += [yval]

                    x_star,y_star = np.array(x_star),np.array(y_star)
                    #if not doncentroid:
                    #    x_star,y_star = cntrd.cntrd(im,x_star,y_star,params.cntrd_fwhm)
                    #    newra,newdec = zip(*w.wcs_pix2world(np.array(zip(xstar_,y_star)),0))

                    mag,magerr,flux,fluxerr,sky,skyerr,badflagx,outstr = \
                        aper.aper(im,x_star,y_star,apr = params.fitrad,verbose=False)

                    hpsf = pyfits.getheader(psffile)
                    magzpt = hpsf['PSFMAG']

                else:

                    hpsf = pyfits.getheader(psffile)
                    magzpt = hpsf['PSFMAG']

                    #self.rdnoise = hdr[params.rdnoise_name]
                    #self.gain = hdr[params.gain_name]



            # begin taking PSF stamps

            #print snparams.psf_model.lower()
            elif snparams.psf_model.lower() == 'psfex':
                print snparams.psf_model.lower()
                hdulist = fits.open(psffile)
                hdulist.info()
                try:
                    psf_fwhm = hdulist[1].header[self.params.hdr_psf_fwhm]
                    #fwhm = open(psffile, 'r').readline().split(self.params.psf_fwhm)[1].split('/')[0].split('=')[1]
                    #print psf_fwhm
                    #raw_input('fwhm')
                except:
                    print 'Could not find pwf_fwhm in fits header'
                    psf_fwhm = np.nan
                print snparams.RA,xsn,snparams
                try:
                    print psffile
                    #raw_input()
                    self.psf, self.psfcenter = self.build_psfex(psffile,xsn,ysn,imfile,psfexworked=psfexworked)
                    self.psf = self.psf/np.sum(self.psf)
                except:
                    self.psf = 0
                    badflag = 1
                    descriptiveflag = 2
                    print 'badflagpsf'*100
            # elif snparams.psf_model.lower() == 'daophot':
            #     self.psf = rdpsf.rdpsf(psffile)[0]/10.**(0.4*(25.-magzpt))
            #     #self.psf = rdpsf.rdpsf(psffile)[0]
            #     #self.psf = self.psf/np.sum(self.psf)
            #     self.psfcenter = None
            else:
                raise exceptions.RuntimeError("Error : PSF_MODEL not recognized!")

            try:
                self.rdnoise = hdr[params.rdnoise_name]
                self.gain = hdr[params.gain_name]
            except:
                print 'could not find readnoise or gain in header setting to default 1'
                self.rdnoise = 1.
                self.gain = 1.
            if self.snparams.survey == 'PS1':
                self.gain = 1.

            print hdr.keys()
            #raw_input()

            if self.snparams.survey == 'DES':
                gaina = hdr['GAINA']
                gainb = hdr['GAINB']
                q1 = [1,4,5,8,9,13,14,15,19,20,21,25,26,27]
                q2 = [32,33,34,39,40,41,45,46,47,51,52,56,57,60]
                q3 = [3,6,7,11,12,16,17,18,22,23,24,29,30,31]
                q4 = [36,37,38,42,43,44,48,49,50,54,55,58,59,62]

                try:
                    ttt = int(self.ccdnum)
                except:
                    self.ccdnum = '0'

                if int(self.ccdnum) in q1:
                    if ysn<1025:
                        self.gain = gaina
                    else:
                        self.gain = gainb
                elif int(self.ccdnum) in q2:
                    if ysn>1024:
                        self.gain = gaina
                    else:
                        self.gain = gainb
                elif int(self.ccdnum) in q3:
                    if ysn < 1025:
                        self.gain = gaina
                    else:
                        self.gain = gainb
                elif int(self.ccdnum) in q4:
                    if ysn > 1024:
                        self.gain = gaina
                    else:
                        self.gain = gainb
                else:
                    if xsn < 2049:
                        if ysn > 1024:
                            self.gain = gaina
                        else:
                            self.gain = gainb
                    else:
                        if ysn < 1025:
                            self.gain = gaina
                        else:
                            self.gain = gainb

            #print self.gain
            #raw_input('gain')
            if self.snparams.survey == 'DES':
                im += 10000.
            mjdoff = 0.
            mjdslopeinteroff = 0.
            gogo = True
            print '1794',nozpt
            if not nozpt:
                if True:
                #try:
                    if doglobalstar:
                        if self.snparams.survey == 'PS1':

                            name = longimfile.split('/')[-1][:-8]
                            #print name
                            zpt_file = os.path.join(self.outdir,'stardata',filt, name + '_zptstardata.npz')

                        elif dogalsimpixfit:
                            zpt_file = os.path.join(longimfile.split('.')[-2] + '_' + str(filt) + 'band_dillonzptinfo_galsimglobalstar.npz')
                        else:
                            zpt_file = os.path.join(longimfile.split('.')[-2] + '_'+str(filt)+'band_dillonzptinfo_globalstar.npz')
                    else:
                        if self.snparams.survey == 'PS1':

                            name = longimfile.split('/')[-1][:-8]
                            #print name
                            zpt_file = os.path.join(self.outdir,'stardata',filt, name + '_zptstardata.npz')
                        else:
                            zpt_file = os.path.join(longimfile.split('.')[-2] + '_'+str(filt)+'band_dillonzptinfo.npz')
                    print zpt_file
                    self.ishead = True
                    if self.fermigrid and self.worker:
                        if self.ishead:
                            os.popen('cp '+zpt_file+' .')
                            zpt_file = zpt_file.split('/')[-1]
                        else:
                            ls = os.popen('ifdh ls ' + zpt_file).read()
                            print ls
                            if len(ls) > 0:
                                os.popen('ifdh cp --force=xrootd '+zpt_file+' .')
                                zpt_file = zpt_file.split('/')[-1]
                            else:
                                print '1812 setting nozpt to True'
                                nozpt = True
                                gogo = False

                    if gogo:
                        try:
                        #if True:
                            print zpt_file

                            zptdata = np.load(zpt_file.replace('+fakeSN','')) #load previous zpt information

                            print zptdata.keys()
                            zpt = zptdata['fit_zpt']
                            zpterr = zptdata['fit_zpt_std']
                            #rmsaddin = zptdata['rmsaddin']
                            rmsaddin=1
                            mjdoff = zptdata['mjdoff']
                            mjdslopeinteroff = zptdata['mjdslopeinteroff']
                            if zptdata['fit_zpt'] == 0.:
                                print 'loaded bad zpt, continuing...'
                                nozpt = True
                                #continue
                            #print 'thisworked'
                        except:
                            print '2099 setting nozpt to True'
                            nozpt = True
                # except:
                #     print('Warning : IMAGE_ZPT field does not exist!  Calculating')
                #     nozpt = True
                #sys.exit()
            if nozpt:
                if not wehavestarcat:
                    starcat = txtobj(self.starcatfile, useloadtxt=True)
                    print 'done reading in starcatfile'
                    wehavestarcat = True
                    starcat.bigra = np.array(starcat.__dict__['RA'][1:], dtype='float')
                    starcat.bigdec = np.array(starcat.__dict__['DEC'][1:], dtype='float')
                    starcat.bigmag = np.array(starcat.__dict__['MAG_PSF_MEAN_%s' % filt.upper()][1:], dtype='float')
                    starcat.bigid = np.array(starcat.__dict__['MATCH_OBJECT_ID'][1:], dtype='float')
                    starcat.id = starcat.bigid
                    starcat.ra = starcat.bigra
                    starcat.dec = starcat.bigdec
                    starcat.mag = starcat.bigmag

                #self.rdnoise = hdr[params.rdnoise_name]
                #self.gain =  hdr[params.gain_name]

                # cols = np.where((starcat.bigra > ra_low) &
                #                 (starcat.bigra < ra_high) &
                #                 (starcat.bigdec > dec_low) &
                #                 (starcat.bigdec < dec_high))[0]

                cols = (starglobalras > ra_low) & (starglobalras < ra_high) & (starglobaldecs > dec_low) & (starglobaldecs < dec_high)
                tras = starglobalras
                tdecs = starglobaldecs
                tids = starglobalids
                mag_star = starglobalmags

                cols = (starcat.ra > ra_low) & (starcat.ra < ra_high) & (starcat.dec > dec_low) & (starcat.dec < dec_high)
                tras = starcat.ra
                tdecs = starcat.dec
                tids = starcat.id
                mag_star = starcat.mag



                if not len(cols):
                    #print 'Error: no stars in image!'
                    #continue
                    #print starcat.ra
                    #print 'ralow,rahi',ra_low,ra_high
                    raise exceptions.RuntimeError("Error : No stars in image!!")
                # try:
                #     if band.lower() == 'g':
                #         mag_star = starcat.mag_g[cols]
                #     elif band.lower() == 'r':
                #         mag_star = starcat.mag_r[cols]
                #     elif band.lower() == 'i':
                #         mag_star = starcat.mag_i[cols]
                #     elif band.lower() == 'z':
                #         mag_star = starcat.mag_z[cols]
                #     else:
                #         raise Exception("Throwing all instances where mag_%band fails to mag. Should not appear to user.")
                # except:
                #     mag_star = starcat.bigmag[cols]



                if wcsworked:
                    # x_starold,y_starold = [],[]
                    x_star,y_star = [],[]
                    # coords = zip(*w.wcs_world2pix(np.array(zip(starcat.ra[cols],starcat.dec[cols])),0))
                    # tras = []
                    # tdecs = []
                    # for ira,idec in zip(starcat.ra[cols],starcat.dec[cols]):
                    #     tras.append(ira)
                    #     tdecs.append(idec)
                    #
                    # for xval,yval in zip(*coords):
                    #     x_starold += [xval]
                    #     y_starold += [yval]
                    # #doglobalstar = False
                    # if doglobalstar:
                    #     print 'doing globalstars'
                    #     #raw_input()
                    #     tras = []
                    #     tdecs = []
                    #     try:
                    #         starcat.objid += 0.
                    #     except:
                    #         starcat.objid = np.arange(len(starcat.mag))
                    # for id,r,d in zip(starglobalids,starglobalras,starglobaldecs):
                    #         # tra = starglobalras[starglobalids == ide]
                    #         # tdec = starglobaldecs[starglobalids == ide]
                    #         # tras.extend(tra)
                    #         # tdecs.extend(tdec)
                    #     coords = zip(*w.wcs_world2pix(np.array(zip(r,d)),0))
                    #     for xval,yval in zip(*coords):
                    #         x_star += [xval]
                    #         y_star += [yval]

                    x_star, y_star = zip(*w.wcs_world2pix(np.array(zip(starglobalras,starglobaldecs)),0))
                    x_star, y_star = zip(*w.wcs_world2pix(np.array(zip(starcat.ra,starcat.dec)),0))


                    # else:
                    #     x_star = x_starold
                    #     y_star = y_starold
                    #     tra = starcat.ra[cols]
                    #     tdec = starcat.dec[cols]

                else:
                    coords = wcsinfo.tran([starcat.ra[cols]/radtodeg,starcat.dec[cols]/radtodeg],False)

                #for xi,yi,xo,yo in zip(x_star,y_star,x_starold,y_starold):
                #    print xi,xo,yi,yo
                #print 'printed globalstar locations'
                #sys.exit()

                x_star1,y_star1 = np.array(x_star),np.array(y_star)
                print x_star1
                badflagarr = (x_star1 < 0) | (y_star1 < 0) | (x_star1 - im.shape[1] > 0) | (y_star1 - im.shape[0] > 0)
                x_star1[badflagarr] = 100.
                y_star1[badflagarr] = 100.
                x_star1 = x_star1[~badflagarr]
                y_star1 = y_star1[~badflagarr]
                tras = tras[~badflagarr]
                tdecs = tdecs[~badflagarr]
                tids = tids[~badflagarr]
                mag_star = mag_star[~badflagarr]
                #print x_star1[~badflagarr]

                #raw_input()
                #for xx in x_star1:
                #    print xx
                #print im.shape
                #try:
                mag,magerr,flux,fluxerr,sky,skyerr,badflagx,outstr = \
                    aper.aper(im,x_star1,y_star1,apr = 30,verbose=False,ignoreneg=True)
                #except:
                #    print 'failed aper'
                #    badflag = 1
                #badflagx[badflagarr] = 1
                #print sky
                #raw_input()
                #I REMOVED CENTROIDING BECAUSE WE NOW FIND A GLOBAL RA AND DEC FOR THE STAR SIMILARLY TO THE SN
                #newx_star,newy_star = cntrd.cntrd(im,x_star1,y_star1,params.cntrd_fwhm)
                # newx_star,newy_star = x_star1,y_star1
                # if wcsworked:
                #     newra,newdec = zip(*w.wcs_pix2world(np.array(zip(newx_star,newy_star)),0))
                # else:
                #     ccc = wcsinfo.tran([x_star,y_star])
                #     newra,newdec = ccc[0]*radtodeg,ccc[1]*radtodeg

                # if self.dobigstarcat:
                #     catra,catdec,mag_starwrong = self.getProperCatRaDec(starcat.ra[cols],starcat.dec[cols])
                #     #print 'got proper cat ra and dec'
                #     #raw_input()
                # else:
                # catra,catdec = starcat.ra[cols],starcat.dec[cols]
                # deltara = catra - newra
                # deltadec = catdec - newdec
                # deltamjd = copy(deltara)*0. + snparams.mjd[j]
                #
                # jjj = abs(deltara) < 1
                # deltaram = deltara[jjj]
                # ram = np.array(newra)[jjj]
                # mm, s, iii = meanclip.meanclip( deltaram, clipsig = 3., maxiter = 8, returnSubs=True)
                # bra = np.polyfit(ram[iii], deltaram[iii], 0)
                # mra,cra = np.polyfit(ram[iii], deltaram[iii], 1)
                #
                # jjj = abs(deltadec) < 1
                # deltadecm = deltadec[jjj]
                # decm = np.array(newdec)[jjj]
                # mm, s, iii = meanclip.meanclip( deltadecm, clipsig = 3., maxiter = 8, returnSubs=True)
                # bdec = np.polyfit(decm[iii], deltadecm[iii], 0)
                # mdec,cdec = np.polyfit(decm[iii], deltadecm[iii], 1)
                #
                #
                # mjdoff = [bra[0],bdec[0]]
                # mjdslopeinteroff = [[mra,cra],[mdec,cdec]]
                #
                # self.deltaras.extend(deltara)
                # self.deltadecs.extend(deltadec)
                # self.deltamjds.extend(deltamjd)
                # self.ras.extend(catra)
                # self.decs.extend(catdec)
                # self.x_stars.extend(x_star)
                # self.y_stars.extend(y_star)
                #
                # self.airmasses.extend(starcat.ra[cols]*0. + round(snparams.airmass,2))

                self.psf = self.psf/np.sum(self.psf)
                print badflag
                #print badflag.shape
                #raw_input('we are before the zpt calc')
                skipactualzeropoint = False
                #sys.exit()
                #print 'nozpt',nozpt
                #sys.exit()
                print '1854',nozpt
                print 'exiting gracefully'
                #sys.exit()

                #print 'comparo',len(starglobalras),len(tras),len(starcat.ra)
                #raw_input()

                if not nozpt:
                    skipactualzeropoint = True
                if not skipactualzeropoint:
                    #x_star1, y_star1 = cntrd.cntrd(im, x_star1, y_star1, params.cntrd_fwhm)
                    #print x_star1[0]-x_star1n[0]
                    #raw_input('stopped')
                    #sys.exit()
                    zpttime = time.time()
                    # print x_star1.shape,tras.shape,mag_star.shape,mag.shape,sky.shape,badflagx.shape
                    # print min(mag_star),np.median(mag_star),max(mag_star)
                    # raw_input()

                    dosextractor = True
                    if dosextractor:
                        print imfile
                        print snparams.snfile.split('/')[-1].replace('.dat','')

                        sexsky, sexrms, bkgrnd, bkgrndrms = runsextractor.getsky_and_skyerr(imfile,weightsfile, im, 500 + 50,
                                                                         500 - 50,
                                                                         500 + 50, 500 - 50,
                                                                         snparams.survey, bigreturn=True,
                                                            index=snparams.snfile.split('/')[-1].replace('.dat',''))

                        print sexsky, sexrms

                        #dt.save_fits_image(im-sexsky,'testsky.fits',go=True)
                        # raw_input('sextractor')
                        #bkgrnd = pf.getdata(imfile + '.background')
                        # print bkgrnd.shape
                        #sexsky = np.mean(bkgrnd[ylow:yhi, xlow:xhi].ravel()) * scalefactor
                        #bkgrndrms = pf.getdata(imfile + '.background_rms')
                        #sexrms = (np.mean(bkgrndrms[ylow:yhi, xlow:xhi].ravel() ** .5) * scalefactor) ** 2.

                    else:
                        bkgrnd = None
                        bkgrndrms = None

                    #print sexsky,sexrms,sky,skyerr
                    #raw_input('sss')
                    zpt,zpterr,zpt_file, rmsaddin, thisra,thisdec, thisids = self.getzpt(x_star1,y_star1,tras,tdecs,tids,mag,sky,skyerr,snparams.mjd[j],
                                         badflagx,mag_star,im,weights,mask,maskfile,weightsfile,psffile,imfile,w,snparams,params.substamp,mjdoff,mjdslopeinteroff,j,
                                         longimfile,bkgrnd,bkgrndrms,psf=self.psf,mjd=str(float(snparams.mjd[j])))
                    print 'zpttime',time.time()-zpttime

                    if zpt == 0:
                        badflag = 1
                        descriptiveflag = 4



            else:
                if doglobalstar:
                    zpt_file = imfile.split('.')[-2].replace('+fakeSN','') + '_'+str(filt)+'band_dillonzptinfo_globalstar.npz'
                else:
                    zpt_file = imfile.split('.')[-2].replace('+fakeSN','') + '_'+str(filt)+'band_dillonzptinfo.npz'

                zptdata = np.load(zpt_file) #load previous zpt information
                zpt = zptdata['fit_zpt']
                zpterr = zptdata['fit_zpt_std']
                mjdoff = zptdata['mjdoff']
                mjdslopeinteroff = zptdata['mjdslopeinteroff']
                rmsaddin = zptdata['rmsaddin']
                thisra = zptdata['thisra']
                print 'getting rasssssss'
                thisdec = zptdata['thisdec']
                thisids = zptdata['thisids']

            try:
                #if True:
                if len(thisra) < params.minzptstars:
                    print 'COULD NOT GET GOOD FIT OF ZEROPOINT... N stars is too small. Nstars=',len(thisra)
                    scalefactor = 1.
                    badflag = 1
                    descriptiveflag = 512
                if zpterr / np.sqrt(float(len(thisra))) > 0.01:
                    badflag = 1
                    print 'COULD NOT GET GOOD FIT OF ZEROPOINT... SCATTER/SQRT(N) LARGER THAN .01 MAGS'
                    scalefactor = 1.
                    badflag = 1
                    descriptiveflag = 1024
            except:
                badlfag = 1
                print 'COULD NOT GET GOOD FIT OF ZEROPOINT... N stars is too small'
                scalefactor = 1.
                descriptiveflag = 512

            dotestoff = False
            if zpt == 0:
                print 'zerpoint badflag'
                badflag = 1
                descriptiveflag = 4

            if dotestoff:
                self.teststarpos(self.rickfakestarfile,w,zpt,sky,skyerr,im,weights,mask,psffile,imfile,snparams,params.substamp,snparams.zp[j],psf=self.psf)


            if not ('firstzpt' in locals()): firstzpt = 31. ####firstzpt = zpt
            if self.usediffimzpt:
                scalefactor = 10**(-.4*(snparams.zp[j]-firstzpt))
            else:
                if zpt != 0.0 and np.min(self.psf) > -10000 and np.isfinite(zpt):
                    scalefactor = 10.**(-0.4*(zpt-firstzpt))
                if zpt == 0.:
                    #aw_input('zpt badflag')
                    badflag = 1
                    scalefactor = 1.
                    descriptiveflag = 4

            print 'zpt',zpt,'nominal',firstzpt
            print 'scalefactor',scalefactor
            #print 'before',np.max(im)
            #im *= self.gain


            #SCALING THE iMAGE HEREEEEE
            im *= scalefactor


            #print 'after',np.max(im)
            #raw_input()

            #if self.useweights:
            #    #raw_input('USING WEIGHTS SETTIN')
            #    im[np.where(mask != 0)] =-999999.0

            badflagd = 0
            if dailyoff:
                print 'dailyoff'
                raw_input()
                if mjdoff[0] > .5:
                    raw_input('mjdoff badflag')
                    badflagd = 1
                if mjdoff[1] > .5:
                    raw_input('mjdoff badflag')
                    badflagd = 1

                xsn,ysn = zip(*w.wcs_world2pix(np.array([[snparams.RA+mjdoff[0],snparams.DECL+mjdoff[1]]]), 0))
                xsn = xsn[0]
                ysn = ysn[0]
                opsf = copy(self.psf)
                opsfcenter = copy(self.psfcenter)

                self.psf, self.psfcenter= self.build_psfex(psffile,xsn,ysn,imfile,psfexworked=psfexworked)
                self.psf = self.psf/np.sum(self.psf)


            zptmin = float(np.array(self.params.zptmin)[np.array(self.params.cutinfobands,dtype='str') == band][0])
            skymin = float(np.array(self.params.skymin)[np.array(self.params.cutinfobands,dtype='str') == band][0])
            skymax = float(np.array(self.params.skymax)[np.array(self.params.cutinfobands,dtype='str') == band][0])
            skyerrmax = float(np.array(self.params.skyerrmax)[np.array(self.params.cutinfobands,dtype='str') == band][0])

            print 'zptmin',zptmin
            if zpt < zptmin:
                badflag = 1
                descriptiveflag = 4
                print ('zpt '+str(zpt)+' less than zptmin '+str(zptmin)+'\n')*10


            if xsn > 25 and ysn > 25 and xsn < snparams.nxpix-25 and ysn < snparams.nypix-25 and np.isfinite(scalefactor):
                index += 1
                radius1 = 5
                radius2 = 8 
                fff = float(snparams.psf[j])
                skyrad=[radius1*fff,radius2*fff]
                np.set_printoptions(threshold=50000)
                #for i in im[:,round(ysn)]:
                #    print i
                magsn,magerrsn,fluxsn,fluxerrsn,skysn,skyerrsn,badflagaper,outstr = aper.aper(im,xsn,ysn,apr = 60.,skyisempty=True,verbose=False)#,skyrad=skyrad)
                print imfile
                print 'hhh2',skysn,skyerrsn,badflag
                #raw_input('checking sky and skyerr')
                #magsn,magerrsn,fluxsn,fluxerrsn,skysn,skyerrsn2,badflag,outstr = aper.aper(im,xsn,ysn,apr = 15.,verbose=False)#,skyrad=skyrad)
                #print skyerrsn,skyerrsn2
                #raw_input()

                #raw_input('skysn'+str(skysn))
                mygain = ((1/skyerrsn**2)*skysn)


                if badflagd == 1:
                    #raw_input('badflagd')
                    badflag = 1

                # if np.sum(mask[ysn-params.fitrad:ysn+params.fitrad+1,xsn-params.fitrad:xsn+params.fitrad+1]) != 0:
                #     if self.useweights:
                #         badflag = 1
                #         #raw_input('mask badflag')
                #         print 'mask badflag'



                #if skysn < -1e5:
                #    badflag = 1
                #    #raw_input('skysn badflag')
                #    #print 'skysn badflag'

                if not badflag:
                    stampsize = 40

                    if ysn-stampsize < 0:
                        ylow = int(0)
                    else:
                        ylow = int(ysn-stampsize)
                    if ysn+stampsize>im.shape[0]:
                        yhi = int(im.shape[0]-1)
                    else:
                        yhi = int(ysn+stampsize)
                    if xsn-stampsize < 0:
                        xlow = int(0)
                    else:
                        xlow = int(xsn-stampsize)
                    if xsn+stampsize>im.shape[1]:
                        xhi = int(im.shape[1]-1)
                    else:
                        xhi = int(xsn+stampsize)

                    mean,st,vals = sigma_clip.meanclip(im[ylow:yhi,xlow:xhi],clipsig = 4, maxiter = 8)

                    # plt.clf()
                    # plt.hist(im[ylow:yhi,xlow:xhi].ravel(),bins=np.arange(5405,5600,10))
                    #
                    # plt.savefig('testimvar.png')

                    #skysig=1.48*np.median(abs(vals-np.median(vals)))
                    skysig = st
                    mysky = mean
                    #mysky = np.median(vals)
                    #print mysky
                    #raw_input('mysky above')
                    mygain = (np.sqrt(mysky)/(skysig))**2
                    mygainsn =  (np.sqrt(skysn)/(skyerrsn))**2
                    #print mygain,mygainsn,hdr['GAINA'],hdr['GAINB']
                    print imfile
                    dosextractor = True
                    if dosextractor:
                        sexsky,sexrms,bkgrnd,bkgrndrms = runsextractor.getsky_and_skyerr(imfile,weightsfile,im,xlow-30+stampsize,xhi+30-stampsize,
                                                                        ylow-30+stampsize,yhi+30-stampsize,snparams.survey,bigreturn=True,index=snparams.snfile.split('/')[-1].replace('.dat',''))

                        print sexsky,sexrms
                        #raw_input('sextractor')
                        #bkgrnd = pf.getdata(imfile + '.background')
                        #print bkgrnd.shape
                        sexsky = np.mean(bkgrnd[ylow:yhi,xlow:xhi].ravel()) #* scalefactor
                        #bkgrndrms = pf.getdata(imfile + '.background_rms')
                        sexrms = np.mean(bkgrndrms[ylow:yhi,xlow:xhi].ravel())#* scalefactor
                        #sexsky *= scalefactor
                        #sexrms *= scalefactor
                        print snparams.mjd[i]
                        print 'sextractor sky ',sexsky,'sextractor rms', sexrms#scalefactor
                        print 'mysky',mysky,'myskysig',skysig
                        temp, mysexskysig, vals = sigma_clip.meanclip(im[ylow:yhi, xlow:xhi]-bkgrnd[ylow:yhi,xlow:xhi], clipsig=4, maxiter=8)
                        print 'mysexskysig',mysexskysig
                        #print 'skysn',skysn,'skyerrsn',skyerrsn
                        #raw_input('youyo')

                    skyvals = im[ylow:yhi,xlow:xhi].ravel()
                    #print im.shape
                    #raw_input('xlow'+str(xlow)+' xhi'+str(xhi)+' ylow'+str(ylow)+' yhi'+str(yhi))

                    pk = pkfit_norecent_noise_smp.pkfit_class(im,self.psf,self.psfcenter,self.rdnoise,self.gain,weights,mask)
                    #pk = pkfit_norecent_noise_smp.pkfit_class(im,self.gauss,self.psf,self.rdnoise,self.gain,noise,mask)
                    #try:
                    if self.snparams.survey == 'PS1':
                        # save_fits_image(mask,'test/fullmask.fits')
                        # print 'maskfile',maskfile
                        # raw_input('printed maskfile')
                        print 'getting multiple image stamps'
                        #from time import time

                        #for i in range(1000):
                        if scalefactor == 1.:
                            continue
                        scale, errmag, chisq, dms, good, image_stamp, psf_stamp, skysig, fitrad, skyysn, psfmag, msk = \
                            chkpsf.fit(imfile.split('.fits')[0], xpos=xsn, ypos=ysn, radius=params.substamp/2.-1.,
                                       returnstamps=True, maskfile=maskfile,mysky=sexsky/scalefactor,mysig=sexrms/scalefactor)


                        #print 'psfmag', psfmag
                        print 'modelshape',image_stamp.shape

                        #save_fits_image(psf_stamp,'test/psf.fits')
                        #raw_input('saved psf stamp')
                        psf_stamp = psf_stamp / 10 ** (-0.4 * (psfmag - 25))
                        noise_stamp = copy(image_stamp)*0+1
                        noise_stamp = noise_stamp*msk
                        if not good:
                            badflag = 1
                            #raw_input('badflagggggooood')
                        #save_fits_image(msk,'test/mask.fits')
                        #print scalefactor
                        #raw_input('scallllllllll')
                        image_stamp *= scalefactor
                        skysig *= scalefactor
                        skyysn  *= scalefactor

                        mask = msk
                        #raw_input('saved mask')

                    else:
                        #try:
                        if True:
                            print 'about to pkfit'
                            #sys.exit()
                            psf_stamp = self.psf
                            print xsn,ysn,im.shape,skysn,skyerrsn
                            if self.usefake:
                                scale = 10**(.4*(31.-float(snparams.fake_truemag[j])))
                            else:
                                scale = float(snparams.flux[i])

                            #errmag =
                            # errmag,chi,niter,scale,iylo,iyhi,ixlo,ixhi,image_stamp,noise_stamp,mask_stamp,pkpsf_stamp = \
                            #     pk.pkfit_norecent_noise_smp(1,xsn,ysn,skysn,skyerrsn,params.fitrad,returnStamps=True,stampsize=params.substamp)

                            # image_stamp = im[np.floor(ysn) - (params.substamp - 1) / 2:np.floor(ysn) + (
                            #                                                                            params.substamp - 1) / 2 + 1,
                            #               np.floor(xsn) - (params.substamp - 1) / 2:np.floor(xsn) + (
                            #                                                                         params.substamp - 1) / 2 + 1]

                            image_stamp = im[int(self.psfcenter[1] - params.substamp/2.):int(self.psfcenter[1] + params.substamp/2.),
                                          int(self.psfcenter[0] - params.substamp/2.):int(self.psfcenter[0] + params.substamp/2.)]
                            noise_stamp = (np.sqrt(weights[int(self.psfcenter[1] - params.substamp/2.):int(self.psfcenter[1] + params.substamp/2.),
                                          int(self.psfcenter[0] - params.substamp/2.):int(self.psfcenter[0] + params.substamp/2.)])*scalefactor)**2

                            #print weightsfile
                            #print self.psfcenter
                            #raw_input()
                            mask_stamp = mask[int(self.psfcenter[1] - params.substamp/2.):int(self.psfcenter[1] + params.substamp/2.),
                                          int(self.psfcenter[0] - params.substamp/2.):int(self.psfcenter[0] + params.substamp/2.)]
                            #save_fits_image(mask_stamp,weightsfile.split('.fits')[0]+'_stamp.fits')
                            #raw_input()
                            mask = mask_stamp
                            #errmag, chi, niter, scale, iylo, iyhi, ixlo, ixhi, image_stamp, noise_stamp, mask_stamp, pkpsf_stamp = 1,1,1,0,100,200,100,200,image_stamp,noise_stamp,image_stamp*0+1,image_stamp*0 + 1


                            # self.dobackgroundstamp = True
                            # if self.dobackgroundstamp:
                            #     bkgrnd = pf.getdata(imfile+'.background')
                            #     bkg_stamp = bkgrnd[self.psfcenter[1] - params.substamp/2.:self.psfcenter[1] + params.substamp/2.,
                            #               self.psfcenter[0] - params.substamp/2.:self.psfcenter[0] + params.substamp/2.]*scalefactor
                            #     # bkgrnd = pf.getdata(imfile + '.background_rms')*scalefactor
                            #     # noise_stamp = 1/bkgrnd[self.psfcenter[1] - params.substamp/2.:self.psfcenter[1] + params.substamp/2.,
                            #     #           self.psfcenter[0] - params.substamp/2.:self.psfcenter[0] + params.substamp/2.]**2.
                            #
                            #     bkg_stamp = bkg_stamp * 0. + skysn

                            msk = copy(image_stamp)
                            msk[msk != 0.] = 1

                        # except ValueError:
                        #     print params.fra,params.fdec,xsn,ysn,im.shape
                        #     print 'SN too close to edge of CCD!'
                        #     badflag = 1
                            #raise ValueError('SN too close to edge of CCD!')




                    # model = np.append(np.zeros(len(image_stamp.ravel())),scale)
                    # newsub = int(image_stamp.shape[0])
                    # stdev = np.zeros(len(model))
                    # stdev[-1] = np.sqrt(model[-1])

                #else:
                #    print 'badflaggggg'

                #if snparams.psf_model.lower() == 'psfex':
                fwhm = float(snparams.psf[j])
                if snparams.psf_unit.lower() == 'arcsec':
                    fwhm_arcsec = fwhm
                elif snparams.psf_unit.lower().startswith('sigma-pix') or snparams.psf_unit.lower().startswith('pix'):
                    fwhm_arcsec = fwhm*2.235*snparams.platescale

                else:
                    raise exceptions.RuntimeError('Error : FWHM units not recognized!!')

                
                if not badflag:
                    if not np.isfinite(skysig):
                        print 'infinite skysig'
                        #raw_input('skysig badflag')
                        badflag = 1
                        descriptiveflag = 8
                    if sexrms < 0.:
                        print 'skysig less than one'
                        #raw_input('skysig1 badflag')
                        badflag = 1
                        descriptiveflag = 32

                badflags.append(badflag)
                if badflag:
                    print 'badflaggg'*100


                if sexsky/scalefactor - 10000. > skymax:
                    badflag = 1
                    descriptiveflag = 8
                    print ('sky '+str(sexsky/scalefactor- 10000.)+' greater than skymax '+str(skymax)+'\n')*10
                elif sexsky/scalefactor - 10000. < skymin:
                    badflag = 1
                    descriptiveflag = 16
                    print ('sky '+str(sexsky/scalefactor- 10000.)+' less than skymin '+str(skymin)+'\n')*10
                elif sexrms/scalefactor > skyerrmax:
                    badflag = 1
                    descriptiveflag = 32
                    print ('skyerr '+str(sexrms/scalefactor)+' greater than skyerrmax '+str(skyerrmax)+'\n')*10

                if not badflag:
                    if fwhm_arcsec < params.fwhm_max:
                        if np.min(im[int(ysn-2):int(ysn+3),int(xsn-2):int(xsn+3)]) != np.max(im[int(ysn-2):int(ysn+3),int(xsn-2):int(xsn+3)]):
                            #if len(np.where(mask[ysn-25:ysn+26,xsn-25:xsn+26] != 0)[0]) < params.max_masknum
                                if np.max(psf_stamp[int(params.substamp/2+1-3):int(params.substamp/2+1+4),int(params.substamp/2+1-3):int(params.substamp/2+1+4)]) == np.max(psf_stamp[:,:]):
                                    #i = j


                                    noise_stamp[image_stamp > 500000.] = 0.

                                    print 'here'*100

                                    #skysn *= scalefactor
                                    #skyerrsn *= scalefactor

                                    #for iiii in noise_stamp.ravel():
                                    #    print iiii
                                    #print 1/(skyerrsn*scalefactor)**2



                                    if self.snparams.survey == 'PS1':
                                        noise_stamp[noise_stamp > 0.] = 1
                                        noise_stamp[noise_stamp == 0.] = 0
                                        #print skyerrsn,skysig
                                        #print np.max(noise_stamp.ravel())
                                        #raw_input()
                                        print self.snparams.skysig[j] * scalefactor,skyerrsn
                                        #raw_input()
                                        smp_noise[i, :, :] = noise_stamp * 1. / (self.snparams.skysig[j] * scalefactor) ** 2 * mask
                                        #smp_noise[i,:,:] = noise_stamp*1./(skyerrsn)**2 * mask
                                        dt.save_fits_image(image_stamp,'im.fits',go=True)

                                        smp_dict['sky'][i] = skyysn
                                        smp_dict['skyerr'][i] = self.snparams.skysig[j] * scalefactor

                                        print skysn,skyerrsn
                                        #print np.max(smp_noise[i,:,:].ravel())
                                        #raw_input()

                                    else:
                                        noise_stamp[noise_stamp > 0.] = 1
                                        noise_stamp[noise_stamp <= 0.] = 0
                                        mask *= noise_stamp
                                        print smp_noise.shape
                                        print noise_stamp.shape
                                        smp_noise[i,:,:] = noise_stamp*0.+1./(mysexskysig**2) * mask


                                        if inflateskyerrworked:
                                            infl1 = self.inflerr(sexsky/scalefactor -10000, inflateskyerr.sky, inflateskyerr.inflatebysky, filt)
                                            infl2 = self.inflerr(mysexskysig/scalefactor, inflateskyerr.skyerr, inflateskyerr.inflatebyskyerr, filt)
                                            if infl1 == 1:
                                                if infl2 == 1:
                                                    print 'here1'
                                                    infltot = 1.
                                                else:
                                                    print 'here2'
                                                    infltot = np.sqrt((infl1-1.)**2 + ( infl2 - 1. )**2) + 1.
                                            else:
                                                print 'here3'
                                                infltot = np.sqrt((infl1 - 1.) ** 2 + (infl2 - 1.) ** 2) + 1.

                                            print 'inflate skyerr!'*100
                                            print 'sky is ',sexsky/scalefactor -10000,'skyerr is ',mysexskysig/scalefactor,'inflating by ',infltot
                                            mysexskysig *= infltot

                                        # smp_dict['sky'][i] = mysky
                                        smp_dict['sky'][i] = sexsky
                                        # smp_dict['skyerr'][i] = skysig
                                        smp_dict['skyerr'][i] = mysexskysig

                                            #smp_noise[i,:,:] = np.sqrt(weights[int(self.psfcenter[1] - params.substamp/2.):int(self.psfcenter[1] + params.substamp/2.),
                                        #  int(self.psfcenter[0] - params.substamp/2.):int(self.psfcenter[0] + params.substamp/2.)]*scalefactor)**2
                                        #smp_noise[i,:,:] = noise_stamp*1./(skyerrsn)**2 * mask
                                        #smp_noise[i,:,:] = noise_stamp*1./(sexrms)**2 * mask

                                    #if round(float(snparams.mjd[j])) == 57011:
                                    #    raw_input()

                                        #if self.dobackgroundstamp:
                                    #    smp_bkg[i,:,:] = bkg_stamp

                                        #print noise_stamp
                                        #raw_input()
                                        #smp_noise[i, :, :] = noise_stamp  # using noise stamp for v5

                                    #print 'image-stamp',image_stamp.shape
                                    #print 'smp_im',smp_im[i,:,:].shape,i

                                    # if not self.oldformat:
                                    #     smp_im[i, :, :] = image_stamp + sexrms**2
                                    # else:
                                    smp_im[i,:,:] = image_stamp
                                    #if float(snparams.fake_truemag[j]) < 24:
                                    if '/pnfs/des/persistent/smp/v62/20131009_SN-S2/r_12/' in longimfile:
                                        print xsn,ysn
                                        print sexsky, sexsky/(10**(.4*(31.-zpt)))
                                        print longimfile
                                        print snparams.fake_truemag[j], zpt
                                        save_fits_image(image_stamp,str(self.snparams.mjd[j])+'.fits')
                                        raw_input()
                                    #save_fits_image(psf_stamp,'test/cpsf.fits')
                                    #raw_input('savedpsf')
                                    if not self.snparams.survey == 'PS1':
                                        smp_psf[i,:,:] = psf_stamp/np.sum(psf_stamp)
                                    else:
                                        smp_psf[i, :, :] = psf_stamp
                                    #save_fits_image(psf_stamp,'testpsf.fits')
                                    c = 20
                                    psa = self.snparams.platescale

                                    smp_starra[i,:len(thisra)] = thisra
                                    smp_stardec[i,:len(thisdec)] = thisdec
                                    smp_starind[i,:len(thisids)] = thisids

                                    smp_dict['scale'][i] = scale
                                    #smp_dict['scale_err'][i] = errmag
                                    # smp_dict['sky'][i] = skysn
                                    # #raw_input('sssssss')
                                    # dosextractor=False
                                    # if dosextractor:
                                    #     smp_dict['sky'][i] = sexsky
                                    #     smp_dict['skyerr'][i] = sexrms
                                    # else:
                                    #     print 'skysn',skysn
                                    #     print 'mysky',mysky
                                    #     #raw_input()


                                    #print sexsky,sexrms,np.mean(image_stamp.ravel()),np.std(image_stamp.ravel())
                                    #raw_input('sss')
                                    smp_dict['flag'][i] = 0
                                    #print smp_dict['flag'][i]

                                    smp_dict['descriptiveflag'][i] = descriptiveflag
                                    #CHECK FOR DIFFIM FLAGS
                                    #if self.snparams.photflag[j] == True:

                                    #    print 'photometry flag!'
                                    #    smp_dict['flag'][i] = 1
                                    #    smp_dict['scale'][i] = np.nan
                                    #    smp_dict['scale_err'][i] = np.nan
                                    #    #raw_input('photometry flag')
                                    #raw_input('flags above')



                                    fileroot = imfile.split('.fits')[0]
                                    smp_dict['fileroots'][i] = fileroot
                                    smp_dict['expnum'][i] = self.expnum
                                    #if self.snparams.survey == 'PS1':
                                    #    smp_dict['fullims'].append('%s.fits' % fileroot)
                                    #    smp_dict['impsfs'].append('%s.dao.psf.fits' % fileroot)
                                    #    tmpp, hp = rdpsf.rdpsf('%s.dao.psf.fits' % fileroot)
                                    #    smp_dict['hpsfs'][i] =

                                    smp_dict['zpt'][i] = zpt
                                    smp_dict['zpterr'][i] = zpterr
                                    #print smp_dict['zpterr']
                                    #raw_input()
                                    smp_dict['mjd'][i] = float(snparams.mjd[j])
                                    #smp_dict['mjdoff'].append( mjdoff )
                                    #smp_dict['mjdslopeinteroff'].append(mjdslopeinteroff)
                                    smp_dict['image_scalefactor'][i] = scalefactor
                                    smp_dict['snx'][i] = xsn
                                    smp_dict['sny'][i] = ysn
                                    smp_dict['scalefactor'][i] = scalefactor
                                    smp_dict['snra'][i] = snparams.RA
                                    smp_dict['sndec'][i] = snparams.DECL
                                    smp_dict['skysig'][i] = mysexskysig
                                    smp_dict['sexrms'][i] = sexrms
                                    smp_dict['sexsky'][i] = sexsky
                                    smp_dict['imwcs'].append(w)
                                    msk = copy(image_stamp)
                                    msk[msk!=0.] = 1
                                    smp_mask[i,:,:] = mask
                                    #smp_dict['mask'][i] = msk
                                    smp_dict['fwhm_arcsec'][i] = fwhm_arcsec
                                    #print fwhm_arcsec
                                    #raw_input('fwhm_arcsec')
                                    smp_dict['image_filename'][i] = longimfile
                                    smp_dict['zpt_file'][i] = os.path.join('/'.join(imfile.split('/')[:-1]), zpt_file)
                                    smp_dict['psf_filename'][i] = longpsffile
                                    if self.snparams.survey == 'PS1':
                                        self.psfcenter = [round(xsn),round(ysn)]
                                    smp_dict['psfcenter'].append(self.psfcenter)
                                    #print xsn,ysn,psffile,self.psfcenter
                                    #raw_input('testing new dict params')
                                    #smp_dict['psf_fwhm'][i] = psf_fwhm
                                    smp_dict['fakepsf'][i] = snparams.psf[j]
                                    if self.useweights:
                                        smp_dict['weight_filename'][i] = longnoisefile
                                    else:
                                        smp_dict['weight_filename'][i] = longnoisefile[0]

                                    if self.usefake:
                                        smp_dict['fakemag'][i] = snparams.fake_truemag[j]
                                    smp_dict['fakezpt'][i] = snparams.zp[j]
                                    smp_dict['diffim_flux'][i] = snparams.flux[j]
                                    smp_dict['diffim_fluxerr'][i] = snparams.fluxerr[j]
                                    smp_dict['id_obs'][i] = snparams.id_obs[j]
                                    smp_dict['id_coadd'][i] = snparams.id_coadd[j]

                                    if self.snparams.survey == 'DES':
                                        smp_dict['psf_centerx'][i] = int(self.psfcenter[1])
                                        smp_dict['psf_centery'][i] = int(self.psfcenter[0])

                                    smp_dict['gain'][i] = self.gain
                                    if self.snparams.survey == 'PS1':
                                        smp_dict['gain'][i] = 1.

                                    smp_dict['rmsaddin'][i] = rmsaddin
                                    fs = snparams.flux
                                    brightlimit = fs[np.argsort(fs)][::-1][:15]
                                    brightlimit = brightlimit[-1]
                                    if snparams.flux[i] > brightlimit:
                                        smp_dict['notbrightflag'][i] = 0
                                    else:
                                        smp_dict['notbrightflag'][i] = 1

                                    print 'zpts',snparams.zp[i], zpt
                                    #raw_input()

                                    #START HERE TOMORROW

                                    # if not nozpt:
                                    #     fname = self.checkstarfile.split('.')[0]+'_deltaradec.npz'
                                    #     self.deltastarsfile = fname
                                    #     df = np.load(self.deltastarsfile)
                                    #     self.usedeltaras = np.array(df['deltaras'])
                                    #     self.usedeltadecs = np.array(df['deltadecs'])
                                    #     self.usemjds = np.array(df['mjds'])
                                    #     self.useras = np.array(df['ras'])
                                    #     self.usedecs = np.array(df['decs'])
                                    #     self.usexstar = np.array(df['x_star'])
                                    #     self.useystar = np.array(df['y_star'])
                                    # else:
                                    #     self.usedeltaras = np.array(copy(self.deltaras))
                                    #     self.usedeltadecs = np.array(copy(self.deltadecs))
                                    #     self.usemjds = np.array(copy(self.deltamjds))
                                    #     self.useras = np.array(copy(self.ras))
                                    #     self.usedecs = np.array(copy(self.decs))
                                    #     self.usexstar = np.array(copy(self.x_stars))
                                    #     self.useystar = np.array(copy(self.y_stars))
                                    #
                                    # print 'nearbystarcalc'
                                    # srad = params.nearby_stars_pixel_rad
                                    # rad = ((self.usexstar-xsn)**2+(self.useystar-ysn)**2)**.5
                                    #
                                    # nearbystars_onthisCCDepoch_indices = [(rad < srad) & (self.usemjds == float(snparams.mjd[j]))]
                                    # nearby_xstar = self.usexstar[nearbystars_onthisCCDepoch_indices]
                                    # nearby_ystar = self.useystar[nearbystars_onthisCCDepoch_indices]
                                    # print xsn
                                    # print nearby_xstar
                                    # print self.usedeltaras[nearbystars_onthisCCDepoch_indices]
                                    # print np.mean(self.usedeltaras[nearbystars_onthisCCDepoch_indices])
                                    # print 'hshshshshshshshs'
                                    #raw_input()

                                    #NOW FIND NEARbY STARS AND GRAB OFFSETS AND APPLY TO SN (MAJE A NEW SMP_DICT ELEMENT)
                                        

                                    if filt == 'u':
                                        #smp_dict['hostgal_mag'][i] = snparams.fake_hostmag_u
                                        smp_dict['hostgal_sbmag'][i] = -99
                                    if filt == 'g':
                                        #smp_dict['hostgal_mag'][i] = snparams.fake_hostmag_g
                                        smp_dict['hostgal_sbmag'][i] = 27.5 - 2.5*np.log10(float(snparams.hostgal_sb_fluxcal[0]))
                                    if filt == 'r':
                                        #smp_dict['hostgal_mag'][i] = snparams.fake_hostmag_r
                                        smp_dict['hostgal_sbmag'][i] = 27.5 - 2.5*np.log10(float(snparams.hostgal_sb_fluxcal[1]))
                                    if filt == 'i':
                                        #smp_dict['hostgal_mag'][i] = snparams.fake_hostmag_i
                                        smp_dict['hostgal_sbmag'][i] = 27.5 - 2.5*np.log10(float(snparams.hostgal_sb_fluxcal[2]))
                                    if filt == 'z':
                                        #smp_dict['hostgal_mag'][i] = snparams.fake_hostmag_z
                                        smp_dict['hostgal_sbmag'][i] = 27.5 - 2.5*np.log10(float(snparams.hostgal_sb_fluxcal[3]))
                                    if filt== 'y':
                                        #smp_dict['hostgal_mag'][i] = snparams.fake_hostmag_y
                                        smp_dict['hostgal_sbmag'][i] = -99

                                    #raw_input()

                                    if smp_dict['mjd'][i] < snparams.peakmjd - params.mjdminus or \
                                        smp_dict['mjd'][i] > snparams.peakmjd + params.mjdplus:
                                        smp_dict['mjd_flag'][i] = 1


                                    #print 'fakemag',snparams.fake_truemag[j],'mjdflag',smp_dict['mjd_flag'][i],'scale',scale
                                    #raw_input()

                                    #print   snparams.peakmjd - params.mjdminus ,smp_dict['mjd'][i],snparams.peakmjd + params.mjdplus
                                    #print smp_dict['flag'][i],'flagggggg'
                                    #raw_input()

                                    if self.dogalfit:
                                        if smp_dict['mjd'][i] > snparams.peakmjd + params.mjdplus:
                                            smp_dict['fitflag'][i] = 0
                                        if smp_dict['mjd'][i] < snparams.peakmjd - params.mjdminus:
                                            smp_dict['fitflag'][i] = 0

                                    #sys.exit()
                                    # if self.fermigrid and self.worker:
                                    #     #print 'cleaning up copied files'
                                    #     os.popen('rm '+imfile).read()
                                    #     os.popen('rm '+imfile+'.fz').read()
                                    #     os.popen('rm '+noisefile).read()
                                    #     os.popen('rm '+psffile).read()


                                    #i += 1
                                    print 'epochtime',time.time()-epochtime, float(snparams.mjd[j]),imfile
                                else:
                                    print 'min max psf equal'*100
                                    smp_dict['imwcs'].append(np.nan)
                                    smp_dict['psfcenter'].append(np.nan)
                                    smp_dict['sky'][i] = sexsky  # smp_dict['skyerr'][i] = skysig
                                    smp_dict['skyerr'][i] = mysexskysig
                                    smp_dict['mjd'][i] = float(snparams.mjd[j])
                                    smp_dict['fwhm_arcsec'][i] = fwhm_arcsec
                                    smp_dict['zpt'][i] = zpt
                                    smp_dict['zpterr'][i] = zpterr
                                    smp_dict['expnum'][i] = self.expnum
                                    smp_dict['image_filename'][i] = imfile
                                    smp_dict['zpt_file'][i] = 'na'
                                    smp_dict['psf_filename'][i] = psffile
                                    smp_dict['weight_filename'][i] = weightsfile
                                    smp_dict['id_obs'][i] = snparams.id_obs[j]
                                    smp_dict['id_coadd'][i] = snparams.id_coadd[j]

                        else:
                            print 'min max image equal'*100
                            if descriptiveflag == 0: descriptiveflag=128
                            smp_dict['descriptiveflag'][i] = descriptiveflag
                            smp_dict['imwcs'].append(np.nan)
                            smp_dict['psfcenter'].append(np.nan)
                            smp_dict['sky'][i] = sexsky  # smp_dict['skyerr'][i] = skysig
                            smp_dict['skyerr'][i] = mysexskysig
                            smp_dict['mjd'][i] = float(snparams.mjd[j])
                            smp_dict['fwhm_arcsec'][i] = fwhm_arcsec
                            smp_dict['zpt'][i] = zpt
                            smp_dict['zpterr'][i] = zpterr
                            smp_dict['expnum'][i] = self.expnum
                            smp_dict['image_filename'][i] = imfile
                            smp_dict['zpt_file'][i] = 'na'
                            smp_dict['psf_filename'][i] = psffile
                            smp_dict['weight_filename'][i] = weightsfile

                            smp_dict['id_obs'][i] = snparams.id_obs[j]
                            smp_dict['id_coadd'][i] = snparams.id_coadd[j]

                    else:
                        print 'failed fwhm'*100
                        print 'fwhm was',fwhm_arcsec,' and fwhm cut was',params.fwhm_max
                        if descriptiveflag == 0: descriptiveflag = 256
                        smp_dict['descriptiveflag'][i] = descriptiveflag

                        smp_dict['imwcs'].append(np.nan)
                        smp_dict['psfcenter'].append(np.nan)
                        smp_dict['sky'][i] = sexsky  # smp_dict['skyerr'][i] = skysig
                        smp_dict['skyerr'][i] = mysexskysig
                        smp_dict['mjd'][i] = float(snparams.mjd[j])
                        smp_dict['fwhm_arcsec'][i] = fwhm_arcsec
                        smp_dict['zpt'][i] = zpt
                        smp_dict['zpterr'][i] = zpterr
                        smp_dict['expnum'][i] = self.expnum
                        smp_dict['image_filename'][i] = imfile
                        smp_dict['zpt_file'][i] = 'na'
                        smp_dict['psf_filename'][i] = psffile
                        smp_dict['weight_filename'][i] = weightsfile

                        smp_dict['id_obs'][i] = snparams.id_obs[j]
                        smp_dict['id_coadd'][i] = snparams.id_coadd[j]
                        #raw_input()

                else:
                    print 'badflag'*100
                    print 'descriptive flag',descriptiveflag
                    smp_dict['descriptiveflag'][i] = descriptiveflag
                    smp_dict['imwcs'].append(np.nan)
                    smp_dict['psfcenter'].append(np.nan)
                    smp_dict['sky'][i] = sexsky  # smp_dict['skyerr'][i] = skysig
                    smp_dict['skyerr'][i] = mysexskysig
                    smp_dict['mjd'][i] = float(snparams.mjd[j])
                    smp_dict['fwhm_arcsec'][i] = fwhm_arcsec
                    smp_dict['zpt'][i] = zpt
                    smp_dict['zpterr'][i] = zpterr
                    smp_dict['expnum'][i] = self.expnum
                    smp_dict['image_filename'][i] = imfile
                    smp_dict['zpt_file'][i] = 'na'
                    smp_dict['psf_filename'][i] = psffile
                    smp_dict['weight_filename'][i] = weightsfile

                    smp_dict['id_obs'][i] = snparams.id_obs[j]
                    smp_dict['id_coadd'][i] = snparams.id_coadd[j]
            else:
                badflag = 1
                if descriptiveflag == 0: descriptiveflag = 64
                print 'xsn > 25 and ysn > 25 and xsn < snparams.nxpix-25 and ysn < snparams.nypix-25\n'*20
                smp_dict['descriptiveflag'][i] = descriptiveflag
                smp_dict['imwcs'].append(np.nan)
                smp_dict['psfcenter'].append(np.nan)
                smp_dict['sky'][i] = sexsky  # smp_dict['skyerr'][i] = skysig
                smp_dict['skyerr'][i] = mysexskysig
                smp_dict['mjd'][i] = float(snparams.mjd[j])
                smp_dict['fwhm_arcsec'][i] = fwhm_arcsec
                smp_dict['zpt'][i] = zpt
                smp_dict['zpterr'][i] = zpterr
                smp_dict['expnum'][i] = self.expnum
                smp_dict['image_filename'][i] = imfile
                smp_dict['zpt_file'][i] = 'na'
                smp_dict['psf_filename'][i] = psffile
                smp_dict['weight_filename'][i] = weightsfile

                smp_dict['id_obs'][i] = snparams.id_obs[j]
                smp_dict['id_coadd'][i] = snparams.id_coadd[j]
        #sys.exit()


        #now loop over images again and get nightly offsets...

        print np.where([smp_dict['mjd_flag'] == 1])
        print len(smp_dict['mjd_flag'][np.where(smp_dict['mjd_flag'] == 1)])
        #raw_input()
        if len(smp_dict['mjd_flag'][np.where(smp_dict['mjd_flag'] == 1)]) < 3:
            raise ValueError(
                "Not enough epochs without SN flux ( > 1 )")

        meanstarras = {}
        meanstardecs = {}
        for ind in np.unique(smp_starind):
            ww = np.where(smp_starind == ind)
            meanstarras[ind] = np.mean(smp_starra[ww].ravel())
            meanstardecs[ind] = np.mean(smp_stardec[ww].ravel())
        print np.unique(smp_starind)

        for im,fl,k in zip(smp_dict['image_filename'],smp_dict['flag'],range(len(smp_dict['image_filename']))):

            if fl != 1:
                nightlyoffra = []
                nightlyoffdec = []
                nightlydist = []
                for ind in np.unique(smp_starind):
                    if ind != 0:
                        #print ind
                        ww = smp_starind[k,:] == ind
                        if len(ww) > 0:
                            try:
                                float(meanstarras[ind] - smp_starra[k, ww])
                                nightlyoffra.append(float(meanstarras[ind] - smp_starra[k,ww]))
                                nightlyoffdec.append(float(meanstardecs[ind] - smp_stardec[k, ww]))
                                nightlydist.append(((smp_dict['snra'][k] - smp_starra[k,ww])**2+(smp_dict['sndec'][k] - smp_stardec[k,ww])**2)[0]**.5)
                                #print ind,'worked'
                            except:
                                #print ind, 'didnt work'
                                pass

                nightlyoffra = np.array(nightlyoffra)
                nightlyoffdec = np.array(nightlyoffdec)
                nightlydist = np.array(nightlydist)
                #for nd in nightlydist:
                #    print nd
                #raw_input('nightly distnaces')

                #FINDING THE CLOSEST 25 STARS TO CALCULATE OFFSETS
                goodindices = []

                for e,j in enumerate(np.argsort(nightlydist)):
                    #print int(j)
                    #if e  < 40:
                    goodindices.append(int(j))

                goodindices = np.array(goodindices,dtype='int')

                try:
                    goodindices = goodindices[:30]
                    smp_dict['raoff'][k] = np.mean(nightlyoffra[goodindices])
                    smp_dict['decoff'][k] = np.mean(nightlyoffdec[goodindices])
                except:
                    print 'WARNING: Not enough nearyby stars to compute nightly offset... \nskipping', im
                    smp_dict['flag'][k] = 1

                # if len(goodindices) < 10:
                #     print 'WARNING: Not enough nearyby stars to compute nightly offset... \nskipping',im
                #     smp_dict['flag'][k] = 1
                # else:
                #     smp_dict['raoff'][k] = np.mean(nightlyoffra[goodindices])
                #     smp_dict['decoff'][k] = np.mean(nightlyoffdec[goodindices])
                #     #print 'radec off', np.mean(nightlyoffra[goodindices]),np.mean(nightlyoffdec[goodindices])




                    w = wcs.WCS(im)



                    xsn, ysn = zip(*w.wcs_world2pix(np.array([[snparams.RA , snparams.DECL]]), 0))
                    xsn = xsn[0]
                    ysn = ysn[0]

                    xsno, ysno = zip(*w.wcs_world2pix(np.array([[snparams.RA + smp_dict['raoff'][k], snparams.DECL + smp_dict['decoff'][k]]]), 0))
                    xsno = xsno[0]
                    ysno = ysno[0]

                    smp_dict['yoff'][k] = xsno-xsn #+ smp_dict['snx'][k] - round(smp_dict['snx'][k])
                    smp_dict['xoff'][k] = ysno-ysn #+ smp_dict['sny'][k] - round(smp_dict['sny'][k])

                    if smp_dict['yoff'][k] > .5:
                        print 'NIGHTLY OFFSET IS TOO LARGE'
                        smp_dict['flag'][k] = 1
                    if smp_dict['xoff'][k] > .5:
                        print 'NIGHTLY OFFSET IS TOO LARGE'
                        smp_dict['flag'][k] = 1

                    # x_psf,y_psf = cntrd.cntrd(smp_psf[k,:,:],15,15,2.)
                    # print x_psf,y_psf,-smp_dict['snx'][k]+round(smp_dict['snx'][k]),\
                    #     -smp_dict['sny'][k]+round(smp_dict['sny'][k])
                    # #
                    # print 'pix off',smp_dict['xoff'][k],smp_dict['yoff'][k]
                    # raw_input()


        #print 'dillscale',smp_dict['scale']
        if self.fermilog:
            self.tmpwriter.appendfile('\nDone with zeropoints\n ', self.fermilogfile)
        if self.fermigrid and self.worker:
            try:
                os.environ('PROCESS')
                print os.popen('rm -r data')
            except:
                print 'this is not really a worker'
            #if self.fermilog:
            #    self.tmpwriter.appendfile(os.popen('ls -ltr data/').read(), self.fermilogfile)

        #print os.popen('ls -ltr').read()
        #sys.exit()
        if mergeno == 0:
            zeroArray = np.zeros(smp_noise.shape)
            largeArray = zeroArray + 1E10

            smp_psfWeight = smp_psf
            smp_psf = np.fmax(smp_psf,zeroArray)
            smp_im = np.fmax(smp_im,zeroArray)
        mergectr = 0
        '''
        while mergectr < mergeno:
            print "Matrix Merger {0}".format(mergectr + 1)
            rem = -1.0 * (smp_noise.shape[1] % 2)
            if np.abs(rem) != 0:
                zeroArray = np.zeros(smp_noise[:,:rem:2,:rem:2].shape)
                largeArray = zeroArray + 1E10
                smp_noise = (np.fmin(smp_noise[:,:rem:2,:rem:2],largeArray) + np.fmin(smp_noise[:,1:rem:2,1:rem:2],largeArray) + np.fmin(smp_noise[:,:rem:2,1:rem:2],largeArray) + np.fmin(smp_noise[:,1:rem:2,:rem:2],largeArray))/4.0
                smp_psfWeight = (np.fmin(smp_psf[:,:rem:2,:rem:2],largeArray) + np.fmin(smp_psf[:,1:rem:2,1:rem:2],largeArray) + np.fmin(smp_psf[:,:rem:2,1:rem:2],largeArray) + np.fmin(smp_psf[:,1:rem:2,:rem:2],largeArray))/4.0
                smp_psf = (np.fmax(smp_psf[:,:rem:2,:rem:2],zeroArray) + np.fmax(smp_psf[:,1:rem:2,1:rem:2],zeroArray) + np.fmax(smp_psf[:,1:rem:2,:rem:2],zeroArray) + np.fmax(smp_psf[:,:rem:2,1:rem:2],zeroArray))/4.0
                smp_im = (np.fmax(smp_im[:,:rem:2,:rem:2],zeroArray) + np.fmax(smp_im[:,1:rem:2,1:rem:2],zeroArray) + np.fmax(smp_im[:,1:rem:2,:rem:2],zeroArray) + np.fmax(smp_im[:,:rem:2,1:rem:2],zeroArray))/4.0
                params.substamp+=rem
                params.substamp/=2.0
                mergectr+=1
            else:
                zeroArray = np.zeros(smp_noise[:,::2,::2].shape)
                largeArray = zeroArray + 1E10
                smp_noise = (np.fmin(smp_noise[:,::2,::2],largeArray) + np.fmin(smp_noise[:,1::2,1::2],largeArray) + np.fmin(smp_noise[:,1::2,::2],largeArray) + np.fmin(smp_noise[:,::2,1::2],largeArray))/4.0
                smp_psfWeight = (np.fmin(smp_psf[:,::2,::2],largeArray) + np.fmin(smp_psf[:,1::2,1::2],largeArray) + np.fmin(smp_psf[:,1::2,::2],largeArray) + np.fmin(smp_psf[:,::2,1::2],largeArray))/4.0
                smp_psf = (np.fmax(smp_psf[:,::2,::2],zeroArray) + np.fmax(smp_psf[:,1::2,1::2],zeroArray) + np.fmax(smp_psf[:,1::2,::2],zeroArray) + np.fmax(smp_psf[:,::2,1::2],zeroArray))/4.0
                smp_im = (np.fmax(smp_im[:,::2,::2],zeroArray) + np.fmax(smp_im[:,1::2,1::2],zeroArray) + np.fmax(smp_im[:,1::2,::2],zeroArray) + np.fmax(smp_im[:,::2,1::2],zeroArray))/4.0
                params.substamp/=2.0
                mergectr+=1
        '''

        badnoisecols = np.where(smp_noise < 1e-9)
        if self.snparams.survey == 'PS1':
            smp_noise[badnoisecols] = 0.
        badpsfcols = np.where(smp_psf < 0)
        if self.snparams.survey == 'PS1':
            smp_noise[badpsfcols] = 0.0
        smp_psf[badpsfcols] = 0.0


        infinitecols = np.where((smp_im == 0) | (np.isfinite(smp_im) == 0) | (np.isfinite(smp_noise) == 0))
        if self.snparams.survey == 'PS1':
            smp_noise[infinitecols] = 0.0
        smp_im[infinitecols] = 0
        mpparams = np.concatenate((np.zeros(int(float(params.substamp)**2.)),smp_dict['scale'],smp_dict['sky']))

        mpdict = [{'value':'','step':0,
                  'relstep':0,'fixed':0, 'xtol': 1E-15} for i in range(len(mpparams))]

        mpparams[:params.substamp**2] = np.fmax((np.nanmax(smp_im, axis=0)/np.nanmax(smp_psfWeight, axis =0)),np.zeros(smp_im[0].shape)).flatten()


        for i in range(len(mpparams)):
            thisparam = mpparams[i]
            if thisparam == thisparam and thisparam < 1E305 and i >= params.substamp**2:
                mpdict[i]['value'] = thisparam
                if i >= (params.substamp**2 + len(smp_dict['mjd'])):
                            mpdict[i]['fixed'] = 1
            else:
                mpdict[i]['value'] = 0.0
                mpdict[i]['fixed'] = 1

        for col in range(int(params.substamp)**2+len(smp_dict['scale'])):
            mpdict[col]['step']=np.sqrt(np.max(smp_dict['scale']))

        for col in np.where((smp_dict['flag'] == 1))[0]+int(params.substamp)**2:
            print 'flagged '+str(col)
            mpdict[col]['fixed'] = 1
            mpdict[col]['value'] = 0


        if verbose: print('Creating Initial Scene Model')

        

        arg = -1
        mn = 999999999
        for fwhma in smp_dict['fwhm_arcsec']:
            arg += 1
            if smp_dict['flag'][arg] == 0:
                if smp_dict['mjd_flag'][arg] == 1:
                    if fwhma < mn:
                        if fwhma > 0.:
                            if np.any(smp_im[arg,:,:]):
                                mn = fwhma
                                usearg = arg


        for i in np.arange(len(smp_dict['sky'])):
            if np.max(smp_psf[i,:,:]) == np.min(smp_psf[i,:,:]):
                #save_fits_image(smp_psf[i,:,:],'test/culpritpsf.fits')
                #print 'hererererere psffsfsffsfsf',smp_dict['mjd'][i]
                print 'MJD:',smp_dict['mjd'][i],'BAD PSF... SETTING FLAG TO 1'
                #raw_input()
                smp_dict['flag'][i] = 1
                #smp_dict['mjd'][i]

        try:
            ww = usearg
        except:
            ww = 0


        model, modelpeak = self.create_model(smp_im[ww,:,:]-smp_dict['sky'][ww],smp_psf[ww,:,:],smp_dict['scale'])


        stdev = np.sqrt(copy(model))
        newsub = int(smp_im[ww].shape[0])

        modelvec = np.zeros(len(smp_dict['scale']))
        modelstd = np.zeros(len(smp_dict['scale']))

        for i,scale in enumerate(smp_dict['scale']):
            if i in np.where((smp_dict['mjd_flag'] == 1) | (smp_dict['flag'] == 1))[0]:
                model[newsub**2+i] = 0
                modelvec[i] = 0.
                print i,'modelveczero'
            else:
                #if scale < 25.:
                #    model[newsub**2+i] = 0.
                #    modelvec[i] = 0.
                #else:
                modelvec[i] = scale

        #SET GALAXY MODEL HERE!!!! TO IMAGE - SKY
        galmodel = smp_im[ww,:,:]-smp_dict['sky'][ww]

        pkyerr = -2.5*np.log10(smp_dict['mcmc_scale']) + 2.5*np.log10(smp_dict['mcmc_scale'] + smp_dict['mcmc_scale_err'])
        print 'oldoutdir',oldoutdir

        outfolder = oldoutdir
        out = os.path.join(oldoutdir,'SNe/'+snparams.snfile.split('/')[-1].split('.')[0] + '/'+filt+'/')
        outimages = os.path.join(oldoutdir,'SNe/'+snparams.snfile.split('/')[-1].split('.')[0] + '/'+filt+'/','image_stamps')
        print 'outimages',outimages

        if self.fermigrid and self.worker:
            os.popen('ifdh mkdir '+outfolder+'/SNe')
            os.popen('ifdh mkdir '+outfolder+'/SNe/'+snparams.snfile.split('/')[-1].split('.')[0])
            os.popen('ifdh mkdir '+out)
            os.popen('ifdh mkdir '+outimages)
            #raw_input('made directories')
        else:
            if not os.path.exists(out):
                os.makedirs(out)
            if not os.path.exists(outimages):
                os.makedirs(outimages)

        
        #pix_model, newsubstamp = self.pixelate(model,self.params.pixelation_factor,substamp=params.substamp)
        '''pix_data, newsubstamp = self.pixelate(smp_im,self.params.pixelation_factor,ndim=len(smp_im))
        pix_noise, newsubstamp = self.pixelate(smp_noise,self.params.pixelation_factor,ndim=len(smp_noise))
        pix_psfs, newsubstamp = self.pixelate(smp_psf,self.params.pixelation_factor,ndim=len(smp_psf))
        pix_mask, newsubstamp = self.pixelate(smp_dict['mask'][0],self.params.pixelation_factor)

        params.substamp = newsubstamp #setting updated pixelated substamp
        '''

        self.tmpwriter.savefits(model[:params.substamp**2].reshape(params.substamp,params.substamp),os.path.join(outimages,'initialmodel.fits'))

        imodel = copy(model)

        fakemag = smp_dict['fakemag']
        zpt = smp_dict['zpt']
        fakezpt = smp_dict['fakezpt']


        
        fake_flux = 10.**(.4*(fakezpt-fakemag))
        nfake_flux = 10.**(.4*(27.5-fakemag))*10.**(.4*(31.-27.5))
        fake_flux_adjusted = fake_flux*10**(.4*(fakezpt-zpt))

        diffim_mag = 27.5-2.5*np.log10(smp_dict['diffim_flux'])
        diffim_mag_err = -2.5*np.log10(smp_dict['diffim_flux'])+2.5*np.log10(smp_dict['diffim_flux']+smp_dict['diffim_fluxerr'])

        diffim_flux = 10.**(.4*(31.-diffim_mag))
        diffim_fluxerr = 10.**(.4*(31.-(diffim_mag-diffim_mag_err)))-10.**(.4*(31.-(diffim_mag+diffim_mag_err)))

        scaled_diffim_flux = smp_dict['diffim_flux']*10.**(.4*(31.-27.5))
        scaled_diffim_fluxerr = smp_dict['diffim_fluxerr']*10.**(.4*(31.-27.5))

        if self.snparams.survey == 'PS1':
            scaled_diffim_flux = smp_dict['diffim_flux'] * 10. ** (.4 * (31. - smp_dict['fakezpt']))
            scaled_diffim_fluxerr = smp_dict['diffim_fluxerr'] * 10. ** (.4 * (31. - smp_dict['fakezpt']))

        self.scaled_diffim_flux = scaled_diffim_flux
        self.scaled_diffim_fluxerr = scaled_diffim_fluxerr

        filename = snparams.snfile.split('/')[-1].split('.')[0] +'_'+ filt

        npoutdir = os.path.join(oldoutdir,'np_data/'+filt+'/')

        galaxyoutdir = os.path.join(galaxyfoldername,'np_data/'+filt+'/')

        if self.fermigrid and self.worker:
            os.system('ifdh mkdir '+os.path.join(galaxyfoldername,'np_data'))
            os.system('ifdh mkdir '+os.path.join(oldoutdir,'np_data'))
            os.system('ifdh mkdir '+npoutdir)
        else:
            if not os.path.exists(npoutdir):
                try:
                    os.makedirs(outdir)
                except:
                    pass

        maxiter = 1200
        print os.path.join(outdir,filename+'_mcmc_input.npz')

        filename = snparams.snfile.split('/')[-1].split('.')[0] +'_'+ filt
        #lightcurves = os.path.join(oldoutdir,'/lightcurves/'+filt+'/')
        #if not os.path.exists(lightcurves):
        #    os.makedirs(lightcurves)

        for i in range(len(smp_dict['mjd'])):
            print smp_dict['mjd'][i]

        print 'mjds'*100


        print smp_dict['image_filename'][-1]        
        print 'MJD','\t','BAND','\t','FIT_ZPT','\t','FAKE_ZPT','\t','PSF','\t','SKY','\t','SexSky','\t','Skyerr','\t','Skysig','\t','SexRMS','\t','IMAGE_FILENAME','\t','GAIN','\t','XOFF','\t','YOFF'
        psfs = []
        for i,scale in enumerate(smp_dict['scale']):
            if i in np.where((smp_dict['flag'] == 1))[0]:
                psfs.append(-9.)
                print smp_dict['mjd'][i],'\t','FLAGGED'
            else:
                fitzpt = smp_dict['zpt'][i]
                fakezpt = smp_dict['fakezpt'][i]
                psfs.append(smp_dict['fwhm_arcsec'][i])
                try:
                    print smp_dict['mjd'][i],'\t',filt,round(fitzpt,2),'\t','\t',round(fakezpt,2),'\t',\
                    round(smp_dict['fwhm_arcsec'],2),round(smp_dict['sky'][i],2),round(smp_dict['sexsky'][i],2),\
                    round(smp_dict['skyerr'][i],2),round(smp_dict['skysig'][i],2),round(smp_dict['sexrms'][i],2),\
                    smp_dict['image_filename'][i],smp_dict['gain'][i],smp_dict['xoff'][i],smp_dict['yoff'][i]
                except:
                    print '999'
                #if abs(fitzpt - fakezpt) > .025:
                #    smp_dict['fitflag'][i] = 1
                #if smp_dict['mjd'][i] < snparams.peakmjd +100:
                #    smp_dict['flag'][i] = 1
                if smp_dict['skyerr'][i] < .001:
                    smp_dict['flag'][i] = 1
                if smp_dict['skyerr'][i] > 1000.:
                    smp_dict['flag'][i] = 1

        print 'MJD', '\t', 'BAND', '\t', 'CCD', '\t', 'FIT_ZPT', '\t', 'EXPNUM', '\t', 'SKY', '\t', 'Skyerr', '\t', 'GAIN', ''
        psfs = []
        for i, scale in enumerate(smp_dict['scale']):
            if i in np.where((smp_dict['flag'] == 1))[0]:
                print smp_dict['mjd'][i],'\t','FLAGGED'
            else:
                fitzpt = smp_dict['zpt'][i]
                fakezpt = smp_dict['fakezpt'][i]
                try:
                    print smp_dict['mjd'][i], filt, self.ccdnum ,round(fitzpt, 2), smp_dict['expnum'][i],\
                    round(smp_dict['sky'][i]/(10**(.4*(31-fitzpt))), 2), \
                    round(smp_dict['skyerr'][i]/(10**(.4*(31-fitzpt))),2),  smp_dict['gain'][i]
                except:
                    print '999'
        #raw_input()
        zptnpz = os.path.join(npoutdir,filename+'_imagezpts.npz')
        self.tmpwriter.savez(zptnpz
                , mjd = smp_dict['mjd']
                , band = filt
                , fitzpt = smp_dict['zpt']
                , fakezpt = smp_dict['fakezpt']
                , psf = psfs
                , imfile = smp_dict['image_filename']
                , flags = smp_dict['flag']
                )


        mjdvec = []
        fluxvec = []
        fluxerrvec = []
        fakefluxvec = []

        tstart = time.time()

        # nm = self.checkstarfile.split('.')[0].split('/')[-1] + '_deltaradec.npz'
        # fname = os.path.join(outdir, foldername, 'np_data', filt, nm)
        # self.deltastarsfile = fname
        
        # if nozpt:
        #     self.tmpwriter.savez(self.deltastarsfile,deltaras=self.deltaras,deltadecs=self.deltadecs,mjds=self.deltamjds,ras=self.ras,decs=self.decs,airmasses=self.airmasses,x_star=self.x_stars,y_star=self.y_stars)
        #     print self.deltastarsfile,'SAVED'
        # else:
        #     dsf  = np.load(self.deltastarsfile)
        #     self.deltaras = dsf['deltaras']
        #     self.deltadecs = dsf['deltadecs']
        #     self.deltamjds = dsf['mjds']
        #     self.ras = dsf['ras']
        #     self.decs = dsf['decs']
        #     self.airmasses = dsf['airmasses']
        #fname = self.checkstarfile.split('.')[0]+'_20magfakes.npz'
        #self.moneyfile = fname
        #if nozpt:
        #    self.tmpwriter.savez(self.moneyfile,flux=self.fakestarfluxes,fluxerr=self.fakestarfluxerrs,zpt=self.fakestarzpts)
        #    print fname,'SAVED'

        #self.plotcheckstars()
        #self.plotallfake20staroffsets()


        for i in range(len(smp_dict['sky'])):
            print smp_dict['mjd'][i],len(smp_noise[i][smp_noise[i]*smp_mask[i] > 0.].ravel()),
            if len(smp_noise[i][smp_noise[i]*smp_mask[i] > 0.].ravel()) < 500.:
                smp_dict['flag'][i] = 1
                print '<--FLAGGED',

        print 'skyerr',smp_dict['skyerr']
        print 'flag',smp_dict['flag']
        print 'mjdflag',smp_dict['mjd_flag']
        print smp_dict['mjd']
        #print os.path.join(outdir,filename+'_mcmc_input.npz')

        print 'mjdslopeinteroff',smp_dict['mjdslopeinteroff']

        print os.path.join(outdir,filename+'_mcmc_input.npz')
        print 'idobs',smp_dict['id_obs']
        print 'idcoadd',smp_dict['id_coadd']

        try:
            print snparams.RA
        except:
            print sys.exc_info()
            if not os.path.exists(isdonedir):
                os.mkdir(isdonedir)
            os.system('touch '+os.path.join(isdonedir,snparams.snfile.split('/')[-1].split('.')[0] + '.done'))
            sys.exit()

        self.tmpwriter.savez( os.path.join(npoutdir,filename+'_mcmc_input.npz'),
                galmodel = galmodel
                , modelvec = modelvec*0.
                , galstd = np.sqrt(galmodel*0. +1.)*2.
                , modelstd = np.sqrt(modelvec*0.)
                , data = smp_im
                , psfs = smp_psf
                , weights = smp_noise
                , substamp = params.substamp
                , Nimage = len(smp_dict['sky'])
                , maxiter = 100000
                , mask = None
                , sky=smp_dict['sky']
                , mjd=smp_dict['mjd']
                , gewekenum=9999999
                , skyerr=smp_dict['skyerr']
                , useskyerr = True
                , flags = smp_dict['flag']
                , fitflags = smp_dict['fitflag']
                , psf_shift_std = .0005
                , shiftpsf = False
                , hostgal_mag = smp_dict['hostgal_mag']
                , hostgal_sbmag = smp_dict['hostgal_sbmag']
                , fileappend = ''
                , stop = False
                , skyerr_radius = 16.
                , outpath = outimages
                , compressionfactor = 10
                , fix_gal_model = None
                , pixelate_model = 1.
                , total_skyerr = smp_dict['total_skyerr']
                , skysig = smp_dict['skysig']
                , mjdflag = smp_dict['mjd_flag'],
                fakefluxadj=fake_flux_adjusted,
                fakeflux=nfake_flux,
                fakemag=smp_dict['fakemag'],
                fakezpt=smp_dict['fakezpt'],
                zpterr = smp_dict['zpterr'],
                peakmjd =snparams.peakmjd,
                snra = snparams.RA,
                sndec = snparams.DECL,
                xsn = xsn,
                ysn = ysn,
                zpt = smp_dict['zpt'],
                diffim_flux=scaled_diffim_flux,
                diffim_fluxerr=scaled_diffim_fluxerr,
                fwhm = smp_dict['fwhm_arcsec']*2.355,
                snname = filename,
                npzloc = outdir,
                lcout = os.path.join(self.lcfilepath+'/'+filename),
                model_data_index = ww,
                mjdoff = smp_dict['mjdoff'],
                mjdslopeinteroff = smp_dict['mjdslopeinteroff'],
                star_offset_file = star_offset_file,
                zpt_files = smp_dict['zpt_file'],
                starglobalids = starglobalids,
                globalraoffsets = offsetra,
                globaldecoffsets = offsetdec,
                id_obs = smp_dict['id_obs'],
                id_coadd = smp_dict['id_coadd']
                )
        
        try:
            self.tmpwriter.savez(os.path.join(npoutdir,filename+'_smpDict.npz'),**smp_dict)
        except:
            print 'could not save ',os.path.join(npoutdir,filename+'_smpDict.npz')


        self.outdir = outdir
        self.galaxyoutdir = galaxyoutdir
        self.filename = filename
        #self.foldername = foldername
        self.filt = filt
        self.smp_dict = smp_dict
        self.smp_im = smp_im
        self.smp_psf = smp_psf
        self.smp_noise = smp_noise




        if self.dosnfit:
            for passv in [0,1]:
                print 'zpterr', smp_dict['zpterr']

                if not self.dogalfit: continue

                minsky = np.argmin(smp_dict['sky'])
                galmodel_params = smp_im[minsky]-smp_dict['sky'][minsky]
                modelvec = scaled_diffim_flux

                xoff = 0.
                yoff = 0.

                galmodel = galmodel_params
                galmodel = smp_im[ww, :, :] - smp_dict['sky'][ww]
                import scipy.ndimage

                print 'galmodel',galmodel

                galstd = np.sqrt(abs(galmodel))/params.galmodel_div
                print 'galstd',galstd

                modelstd = np.sqrt(abs(smp_dict['scale']))/params.flux_std_div
                print 'mjd_flag',smp_dict['mjd_flag']
                print 'modelstd before',modelstd
                print smp_dict['scale']
                modelstd[(modelstd < 10.) & (modelstd > 0.)] = 10.
                modelstd[smp_dict['mjd_flag'] == 1] = 0.
                modelvec[smp_dict['mjd_flag'] == 1] = 0.
                modelvec[smp_dict['flag'] == 1] = 0.


                print 'modelstd after',modelstd
                tstart = time.time()
                st = time.time()

                passflags = copy(smp_dict['flag'])
                galstd = galstd * 0 + .9
                numiter = self.params.sn_plus_galmodel_steps
                shiftstd = self.params.sn_shift_std
                burnin = .4
                if passv == 0:
                    dontfitflags = np.zeros(len(smp_dict['sky']))
                    dontfitflags[smp_dict['mjd_flag'] == 0 ] = 1
                    dontfitflags[smp_dict['flag'] == 1] = 1
                    if len(dontfitflags[dontfitflags == 0])>10.:
                        aw = np.argwhere(dontfitflags == 0)
                        #print 'aw', aw
                        dontfitflags[aw[:-10]] = 1

                    dontfitflags[smp_dict['flag'] == 1] = 0
                    passflags += dontfitflags
                    numiter = 5000
                    shiftstd = 0.
                    galstd = galstd * 0 + 2.
                    burnin = .8
                #print 'smp_dict[flag]',smp_dict['flag']
                print 'passflags',passflags
                #print 'dontfitflags',dontfitflags

                if passv == 1: self.continu = True
                try:
                    if self.continu:
                        galmodel = np.load(self.lcfilepath+'/'+snparams.snfile.split('/')[-1].split('.')[0] + '_' + self.filt + '.npz')['galmodel_params']
                        modelvec = np.load(self.lcfilepath+'/'+snparams.snfile.split('/')[-1].split('.')[0] + '_' + self.filt + '.npz')['modelvec']
                except:
                    print 'No previous fit'

                plt.clf()
                plt.close()
                import gc
                collected = gc.collect()

                stdoutfile = os.path.join(self.lcfilepath,
                                          snparams.snfile.split('/')[-1].split('.')[0] + '_' + self.filt + '.mcmcout')

                import mcmc_geweke as mcmc3
                aaa = mcmc3.metropolis_hastings(
                        galmodel = galmodel
                        , modelvec = modelvec
                        , galstd = galstd
                        , modelstd = modelstd*2.5
                        , data = smp_im
                        , psfs = smp_psf
                        , weights = smp_noise
                        , substamp = params.substamp
                        , Nimage = len(smp_dict['sky'])
                        , maxiter = numiter
                        , mask = smp_mask
                        , sky=smp_dict['sky']
                        , mjd=smp_dict['mjd']
                        , gewekenum=-9999999*(passv-1) + 10000 #cannot be less than 100 after compression
                        , skyerr=smp_dict['skyerr']
                        , useskyerr = False
                        , usesimerr = False
                        , flags = passflags
                        , fitflags = smp_dict['fitflag']*0.
                        , psf_shift_std = shiftstd
                        , xoff = 0.
                        , yoff = 0.#.06
                        , shiftpsf = shiftstd > 0.
                        , fileappend = ''
                        , stop = False
                        , skyerr_radius = 12.
                        , outpath = outimages
                        , compressionfactor = 10
                        , fix_gal_model = False
                        , pixelate_model = None
                        , burnin = burnin
                        , lcout = os.path.join(self.lcfilepath,filename)
                        , chainsnpz = os.path.join(npoutdir,filename+'_withSn.npz')
                        , mjdoff = smp_dict['mjdoff']
                        , dontsavegalaxy=True
                        , log = self.fermilogfile
                        , isfermigrid=self.fermigrid
                        , isworker=self.worker
                        , x=smp_dict['snx']
                        , y=smp_dict['sny']
                        , psffile=smp_dict['psf_filename']
                        , psfcenter=smp_dict['psfcenter']
                        , model_errors=True
                        , survey = self.snparams.survey
                        #, fullims = smp_dict['fullims']
                        #, impsfs = smp_dict['impsfs']
                        #, hpsfs = smp_dict['hpsfs']
                        , fileroots = smp_dict['fileroots']
                        , scalefactor = smp_dict['scalefactor']
                        , gain=smp_dict['gain']
                        , dobkg = False
                        , bkg = smp_bkg
                        , sigmazpt = smp_dict['zpterr']
                        , fakemag = smp_dict['fakemag']
                        , fitzpt = smp_dict['zpt']
                        , fakezpt = smp_dict['fakezpt']
                        , datafilenames = smp_dict['image_filename']
                        , nightlyoffx = smp_dict['xoff']
                        , nightlyoffy = smp_dict['yoff']
                        , sstime = sstime
                        , stdoutfile = stdoutfile

                        )
                modelveco = copy(modelvec)
                #if self.fermilog:
                #    self.tmpwriter.appendfile('DONE... saving snfit\n', self.fermilogfile)
                #sys.exit()
                modelvec, modelvec_uncertainty, galmodel_params, galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, sims, xhistory,yhistory,accepted_history,pix_stamp,chisqhist,redchisqhist,stamps,chisqs,ndof,gewekediag,covar,corr  = aaa.get_params()
                fitprob = dt.fitprobfromchisq(chisqs,ndof)
                print fitprob
                print 'TOTAL SMP SN TIME ',time.time()-tstart
                print os.path.join(npoutdir,filename+'_withSn.npz')
                npzoutfile = os.path.join(self.lcfilepath,
                                             snparams.snfile.split('/')[-1].split('.')[0] + '_' + self.filt + '.npz')
                np.savez(npzoutfile, modelvec=modelvec,
                                 modelvec_uncertainty=modelvec_uncertainty, galmodel_params=galmodel_params,
                                 galmodel_uncertainty=galmodel_uncertainty, modelvec_nphistory=modelvec_nphistory,
                                 galmodel_nphistory=galmodel_nphistory, sims=sims, data=smp_im,
                                 xhistory=xhistory,yhistory=yhistory,stamps=stamps,chisqs=chisqs,weights=smp_noise,
                                 flags=smp_dict['flag'], sky=smp_dict['sky'], mjd=smp_dict['mjd'], skyerr=smp_dict['skyerr'],
                                 covar=covar,corr=corr,
                                 accepted_history=accepted_history, chisqhist=chisqhist,gewekediag=gewekediag)

        #self.dogalsimfit = True
        if self.dogalsimfit:

            # if not self.dogalfit:
            #     chains = np.load(os.path.join(galaxyoutdir,filename+'_nosn.npz'))
            #     galmodel_params = chains['galmodel_params']
            #     galmodel_uncertainty = chains['galmodel_uncertainty']
            #


            galstd = np.sqrt(abs(galmodel)) / params.galmodel_div

            galmodel = smp_im[ww, :, :] - smp_dict['sky'][ww]

            brighttest = False
            if brighttest:
                modelvec[:15] = 0
                modelstd[:15] = 0
                smp_dict['flag'][15:] = 1
                smp_dict['flag'][np.array(smp_dict['mjd'], dtype='int') == 56547] = 1

            # minsky = np.argmin(smp_dict['sky'])
            # galmodel_params = smp_im[minsky] - smp_dict['sky'][minsky]
            # modelvec = scaled_diffim_flux
            #
            # galmodel = galmodel_params
            # galstd = np.sqrt(abs(galmodel))/5.

            tstart = time.time()

            # fixmodels = np.array(fixmodelvec*((smp_dict['mjd_flag']+1)%2))
            # fixmodelvec = []
            # for ff in fixmodels:
            #     fixmodelvec.append(float(ff))
            # fixmodelvec = np.array(fixmodelvec)
            # print fixmodelvec
            # modelstd = np.sqrt(abs(fixmodelvec))/5.


            modelvec = scaled_diffim_flux
            modelstd = abs(scaled_diffim_fluxerr)/params.flux_std_div
            modelvec[smp_dict['mjd_flag'] == 1] = 0
            modelstd[smp_dict['mjd_flag'] == 1] = 0


            if self.continu:
               galmodel = np.load(self.lcfilepath+'/'+snparams.snfile.split('/')[-1].split('.')[0] + '_' + self.filt + '.npz')['galmodel_params']/10.
               modelvec = np.load(self.lcfilepath+'/'+snparams.snfile.split('/')[-1].split('.')[0] + '_' + self.filt + '.npz')['modelvec']

            print modelstd
            print 'galmodelshape', galmodel.shape
            import mcmcgalsimsmall as mcmcgalsim
            #smp_dict['flag'][:-15] = 1
            aaa = mcmcgalsim.metropolis_hastings(
                    galmodel = galmodel
                    , modelvec = modelvec
                    , galstd = galstd*0+8.
                    , modelstd = modelstd
                    , data = smp_im
                    , psfs = smp_psf
                    , weights = smp_noise
                    , substamp = params.substamp
                    , Nimage = len(smp_dict['sky'])
                    , maxiter = self.params.sn_plus_galmodel_steps_galsim
                    , mask = smp_mask
                    , sky=smp_dict['sky']
                    , mjd=smp_dict['mjd']
                    , gewekenum=9999999
                    , skyerr=smp_dict['skyerr']
                    , useskyerr = True
                    , flags = smp_dict['flag']
                    , fitflags = smp_dict['fitflag']*0.
                    , psf_shift_std = self.params.sn_shift_std
                    , shiftpsf = self.params.sn_shift_std > 0.
                    , fileappend = ''
                    , stop = False
                    , skysig = smp_dict['skysig']
                    , skyerr_radius = 16.
                    , outpath = outimages
                    , compressionfactor = 25
                    , fix_gal_model = None
                    , pixelate_model = 1.
                    , imagefiles = smp_dict['image_filename']
                    , psffiles = smp_dict['psf_filename']
                    , weightfiles = smp_dict['weight_filename']
                    , snra = snparams.RA
                    , sndec = snparams.DECL
                    #, burnin = 100
                    , model_pixel_scale = .27
                    , lcout = os.path.join(self.lcfilepath,filename)
                    , chainsnpz = os.path.join(npoutdir, filename + '_galsim.npz')
                    , isfermigrid = self.fermigrid
                    , psfcenterx = smp_dict['psf_centerx']
                    , psfcentery = smp_dict['psf_centery']
                    , gain = smp_dict['gain']
                    )

            modelveco = copy(modelvec)
            
            #modelvec, modelvec_uncertainty, galmodel_params, galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, sims, xhistory,yhistory,accepted_history,pix_stamp,chisqhist  = aaa.get_params()
            modelvec, modelvec_uncertainty, galmodel_params, galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, sims, xhistory, yhistory, accepted_history, pix_stamp, chisqhist, redchisqhist, stamps, chisqs = aaa.get_params()

            print 'TOTAL SMP SN TIME ',time.time()-tstart

            xoff = np.mean(xhistory[int(3 * len(xhistory) / 4.):])
            yoff = np.mean(yhistory[int(3 * len(yhistory) / 4.):])


            self.tmpwriter.savez(os.path.join(outdir,filename+'_Galsim.npz'),modelvec=modelvec, modelvec_uncertainty=modelvec_uncertainty, galmodel_params=galmodel_params, galmodel_uncertainty=galmodel_uncertainty, modelvec_nphistory=modelvec_nphistory, galmodel_nphistory=galmodel_nphistory, sims=sims,data=smp_im,accepted_history=accepted_history,chisqhist=chisqhist)
            print os.path.join(outdir,filename+'_Galsim.npz')
        
        if self.dogalsimpixfit:    
            if not self.dogalfit:
                chains = np.load(os.path.join(galaxyoutdir,filename+'_nosn.npz'))
                galmodel_params = chains['galmodel_params']
                galmodel_uncertainty = chains['galmodel_uncertainty']
            if not pixstart == None:
                usedir = os.path.join(outdir,pixstart+'/np_data/'+filt+'/')
                chains = np.load(os.path.join(usedir,filename+'_nosn.npz'))
                modelvec = chains['modelvec']
                modelstd = scaled_diffim_fluxerr/5.
                galmodel = chains['galmodel_params']
            else:
                modelvec = scaled_diffim_flux
                modelstd = scaled_diffim_fluxerr/5.
                galmodel = galmodel_params
            galstd = np.sqrt(abs(galmodel))/10.

            tstart = time.time()

            # fixmodels = np.array(fixmodelvec*((smp_dict['mjd_flag']+1)%2))
            # fixmodelvec = []
            # for ff in fixmodels:
            #     fixmodelvec.append(float(ff))
            # fixmodelvec = np.array(fixmodelvec)
            # print fixmodelvec
            # modelstd = np.sqrt(abs(fixmodelvec))/5.


            modelvec = scaled_diffim_flux
            modelstd = abs(scaled_diffim_fluxerr)/5.
            if not self.floatallepochs:
                modelvec[smp_dict['mjd_flag'] == 1] = 0
                modelstd[smp_dict['mjd_flag'] == 1] = 0

            print 'modelvec',modelvec
            print 'modelstd',modelstd
            print 'galmodelshape', galmodel.shape
            fixgalzero = False
            if self.fixgalzero:
                galmodel = galmodel*0.+.1
                galstd = galstd*0.
                fixgal = True
            else:
                fixgal = False

            extraflag = smp_dict['fitflag'] * 0.

            print 'image shapes',smp_im.shape
            aaa = mcmc3galsimpixshift.metropolis_hastings(
                    galmodel = galmodel
                    , modelvec = modelvec
                    , galstd = galstd
                    , modelstd = modelstd
                    , data = smp_im
                    , psfs = smp_psf
                    , weights = smp_noise
                    , substamp = params.substamp
                    , Nimage = len(smp_dict['sky'])
                    , maxiter = 500000
                    , mask = None
                    , sky=smp_dict['sky']
                    , mjd=smp_dict['mjd']
                    , gewekenum=9999999
                    , skyerr=smp_dict['skyerr']
                    , useskyerr = True
                    , flags = smp_dict['flag']
                    , fitflags = extraflag
                    , psf_shift_std = .00008
                    , shiftpsf = True
                    , fileappend = ''
                    , stop = False
                    , skysig = smp_dict['skysig']
                    , skyerr_radius = 16.
                    , outpath = outimages
                    , compressionfactor = 100
                    , fix_gal_model = fixgal
                    , pixelate_model = 1.
                    , imagefiles = smp_dict['image_filename']
                    , psffiles = smp_dict['psf_filename']
                    , weightfiles = smp_dict['weight_filename']
                    , snra = snparams.RA
                    , sndec = snparams.DECL
                    , model_pixel_scale = .27
                    , lcout = lightcurves+filename
                    , chainsnpz = os.path.join(outdir,filename+'_withSnAndGalsimPix.npz')
                    , platescale = .27
                    , snraoff = 0.
                    , sndecoff = 0.
                    )

            modelveco = copy(modelvec)
            
            modelvec, modelvec_uncertainty, galmodel_params, galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, sims, xhistory,yhistory,accepted_history,pix_stamp,chisqhist,rahistory,dechistory  = aaa.get_params(dosave=True)

            print 'TOTAL SMP SN TIME ',time.time()-tstart

            print os.path.join(outdir,filename+'_withSnAndGalsimPix.npz')

        self.outdir = outdir
        self.filename = filename
        self.filt = filt
        self.smp_dict = smp_dict
        self.smp_im = smp_im
        self.smp_psf = smp_psf
        self.smp_noise = smp_noise

        image_stampf,sim_stampf,galmodel_stampf,weight_stampf,psf_stampf,chisq_stampf = stamps[0],stamps[1],stamps[2],stamps[3],stamps[4],stamps[5]
        #print 'lcfilepath',self.lcfilepath

        #print snparams.snfile
        #print snparams.snfile.split('/')[-1]
        #print snparams.snfile.split('/')[-1].split('.')[0]
        #print 'Chi Squares for Each Epoch'
        #print chisqs
        smplightcurvefile = os.path.join(self.lcfilepath,
                                         snparams.snfile.split('/')[-1].split('.')[0] + '_' + self.filt + '.smp')

        smpcovarfile = os.path.join(self.lcfilepath,
                                         snparams.snfile.split('/')[-1].split('.')[0] + '_' + self.filt + '_covar.npz')


        np.savez(smpcovarfile,covar=covar,corrcoef=corr)


        #os.popen('ifdh mkdir '+self.lcfilepath).read()
        if self.fermigrid and self.worker:
            tmp = 'lc.txt'
        else:
            tmp = smplightcurvefile
        #print len(gewekediag),len(smp_dict['snx'])
        #print 'lencheck'*100
        fout = open(tmp, 'w')
        print >> fout, '# MJD DPMJD ID_OBS ID_COADD BAND ZPT ZPTERR FLUX FLUXERR FAKEMAG FAKEZPT DIFFIM_FLUX DIFFIM_FLUXERR ' \
                       'XPOS YPOS XOFF YOFF RA DEC CHI2 NDOF ' \
                       'SMP_FLAG MJD_FLAG SKY SKYERR RMSADDIN GEWKEDIAG ' \
                       'IMAGE_FILE PSF_FILE WEIGHT_FILE ZPTFILE FITGALMODEL_STAMP ' \
                       'IMAGE_STAMP PSF_STAMP WEIGHT_STAMP SIM_STAMP CHISQ_STAMP'
        for i in range(len(smp_dict['snx'])):
            print >> fout, '%.5f %.5f %i %i %s %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %s %i %i %.5f ' \
                           '%.5f %.5f %.5f %s %s %s %s %s %s %s %s %s %s' % (
                                smp_dict['mjd'][i],float(smp_dict['mjd'][i])-snparams.peakmjd, smp_dict['id_obs'][i],smp_dict['id_coadd'][i], self.filt,
                                smp_dict['zpt'][i], smp_dict['zpterr'][i],
                                modelvec[i], modelvec_uncertainty[i], smp_dict['fakemag'][i], smp_dict['fakezpt'][i],
                                smp_dict['diffim_flux'][i],smp_dict['diffim_fluxerr'][i],
                                smp_dict['snx'][i], smp_dict['sny'][i],xoff,yoff,
                                smp_dict['snra'][i], smp_dict['sndec'][i],
                                chisqs[i], -999, smp_dict['flag'][i],smp_dict['mjd_flag'][i],smp_dict['descriptiveflag'][i],
                                smp_dict['sky'][i], smp_dict['skyerr'][i], smp_dict['rmsaddin'][i],gewekediag[i],
                                smp_dict['image_filename'][i], smp_dict['psf_filename'][i],
                                smp_dict['weight_filename'][i], smp_dict['zpt_file'][i],
                                galmodel_stampf[i],
                                image_stampf[i],psf_stampf[i],weight_stampf[i],sim_stampf[i],chisq_stampf[i])
        fout.close()
        if self.fermigrid and self.worker:
            ls = os.popen('ifdh ls ' + smplightcurvefile).read()
            if len(ls) > 0:
                os.popen('ifdh rm ' + smplightcurvefile)
            os.popen('ifdh cp --force=xrootd lc.txt '+smplightcurvefile)
        if self.fermilog:
            self.tmpwriter.appendfile('SMP Successful\n', self.fermilogfile)
        print('SMP was successful!!!')
        print('See stamps/mcmc_chains in',self.outdir)
        print('See lightcurve file',smplightcurvefile)
        if self.snparams.survey == 'PS1':
            if not os.path.exists(self.isdonedir):
                os.mkdirs(self.isdonedir)
            os.system('touch '+os.path.join(self.isdonedir,snparams.snfile.split('/')[-1].split('.')[0] + '.done'))
        #sys.exit()
        return

    def inflerr(self, value, x, y, band):
        index = 0
        if band == 'g':
            index = 0
        elif band == 'r':
            index = 1
        elif band == 'i':
            index = 2
        elif band == 'z':
            index = 3
        scaling = 1.

        scaling = np.interp(value,x[index],y[index])

        return scaling

    def closest_node(self,ra,dec):
        tra = self.bigcatalogras*0. + ra
        tdec = self.bigcatalogdecs*0. + dec

        dist_2 = ((self.bigcatalogras-tra)**2+(self.bigcatalogdecs-tdec)**2)**.5
        return np.argmin(dist_2)


    def getProperCatRaDec(self,ra,dec):
        properra = np.zeros(len(ra))
        properdec = np.zeros(len(dec))
        propermag = np.zeros(len(ra))
        #self.bigcatalog
        for i in np.arange(0,len(ra)):
            j = self.closest_node(ra[i],dec[i])
            properra[i] = self.bigcatalogras[j]
            properdec[i] = self.bigcatalogdecs[j]
            propermag[i] = self.bigcatalogmags[j]
        return properra,properdec,propermag
        


    def plotcheckstars(self):
        outpath = os.path.join(self.outfile,foldername+'/standards/'+filt+'/')
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        
        print self.checkstarfile
        skipindivstars = True
        if not skipindivstars:
            ExposureNum,mjd,ra,dec,xstar,ystar,catzpt,mpfitzpt,mpfitzpterr,fitflux,fitfluxerr,catmag = np.loadtxt(self.checkstarfile,delimiter='\t',skiprows=1,unpack=True)
        
            plt.clf()
            fig, axs = plt.subplots(2,2,figsize=(13,8))
            axes = axs.ravel()
            tra = np.unique(ra)
            ctr = -1
            for r in tra:
                ctr += 1
                ww = ra == r
                umjd = mjd[ww]
                ura = ra[ww]
                udec = dec[ww]
                uzpt = mpfitzpt[ww]
                uzpterr = mpfitzpterr[ww]
                uflux = fitflux[ww]
                ufluxerr = fitfluxerr[ww]
                ucatmag = catmag[ww]
                uxstar = xstar[ww]
                uystar = ystar[ww]
                #mm = np.argsort(umjd)
                #print uzpt,uflux,ucatmag
                #raw_input()
                try:
                    #axes[ctr].plot(umjd,uzpt-2.5*np.log10(uflux)-ucatmag)
                    axes[ctr].errorbar(umjd,uzpt-2.5*np.log10(uflux)-ucatmag,yerr=-2.5*np.log10(uflux)+2.5*np.log10(uflux+ufluxerr),color='black',alpha=.7,fmt='o',label='RA: '+str(round(ura[0],3))+' DEC: '+str(round(udec[0],3))+'\nXstar: '+str(round(np.mean(uxstar)))+' Ystar: '+str(round(np.mean(uystar))))
                    axes[ctr].legend()
                    axes[ctr].plot(umjd,umjd*0.,color='black')
                    axes[ctr].set_xlabel('MJD')
                    axes[ctr].set_ylabel('Fit Star Mag - Catalog Mag')
                except:
                    continue
            #plt.plot(umjd,ucatmag,color='black')
                
            plt.xlabel('MJD')
            plt.ylabel('Fit Star Mag - CAtalog Mag')
            fname = outpath+self.snfn+'_standards.png'
        #fff = open('standardplotslist.txt','a')
        #fff.write(fname+'\n')
        #fff.close()
            plt.tight_layout()
            self.savefig(fname)
            #print fname
        
        dta = 3600.

        a = np.load(self.deltastarsfile)
        print a.keys()
        deltamjds = a['mjds']
        deltaras = a['deltaras']*dta
        deltadecs = a['deltadecs']*dta
        ras = a['ras']
        decs = a['decs']
        airm = a['airmasses']


        

        dra = np.array(deltaras)
        print float(len(dra[abs(dra) > .0005*dta]))/float(len(dra))
        #print 'heahehehahehe'
        #raw_input()
        plt.clf()
        fig, axs = plt.subplots(2,1,figsize=(13,8))
        axes = axs.ravel()
        
        axes[0].scatter(deltamjds,deltaras,color='black',alpha=.1)
        axes[0].plot([min(deltamjds),max(deltamjds)],[0,0],color='black')
        axes[0].set_ylim(-.0002*dta,.0002*dta)
        axes[0].set_xlabel('MJD')
        axes[0].set_ylabel('Delta RA (arcsec)')
        axes[1].scatter(deltamjds,deltadecs,color='black',alpha=.1)
        axes[1].plot([min(deltamjds),max(deltamjds)],[0,0],color='black')
        axes[1].set_ylim(-.0002*dta,.0002*dta)
        axes[1].set_xlabel('MJD')
        axes[1].set_ylabel('Delta DEC (arcsec)')

        fname = outpath+self.snfn+'_deltaradec.png'
        plt.tight_layout()
        self.savefig(fname)
        #print fname

        print min(deltaras),max(deltaras)

        plt.clf()
        fig, axs = plt.subplots(2,1,figsize=(13,13))
        axes = axs.ravel()
        jjj = abs(deltaras) < 1
        deltaras = deltaras[jjj]
        ras = ras[jjj]
        mm, s, iii = meanclip.meanclip( deltaras, clipsig = 3., maxiter = 8, returnSubs=True)
        print iii
        print len(deltaras[iii]),len(deltaras)
        m = 0.
        b = np.polyfit(ras[iii], deltaras[iii], 0)
        a,c = np.polyfit(ras[iii], deltaras[iii], 1)

        print 'rafit slope inter mean',a,c,b,mm
        axes[0].scatter(ras[iii],deltaras[iii],color='black',alpha=.1)
        axes[0].plot([min(ras),max(ras)],[0,0],color='black')
        axes[0].plot([min(ras),max(ras)],[m*min(ras)+b,m*max(ras)+b],color='green',linewidth=3,alpha=.4)
        axes[0].plot([min(ras),max(ras)],[a*min(ras)+c,a*max(ras)+c],color='red',linewidth=3,alpha=.4)
        axes[0].set_xlabel('RA')
        axes[0].set_ylim(-.0002*dta,.0002*dta)
        axes[0].set_ylabel('Delta RA (arcsec)')


        jjj = abs(deltadecs) < 1.
        deltadecs = deltadecs[jjj]
        decs = decs[jjj]

        m, s, iii = meanclip.meanclip( deltadecs, clipsig = 3, maxiter = 8, returnSubs=True)
        m = 0.
        b = np.polyfit(decs[iii], deltadecs[iii], 0)
        a,c = np.polyfit(decs[iii], deltadecs[iii], 1)

        print 'decfit slope inter mean',a,c,b

        axes[1].scatter(decs[iii],deltadecs[iii],color='black',alpha=.1)
        axes[1].plot([min(decs),max(decs)],[0,0],color='black')
        axes[1].plot([min(decs),max(decs)],[m*min(decs)+b,m*max(decs)+b],color='green',linewidth=3,alpha=.4)
        axes[1].plot([min(decs),max(decs)],[a*min(decs)+c,a*max(decs)+c],color='red',linewidth=3,alpha=.4)
        axes[1].set_ylim(-.0002*dta,.0002*dta)
        axes[1].set_xlabel('DEC')
        axes[1].set_ylabel('Delta DEC (arcsec)')

        fname = outpath+self.snfn+'_deltaradecvsRADEC.png'
        plt.tight_layout()
        self.savefig(fname)
        #print fname
        #raw_input()
        plt.clf()
        fig, axs = plt.subplots(2,1,figsize=(13,13))
        axes = axs.ravel()

        axes[0].hist(deltaras,color='black',bins=np.arange(-.000205*dta,.0002*dta,.00001*dta))
        axes[0].set_ylabel('Count')
        axes[0].set_xlabel('Delta RA (arcsec)')
        axes[1].hist(deltadecs,color='black',bins=np.arange(-.000205*dta,.0002*dta,.00001*dta))
        axes[1].set_ylabel('Count')
        axes[1].set_xlabel('Delta DEC (arcsec)')

        fname = outpath+self.snfn+'_deltaHist.png'
        plt.tight_layout()
        self.savefig(fname)
        #print fname
        
        plt.clf()
        fig, axs = plt.subplots(5,5,figsize=(35,28))
        axes = axs.ravel()
        for m,i in zip(np.unique(deltamjds),np.arange(len(np.unique(deltamjds)))):
            ww = (deltamjds == m)
            try:
                axes[i].scatter(ras[ww],deltaras[ww],color='black',alpha=.6,label='Airmass: '+str(airm[ww][0])+'\nMJD: '+str(deltamjds[ww][0]))
            except:
                continue
            axes[i].plot([min(ras),max(ras)],[0,0],color='black')
            axes[i].set_ylim(-.0002*dta,.0002*dta)
            axes[i].set_xlabel('RA')
            axes[i].set_ylabel('Delta RA (arcsec)')
            axes[i].legend()
        fname = outpath+self.snfn+'_deltara_singlemjd.png'
        plt.tight_layout()
        self.savefig(fname)
        #print fname

        plt.clf()
        fig, axs = plt.subplots(5,5,figsize=(35,28))
        axes = axs.ravel()
        for m,i in zip(np.unique(deltamjds),np.arange(len(np.unique(deltamjds)))):
            ww = (deltamjds == m)
            try:
                axes[i].scatter(decs[ww],deltadecs[ww],color='black',alpha=.6,label='Airmass: '+str(airm[ww][0])+'\nMJD: '+str(deltamjds[ww][0]))
            except:
                continue
            axes[i].plot([min(decs),max(decs)],[0,0],color='black')
            axes[i].set_ylim(-.0002*dta,.0002*dta)
            axes[i].set_xlabel('DEC')
            axes[i].set_ylabel('Delta DEC (arcsec)')
            axes[i].legend()
        fname = outpath+self.snfn+'_deltadec_singlemjd.png'
        plt.tight_layout()
        self.savefig(fname)
        #print fname

        plt.clf()
        fig, axs = plt.subplots(1,2,figsize=(12,12))
        axes = axs.ravel()
        axes[0].scatter(airm,deltaras,color='black',alpha=.1)
        axes[0].set_ylim(-.0002*dta,.0002*dta)
        axes[0].set_xlabel('Airmass')
        axes[0].set_ylabel('Delta RA (arcsec)')
        ax,ay,aystd = bindata(airm,deltaras,np.arange(-.0002*dta,.0002*dta,.00005*dta))
        axes[0].errorbar(ax,ay,aystd,markersize=10,color='lime',fmt='o')



        axes[1].scatter(airm,deltadecs,color='black',alpha=.1)
        #axes[1].plot([min(deltamjds),max(deltamjds)],[0,0],color='black')
        axes[1].set_ylim(-.0002*dta,.0002*dta)
        axes[1].set_xlabel('Airmass')
        axes[1].set_ylabel('Delta DEC (arcsec)')
        ax,ay,aystd = bindata(airm,deltadecs,np.arange(-.0002*dta,.0002*dta,.00005*dta))
        axes[1].errorbar(ax,ay,aystd,markersize=10,color='lime',fmt='o')

        fname = outpath+self.snfn+'_airmass_vs_deltaradec.png'
        plt.tight_layout()
        self.savefig(fname)
        #print fname

        #print mjd,ra,dec
        #print 'heyyyy'
        #raw_input()
    def afterfit(self,snparams,params,donesn=False):
        print os.path.join(self.outdir,self.filename+'_nosn.npz')
        wosn = np.load(os.path.join(self.galaxyoutdir,self.filename+'_nosn.npz'))

        if not donesn:
            galmodel_wsn = wosn['galmodel_params']*0.
            float_gal_scale = wosn['modelvec']*0.
            float_gal_std = wosn['modelvec_uncertainty']*0.
            ggalmodel_wsn = wosn['galmodel_params']*0.
            gfloat_gal_scale = wosn['modelvec']*0.
            gfloat_gal_std = wosn['modelvec_uncertainty']*0.
        else:
            wsn = np.load(os.path.join(self.outdir,self.filename+'_withSn.npz'))
            galmodel_wsn = wsn['galmodel_params']
            float_gal_scale = wsn['modelvec']
            float_gal_std = wsn['modelvec_uncertainty']
            try:
                gwsn = np.load(os.path.join(self.outdir,self.filename+'_withSnAndGalsim.npz'))
                ggalmodel_wsn = wsn['galmodel_params']
                gfloat_gal_scale = wsn['modelvec']
                gfloat_gal_std = wsn['modelvec_uncertainty']
                print gfloat_gal_scale
                print gfloat_gal_std
                print 'sssss'
            except:
                ggalmodel_wsn = galmodel_wsn*0.
                gfloat_gal_scale = float_gal_scale*0.
                gfloat_gal_std = float_gal_std*0.
            #raw_input()

        galmodel_wosn = wosn['galmodel_params']

        #save_fits_image(galmodel_wsn,os.path.join('/global/cscratch1/sd/dbrout/',self.foldername,'SNe/'+snparams.snfile.split('/')[-1].split('.')[0]+'/'+self.filt+'/image_stamps/finalgalmodel_withSn.fits'))
        #save_fits_image(ggalmodel_wsn,os.path.join('/global/cscratch1/sd/dbrout/',self.foldername,'SNe/'+snparams.snfile.split('/')[-1].split('.')[0]+'/'+self.filt+'/image_stamps/finalgalmodel_withSnAndGalsim.fits'))
        #save_fits_image(galmodel_wosn,os.path.join('/global/cscratch1/sd/dbrout/',self.foldername,'SNe/'+snparams.snfile.split('/')[-1].split('.')[0]+'/'+self.filt+'/image_stamps/finalgalmodel_nosn.fits'))

        skyerr_radius = 16.

        fitrad = np.zeros([params.substamp,params.substamp])
        for x in np.arange(params.substamp):
            for y in np.arange(params.substamp):
                if np.sqrt((params.substamp/2. - x)**2 + (params.substamp/2. - y)**2) < skyerr_radius:
                    fitrad[int(x),int(y)] = 1.

        final_mcmc_fixfluxes = []
        final_mcmc_fixstd = []
        final_results_chisq = []
        final_results_dms = []
        final_mjd = []
        final_fakemag = []
        final_diffim_flux = []
        final_diffim_fluxerr = []
        final_mcmc_floatfluxes = []
        final_mcmc_floatstd = []
        gfinal_mcmc_floatfluxes = []
        gfinal_mcmc_floatstd = []
        for mjd,im,sky,psf,weight,flag,fakemag,diffim_flux,diffim_fluxerr,floatscale,floatstd,gfloatscale,gfloatstd in zip(self.smp_dict['mjd'],self.smp_im,self.smp_dict['sky'],self.smp_psf,self.smp_noise,self.smp_dict['flag'],self.smp_dict['fakemag'],self.scaled_diffim_flux,self.scaled_diffim_fluxerr,float_gal_scale,float_gal_std,gfloat_gal_scale,gfloat_gal_std):
            if flag == 1:
                final_mcmc_fixfluxes.append(np.nan)
                final_mcmc_fixstd.append(np.nan)
                final_results_chisq.append(np.nan)
                final_results_dms.append(np.nan)
                final_mjd.append(mjd)
                final_fakemag.append(np.nan)
                final_diffim_flux.append(np.nan)
                final_diffim_fluxerr.append(np.nan)
                final_mcmc_floatfluxes.append(np.nan)
                final_mcmc_floatstd.append(np.nan)
                gfinal_mcmc_floatfluxes.append(np.nan)
                gfinal_mcmc_floatstd.append(np.nan)
            else:
                SNscale,SNscale_std,chisq,dms = self.getfluxsmp(im,psf,sky,weight,fitrad,galmodel_wosn,mjd,guess_scale=None)
                final_mcmc_fixfluxes.append(SNscale)
                final_mcmc_fixstd.append(SNscale_std)
                final_results_chisq.append(chisq)
                final_results_dms.append(dms)
                final_mjd.append(mjd)
                final_fakemag.append(fakemag)
                final_diffim_flux.append(diffim_flux)
                final_diffim_fluxerr.append(diffim_fluxerr)
                final_mcmc_floatfluxes.append(floatscale)
                final_mcmc_floatstd.append(floatstd)
                gfinal_mcmc_floatfluxes.append(gfloatscale)
                gfinal_mcmc_floatstd.append(gfloatstd)

        print 'Final Fixed Gal fluxes'
        print final_mcmc_fixfluxes
        print 'Final Float Gal fluxes'
        print final_mcmc_floatfluxes

        self.tmpwriter.savez(os.path.join(self.outdir,snparams.snfile.split('/')[-1].split('.')[0]+'_'+self.filt+'_finalresults.npz'),
            fixedgal_scale = np.array(final_mcmc_fixfluxes),
            fixedgal_std = np.array(final_mcmc_fixstd),
            floatgal_scale = np.array(final_mcmc_floatfluxes),
            floatgal_std = np.array(final_mcmc_floatstd),
            gfloatgal_scale = np.array(gfinal_mcmc_floatfluxes),
            gfloatgal_std = np.array(gfinal_mcmc_floatstd),
            zpt = self.smp_dict['zpt'],
            fakezpt = self.smp_dict['fakezpt'],
            chisq = np.array(final_results_chisq),
            dms = np.array(final_results_dms),
            mjd = np.array(final_mjd),
            fakemag = np.array(final_fakemag),
            fakeflux = 10.**(.4*(31.-np.array(final_fakemag))),
            diffim_flux = np.array(final_diffim_flux),
            diffim_fluxerr = np.array(final_diffim_fluxerr),
            peakmjd = snparams.peakmjd
            )

        if not donesn:
            self.tmpwriter.savez(os.path.join(self.outdir,snparams.snfile.split('/')[-1].split('.')[0]+'_'+self.filt+'_nosnresults.npz'),
                fixedgal_scale = np.array(final_mcmc_fixfluxes),
                fixedgal_std = np.array(final_mcmc_fixstd)
            )
        print 'SMP3 Completed Successfully'
        try:
            print np.array(final_mcmc_fixfluxes)[:,0]
            rr = np.array(final_mcmc_fixfluxes)[:,0]
        except:
            rr = np.array(final_mcmc_fixfluxes)





        '''
        mcmc_result = mcmc.metropolis_hastings(galmodel=galmodel,
                                               modelvec=modelvec,
                                               galstd = galstd,
                                               modelstd = modelstd,
                                               psfs=smp_psf,
                                               data=smp_im,
                                               weights=smp_noise,
                                               substamp=params.substamp,
                                               flags=smp_dict['flag'],
                                               mjdflag=smp_dict['mjd_flag'],
                                               fakeflux=fake_flux_adjusted,
                                               diffim_flux=diffim_flux,
                                               diffim_fluxerr=diffim_fluxerr,
                                               sky=smp_dict['sky'],
                                               model_errors=False,
                                               mask=smp_dict['mask'][0],
                                               Nimage=len(smp_dict['sky']),
                                               #Nimage=1,
                                               maxiter=maxiter,
                                               mjd=smp_dict['mjd'],
                                               gewekenum=1000000,
                                               skyerr=smp_dict['skyerr'],
                                               useskyerr=True,
                                               skyerr_radius=self.params.skyerr_radius,
                                               psf_shift_std = .01,
                                               shiftpsf=False,
                                               stop=False,
                                               outpath=outimages,
                                               compressionfactor=10,
                                               pixelate_model=2.5,
                                               npzfile = os.path.join(outdir,filename+'.npz'),
                                               snname = filename,
                                               npzloc = outdir,
                                               dorun = False,
                                               lcout=lightcurves+filename
                                               #gal_model=pyfits.getdata(self.gal_model)
                                               )
        print 'hereerreere'
        model, uncertainty, modelhistory, sims, xhistory, yhistory,acceptedrate,substamp = mcmc_result.get_params()
        del mcmc_result
        #model, theta_std, vals, sims,[],[],np.mean(self.sampler.acceptance_fraction),self.pix_stamp
        print 'bak inside'
        save_fits_image(galmodel,os.path.join(outimages,'finalmodel.fits'))


        #plt.plot(history[:,-10])
        #plt.show()

        fakemag = smp_dict['fakemag']
        zpt = smp_dict['zpt']
        fakezpt = smp_dict['fakezpt']

        fake_flux = 10.**(.4*(31.-fakemag))
        fake_flux_adjusted = fake_flux*10**(.4*(fakezpt-zpt))
        '''
        if not self.gal_model is None:
            import scipy.signal as sci
            #stamp_path = '/global/scratch2/sd/dbrout/smp_y1y2_deep4/SNe/des_fake_00235910/r/image_stamps/'
            #stamp_path = os.path.join(outfile,foldername+'/SNe/'+snparams.snfile.split('/')[-1].split('.')[0]+'/'+filt+'/image_stamps/') 
            stamp_path = outimages
            stamp_path2 =  outimages
            iii = -1
            fig, ax = plt.subplots(6, 7,figsize=(40, 25))
            axs = ax.ravel()
            axi = -1
            for mjd in smp_dict['mjd']:
                iii += 1
                if smp_dict['flag'][iii] == 1:
                    continue
                axi += 1
                data = os.path.join(stamp_path,'MDJ'+str(mjd)+'_fluxdata.fits')
                psf = os.path.join(stamp_path2,'MDJ'+str(mjd)+'_psf.fits')
                skyerr = os.path.join(stamp_path2,'MDJ'+str(mjd)+'_skyerr.fits')
                sky = float(open(os.path.join(stamp_path2,'MDJ'+str(mjd)+'_skyval.txt'),'r').read())
                galaxy_model = stamp_path+'finalmodel.fits'

                flux_range = np.arange(-500,1000,1)
                gal_stamp = pyfits.getdata(galaxy_model)
                try:
                    data_stamp = pyfits.getdata(data)
                except:
                    axi -= 1
                    continue
                psf_stamp = pyfits.getdata(psf)
                skyerr_stamp = pyfits.getdata(skyerr)
                gal_conv = sci.convolve2d(gal_stamp,psf_stamp,mode='same')

                chisqvec = copy(flux_range)*0.
                for i,flux in enumerate(flux_range):
                    sim = gal_conv + flux*psf_stamp + sky
                    chisqvec[i] = np.sum(((sim - data_stamp)**2 / skyerr_stamp**2).ravel())
                minflux = flux_range[np.argmin(chisqvec)]
                axs[axi].plot(flux_range,chisqvec,'r',label=str(mjd))
                axs[axi].axhline(min(chisqvec)+2.3, color='b', linestyle='dashed', linewidth=2,label='Chisq + 2.3')
                axs[axi].axvline(fake_flux[iii], color='black', linestyle=':', linewidth=2,label='Fake Flux Zpt Adjusted')
                axs[axi].legend(prop={'size':6})
        
        filename = snparams.snfile.split('/')[-1].split('.')[0] +'_'+ filt

        lightcurves = os.path.join(outdir,foldername+'/lightcurves/'+filt+'/')
        self.savefig(lightcurves+filename+'_chisqlike.pdf')
        
        #print lightcurves+filename+'_chisqlike.pdf'
        '''
        #raw_input()
        

        filename = snparams.snfile.split('/')[-1].split('.')[0] +'_'+ filt
        lightcurves = os.path.join(outfile,foldername+'/lightcurves/'+filt+'/')    
        if not os.path.exists(lightcurves):
            os.makedirs(lightcurves)

        peakmjd = snparams.peakmjd
        fakemag = smp_dict['fakemag']
        zpt = smp_dict['zpt']
        fakezpt = smp_dict['fakezpt']
        zpterr = smp_dict['zpterr']

        #GET DIFFIM MAGNITUDE HERE!
        diffim_mag = 27.5-2.5*np.log10(smp_dict['diffim_flux'])
        diffim_mag_err = -2.5*np.log10(smp_dict['diffim_flux'])+2.5*np.log10(smp_dict['diffim_flux']+smp_dict['diffim_fluxerr'])

        #mcmc_mag = 31.-2.5*np.log10(model[int(newsub)**2:int(newsub)**2+len(smp_psf)])
        plt.figure(1)
        

        #yerr = -2.5*np.log10(model[int(newsub)**2:int(newsub)**2+len(smp_psf)]) + 2.5*np.log10(model[int(newsub)**2:int(newsub)**2+len(smp_psf)] + uncertainty[int(newsub)**2:int(newsub)**2+len(smp_psf)])
        #mcmc_mag_err = yerr

        #half = int(round(len(modelhistory[:,0])/2,0))
        '''
        snepochs = modelhistory[:,substamp**2:]
        mcmc_flux = np.median(snepochs,axis=0)
        mad = np.median(abs(snepochs - np.median(snepochs,axis=0)),axis=0)
        madstd = 1.48*mad
        mcmc_fluxerr = madstd

        mcmc_mag = 31.-2.5*np.log10(mcmc_flux)
        mcmc_mag_err = abs(-2.5*np.log10(mcmc_flux+mcmc_fluxerr)+2.5*np.log10(mcmc_flux))

        diffim_flux = 10.**(.4*(31.-diffim_mag))
        diffim_fluxerr = 10.**(.4*(31.-(diffim_mag-diffim_mag_err)))-10.**(.4*(31.-(diffim_mag+diffim_mag_err)))
        fake_flux = 10.**(.4*(31.-fakemag))
        fake_flux_adjusted = fake_flux*10**(.4*(fakezpt-zpt))
        
        f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False)
        ax1.errorbar(smp_dict['mjd'],mcmc_mag,yerr=mcmc_mag_err,fmt='o',color='black',label='Scene Model Fit Mag',alpha=.7)
        ax1.scatter(smp_dict['mjd'],fakemag,color='blue',label='Fake True Mag',alpha=.7)
        ax1.errorbar(smp_dict['mjd'],diffim_mag,yerr=diffim_mag_err,fmt='o',color='green',label='Diffim Mag',alpha=.7)
        ax1.set_ylabel('Mag')
        ax1.set_ylim(min(fakemag[fakemag > 10.]) - 2.,max(fakemag[fakemag < 50.])+2.)
        ax1.set_xlim(peakmjd-30.,peakmjd+100)
        ax1.legend(fontsize=8)
        ax1.invert_yaxis()
        ax1.set_title(filename.split('.')[0])
        ax2.errorbar(smp_dict['mjd'],mcmc_mag-fakemag,yerr=mcmc_mag_err,fmt='o',color='black',alpha=.7)
        ax2.errorbar(smp_dict['mjd'],diffim_mag-fakemag,yerr=diffim_mag_err,fmt='o',color='green',label='Diffim Mag',alpha=.7)
        ax2.plot([peakmjd-30.,peakmjd+100],[0,0],color='black')
        ax2.set_xlim(peakmjd-30.,peakmjd+100)
        ax2.set_ylim(-.5,.5)
        ax2.set_xlabel('MJD')
        ax2.set_ylabel('Fit Mag - Fake Mag')
        # Fine-tune figure; make subplots close to each other and hide x ticks for
        # all but bottom plot.
        ax3.scatter(smp_dict['mjd'],(mcmc_flux-fake_flux)/mcmc_fluxerr,color='black',alpha=.7)
        ax3.scatter(smp_dict['mjd'],(diffim_flux-fake_flux)/diffim_fluxerr,color='green',label='Diffim Mag',alpha=.7)
        ax3.plot([peakmjd-30.,peakmjd+100],[0,0],color='black')
        ax3.set_xlim(peakmjd-30.,peakmjd+100)
        ax3.set_ylim(-3.,3.)
        ax3.set_xlabel('MJD')
        ax3.set_ylabel('(Fit Flux - Fake Flux)\n (Fit Flux Err)')
        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    
        self.savefig(lightcurves+filename+'_lightcurve.pdf')
        #print lightcurves+filename+'_lightcurve.pdf'
        '''
        # Plot the three kernel density estimates
        #numepochs = len(history[0,params.substamp**2:])
        #half = int(round(len(history[:,params.substamp**2+1])/2,0))
        fig, ax = plt.subplots(5, 5,figsize=(40, 25))
        #fig.subplots_adjust(wspace=0)

        ############# PLOT SAMPLE HISTOGRAMS/KDE #############
        #from statsmodels.nonparametric.kde import KDEUnivariate
        from scipy.stats import gaussian_kde
        axs = ax.ravel()
        axi = -1
        #half = int(round(len(modelhistory[:,0])/2,0))
        model = model[substamp**2:]
        numepochs = len(model)
        mjd = smp_dict['mjd']
        for epoch in np.arange(numepochs):
            thisepoch = modelhistory[:,substamp**2+epoch]
            try:
                axie = axi + 1
                my_pdf = gaussian_kde(thisepoch,)
                x_grid = np.linspace(min(thisepoch),max(thisepoch),num=1000)
                axs[axie].plot(x_grid,my_pdf(x_grid),'r',label=str(mjd[epoch]))    
                #axs[axie].hist(thisepoch,10,normed=True,label=str(mjd[epoch]))
                axi += 1
            except:
                print 'exception'
                continue
            mad = np.median(abs(thisepoch - np.median(thisepoch)))
            madstd = 1.48*mad
            axs[axi].axvline(np.mean(thisepoch), color='b', linestyle='dashed', linewidth=2,label='Mean')
            axs[axi].axvline(np.median(thisepoch), color='red', linestyle='dashed', linewidth=2,label='Median')
            axs[axi].axvline(np.median(thisepoch) + madstd,color='green',linestyle='-.',linewidth=2,label='Std')
            axs[axi].axvline(np.median(thisepoch) - madstd,color='green',linestyle='-.',linewidth=2,label='Std')
            axs[axi].axvline(fake_flux_adjusted[epoch], color='black', linestyle=':', linewidth=2,label='Fake Flux Zpt Adjusted')
            axs[axi].legend(prop={'size':6})
            axs[axi].set_ylabel('$\log\,L$')
            axs[axi].set_xlabel('Param Value')
            
        plt.savefig(lightcurves+filename+'_loglike.pdf')
        print lightcurves+filename+'_loglike.pdf'


        '''
        keeptrying = True
        fx = 5
        fy = 5
        axi = -1
        while keeptrying:
            fy += 1
            fx += 1
            fig, ax = plt.subplots(fx, fy,figsize=(60, 30))
            from scipy.stats.kde import gaussian_kde
            axs = ax.ravel()
            for epoch in np.arange(numepochs):
                if smp_dict['mjd_flag'][epoch] == 1 or smp_dict['flag'][epoch] == 1:
                    #thisepoch = history[half:,params.substamp**2+epoch]
                    #axs[epoch].hist(thisepoch, 25, normed=1, histtype='stepfilled',alpha=.5,color='green',label=str(smp_dict['mjd'][epoch]))
                    #axs[epoch].legend()
                    continue
                else:
                    thisepoch = history[half:,params.substamp**2+epoch]
                    my_pdf = gaussian_kde(thisepoch,)
                    x_grid = np.linspace(min(thisepoch),max(thisepoch),num=1000)
                    #axs[epoch].hist(thisepoch, 25, normed=1, histtype='stepfilled',alpha=.5,color='green',label=str(smp_dict['mjd'][epoch]))
                    try:
                        axs[axi].plot(x_grid,my_pdf(x_grid),'r',label=str(smp_dict['mjd'][epoch]))
                        keeptrying = False
                        mad = np.median(abs(thisepoch - np.median(thisepoch)))
                        madstd = 1.48*mad
                        axi += 1
                        axs[axi].axvline(np.mean(thisepoch), color='b', linestyle='dashed', linewidth=2,label='Mean')
                        axs[axi].axvline(np.median(thisepoch), color='red', linestyle='dashed', linewidth=2,label='Median')
                        axs[axi].axvline(np.median(thisepoch) + madstd,color='green',linestyle='-.',linewidth=2,label='Std')
                        axs[axi].axvline(np.median(thisepoch) - madstd,color='green',linestyle='-.',linewidth=2,label='Std')
                        axs[axi].axvline(fake_flux_adjusted[epoch], color='black', linestyle=':', linewidth=2,label='Fake Flux Zpt Adjusted')
                        axs[axi].legend(prop={'size':6})
                        axs[axi].set_ylabel('$\log\,L$')
                        axs[axi].set_xlabel('Param Value')
                        #axs[epoch].set_ylim(0,.2)
                    except:
                        keeptrying = True
                        continue


        plt.savefig(lightcurves+filename+'_loglike.pdf')
        
        print lightcurves+filename+'_loglike.pdf'
        '''
        '''
        plt.figure(10)
        plt.plot(np.arange(len(xhistory)),xhistory)
        plt.plot(np.arange(len(yhistory)),yhistory)
        plt.ylabel('Pixel offset RA and DEC')
        plt.xlabel('MCMC Iteration')
        plt.savefig('pix_shift_histories2.pdf')

        #print history[:,-1]
        #raw_input()
        '''

        '''
        plt.figure(2)
        plt.plot([0,0],[0,0])
        plt.savefig('clear.png')

        plt.figure(20)
        for h in np.arange(len(smp_dict['mjd'])):
            plt.plot(np.arange(len(history[:,-1*(h+1)])),history[:,-1*(h+1)])
            #plt.plot(np.arange(len(history[:,h])),history[:,h])
        plt.ylabel('Counts')
        plt.xlabel('MCMC Iteration')
        plt.savefig(out+'flux_histories.pdf')
        '''
        '''
        plt.figure(5)
        for h in np.arange(0,params.substamp**2,15):
            plt.plot(np.arange(len(history[-50000:,h])),history[-50000:,h])
            #plt.plot(np.arange(len(history[:,h])),history[:,h])
        plt.ylabel('Counts')
        plt.xlabel('MCMC Iteration')
        plt.savefig(out+'galflux_histories.png')
        '''

        if filt == 'u':
            hostgal_bandmag = snparams.fake_hostmag_u
        if filt == 'g':
            hostgal_bandmag = snparams.fake_hostmag_g
        if filt == 'r':
            hostgal_bandmag = snparams.fake_hostmag_r
        if filt == 'i':
            hostgal_bandmag = snparams.fake_hostmag_i
        if filt == 'z':
            hostgal_bandmag = snparams.fake_hostmag_z
        if filt== 'y':
            hostgal_bandmag = snparams.fake_hostmag_y

        outdir = os.path.join(outfile,foldername+'/np_data/'+filt+'/')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        plots = os.path.join(outfile,foldername+'/plots/'+filt+'/')
        if not os.path.exists(plots):
            os.makedirs(plots)



        fileout = os.path.join(outdir,filename+'_smp.npz')
        self.tmpwriter.savez( fileout
                  ,diffim_mag = diffim_mag
                  ,diffim_mag_err = diffim_mag_err
                  ,diffim_flux = diffim_flux
                  ,diffim_fluxerr = diffim_fluxerr
                  ,mjd = smp_dict['mjd']
                  ,mcmc_mag = mcmc_mag
                  ,mcmc_mag_err = mcmc_mag_err
                  ,mcmc_flux = mcmc_flux
                  ,mcmc_fluxerr = mcmc_fluxerr
                  ,skyerr = smp_dict['skyerr']
                  ,mjd_flag = smp_dict['mjd_flag']
                  ,fakemag = smp_dict['fakemag']
                  ,zpterr = smp_dict['zpterr']
                  ,zpt = smp_dict['zpt']
                  ,fakezpt = smp_dict['fakezpt']
                  ,peakmjd=snparams.peakmjd
                  ,hostgal_bandmag = hostgal_bandmag
                  ,hostgal_sb_fluxcal= snparams.hostgal_sb_fluxcal
                  ,psf = smp_dict['psf']
                  ,fwhm = smp_dict['fwhm_arcsec']*2.355
                  #,model = model
                  #,uncertainty = uncertainty
                  #,paramhistory = history
                  ,substamp = substamp
                  ,acceptedrate = acceptedrate
                  )

        fileout = os.path.join(outdir,filename+'_fluxhistory.npz')
        self.tmpwriter.savez( fileout
                ,history = modelhistory
                ,substamp = substamp)

        #fileout = os.path.join(outdir,filename+'_galfluxhistory.npz')
        #self.tmpwriter.savez( fileout
        #        ,history = galmodelhistory)
        
        sindices = []
        for i in np.arange(len(smp_dict['mjd'])):
            if smp_dict['mjd_flag'][i] == 1 or smp_dict['flag'][i] == 1:
                continue
            else:
                if max(history[half:,substamp**2+i])-min(history[half:,substamp**2+i]) == 0.:
                    continue
                else:
                    sindices.append(int(i))
        
        
        import triangle
        shistory = history[half:,params.substamp**2:]
        parray = np.zeros((len(shistory[:,0]),len(sindices)))
        thismjds = []
        for i in np.arange(len(sindices)):
            parray[:,i] = shistory[:,sindices[i]]
            thismjds.append(str(smp_dict['mjd'][sindices[i]]))
        #p = np.take(shistory,np.arange(half-1,dtype=int),sindices)
        #pmjds = np.take(np.asarray(smp_dict['mjd']),sindices)
        figure = triangle.corner(parray,labels=thismjds,
                                show_titles=True, title_args={"fontsize": 12})
        figure.gca().annotate("A Title", xy=(0.5, 1.0), xycoords="figure fraction",
                              xytext=(0, -5), textcoords="offset points",
                              ha="center", va="top")
        figure.savefig(lightcurves+filename+'_covar.pdf')
        
        print 'DONEEEE'
        sys.exit()
        print 'mcmc worked!'
        print "first_result"
        print first_result
        
        for i in range(len(first_result.params)):
            mpdict[i]['value'] = first_result.params[i]
        if verbose: print('Creating Final Scene Model')
        second_result = mpfit(scene,parinfo=mpdict,functkw=mpargs, debug = True, quiet=False)
        print "second_result"
        print second_result

        chi2 = scene_check(second_result.params,x=smp_psf,y=smp_im,err=smp_noise,params=params)
        # write the results to file
        fout = open(outfile,'w')
        print >> fout, '# MJD ZPT Flux Fluxerr Mag Magerr pkflux pkfluxerr xpos ypos chi2 mjd_flag flux_firstiter fluxerr_firstiter mag_firstiter magerr_firstiter'
        for i in range(len(smp_dict['snx'])):
            print "first result error"
            print  first_result.perror[params.substamp**2.+i]
            print "type of first result error"
            print type(first_result.perror[params.substamp**2.+i])
            print >> fout, '%.1f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.2f %.2f %.2f %i %.3f %.3f %.3f %.3f'%(smp_dict['mjd'][i],smp_dict['zpt'][i],
                                                                                       second_result.params[params.substamp**2.+i],
                                                                                       second_result.perror[params.substamp**2.+i],
                                                                                       (-2.5*np.log10(second_result.params[params.substamp**2.+i]) + 31) ,
                                                                                       0.000,
                                                                                       smp_dict['scale'][i],smp_dict['scale_err'][i],
                                                                                       smp_dict['snx'][i],smp_dict['sny'][i],chi2[i],
                                                                                       smp_dict['mjd_flag'][i],
                                                                                       first_result.params[params.substamp**2.+i],
                                                                                       first_result.perror[params.substamp**2.+i],
                                                                                       (-2.5*np.log10(first_result.params[params.substamp**2.+i]) + 31) ,
                                                                                       0.000)
        fout.close()
        #self.big_zpt_plot()
        print('SMP was successful!!!')
    '''    

    def savefig(self,fname):
        tempfile = os.path.join(self.tmpwriter.tmpdir,'tmp_'+self.tmpwriter.tmp_index+'.png')
        plt.savefig(tempfile)
        if os.path.isfile(fname):
            os.remove(fname)
        os.system('mv '+tempfile+' '+fname)
        print 'saved',fname

    def getfluxsmp(self,im,psf,sky,weight,fitrad,gal,mjd,skyerr,gain,guess_scale=100000,index='',mypsf=None,imfile=None,x=None,y=None,pdf_pages=None,dosimultaneous=True):
        #print 'inside getfluxsmp'
        chisqvec = []
        fluxvec = []
        bad = False
        #galconv = scipy.signal.fftconvolve(gal,psf,mode='same')
        galconv = 0.
        #substamp = galconv.shape[0]
        #Make a mask with radius
        # fitrad = np.zeros([substamp,substamp])
        psf = psf/np.sum(psf.ravel())

        # fluxerrls = cov
        # # print fluxls,cov
        #
        #


        # def f(x):
        #     return (x*psf.ravel()-im.ravel()+sky.ravel())/skyerr
        #
        # m = Minuit(f,x=fluxls)
        # m.migrad()
        # #print m.values['x'],flux,fluxls  # {'x': 2,'y': 3,'z': 4}
        # #print m.errors['x'],errmag,cov  # {'x': 1,'y': 1,'z': 1}
        # fluxminuit = m.values['x']
        # fluxerrminuit = m.errors['x']

        def f(scale,x):
            #scale = prms['scale']
            #power = prms['pow']
            return scale*psf.ravel()
        def fastier(x):
            scale = x
            return 1000. + scale * np.sum((1./skyerr**2)*psf*psf) - np.sum((1./skyerr**2)*psf*(im-sky))

        # params = Parameters()
        # params.add('scale', value=guess_scale, min=1.)
        # #params.add('pow', value=.5, vary=False)
        #
        # fitter = Minimizer(f, params)
        # try:
        #     v = fitter.minimize(method='brute')
        # except:
        #     print 'FAILED'*5
        #     return 1, 1, 1, 1, 1, True
        # #print fluxls,v.params['scale'].value
        # #print v.params['scale'].__dict__
        # fluxlm = v.params['scale'].value
        # fluxerrlm = v.params['scale'].stderr
        # guess_scale = fluxlm
        #
        # from scipy.optimize import curve_fit
        # popt, pcov = curve_fit(fastier, [0,0,0], [1000.,1000,1000], p0=guess_scale)

        scales = []

        #for i in range(guess_scale-2000,guess_scale+2000,.5):



        # print v
        # print v.params
        vals = \
            mpfitexpr.mpfitexpr("p[0]*x", psf.ravel()*fitrad.ravel(), (im.ravel() - sky.ravel())*fitrad.ravel(), np.sqrt(skyerr**2), [1], full_output=True)[0]
        try:
            errmag = vals.perror[0]
            fluxmp = vals.params[0]
            fluxerrmp = errmag
        except:
            return 1, 1, 1, 1, 1, True


        # def f(prms):
        #     scale = prms['scale']
        #     #power = prms['pow']
        #     return ((scale * psf.ravel() - im.ravel() + sky.ravel()) / np.sqrt(skyerr**2))*fitrad.ravel()
        #
        # params = Parameters()
        # params.add('scale', value=fluxmp, min=1.)
        # #params.add('pow', value=.5, vary=False)
        #
        # fitter = Minimizer(f, params)
        # #try:
        # v = fitter.minimize(method='leastsq')
        # #except:
        # #    print 'FAILED'*5
        # #    return 1, 1, 1, 1, 1, True
        # # #print fluxls, v.params['scale'].value
        # # #print v.params['scale'].__dict__
        # fluxlm = v.params['scale'].value
        # fluxerrlm = v.params['scale'].stderr

        # from iminuit import Minuit
        # #def f(x):
        # #    return x*psf.ravel()
        #
        # m = Minuit(fastier, x=fluxmp, print_level=1)
        # m.migrad(ncall=100000)
        # m.hesse()
        # print(m.values)  # {'x': 2,'y': 3,'z': 4}
        # print(m.errors)  # {'x': 1,'y': 1,'z': 1}

        # quad_model = Model(f)
        # data = RealData((im.ravel()-sky.ravel())*0.,im.ravel()-sky.ravel(), sy=skyerr)
        # odr = ODR(data, f, beta0=[fluxmp])
        # # Run the regression.
        # out = odr.run()
        # # Use the in-built pprint method to give us results.
        # out.pprint()
        # # print vals

        #fluxls, cov = opti.leastsq(starresid, fluxmp, args=(psf, im, skyerr, fitrad, sky, gain),
        #                           full_output=False)

        # print 'mpfit',fluxmp,errmag,skyerr
        # raw_input()
        #print 'lsq',fluxls,cov
        #print 'lmfit',fluxlm,fluxerrlm
        #print 'astier'#,skyer
        #raw_input()
        # #print len(cov)
        # # raw_input('comparison')
        # sim =  sky + fluxlm * psf
        # gain = 1.
        # weight = 1. / (skyerr ** 2 + max([0,float(fluxlm)]) / gain)
        # mchisq = np.sum((im - sim) ** 2 * weight * fitrad)
        # ndof = len(fitrad[fitrad==1].ravel())
        # guess_scale = fluxlm
        # dosimultaneous = False
        # if dosimultaneous:
        #     # totalarea = 0
        #     # for x in np.arange(substamp):
        #     #     for y in np.arange(substamp):
        #     #         #print np.sqrt((substamp/2. - x)**2 + (substamp/2. - y)**2)
        #     #         if np.sqrt((substamp/2. - x)**2 + (substamp/2. - y)**2) <= 13:
        #     #             totalarea+=1
        #     totalarea = len(fitrad[fitrad > 0])
        #     #print 'totalarea',totalarea
        #     #print 'skyerr2',skyerr**2
        #     guessrange = None
        #     if guess_scale is None:
        #         for i in np.arange(-1000, 10000000, 1000):
        #             sim = galconv + sky + i * psf
        #             #sigtot = np.sqrt((skyerr/4.) + abs(float(i))/4.)
        #             weight = 1./(skyerr**2 + abs(float(fluxlm))/gain+ 1.) #holtzman
        #             chisqvec.append(np.sum((im - sim) ** 2 * weight * fitrad))
        #             fluxvec.append(i)
        #             #print 'sigtot',sigtot,'weight',weight,'chisqvec',chisqvec[-1]
        #
        #         fluxvec = np.array(fluxvec)
        #         chisqvec = np.array(chisqvec)
        #         #print 'argmin guesscale',np.argmin(chisqvec)
        #         guess_scale = fluxvec[np.argmin(chisqvec)]
        #         guessrange = 1000
        #
        #     #print guess_scale
        #     #raw_input()
        #     chisqvec = []
        #     fluxvec = []
        #     if guessrange is None:
        #         guessrange = .15 * abs(guess_scale)
        #     try:
        #         guess_scale_step = min([abs(guess_scale) / 1000., 1.])
        #         for i in np.arange(guess_scale - guessrange, guess_scale + guessrange, guess_scale_step):
        #             sim = galconv + sky + i * psf
        #             #sigtot = np.sqrt(skyerr ** 2 + abs(float(i)) / 4.)
        #             weight = 1./(skyerr**2 + abs(float(fluxlm))/gain ) #holtzman
        #             #weight = 1./((skyerr/4.) + abs(float(i))/4. + 1.)#first time around
        #             chisqvec.append(np.sum((im - sim) ** 2 * weight * fitrad))
        #             fluxvec.append(i)
        #     except:
        #         bad=True
        #
        # else:
        #
        #     # guessrange = None
        #     # if guess_scale is None:
        #     #     for i in np.arange(-100000,2000000,5000):
        #     #         sim = galconv + sky + i*psf
        #     #         chisqvec.append(np.sum((im-sim)**2*weight*fitrad))
        #     #         #print i,weight,chisqvec[-1]
        #     #
        #     #         fluxvec.append(i)
        #     #     fluxvec = np.array(fluxvec)
        #     #     chisqvec = np.array(chisqvec)
        #     #     guess_scale = fluxvec[np.argmin(chisqvec)]
        #     #     guessrange = 5000
        #
        #     chisqvec = []
        #     fluxvec = []
        #     # if guessrange is None:
        #     #     guessrange = .1*guess_scale
        #     guessrange = 2000.
        #     guess_scale_step = 1.#min([guess_scale/1000.,1.])
        #     weight = 1. / (skyerr ** 2)
        #     guess_scale = fluxmp
        #     try:
        #         for i in np.arange(guess_scale-guessrange,guess_scale+guessrange,guess_scale_step):
        #             sim = sky + i*psf
        #             chisqvec.append(np.sum((im-sim)**2*weight*fitrad))
        #             fluxvec.append(i)
        #     except:
        #         bad = True
        #
        # if not bad:
        #     ii = fitrad.ravel()
        #     i = ii[ii != 0]
        #
        #     ndof = len(i)
        #
        #     fluxvec = np.array(fluxvec)
        #     chisqvec = np.array(chisqvec)
        #     try:
        #         hh = chisqvec*0 + min(chisqvec)
        #     except:
        #         bad = True
        #
        # if not bad:
        #     mchisq = min(chisqvec)
        #     idx = np.isclose(chisqvec, hh, atol=1.0)
        #
        #     argm =  np.argmin(chisqvec)
        #
        #     try:
        #         sim = galconv + sky + fluxvec[argm]*psf
        #     except:
        #         bad = True
        #
        # from scipy.optimize import curve_fit
        # from lmfit import Model
        # def f(scale):
        #     #scale = prms['scale']
        #     #power = prms['pow']
        #     return (scale * psf.ravel() - im.ravel() + sky.ravel()) / (skyerr + fluxlm**.5)
        # gmodel = Model(f)
        # params = gmodel.make_params(scale=fluxlm)
        # result = gmodel.fit(f,params)
        # print(result.fit_report())
        # plt.clf()
        # plt.plot(fluxvec,chisqvec)
        # plt.xlim(fluxmp-10000.,fluxmp+10000.)
        # #plt.ylim(chisqvec[argm]-1.,chisqvec[argm]+10.)
        # plt.savefig('fluxtest.png')
        #print 'mychisq',fluxvec[argm], fluxvec[argm] - fluxvec[idx][0]
        print 'mpfit',fluxmp,fluxerrmp
        #print 'lmfit',fluxlm,fluxerrlm
        #print 'minchisq',chisqvec[argm]
        #print 'minchisq/ndof',chisqvec[argm]/len(fitrad[fitrad == 1].ravel())
        #raw_input()
        #print 'lmfit',fluxlm,fluxerrlm
        # print result.__dict__

        if not bad:
            if not pdf_pages is None:
                fig = plt.figure(figsize=(20,10))
                axim = plt.subplot(141)
                axpsf = plt.subplot(142)
                axdiff = plt.subplot(143)
                axchi = plt.subplot(144)
                for ax, title in zip([axim, axpsf, axdiff, axchi], ['image', 'model','resid', 'chisq']):
                    ax.set_title(title)
                ax = axim.imshow(im*fitrad, cmap='gray', interpolation='nearest')
                cbar = fig.colorbar(ax,ax=axim)
                ax = axpsf.imshow(sim*fitrad, cmap='gray', interpolation='nearest')
                cbar = fig.colorbar(ax,ax=axpsf)
                ax = axdiff.imshow((im - sim)*fitrad, cmap='gray', interpolation='nearest',vmin=-100, vmax=10)
                cbar = fig.colorbar(ax,ax=axdiff)
                ax = axchi.imshow((im - sim)**2*weight*fitrad, cmap='gray', interpolation='nearest',vmin=0, vmax=10.)
                cbar = fig.colorbar(ax,ax=axchi)
                # plt.imshow((subim-scaledpsf)/imhdr['SKYSIG'],cmap='gray',interpolation='nearest')
                # plt.colorbar()
                plt.title(title)
                pdf_pages.savefig(fig)

            #if not mypsf is None:
            #    self.tmpwriter.savefits(mypsf, '/pnfs/des/scratch/pysmp/test/' + str(index) + '_sim.fits')
            # if not imfile is None:
            #     imstamp = pf.getdata(imfile)
            #     print 'imstamp shape', imstamp.shape
            #     print x, y
            #     xo = round(x)
            #     yo = round(y)
            #     imstamp = imstamp[yo - 17:yo + 17 + 1, xo - 17:xo + 17 + 1]
            #     imstamp = imstamp / np.sum(imstamp)
            sim = sky + fluxmp*psf
            sum_data_minus_sim = np.sum(im-sim)
            #sim = galconv + sky + fluxvec[argm]*psf
            mchisq = np.sum((im - sim) ** 2 * 1./(1./weight**2)**.5 * fitrad)
            ndof = len(fitrad[fitrad == 1].ravel())#+10000000
        # if not bad:
        #     if fluxvec[argm] > 8000000:
        #         bad = True
        #         print 'star too bright...'
        #     elif fluxvec[argm] < 100.:
        #         bad = True
        #         print 'star too dim...'
        sim = sky + fluxmp * psf
        print 'chisq',np.sum((im - sim) ** 2 /(skyerr**2+fluxmp)  * fitrad)
        if not bad:
            #return fluxvec[argm], fluxvec[argm] - fluxvec[idx][0], mchisq/ndof, sum_data_minus_sim, np.sum((im - sim) ** 2 * weight * fitrad)/ndof, bad
            return fluxmp, np.sqrt(fluxerrmp**2+fluxmp), np.sum((im - sim) ** 2 /(skyerr**2+fluxmp)  * fitrad), sum_data_minus_sim, np.sum((im - sim) ** 2 /(skyerr**2+fluxmp) * fitrad)/ndof, bad
        else:
            return 1,1,1,1,1,True


    def iterstat(self,d,startMedian=False,sigmaclip=3.0,
             iter=6):
        """Get the sigma-clipped mean of 
        a distribution, d.
        Usage: mean,stdev = iterstat.iterstat
        Input:
        d:           the data
        Optional Inputs:
        sigmaclip:   number of standard deviations to clip
        startMedian: if True, begin with the median of the distribution
        iter:        number of iterations
        """
        
        clip=sigmaclip
        img=d.astype('float64')
        if startMedian:
            md=np.median(img)
        else:
            md=np.mean(img)
        n = float(len(img))
        std = np.sqrt(np.sum((img-md)**2.)/(n-1))

        for ii in range(iter):
            gd=np.where((img < md+clip*std) &
                        (img > md-clip*std))

            md=np.mean(img[gd])
            num = len(img[gd])
            n = float(len(gd[0]))
            std = np.sqrt(np.sum((img[gd]-md)**2.)/(n-1.))

        return(md,std,num)

    def teststarpos(self,fakestarfile,wcs,zpt,sky,skyerr,im,noise,mask,psffile,imfile,snparams,substamp,diffimzp,psf='None'):
        #print zpt
        #raw_input()
        ccds,ras,decs,extra = np.loadtxt(fakestarfile, unpack=True)
        ras = ras[ccds == float(self.ccdnum)][0:100]
        decs = decs[ccds == float(self.ccdnum)][0:100]
        #from matplotlib import colors

        colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']

        coords = zip(*wcs.wcs_world2pix(np.array(zip(ras,decs)),0))
        x_star,y_star = [],[]
        for xval,yval in zip(*coords):
            x_star += [xval]
            y_star += [yval]

        x_star1,y_star1 = np.array(x_star),np.array(y_star)
        mag1,magerr1,flux1,fluxerr1,sky1,skyerr1,badflag1,outstr1 = \
            aper.aper(im,x_star1,y_star1,apr = params.fitrad,verbose=False)
        x_star,y_star = cntrd.cntrd(im,x_star1,y_star1,params.cntrd_fwhm)

        #print len(x_star)
        #raw_input()
        pixel_offsets = np.arange(-.5,.5,.1)
        flux_star = np.zeros((len(x_star),len(pixel_offsets)))
        flux_star_std = np.zeros((len(x_star),len(pixel_offsets)))
        
        for j,po in enumerate(pixel_offsets):
            radius = 10.
            fitrad = np.zeros([substamp,substamp])
            for x in np.arange(substamp):
                for y in np.arange(substamp):
                    if np.sqrt((substamp/2. - x)**2 + (substamp/2. - y)**2) < radius:
                        fitrad[int(x),int(y)] = 1.
            for x,y,s,se,i in zip(x_star,y_star+po,sky,skyerr,range(len(x_star))):
                if x > 51 and y > 51 and x < self.snparams.nxpix-51 and y < self.snparams.nypix-51:
                    if self.stardumppsf:
                        if self.snparams.psf_model.lower() == 'psfex':
                            psf, psfcenter = self.build_psfex(psffile,x,y,imfile,psfexworked=psfexworked)

                            #print psf.shape
                        elif psf == '':
                            raise exceptions.RuntimeError("Error : PSF array is required!")
                    else:
                        psf, psfcenter = self.psf, self.psfcenter
                        print 'psfcenter',psfcenter

                    pk = pkfit_norecent_noise_smp.pkfit_class(im,psf,psfcenter,self.rdnoise,self.gain,noise,mask)
                    try:
                        errmag,chi,niter,scale,iylo,iyhi,ixlo,ixhi,image_stamp,noise_stamp,mask_stamp,psf_stamp = pk.pkfit_norecent_noise_smp(1,x,y,s,se,params.fitrad,returnStamps=True,stampsize=params.substamp)
                        noise_stamp[noise_stamp > 0.] = 1
                        noise_stamp[noise_stamp <= 0.] = 0
                        noise_stamp = noise_stamp*1/(se**2)
                        gal = np.zeros(image_stamp.shape)
                        mjd = 000.
                        cscale,cscale_std,chisq,dms = self.getfluxsmp(image_stamp,psf_stamp,s,noise_stamp,fitrad,gal,mjd,scale)
                        scale = cscale
                        flux_star[i,j] = copy(scale)
                        flux_star_std[i,j] = copy(cscale_std)

                    except:
                        continue
                else:
                    print 'star not within bounds of image',x,y

        if os.path.isfile('fake20staroffsets.txt'):
            b = open('fake20staroffsets.txt','a')
        else:
            b = open('fake20staroffsets.txt','w')
            b.write('ax2\tbx\tc\n')


        plt.clf()
        print flux_star.shape
        print len(colors)
        for j,po in enumerate(pixel_offsets):
            for k in np.arange(4):
                plt.scatter(flux_star[k,j]*0.+po,-2.5*np.log10(flux_star[k,j])+diffimzp,alpha=.5)
        from scipy.optimize import curve_fit
        for k in np.arange(4):
            popt, pcov = curve_fit(parabola, pixel_offsets, -2.5*np.log10(flux_star[k,:])+diffimzp)
            perr = np.sqrt(np.diag(pcov))
            if perr[0] < .5:
                plt.plot(pixel_offsets,popt[0]*(pixel_offsets)**2+popt[1]*(pixel_offsets)+popt[2],color='black')
                b.write(str(popt[0])+'\t'+str(popt[1])+'\t'+str(popt[2])+'\n')

        plt.xlabel('pixel offset')
        plt.ylabel('fit star Mag')
        ax = plt.gca()
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        ax.get_yaxis().get_major_formatter().set_scientific(False)

        plt.savefig('teststaroffsety.png')
        print 'saved teststaroffsety.png'

        plt.clf()
        for j,po in enumerate(pixel_offsets):
            plt.scatter(flux_star[:,j]*0.+po,-2.5*np.log10(flux_star[:,j]+flux_star_std[:,j])+2.5*np.log10(flux_star[:,j]),alpha=.5)
        plt.xlabel('pixel offset')
        plt.ylabel('fit star flux error')

        plt.savefig('teststaroffsetyerror.png')
        print 'saved teststaroffsetyerror.png'

        b.close()
        

    def plotallfake20staroffsets(self):
        ays,bs,cs = np.loadtxt('fake20staroffsets.txt',unpack=True,skiprows=1)
        plt.clf()
        pixel_offsets = np.arange(-.5,.5,.05)
        for a,b,c in zip(ays,bs,cs):
            plt.plot(pixel_offsets,a*pixel_offsets**2+b*pixel_offsets+c,alpha=.2)
            
        plt.savefig('allstaroff.png')

        plt.clf()
        pixel_offsets = np.arange(-.5,.5,.05)
        ploty = pixel_offsets*0.
        minoffset = []
        for a,b,c in zip(ays,bs,cs):
            mo = pixel_offsets[np.argmin(a*pixel_offsets**2+b*pixel_offsets+c)]
            if mo < 1.:
                minoffset.append(mo)
            
        plt.hist(minoffset,bins=np.arange(-1.025,1,.05),color='black')
        plt.xlabel('(centroided y-pixel) - (y-pixel of brightest flux)')
        plt.ylabel('counts')
        plt.savefig('allstaroffcombined.png')

        print 'allstaroffcombined.png'

    def fake20magstars(self,fakestarfile,wcs,zpt,sky,skyerr,im,noise,mask,psffile,imfile,snparams,substamp,psf='None'):
        ccds,ras,decs,extra = np.loadtxt(fakestarfile, unpack=True)
        ras = ras[ccds == float(self.ccdnum)]
        decs = decs[ccds == float(self.ccdnum)]
        coords = zip(*wcs.wcs_world2pix(np.array(zip(ras,decs)),0))
        x_star,y_star = [],[]
        for xval,yval in zip(*coords):
            x_star += [xval]
            y_star += [yval]
            
        x_star1,y_star1 = np.array(x_star),np.array(y_star)
        mag1,magerr1,flux1,fluxerr1,sky1,skyerr1,badflag1,outstr1 = \
            aper.aper(im,x_star1,y_star1,apr = params.fitrad)
        x_star,y_star = cntrd.cntrd(im,x_star1,y_star1,params.cntrd_fwhm)

        flux_star = np.array([-999.]*len(x_star))
        flux_star_std = np.array([-999.]*len(x_star))

        radius = 10.
        fitrad = np.zeros([substamp,substamp])
        for x in np.arange(substamp):
            for y in np.arange(substamp):
                if np.sqrt((substamp/2. - x)**2 + (substamp/2. - y)**2) < radius:
                    fitrad[int(x),int(y)] = 1.
        for x,y,s,se,i in zip(x_star,y_star,sky,skyerr,range(len(x_star))):
            if x > 51 and y > 51 and x < self.snparams.nxpix-51 and y < self.snparams.nypix-51:
                if self.stardumppsf:
                    if self.snparams.psf_model.lower() == 'psfex':
                        psf, psfcenter = self.build_psfex(psffile,x,y,imfile,psfexworked=psfexworked)
                        print psf.shape
                    elif psf == '':
                        raise exceptions.RuntimeError("Error : PSF array is required!")
                else:
                    psf, psfcenter = self.psf, self.psfcenter
                    print 'psfcenter',psfcenter

                pk = pkfit_norecent_noise_smp.pkfit_class(im,psf,psfcenter,self.rdnoise,self.gain,noise,mask)
                try:
                    errmag,chi,niter,scale,iylo,iyhi,ixlo,ixhi,image_stamp,noise_stamp,mask_stamp,psf_stamp = pk.pkfit_norecent_noise_smp(1,x,y,s,se,params.fitrad,returnStamps=True,stampsize=params.substamp)
                    noise_stamp[noise_stamp > 0.] = 1
                    noise_stamp[noise_stamp <= 0.] = 0
                    noise_stamp = noise_stamp*1/(se**2)
                    gal = np.zeros(image_stamp.shape)
                    mjd = 000.
                    cscale,cscale_std,chisq,dms = self.getfluxsmp(image_stamp,psf_stamp,s,noise_stamp,fitrad,gal,mjd,scale)
                    scale = cscale
                    flux_star[i] = copy(scale)
                    flux_star_std[i] = copy(cscale_std)
                    
                except:
                    continue
        
        print flux_star
        #print 'moneyflux'
        #raw_input()
        return flux_star,flux_star_std,flux_star*0.+zpt



    def getzpt(self,xstar,ystar,ras, decs,ids,mags,sky,skyerr,thismjd,
                badflag,mag_cat,im,noise,mask,maskfile,weightsfile,psffile,imfile,imwcs,snparams,substamp,
                mjdoff,mjdslopeinteroff,j,longimfile,bkgrnd,bkgrndrms,psf='',mjd=None,
                mpfit_or_mcmc='mpfit',cat_zpt=-999):
        """Measure the zeropoints for the images"""

        #print xstar,ystar
        #raw_input('hhh')


        # xstar, ystar, ras, decs, mags, sky, skyerr, mag_cat = xstar[~badflag], ystar[~badflag], ras[~badflag], \
        #                                                                decs[~badflag], mags[~badflag],\
        #                                                                sky[~badflag], skyerr[~badflag], mag_cat
        xstar, ystar = cntrd.cntrd(im, xstar, ystar, params.cntrd_fwhm)
        # if snparams.survey == 'DES':
        #     xstar += 1.
        #     ystar += 1.
        thisra, thisdec = zip(*imwcs.wcs_pix2world(np.array(zip(xstar, ystar)), 0))
        thisra = np.array(thisra)
        thisdec = np.array(thisdec)
        thisids = np.array(ids)
        print 'Computing zeropoint for',imfile
        print '\n'
        import pkfit_norecent_noise_smp
        counter = 0
        bad = False
        #print 'skies',sky
        #raw_input()
        flux_star = np.array([-999.]*len(xstar))        
        flux_star_std = np.array([-999.]*len(xstar))
        flux_chisq = np.array([-999.]*len(xstar))
        starsky = np.array([-999.]*len(xstar))
        starskyerr = np.array([-999.]*len(xstar))
        psfs = np.zeros((len(xstar),substamp,substamp))
        flux_mychisq = np.array([-999.]*len(xstar))
        flux_dms = np.array([-999.]*len(xstar))
        gsflux = np.array([-999.]*len(xstar))
        gsflux_std = np.array([-999.]*len(xstar))
        gsflux_chisq = np.array([-999.]*len(xstar))
        gsflux_dms = np.array([-999.]*len(xstar))
        isnotcheckstars = np.ones(len(xstar))
        flux_star_mcmc = np.array([-999.]*len(xstar))
        #flux_star_std_mcmc = np.array([-999.]*len(xstar))
        flux_star_mcmc_modelerrors = np.array([-999.]*len(xstar))
        #flux_star_std_mcmc_modelerrors = np.array([-999.]*len(xstar))
        flux_star_mcmc_me_simple = np.array([-999.]*len(xstar))
        #flux_star_std_mcmc_me_simple = np.array([-999.]*len(xstar))
        flux_star_mcmc_me_weighted = np.array([-999.]*len(xstar))
        #flux_star_std_mcmc_me_weighted = np.array([-999.]*len(xstar))
        #mcmc_mag_std = np.array([-999.]*len(xstar))
        mcmc_me_mag_std = np.array([-999.]*len(xstar))

        radius = 9
        cntr = 0
        fitrad = np.zeros([substamp,substamp])
        for x in np.arange(substamp):   
            for y in np.arange(substamp):
                if np.sqrt((substamp/2. - x)**2 + (substamp/2. - y)**2) < radius:
                    fitrad[int(x),int(y)] = 1.
        #print 'xstarrrrrrr222222',len(xstar)
        if self.dogalsimpixfit:
            big_fft_params = galsim.GSParams(maximum_fft_size=2024000)
            full_data_image = galsim.fits.read(imfile)
        #print 'starfits_'+str(thismjd)
        #pdf_pages = PdfPages('starfits_'+str(thismjd)+'.pdf')
        if self.savezptstamps:
            pdf_pagesc = PdfPages('starfits_'+str(thismjd)+'.pdf')
        else:
            pdf_pagesc = None

        # self.dosextractor = True
        # if self.dosextractor:
        #     runsextractor.getsky_and_skyerr(imfile, im, 0, 100,  100, 100, snparams.survey)
        #     bkgrnd = pf.getdata(imfile+'.background')


        #print imfile
        #print thismjd
        #print 'mjdabove'
        #raw_input()
        #for ra,dec,x,y,s,se in zip(ras,decs,xstar,ystar,sky,skyerr):
        #    print ra,dec,x,y,s,se
        #raw_input()
        #sys.exit()
        # if self.snparams.survey == 'PS1':
        #     print mag_cat
        #     raw_input()
        #     mag_cat *= -1.
        prevra = 0
        for x,y,m,s,se,mc,ra,dec,iii,bf,i in zip(xstar,ystar,mags,sky,skyerr,mag_cat,ras,decs,ids,badflag,range(len(xstar))):

            #print 'sky',s,'skyerr',se
            #raw_input()
            #print x,y,mc
            #print self.snparams.nxpix, self.snparams.nypix
            if round(prevra,5) == round(ra,5):
                #print 'same star, so skipping'
                continue
            prevra = ra
            if bf == 1:
                continue
            #print i
            #cntr += 1
            #if i > 150:
            #    continue
            #print 'xstar',xstar
            #raw_input()
            #y -= 2.
            #if i < float(params.numcheckstars):
            #    isnotcheckstars[i] = 0

            if mc > 21.:
                print mc,'star too dim'
                continue
            if mc < 16.5:
                print mc,'star too bright'
                continue
            print i,x,y,mc,s,se
            #print 'nxpix',self.snparams.nxpix
            #print 'nypix',self.snparams.nypix
            #print x,y
            if x > 51 and y > 51 and x < self.snparams.nxpix-51 and y < self.snparams.nypix-51:# and s > 25. and se < 1000.:
                if self.stardumppsf:
                    if self.snparams.psf_model.lower() == 'psfex':
                        try:
                            psf, psfcenter = self.build_psfex(psffile,x,y,imfile,stop=True,psfexworked=psfexworked)
                        except:
                            continue
                        #print 'psfcenter',psfcenter
                        #raw_input()
                        #psf2, psfcenter2 = self.build_psfex(psffile,x+.05,y+.05,imfile,stop=True)
                        #opsf, opsfcenter = self.build_psfex(psffile, np.floor(x) + .2, np.floor(y) + .2, imfile)
                        #ppsf, ppsfcenter = self.build_psfex(psffile, np.floor(x) + .4, np.floor(y) + .4, imfile)
                        #self.tmpwriter.savefits(opsf - ppsf, '/pnfs/des/scratch/pysmp/test/psfsub.fits')
                        #sys.exit()
                        #print psf.shape
                    elif psf == '':
                        raise exceptions.RuntimeError("Error : PSF array is required!")
                else:
                    # pk = pkfit_norecent_noise_smp.pkfit_class(im,self.psf,self.psfcenter,self.rdnoise,self.gain,weights,mask)
                    # #pk = pkfit_norecent_noise_smp.pkfit_class(im,self.gauss,self.psf,self.rdnoise,self.gain,noise,mask)
                    # try:
                    #     errmag,chi,niter,scale,iylo,iyhi,ixlo,ixhi,image_stamp,noise_stamp,mask_stamp,psf_stamp = \
                    #         pk.pkfit_norecent_noise_smp(1,x,y,s,se,self.params.fitrad,returnStamps=True,stampsize=self.params.substamp)
                    # except ValueError:
                    #     raise ValueError('SN too close to edge of CCD!')
                    #if not self.psfcenter is None:
                    #    psf, psfcenter = self.psf, self.psfcenter
                    #else:
                    #print 'xyxyxyxy'
                    psf, psfcenter = self.psf, (x,y)
                    #print psfcenter
                
                counter += 1
                mask = mask*0.
                #print 'index,ra,dec,x,y',i,ra,dec,x,y
                #ppp = 'index,ra,dec '+str(i)+' '+str(ra)+' '+str(dec)
                #if self.fermilog:
                #    self.tmpwriter.appendfile(ppp + '\n', self.fermilogfile)
                if self.snparams.survey == 'PS1':

                    self.dosextractor = False
                    if self.dosextractor:
                        if imfile == '':
                            badflag[i] = 1
                            mag_cat[i] = 99
                            bad = True
                            sexsky = 1000.
                            sexrms = 10.
                        else:
                            print imfile
                            try:
                                sexsky, sexrms = runsextractor.getsky_and_skyerr(imfile,weightsfile,im, x - 100, x + 100, y - 100, y + 100,
                                                                         snparams.survey)
                            except IOError:
                                print 'could not run sextractor'
                                badflag[i] = 1
                                mag_cat[i] = 99
                                bad = True
                                sexsky = 1000.
                                sexrms = 10.
                        #print s, sexsky
                        #print se, sexrms
                        #raw_input('comparison')
                        s = sexsky
                        se = sexrms
                    else:
                        sexsky, sexrms = s, se

                    dosextractor = False
                    if dosextractor:
                        xl = max([x-50,0])
                        xh = min([x + 50, im.shape[1]])
                        yl = max([x - 50, 0])
                        yh = min([y + 50, im.shape[0]])

                        s  = np.mean(bkgrnd[yl:yh,xl:xh].ravel())
                        se = np.mean(bkgrndrms[yl:yh,xl:xh].ravel())

                    scale,errmag,chi,dms,good,image_stamp,simstamp,psf = chkpsf.fit(imfile.split('.fits')[0],xpos=x+1,ypos=y+1,ra=ra,dec=dec,
                                                                 pdf_pages=pdf_pagesc,radius=params.substamp/2.-1.,
                                                                 title=str(ra)+' '+str(dec)+' '+str(i),
                                                                 maskfile=maskfile,mysky=s,mysig=se)

                    if scale == 1:
                        badflag[i] = 1
                        mag_cat[i] = 99
                        print 'bad star'
                        bad = True
                    if not good:
                        #badflag[i] = 1
                        #mag_cat[i] = 99
                        #print 'badflaggg'*10
                        print 'notgood'
                        #print 'star not in fcmp... still including in fit.'
                        #bad = True
                else:
                    #print 'here1'
                    #pk = pkfit_norecent_noise_smp.pkfit_class(im, psf/np.sum(psf), psfcenter, self.rdnoise, self.gain,
                    #                                      noise*0.+1., mask)
                    #pk = pkfit_norecent_noise_smp.pkfit_class(im,psf/np.sum(psf),psfcenter,self.rdnoise,self.gain,noise,mask)
                    #Run for MPFIT
                    #print 'initialized'
                    #try:
                    if True:
                        #try:
                        #print 'here2'
                        #errmag, chi, niter, scale, iylo, iyhi, ixlo, ixhi, image_stamppk, noise_stamp, mask_stamp, psf = \
                        #    pk.pkfit_norecent_noise_smp(1, x, y, s, se, params.fitrad, returnStamps=True,
                        #                                stampsize=params.substamp)
                        #print 'here3',scale
                        # x_star, y_star = cntrd.cntrd(image_stamp, 15, 15, params.cntrd_fwhm * 2.)
                        # xpsf, ypsf = cntrd.cntrd(psf_stamp, 15, 15, params.cntrd_fwhm * 2.)
                        # psf, psfcenter = self.build_psfex(psffile, x+x_star-xpsf, y+y_star-ypsf, imfile, stop=True)
                        # pk = pkfit_norecent_noise_smp.pkfit_class(im, psf, psfcenter, self.rdnoise, self.gain,
                        #                                           noise * 0. + 1., mask)
                        # errmag, chi, niter, scale, iylo, iyhi, ixlo, ixhi, image_stamp, noise_stamp, mask_stamp, psf_stamp = \
                        #     pk.pkfit_norecent_noise_smp(1, x, y, s, se, params.fitrad, returnStamps=True,
                        #                                 stampsize=params.substamp)
                        # xpsf2, ypsf2 = cntrd.cntrd(psf_stamp, 15, 15, params.cntrd_fwhm * 2.)
                        #
                        # print x_star,y_star
                        # print xpsf,ypsf
                        # print xpsf2,ypsf2
                        # print x,y
                        #
                        # print 'rerecentroid'
                        # #raw_input('rerecentroid')

                        #image_stamp = im[np.floor(y+.5) - (params.substamp) / 2:np.floor(y+.5) + (params.substamp ) / 2 ,
                        #               np.floor(x+.5) - (params.substamp ) / 2:np.floor(x+.5) + (params.substamp ) / 2 ]

                        image_stamp = im[int(psfcenter[1]-15):int(psfcenter[1]+15),int(psfcenter[0]-15):int(psfcenter[0]+15)]
                        #gnoise_stamp = noise[psfcenter[1]-15:psfcenter[1]+15,psfcenter[0]-15:psfcenter[0]+15]
                        #noise_stamp = noise[psfcenter[1]-15:psfcenter[1]+15,psfcenter[0]-15:psfcenter[0]+15]
                        #gnoise_stamp = np.ones(image_stamp.shape)/se**2

                        stampsize = 30
                        ysn = psfcenter[0]
                        xsn = psfcenter[1]
                        if ysn - stampsize < 0:
                            ylow = int(0)
                        else:
                            ylow = int(ysn - stampsize)
                        if ysn + stampsize > im.shape[0]:
                            yhi = int(im.shape[0] - 1)
                        else:
                            yhi = int(ysn + stampsize)
                        if xsn - stampsize < 0:
                            xlow = int(0)
                        else:
                            xlow = int(xsn - stampsize)
                        if xsn + stampsize > im.shape[1]:
                            xhi = int(im.shape[1] - 1)
                        else:
                            xhi = int(xsn + stampsize)

                        #sin, sein, vals = sigma_clip.meanclip(im[ylow:yhi, xlow:xhi], clipsig=4, maxiter=8)
                        sin = np.mean(bkgrnd[ylow:yhi, xlow:xhi].ravel())
                        sein = np.mean(bkgrndrms[ylow:yhi, xlow:xhi].ravel())
                        temp, skysig, vals = sigma_clip.meanclip(im[ylow:yhi, xlow:xhi]-bkgrnd[ylow:yhi,xlow:xhi], clipsig=4, maxiter=8)

                        # pm = 100
                        # sin, sein, vals = sigma_clip.meanclip(im[psfcenter[1]-pm:psfcenter[1]+pm,psfcenter[0]-pm:psfcenter[0]+pm],
                        #                                      clipsig=4, maxiter=8)

                        #mysky, skysig, vals = sigma_clip.meanclip(im[ylow:yhi, xlow:xhi], clipsig=4, maxiter=8)

                        #skysig = np.std(vals))
                        #mysky = np.mean(vals)

                        #raw_input()
                        totalarea = len(fitrad[fitrad>0])
                        #print se**2, s/3.8,
                        #raw_input()
                        gnoise_stamp = np.ones(image_stamp.shape)*(se**2)
                        #print 'errrrr',se**2,(se**2/np.sum(psf**2))
                        #print 'sumimresid',np.sum(image_stamp-image_stamppk)

                        # ix, iy = cntrd.cntrd(image_stamp, 15, 15, 3.)
                        # px, py = cntrd.cntrd(psf, 15, 15, 3.)
                        #
                        # #tt = time.time()
                        # #for i in range(100):
                        # #    psf ,psfceter = self.build_psfex(psffile,x-ix+px,y-iy+py,imfile,stop=True)
                        #
                        # #ttt = time.time()
                        # #print 'time per psfdump is '+str((ttt-tt)/100.)
                        # #raw_input()
                        # ix = round(ix,2)
                        # iy = round(iy,2)
                        # px = round(px,2)
                        # py = round(py,2)

                        #print 'scale CHECKEEEEEE', scale, scaleck

                        #raw_input()
                        #self.tmpwriter.savefits(scale*psf_stamp+s-image_stamp,'/pnfs/des/persistent/smp/dtest.fits')
                        #if i == 9:
                        #    sys.exit()
                        # noise_stamp[noise_stamp > 0.] = 1
                        # noise_stamp[noise_stamp <= 0.] = 0
                        # self.dosextractor = False
                        # if self.dosextractor:
                        #     sexsky, sexrms = runsextractor.getsky_and_skyerr(imfile,im, ix-100, ix+100, iy-100, iy+100,snparams.survey)
                        #     print s,sexsky
                        #     print se,sexrms
                        #     #raw_input('comparison')
                        #     s = sexsky
                        #     se = sexrms
                        #     raw_input(  )
                        # else:
                        #     sexsky, sexrms = s,se

                        #usesextractorim = True
                        #if usesextractorim:
                        self.dosextractor = False
                        if self.dosextractor:
                            bkgrndstamp =  bkgrnd[int(psfcenter[1]-15):int(psfcenter[1]+15),int(psfcenter[0]-15):int(psfcenter[0]+15)]*0. + \
                                           np.mean(bkgrnd[int(psfcenter[1]-15):int(psfcenter[1]+15),int(psfcenter[0]-15):int(psfcenter[0]+15)].ravel())
                            sein =np.mean(bkgrndrms[int(psfcenter[1]-15):int(psfcenter[1]+15),int(psfcenter[0]-15):int(psfcenter[0]+15)].ravel())
                        else:
                            bkgrndstamp = im[int(psfcenter[1]-15):int(psfcenter[1]+15),int(psfcenter[0]-15):int(psfcenter[0]+15)]*0. + sin
                            #se_in = se


                        # noise_stamp = noise_stamp*1/(se**2)
                        # noise_stamp = noise_stamp * 1 / (sexrms ** 2)
                        gal = np.zeros(image_stamp.shape)
                        mjd = 000.
                        #oldcscale, cscale_std, chisq, dms = self.getfluxsmp(image_stamp, psf_stamp, s, noise_stamp, fitrad, gal,
                        #                                                    mjd, scale)
                        if self.dogalsimpixfit:
                            fiducial_coord = galsim.CelestialCoord(ra * galsim.degrees, dec * galsim.degrees)
                            stamp_center = full_data_image.wcs.posToImage(fiducial_coord)
                            cx = int(round(stamp_center.x))
                            cy = int(round(stamp_center.y))
                            des_psfex = galsim.des.DES_PSFEx(psffile)
                            thispsf = des_psfex.getPSF(stamp_center)
                            tim = full_data_image[galsim.BoundsI(cx - params.substamp/2, cx + params.substamp/2-1,
                                                                 cy - params.substamp/2, cy + params.substamp/2-1)]
                            galsimpsfworld = tim.wcs.toWorld(thispsf, image_pos=stamp_center)
                            simstamp = full_data_image[ galsim.BoundsI(cx - params.substamp/2, cx + params.substamp/2-1,
                                                                       cy - params.substamp/2, cy + params.substamp/2-1)] * 0.0
                            offset = tim.wcs.toWorld(tim.trueCenter()).project(fiducial_coord)
                            sn = galsim.Gaussian(sigma=1.e-8, flux=1.)
                            sn = sn.shift(offset)
                            conv = galsim.Convolve(sn, galsimpsfworld, gsparams=big_fft_params)
                            conv.drawImage(image=simstamp,method='no_pixel')
                            gpsf = simstamp.array
                            gscale, gscale_std, gchisq, gdms = self.getfluxsmp(image_stamp, gpsf,
                                                                               sexsky, noise_stamp,
                                                                               radius, gal, mjd, scale)
                            gsflux[i] =gscale
                            gsflux_std[i] = gscale_std
                            gsflux_chisq[i]  = gchisq
                            gsflux_dms[i] = gdms
                            #print 'gchisq',gchisq
                            #raw_input()
                        # cscale, cscale_std, chisq, dms = self.getfluxsmp(image_stamp, psf_stamp, sexsky, noise_stamp, params.fitrad,
                        #                                                  gal, mjd, scale,index=i)

                        #print 'running star',x,y
                        #try:
                            # scale, errmag, chi, dms, chinoposs = self.getfluxsmp(image_stamp, psf, sexsky, noise_stamp, fitrad, gal, mjd)
                            # oscale, oerrmag, ochi, odms, ochinoposs = self.getfluxsmp(image_stamp, psf, sexsky, onoise_stamp,
                            #                                                  fitrad, gal, mjd)
                        if True:
                            scale, errmag, chi, dms, chinoposs, bad = self.getfluxsmp(image_stamp, psf, bkgrndstamp, gnoise_stamp,
                                                                             fitrad, gal, mjd, skysig,self.gain,guess_scale=10**(.4*(31.-m)))
                            if not chi > 0.:
                                bad = True
                            #print scale
                            #raw_input('ls fit')
                        # except:
                        #
                        #     print 'could not scale'
                        #     scale, errmag, chi, dms, chinoposs = -999,-999,-999,-999,-999
                        #print 'scale and error',scale,errmag
                        #raw_input()
                        #print scale
                        # print 'DIFFFFF',scale,cscale,(scale-cscale)/scale
                        # schi = np.sum((image_stamp - psf_stamp*scale - sexsky)**2/se**2*fitrad)
                        # cchi = np.sum((image_stamp-psf*cscale-sexsky)**2*noise_stamp*fitrad)
                        # print 'pkfit chisq',schi,'fluxsmp chisq',cchi

                        #for i in range(self.Nimage):
                        if (self.savezptstamps) & (bad == False) & (chi > 0.):
                            fig = plt.figure(figsize=(20, 10))
                            axim = plt.subplot(141)
                            axpsf = plt.subplot(142)
                            axdiff = plt.subplot(143)
                            axchi = plt.subplot(144)


                            #for ax, title in zip([axim, axpsf, axdiff, axchi], ['im', 'mod', 'resid', 'chisq: '+
                            #        str(round(np.sum((image_stamp - s - (psf*scale))**2 * fitrad /se**2)/len(fitrad[fitrad>0].ravel()),2))]):

                            # for ax, title in zip([axim, axpsf, axdiff, axchi],
                            #                      ['im ' + str(ix) + ', ' + str(iy), 'mod ' + str(px) + ", " + str(py),
                            #                       'resid ' + str(round(x, 2)) + ', ' + str(round(y, 2)), 'chisq: ' +
                            #                               str(round(np.sum((image_stamp - s - (
                            #                                   psf * scale)) ** 2 * fitrad / se ** 2) / len(
                            #                                   fitrad[fitrad > 0].ravel()), 2))]):

                                # print ra
                                # if round(ra,5) == round(43.11884695,5):
                                #     print 'found culprit'
                                #     #raw_input()
                                #     ax.set_title(title+' CULPRIT')
                                # else:
                                #ax.set_title(title)
                            axs = axim.imshow(image_stamp * fitrad, cmap='gray', interpolation='nearest',vmin=min(image_stamp.ravel()),vmax=max(image_stamp.ravel()))
                            axim.set_title('Catalog Mag '+str(round(mc,2)))
                            cbar = fig.colorbar(axs, ax=axim)
                            axs = axpsf.imshow(psf/np.sum(psf) * scale * fitrad + s, cmap='gray', interpolation='nearest',vmin=min(image_stamp.ravel()),vmax=max(image_stamp.ravel()))
                            axpsf.set_title('Model')

                            cbar = fig.colorbar(axs, ax=axpsf)
                            axs = axdiff.imshow((image_stamp - s - (psf*scale)) * fitrad, cmap='gray',
                                                interpolation='nearest')
                            axdiff.set_title('Residual')

                            #axs = axdiff.imshow(psfo-psf2,cmap='gray',interpolation='nearest')
                            cbar = fig.colorbar(axs, ax=axdiff)
                            axs = axchi.imshow(
                                (image_stamp - s - (psf*scale))**2 * fitrad /se**2,
                                cmap='gray', interpolation='nearest', vmin=0, vmax=10.)
                            axchi.set_title('Chisq '+str(round(np.sum((image_stamp - s - (psf*scale))**2 * fitrad /se**2))))

                            cbar = fig.colorbar(axs, ax=axchi)
                            # plt.imshow((subim-scaledpsf)/imhdr['SKYSIG'],cmap='gray',interpolation='nearest')
                            # plt.colorbar()
                            title = 'Chisq '+str(round(np.sum((image_stamp - s - (psf*scale))**2 * fitrad /se**2),2))
                            plt.title(title)
                            pdf_pagesc.savefig(fig)

                        # psfx = np.sum(psf,axis=0)
                        # psfy = np.sum(psf,axis=1)
                        # imx = np.sum(image_stamp,axis=0)
                        # imy = np.sum(image_stamp,axis=1)
                        # plt.clf()
                        # plt.plot(np.arange(0,len(psfx)),psfx,label='psfx')
                        # plt.plot(np.arange(0,len(imx)),(imx-sexsky)/np.sum(imx-sexsky),label='imx')
                        # plt.axvline(x-np.round(x) + 15)
                        # plt.savefig('testpsfx.png')
                        # os.system('ifdh cp testpsfx.png /pnfs/des/scratch/pysmp/test/testpsfx'+str(i)+'.png')
                        # plt.clf()
                        # plt.plot(np.arange(0, len(psfy)), psfy, label='psfy')
                        # plt.plot(np.arange(0, len(imy)), (imy-sexsky)/np.sum(imy-sexsky), label='imy')
                        # plt.axvline(y-np.round(y) + 15)
                        # plt.savefig('testpsfy.png')
                        # os.system('ifdh cp testpsfy.png /pnfs/des/scratch/pysmp/test/testpsfy'+str(i)+'.png')
                        #
                        # #cscale, cscale_std, chisq, dms = self.getfluxsmp(image_stamp, psf, sexsky, noise_stamp,
                        # #                                                 params.fitrad,
                        # #                                                 gal, mjd, scale, index=i)
                        #
                        # print 'index',i,'chisq',chisq,'xpix',x,'ypix',y,'xlow',ixlo,'xhi',ixhi,'ylow',iylo,'yhi',iyhi,
                        # # if y-np.floor(y) < 5.:
                        # #     suby = 16
                        # #     addy = 14
                        # # else:
                        # #     suby = 15
                        # #     addy = 15
                        # #
                        # # if x - np.floor(x) < 5.:
                        # #     subx = 16
                        # #     addx = 14
                        # # else:
                        # #     subx = 15
                        # #     addx = 15
                        # suby = 15
                        # addy = 15
                        # subx = 15
                        # addx = 15
                        # mimage_stamp = im[np.round(y)-suby:np.round(y)+addy+1,np.round(x)-subx:np.round(x)+addx+1]
                        # mcscale, mcscale_std, mchisq, mdms = self.getfluxsmp(mimage_stamp, psf_stamp, sexsky, noise_stamp,
                        #                                                  params.fitrad,
                        #                                                  gal, mjd, scale, index=i+1000)
                        #
                        # psfx = np.sum(psf, axis=0)
                        # psfy = np.sum(psf, axis=1)
                        # imx = np.sum(mimage_stamp, axis=0)
                        # imy = np.sum(mimage_stamp, axis=1)
                        # plt.clf()
                        # plt.plot(np.arange(0, len(psfx)), psfx, label='psfx')
                        # plt.plot(np.arange(0, len(imx)), (imx - sexsky) / np.sum(imx - sexsky), label='imx')
                        # plt.axvline(x - np.round(x) + 15)
                        # plt.savefig('testpsfx.png')
                        # os.system('ifdh cp testpsfx.png /pnfs/des/scratch/pysmp/test/mtestpsfx' + str(i) + '.png')
                        # plt.clf()
                        # plt.plot(np.arange(0, len(psfy)), psfy, label='psfy')
                        # plt.plot(np.arange(0, len(imy)), (imy - sexsky) / np.sum(imy - sexsky), label='imy')
                        # plt.axvline(y - np.round(y) + 15)
                        # plt.savefig('testpsfy.png')
                        # os.system('ifdh cp testpsfy.png /pnfs/des/scratch/pysmp/test/mtestpsfy' + str(i) + '.png')
                        # #print 'checking!!!', cscale, oldcscale
                        # print 'DIFFFFFF',scale,cscale,mcscale
                        # self.tmpwriter.savefits(image_stamp-mimage_stamp,'/pnfs/des/scratch/pysmp/test/'+str(i)+'_pk-me.fits')
                        #sys.exit()
                        #scale = cscale
                        #print psfcenter,scale
                        #print 'scaled'
                        #print 'chisq',gchisq,chisq
                        #print 'flux',gscale,cscale
                        #raw_input()
                    #except ValueError:
                    #except:
                    #    print 'skipped star...\n'
                    #    continue
                #print np.median(noise_stamp)
                #mychi = np.sum((image_stamp - psf)**2 * fitrad / s**2)/len(fitrad[fitrad>0.])
                #print 'DONEEEEE',scale,errmag,chi,mychi
                if bad:
                    flux_star[i] = 1
                else:
                    flux_star[i] = scale #write file mag,magerr,pkfitmag,pkfitmagerr and makeplots
                flux_star_std[i] = errmag
                #print 'ccc',scale,errmag,oscale,oerrmag,gscale,gerrmag
                flux_chisq[i] = chi
                starsky[i] = s
                starskyerr[i] = se
                #print scale
                if not bad:
                    psfs[i,:,:] = psf
                    print 'Fit star',i,scale,errmag,s,se, chi
                    print chi
                #raw_input()
                #print flux_chisq[i]
                #raw_input()
                #flux_mychisq[i] = np.sum((image_stamp - s - (psf*scale))**2 * fitrad /se**2) / len(image_stamp.ravel())
                #flux_dms[i] = dms
                # fig = plt.figure()
                # plt.clf()
                #image_stamp[abs(image_stamp) < .1] = s
                # plt.imshow(image_stamp-sexsky-psf_stamp*scale,cmap='gray',interpolation='nearest')
                # pdf_pages.savefig(fig)
                #pdf_pages.savefig()
                #raw_input('saved teststamp.png')
                #scale = scale*.93
                # if self.savezptstamps:
                #     if not self.snparams.survey == 'PS1':
                #         simstamp = s-psf_stamp*scale
                #
                #     dt.save_fits_image(image_stamp-simstamp,os.path.join(self.zptstamps,str(mjd)+'_dms_'+str(i)+'.fits'))
                #     dt.save_fits_image(image_stamp,os.path.join(self.zptstamps,str(mjd)+'_im_'+str(i)+'.fits'))
                #     dt.save_fits_image(simstamp,os.path.join(self.zptstamps,str(mjd)+'_sim_'+str(i)+'.fits'))
                #     #dt.save_fits_image(psf_stamp,os.path.join(self.zptstamps,str(mjd)+'_psf_'+str(i)+'.fits'))
                prevra = ra
        #print 'comparing chi sqaureds'
        #for f, m in zip(flux_chisq,flux_mychisq):
        #    print f, m
        print '-'*100

        #sys.exit()

        if self.savezptstamps:
            print 'star fit stamps saved in ', self.zptstamps
            pdf_pagesc.close()
            # if not self.snparams.survey == 'PS1':
            #     #os.system('ifdh rm '+'/'.join(longimfile.split('/')[:-1]) + '/starfitstamps.pdf')
            #     self.tmpwriter.cp('starfits_' + str(thismjd) + '.pdf',
            #                       '/'.join(longimfile.split('/')[:-1]) + '/starfitstamps.pdf')
            #     self.tmpwriter.cp('starfits_' + str(thismjd) + '.pdf',self.zptstamps)
            #     os.popen('rm starfits_' + str(thismjd) + '.pdf')
            #raw_input()
            if self.fermilog:
                self.tmpwriter.appendfile('saved starfit stamps to pdf\n', self.fermilogfile)
        #self.tmpwriter.cp('mystarfits_'+str(thismjd)+'.pdf','/'.join(longimfile.split('/')[:-1])+'/mystarfitstamps.pdf')

        #raw_input('saved teststamps daophot_resid.pdf')

                #raw_input('saved teststamp.fits')
        #plt.scatter(sky[sky>10],flux_star[sky>10])
        #plt.savefig('testsky.png')
        badflag = badflag.reshape(np.shape(badflag)[0])
        
        #check for only good fits MPFIT
        if not self.dogalsimpixfit:
            goodstarcols = np.where((mag_cat != 0) &
                                (mag_cat < 26) &
                                (flux_star != 1) &
                                (flux_star != -999) &
                                (flux_star < 9e6) &
                                #(flux_chisq < 1.5) &
                                #(flux_chisq > 0) &
                                #(flux_star_mcmc < 1e7) &
                                #(flux_star_mcmc != 0) &
                                #(flux_star_mcmc_modelerrors != 0) &
                                #(flux_star_mcmc_modelerrors < 1e7) &
                                #(flux_star_std_mcmc > 1.0) &
                                #(flux_star_std_mcmc_modelerrors > 1.0) &
                                (np.isfinite(mag_cat)) &
                                (np.isfinite(flux_star)) &
                                (flux_chisq < 500.) &
                                (np.isfinite(flux_chisq)) &
                                (~np.isnan(flux_chisq)) &
                                (flux_star > 0) &
                                (badflag == 0) &
                                (isnotcheckstars == 1))[0]
        else:
            goodstarcols = np.where((mag_cat != 0) &
                                    (mag_cat < 26) &
                                    (gsflux != 1) &
                                    (gsflux < 9e6) &
                                    #(flux_chisq < 1.5) &
                                    #(flux_chisq > 0) &
                                    # (flux_star_mcmc < 1e7) &
                                    # (flux_star_mcmc != 0) &
                                    # (flux_star_mcmc_modelerrors != 0) &
                                    # (flux_star_mcmc_modelerrors < 1e7) &
                                    # (flux_star_std_mcmc > 1.0) &
                                    # (flux_star_std_mcmc_modelerrors > 1.0) &
                                    (np.isfinite(mag_cat)) &
                                    (np.isfinite(flux_star)) &
                                    (flux_chisq < 500.) &
                                    (np.isfinite(flux_chisq)) &
                                    (~np.isnan(flux_chisq)) &
                                    (flux_star > 0) &
                                    (badflag == 0) &
                                    (isnotcheckstars == 1))[0]


        checkstarcols = np.where((mag_cat != 0) &
                                (mag_cat < 26) &
                                (flux_star != 1) &
                                (flux_star < 9e6) &
                                #(flux_chisq < 1.5) &
                                #(flux_chisq > 0) &
                                (np.isfinite(mag_cat)) &
                                (np.isfinite(flux_star)) &
                                (flux_star > 0) &
                                (badflag == 0) &
                                (isnotcheckstars == 0))[0]


        #NEED TO MAKE A PLOT HERE!
        print goodstarcols
        print '-'*100

        if len(goodstarcols) > self.params.minzptstars:

            for fl,fc in zip(mag_cat[goodstarcols],flux_chisq[goodstarcols]):
                print fl,fc
            print '-' * 100

            if not self.dogalsimpixfit:
                fluxcol = flux_star
                #print 'fluxcol  =flux_star'
            else:
                fluxcol = gsflux


            md,std,num = self.iterstat(mag_cat[goodstarcols]+2.5*np.log10(fluxcol[goodstarcols]),
                                       startMedian=True,sigmaclip=3,iter=10)

            ww2 = (mag_cat[goodstarcols] < 19.) & (mag_cat[goodstarcols] > 16.5)
            md2,zptscat,num2 = self.iterstat(mag_cat[goodstarcols][ww2]+2.5*np.log10(fluxcol[goodstarcols][ww2]),
                                       startMedian=True,sigmaclip=3,iter=10)
            print ''
            print ''
            print '-'*100
            print '-'*38,'Done Fitting Zeropoint','-'*38
            print '-'*100
            print 'fitzpt',md,'diffimzpt',snparams.zp[j]
            print 'std',std,'std/sqrt(n)',std/num**.5
            print '-'*100
            print ''
            print ''

            #std = float(std)/float(num**.5)

            med, rmsaddin, num = self.iterstat(float(md) - mag_cat[goodstarcols] - 2.5 * np.log10(fluxcol[goodstarcols]),
                                        startMedian=True, sigmaclip=3, iter=10)
            zptresid = float(md) - mag_cat - 2.5 * np.log10(fluxcol)

            goodstarcols = np.where((mag_cat != 0) &
                                    (mag_cat < 26) &
                                    (gsflux != 1) &
                                    (gsflux < 9e6) &
                                    # (flux_chisq < 1.5) &
                                    # (flux_chisq > 0) &
                                    # (flux_star_mcmc < 1e7) &
                                    # (flux_star_mcmc != 0) &
                                    # (flux_star_mcmc_modelerrors != 0) &
                                    # (flux_star_mcmc_modelerrors < 1e7) &
                                    # (flux_star_std_mcmc > 1.0) &
                                    # (flux_star_std_mcmc_modelerrors > 1.0) &
                                    (np.isfinite(mag_cat)) &
                                    (np.isfinite(flux_star)) &
                                    (abs(zptresid) < 3.*rmsaddin ) &
                                    (flux_star > 0) &
                                    (badflag == 0) &
                                    (isnotcheckstars == 1))[0]
            bad = False
            if self.fermilog:
                self.tmpwriter.appendfile('fitzpt '+str(md)+' +-'+str(std)+' diffimzpt '+str(snparams.zp[j])+'\n', self.fermilogfile)

            #sys.exit()
            #dstd = 1.48*np.median(abs(mag_cat[goodstarcols]+2.5*np.log10(flux_star[goodstarcols])- np.ones(len(flux_star[goodstarcols]))*md))/np.sqrt(len(flux_star[goodstarcols]))
            #print 'reduced std', std
            #print 'dan std',dstd
            #mcmc_md = -999.
            #mcmc_std = -999.

            #mcmc_me_md,mcmc_me_std = self.weighted_avg_and_std(mag_cat[goodstarcols]+2.5*np.log10(flux_star_mcmc_modelerrors[goodstarcols]),1.0/(mcmc_me_mag_std[goodstarcols])**2)

            #zpt_plots_out = mag_compare_out = imfile.split('.')[-2] + '_zptPlots'
            #exposure_num = imfile.split('/')[-1].split('_')[1]
            #print 'writing zeropoints'
            # if nozpt:
            #     fn = self.big_zpt+'.txt'
            #     if os.path.isfile(self.big_zpt+'.txt'):
            #         pass
            #     else:
            #        self.tmpwriter.writefile('Exposure Num\tRA\tDEC\tCat Zpt\tMPFIT Zpt\tMPFIT Zpt Err\tMCMC Zpt\tMCMC Zpt Err\tMCMC Model Errors Zpt\tMCMC Model Errors Zpt Err\tCat Mag\tMP Fit Mag\tMCMC Fit Mag\tMCMC Model Errors Fit Mag\tMCMC Analytical Simple\tMCMC Analytical Weighted\n',fn)
            #     for i in goodstarcols:
            #         self.tmpwriter.appendfile(str(exposure_num)+'\t'+str(ras[i])+'\t'+str(decs[i])+'\t'+str(cat_zpt)+'\t'+str(md)+'\t'+str(std)\
            #             +'\t'+str(mcmc_md)+'\t'+str(mcmc_std)+'\t'+str(mcmc_me_md)+'\t'+str(mcmc_me_std)+'\t'+str(mag_cat[i])\
            #             +'\t'+str(-2.5*np.log10(flux_star[i]))+'\t'+str(-2.5*np.log10(flux_star_mcmc[i]))\
            #             +'\t'+str(-2.5*np.log10(flux_star_mcmc_modelerrors[i]))\
            #             +'\t'+str(-2.5*np.log10(flux_star_mcmc_me_simple[i]))
            #             +'\t'+str(-2.5*np.log10(flux_star_mcmc_me_weighted[i]))
            #             +'\n',fn)
            #     #b.close()
            #
            #     #b = open(self.checkstarfile,'a')
            #     for i in goodstarcols:
            #         self.tmpwriter.appendfile(str(exposure_num)+'\t'+str(thismjd)+'\t'+str(ras[i])+'\t'+str(decs[i])+'\t'+str(xstar[i])
            #                 +'\t'+str(ystar[i])+'\t'+str(cat_zpt)+'\t'+str(md)+'\t'+str(std)+'\t'
            #                 +str(flux_star[i])+'\t'+str(flux_star_std[i])+'\t'+str(flux_chisq[i])+'\t'+str(flux_dms[i])+'\t'
            #                 + str(gsflux[i]) + '\t' + str(gsflux_std[i]) + '\t' + str(gsflux_chisq[i]) + '\t' + str(gsflux_dms[i])
            #                 +'\t'+str(mag_cat[i])+'\n',self.checkstarfile)

                #b.close()
                #print 'checkstarfilea appended'
                #print self.checkstarfile
                #raw_input()
            #hh = mag_cat[goodstarcols]+2.5*np.log10(flux_star[goodstarcols]) - np.ones(len(flux_star[goodstarcols]))*md
            #hh = hh[abs(hh < .25)]
            #print 'plotting zeropoints'
            #
            #plt.clf()
            # plt.hist(mag_cat[goodstarcols]+2.5*np.log10(flux_star[goodstarcols]) - np.ones(len(flux_star[goodstarcols]))*md,bins=np.arange(-.25,.25,.04),label='mean: '+str(np.mean(hh))+' std: '+str(np.std(hh)))
            # plt.xlabel('cat mag + 2.5log10(flux) - zeropoint')
            # plt.ylabel('counts')
            # plt.xlim(-.25,.25)
            # if self.fermigrid:
            #     if self.worker:
            #         if not os.path.exists('./zpts/'):
            #             os.makedirs('./zpts/')
            #         print os.path.join('./zpts/', imfile.split('.fits')[-2].split('/')[-1] + '_'+str(filt)+'band_starfitresids1s.png')
            #         plt.savefig(os.path.join('./zpts/', imfile.split('.fits')[-2].split('/')[-1] + '_'+str(filt)+'band_starfitresids1s.png'))
            #         os.system('ifdh cp -D '+os.path.join('./zpts/', imfile.split('.fits')[-1].split('/')[-1] + '_'+str(filt)+'band_starfitresids1s.png')
            #                   + ' '+self.zptoutpath)
            #     else:
            #         print imfile.split('.fits')
            #         print os.path.join(self.zptoutpath,imfile.split('.fits')[-2].split('/')[-1] + '_'+str(filt)+'band_starfitresids1s.png')
            #         plt.savefig('tmp.png')
            #         if os.path.exists(os.path.join(self.zptoutpath,imfile.split('.fits')[-2].split('/')[-1] + '_'+str(filt)+'band_starfitresids1s.png')):
            #             os.system('yes | rm '+os.path.join(self.zptoutpath,imfile.split('.fits')[-2].split('/')[-1] + '_'+str(filt)+'band_starfitresids1s.png'))
            #         os.system('mv tmp.png '+os.path.join(self.zptoutpath,imfile.split('.fits')[-2].split('/')[-1] + '_'+str(filt)+'band_starfitresids1s.png'))
            # else:
            #     plt.savefig(os.path.join(self.zptoutpath,imfile.split('.fits')[-2].split('/')[-1] + '_'+str(filt)+'band_starfitresids1s.png'))
            #     #plt.savefig(os.path.join(self.zptoutpath,imfile.split('.fits')[-1].split('/')[-1] + '_'+str(filt)+'band_starfitresids1s.png'))
            # #print imfile.split('.')[-2].split('/')[-1] + '_'+str(filt)+'band_starfitresids1s.png'
            # #plt.savefig(imfile.split('.')[-2] + '_'+str(filt)+'band_starfitresids1s.png')
            # #print imfile.split('.')[-2] + '_'+str(filt)+'band_starfitresids1s.png'
            # #r.write(imfile.split('.')[-2] + '_'+str(filt)+'band_starfitresids1s.png\n')
            # #r.close()

            #plt.clf()
            #print 'scatter'
            #print len(mag_cat[goodstarcols])


            doplot = False
            if doplot: plt.clf()
            #plt.scatter(mag_cat[goodstarcols], md-mag_cat[goodstarcols]-2.5*np.log10(flux_star[goodstarcols]))


            from lmfit import Minimizer, Parameters

            def f(p):
                return (p['m']-mag_cat[goodstarcols]-2.5*np.log10(flux_star[goodstarcols]))/(flux_star_std[goodstarcols]/flux_star[goodstarcols])

            lmfparams = Parameters()
            lmfparams.add('m', value=31., min=20.)
            # params.add('pow', value=.5, vary=False)

            fitter = Minimizer(f, lmfparams)
            v = fitter.minimize(method='leastsq')
            #print fluxls, v.params['scale'].value
            #print v.params['m'].__dict__
            mde = v.params['m'].value
            mdeerr = v.params['m'].stderr

            med, rmsaddin, num = self.iterstat(
                float(mde) - mag_cat[goodstarcols] - 2.5 * np.log10(fluxcol[goodstarcols]),
                startMedian=True, sigmaclip=3, iter=10)
            # mde, std, num = self.iterstat(mag_cat[goodstarcols] + 2.5 * np.log10(fluxcol[goodstarcols]),
            #                              startMedian=True, sigmaclip=3, iter=10)


            print '-'*100
            print 'Fit ZPT:',md,'+-',std/np.sqrt(num)
            print 'W Fitzp:',mde,'+-',mdeerr
            print '-'*100

            #doplot = False
            if doplot:
                plt.errorbar(mag_cat[goodstarcols], mde-mag_cat[goodstarcols]-2.5*np.log10(flux_star[goodstarcols]),
                             flux_star_std[goodstarcols]/flux_star[goodstarcols],fmt='o',label='ZPT: '+str(round(mde,3))+' +- '+str(round(mdeerr,3)))
                #print 'plot'
                #plt.plot([min(mag_cat[goodstarcols]),max(mag_cat[goodstarcols])],[min(mag_cat[goodstarcols]),max(mag_cat[goodstarcols])]-md,color='black')
                plt.axhline(0,color='black')
                #plt.plot([min(mag_cat[goodstarcols]),max(mag_cat[goodstarcols])],[min(mag_cat[goodstarcols])-30.734,max(mag_cat[goodstarcols])-30.734],color='red')

                plt.xlabel('cat mag')
                plt.ylabel('zpt-2.5log10(flux)')

                plt.legend()
            #print 'saving'
            #print mag_cat[goodstarcols].shape
            #plt.savefig(imfile.split('.')[-2] + '_'+str(filt)+'band_starfit_zptplot.png')
            #print imfile.split('.')[-2] + '_'+str(filt)+'band_starfit_zptplot.png'
            if self.fermigrid:
                print 'insidefermigrid'

                #print os.path.join(longimfile.split('.fits')[-2].split('/')[-1] + '_'+str(filt)+'band_starfit_zptplot.png')
                plt.savefig(os.path.join(imfile.split('.fits')[-2].split('/')[-1] + '_'+str(filt)+'band_starfit_zptplot.png'))
                os.popen('ifdh rm ' + os.path.join(self.zptoutpath,imfile.split('.fits')[-2].split('/')[-1] + '_' + str(
                    filt) + 'band_starfit_zptplot.png'))
                os.system('ifdh mv ' + os.path.join(imfile.split('.fits')[-2].split('/')[-1] + '_' + str(
                    filt) + 'band_starfit_zptplot.png')
                          + ' ' + os.path.join(self.zptoutpath,imfile.split('.fits')[-2].split('/')[-1] + '_' + str(
                    filt) + 'band_starfit_zptplot.png'))
                os.popen('rm '+os.path.join(imfile.split('.fits')[-2].split('/')[-1] + '_'+str(filt)+'band_starfit_zptplot.png'))
                if self.fermilog:
                    self.tmpwriter.appendfile(
                        os.path.join(imfile.split('.fits')[-2].split('/')[-1] + '_'+str(filt)+'band_starfit_zptplot.png')+'\n',
                        self.fermilogfile)
                # plt.clf()
                # ras = np.array(ras)
                # decs = np.array(decs)
                # plt.scatter(ras[goodstarcols], -2.5 * np.log10(flux_star[goodstarcols]) - mag_cat[goodstarcols] + md)
                # plt.savefig(os.path.join(imfile.split('.fits')[-2].split('/')[-1] + '_' + str(
                #     filt) + 'band_starfit_zptplot_ra.png'))
                # print 'saved', imfile.split('.fits')[-2].split('/')[-1] + '_' + str(
                #     filt) + 'band_starfit_zptplot_ra.png')
                # plt.clf()
                # plt.scatter(decs[goodstarcols], -2.5 * np.log10(flux_star[goodstarcols]) - mag_cat[goodstarcols] + md)
                # plt.savefig(os.path.join('./zpts', imfile.split('.fits')[-2].split('/')[-1] + '_' + str(
                #     filt) + 'band_starfit_zptplot_dec.png'))
                # print 'saved', os.path.join('./zpts', imfile.split('.fits')[-2].split('/')[-1] + '_' + str(
                #     filt) + 'band_starfit_zptplot_dec.png')


            else:
                try:
                    os.mkdir(os.path.join(self.outdir,'stardata'))
                except:
                    pass
                try:
                    os.mkdir(os.path.join(self.outdir,'stardata',filt))
                except:
                    pass
                #print imfile
                #raw_input()
                name = imfile.split('/')[-1][:-8]
                zptplotout = os.path.join(self.outdir,'stardata',filt, name + '_zptplot.png')
                if doplot:
                    plt.savefig(zptplotout)
                    print 'saved',zptplotout
                    plt.clf()
                ras = np.array(ras)
                decs = np.array(decs)
                ids = np.array(ids)
                if doplot:
                    plt.scatter(ras[goodstarcols], -2.5 * np.log10(flux_star[goodstarcols]) - mag_cat[goodstarcols] + md )
                zptplotoutra = os.path.join(self.outdir, 'stardata', filt, name + '_zptplot_ra.png')
                if doplot:
                    plt.savefig(zptplotoutra)
                #print 'saved', os.path.join(self.zptoutpath, imfile.split('.fits')[-2].split('/')[-1] + '_' + str(
                #    filt) + 'band_starfit_zptplot_ra.png')
                if doplot:
                    plt.clf()
                    plt.scatter(decs[goodstarcols], -2.5 * np.log10(flux_star[goodstarcols]) - mag_cat[goodstarcols] + md )
                zptplotoutdec = os.path.join(self.outdir, 'stardata', filt, name + '_zptplot_dec.png')

                if doplot:
                    plt.savefig(zptplotoutdec)
                #print 'saved', os.path.join(self.zptoutpath, imfile.split('.fits')[-2].split('/')[-1] + '_' + str(
                #    filt) + 'band_starfit_zptplot_dec.png')

            rrr = -2.5 * np.log10(flux_star[goodstarcols]) - mag_cat[goodstarcols] + md
            badguys = abs(rrr) > .35

            # fff = open('badguys.reg','w')
            # for x,y in zip(xstar[goodstarcols][badguys],ystar[goodstarcols][badguys]):
            #     fff.write('circle '+str(x)+' '+str(y)+' 3\n')
            # fff.close()
            #print 'wrote badguys.reg'
            #print imfile
            #raw_input()
            #print 'saved properly'
            #raw_input()
            '''print 'mean python', np.mean(hh)
            print 'mean idl ',np.mean(istarmags-istarcats+30.6198)
            print 'std python', np.std(hh)
            print 'std idl ',np.std(istarmags-istarcats+30.6198)
            '''
            #raw_input()
            #Writing mags out to file .zpt in same location as image
            #print 'saving npz'
            name = imfile.split('/')[-1].split('.')[-2]
            if doglobalstar:
                if self.dogalsimpixfit:
                    print name
                    print self.impath
                    mag_compare_out = os.path.join(self.impath,name + '_' + str(filt) + 'band_dillonzptinfo_galsimglobalstar.npz')
                    print mag_compare_out
                else:
                    #mag_compare_out = os.path.join(self.zptoutpath,name + '_'+str(filt)+'band_dillonzptinfo_globalstar.npz')
                    mag_compare_out = os.path.join(self.impath,name + '_'+str(filt)+'band_dillonzptinfo_globalstar.npz')
            else:
                if self.dogalsimpixfit:
                    mag_compare_out = os.path.join(self.impath,name + '_' + str(filt) + 'band_dillonzptinfo_galsim.npz')
                else:
                    mag_compare_out = os.path.join(self.impath,name + '_'+str(filt)+'band_dillonzptinfo.npz')
            #print goodstarcols
            if self.fermilog:
                self.tmpwriter.appendfile(
                    'about to write zpt tempfile\n', self.fermilogfile)

            ras = np.array(ras)
            decs = np.array(decs)
            print 'magerr', -2.5*np.log10(fluxcol[goodstarcols])+2.5*np.log10(fluxcol[goodstarcols]+flux_star_std[goodstarcols])
            if self.fermigrid and self.worker:
                #print fluxcol[goodstarcols]
                #raw_input('testing123' )
                ff = 'temp.npz'
                np.savez( ff
                    #,ra = ras[goodstarcols]
                    #,dec = decs[goodstarcols]
                    ,cat_mag = mag_cat[goodstarcols]
                    ,fit_mag = -2.5*np.log10(fluxcol[goodstarcols])
                    ,fit_mag_err = -2.5*np.log10(fluxcol[goodstarcols])+2.5*np.log10(fluxcol[goodstarcols]+flux_star_std[goodstarcols])
                    ,flux_starmp = fluxcol[goodstarcols]
                    ,flux_star_std = flux_star_std[goodstarcols]
                    ,chisqu=flux_chisq[goodstarcols]
                    #,mcmc_me_fit_mag = -2.5*np.log10(flux_star_mcmc_modelerrors[goodstarcols])
                    #,mcmc_me_fit_mag_std = mcmc_me_mag_std[goodstarcols]
                    , ras=ras[goodstarcols]
                    , decs=decs[goodstarcols]
                    , centroidedras = thisra[goodstarcols]
                    , centroideddecs = thisdec[goodstarcols]
                    , ids = ids[goodstarcols]
                    , rmsaddin = rmsaddin
                    ,fit_zpt = md
                    ,wfit_zpt = med
                    ,wfit_zpt_std=rmsaddin
                    ,fit_zpt_std = std
                    ,sky = starsky[goodstarcols]
                    ,skyerr = starskyerr[goodstarcols]
                    #,mcmc_me_zpt = mcmc_me_md
                    #,mcmc_me_zpt_std = mcmc_me_std
                    ,cat_zpt = cat_zpt
                    ,mjd = thismjd
                    ,psfs = psfs[goodstarcols,:,:]
                    ,mjdoff=mjdoff
                    ,mjdslopeinteroff=mjdslopeinteroff
                    )
                if self.fermilog:
                    self.tmpwriter.appendfile(
                        'about to copy to mag_compare out\n', self.fermilogfile)
                tt = time.time()
                #print 'about to copy mag_compare_out'
                os.popen('ifdh rm '+mag_compare_out)
                print os.popen('ifdh mv '+ff+' '+mag_compare_out).read()


                outtext = []
                for i in range(len(ras[goodstarcols])):
                    outtext.append(str(ids[goodstarcols][i])+'\t'+str(ras[goodstarcols][i])+'\t'+str(decs[goodstarcols][i])+'\t'+str(mag_cat[goodstarcols][i])+'\t'+filt+'\t'+self.expnum+'\n')

                self.tmpwriter.appendfile(outtext, self.laskerstarcat)


                #ttt = time.time()
                #self.tmpwriter.appendfile(
                #    'ifdh took '+str(ttt-tt)+'seconds\n', self.fermilogfile)
                #sys.exit()
                print mag_compare_out
                #raw_input('saved')
            else:
                try:
                    os.mkdir(os.path.join(self.outdir,'stardata'))
                except:
                    pass
                try:
                    os.mkdir(os.path.join(self.outdir,'stardata',filt))
                except:
                    pass
                #print imfile
                #raw_input()
                name = imfile.split('/')[-1][:-5]
                mag_compare_out = os.path.join(self.impath,
                                               name + '_' + str(filt) + 'band_dillonzptinfo_globalstar.npz').replace('+fakeSN','')
                #print 'lengoodstarcols',len(mag_cat[goodstarcols])
                #raw_input()

                fwhm = open(psffile,'r').readline().split('PSF_FWHM')[1].split('/')[0].split('=')[1]
                print 'FWHM:',fwhm

                np.savez(mag_compare_out
                         , cat_magsmp=mag_cat[goodstarcols]
                         , fit_magmy=-2.5 * np.log10(fluxcol[goodstarcols])
                         , fit_mag_err=-2.5 * np.log10(fluxcol[goodstarcols]) + 2.5 * np.log10(fluxcol[goodstarcols] + flux_star_std[goodstarcols])
                         , flux_starmp=fluxcol[goodstarcols]
                         , flux_star_stdmp=flux_star_std[goodstarcols]
                         , chisqu=flux_chisq[goodstarcols]
                         , ras=ras[goodstarcols]
                         , decs=decs[goodstarcols]
                         , centroidedras=thisra[goodstarcols]
                         , centroideddecs=thisdec[goodstarcols]
                         , ids=ids[goodstarcols]
                         , rmsaddin=rmsaddin
                         , fit_zpt=md
                         , fit_zpt_std=std/np.sqrt(num)
                         , wfit_zpt=mde
                         , wfit_zpt_std=mdeerr
                         , sky=starsky[goodstarcols]
                         , skyerr=starskyerr[goodstarcols]
                         , cat_zpt=cat_zpt
                         , mjd=thismjd
                         , psfs=psfs[goodstarcols, :, :]
                         , mjdoff=mjdoff
                         , mjdslopeinteroff=mjdslopeinteroff
                         , thisra = thisra[goodstarcols]
                         , thisdec = thisdec[goodstarcols]
                         , thisids = thisids[goodstarcols]
                         , fwhm = np.float(fwhm)
                         , zptscat = zptscat
                         )
                #raw_input()


                # self.tmpwriter.savez(mag_compare_out
                #                      # ,ra = ras[goodstarcols]
                #                      # ,dec = decs[goodstarcols]
                #                      , cat_mag=mag_cat[goodstarcols]
                #                      , fit_mag=-2.5 * np.log10(fluxcol[goodstarcols])
                #                      , fitflux = fluxcol[goodstarcols]
                #                      # ,mcmc_me_fit_mag = -2.5*np.log10(flux_star_mcmc_modelerrors[goodstarcols])
                #                      # ,mcmc_me_fit_mag_std = mcmc_me_mag_std[goodstarcols]
                #                      , ras = ras[goodstarcols]
                #                      , decs = decs[goodstarcols]
                #                      , fit_zpt=md
                #                      , fit_zpt_std=std
                #                      , sky=starsky[goodstarcols]
                #                      , skyerr=starskyerr[goodstarcols]
                #                      # ,mcmc_me_zpt = mcmc_me_md
                #                      # ,mcmc_me_zpt_std = mcmc_me_std
                #                      , cat_zpt=cat_zpt
                #                      , rmsaddin=rmsaddin
                #                      , mjd=mjd
                #                      , mjdoff=mjdoff
                #                      , mjdslopeinteroff=mjdslopeinteroff
                #                      )
            bad = False
            #raw_input('ZEROPOINTING WAS GOOD')
        else:
            print len(goodstarcols)
            print len(checkstarcols)
            print isnotcheckstars
            print params.numcheckstars
            md = 0
            std = 0
            mag_compare_out = 0
            #raw_input('Error : not enough good stars to compute zeropoint!!!')
            bad = True
            #raise exceptions.RuntimeError('Error : not enough good stars to compute zeropoint!!!')

            name = imfile.split('/')[-1][:-5]
            mag_compare_out = os.path.join(self.impath,
                                           name + '_' + str(filt) + 'band_dillonzptinfo_globalstar.npz').replace(
                                           '+fakeSN', '')

            np.savez(mag_compare_out
                     , cat_magsmp=0
                     , fit_mag=-0
                     , fit_mag_err=0
                     , flux_starh=0
                     , flux_star_std=0
                     , chisqu=0
                     , ras=0
                     , decs=0
                     , centroidedras=0
                     , centroideddecs=0
                     , ids=0
                     , rmsaddin=0
                     , fit_zpt=0
                     , fit_zpt_std=0
                     , sky=0
                     , skyerr=0
                     , cat_zpt=0
                     , mjd=0
                     , psfs=0
                     , mjdoff=0
                     , mjdslopeinteroff=0
                     , thisra=0
                     , thisdec=0
                     , thisids=0
                     , fwhm=0
                     , zptscat = 0
                     )


            print 'Error : not enough good stars to compute zeropoint!!!'*20

        if not bad:
            if std > 0.05:
                #print rmsaddin
                #print 'rmsaddin large'*100.
                bad = True
                md = 0
                std = 0
                mag_compare_out = 0

        if not bad:
            if self.fermilog:
                self.tmpwriter.appendfile(
                    'saved zpt file '+mag_compare_out+'\n', self.fermilogfile)
            #self.tmpwriter.appendfile(
            #    os.popen('ls -ltr ').read(),self.fermilogfile
            #)
            #self.tmpwriter.appendfile(
            #    os.popen('ls -ltr zpts/').read(),self.fermilogfile
            #)

        #if self.verbose:
        print('measured ZPT: %.3f +/- %.3f'%(md,std))
        #sys.exit()

        #if self.fermigrid:
        #    os.system('ifdh cp -D -r ./zpts/ ' + self.zptoutpath)
        #    print 'copied from worker to zpt path',self.zptoutpath
        #sys.exit()
        #raw_input('stopped')
        if bad:
            return 0,0,0,0,0,0,0
        return(mde,mdeerr,mag_compare_out,rmsaddin,thisra[goodstarcols],thisdec[goodstarcols],thisids[goodstarcols])

    def get_fwhm_of_2d_psf(self,psfstamp):

        oned_psf = np.sum(psfstamp,axis=0)
        oned_psft = np.sum(psfstamp,axis=1)

        x = np.arange(len(oned_psf))*self.snparams.platescale
        spline = UnivariateSpline(x, oned_psf-np.max(oned_psf)/2, s=0)
        r1, r2 = spline.roots()

        fwhm1 = r2-r1

        #plt.clf()
        #plt.plot(x,oned_psf,color='b')
        #plt.axvspan(r1, r2, facecolor='b', alpha=0.4)

        x = np.arange(len(oned_psft))*self.snparams.platescale
        spline = UnivariateSpline(x, oned_psft-np.max(oned_psft)/2, s=0)
        r1, r2 = spline.roots()

        #plt.plot(np.arange(len(oned_psft))*.27,oned_psft,color='r')
        #plt.axvspan(r1, r2, facecolor='r', alpha=0.4)
        #plt.xlabel('Arcsec')
        #plt.ylabel('PSF Value')

        fwhm2 = r2-r1

        #plt.savefig('testpsf.png')


        return (fwhm1+fwhm2)/2.

    def sector_mask(self,shape,centre,radius,angle_range=(0,360)):
        """
        Return a boolean mask for a circular sector. The start/stop angles in  
        `angle_range` should be given in clockwise order.
        """

        x,y = np.ogrid[:shape[0],:shape[1]]
        cx,cy = centre
        tmin,tmax = np.deg2rad(angle_range)

        # ensure stop angle > start angle
        if tmax < tmin:
                tmax += 2*np.pi

        # convert cartesian --> polar coordinates
        r2 = (x-cx)*(x-cx) + (y-cy)*(y-cy)
        theta = np.arctan2(x-cx,y-cy) - tmin

        # wrap angles between 0 and 2*pi
        theta %= (2*np.pi)

        # circular mask
        circmask = r2 <= radius*radius

        # angular mask
        anglemask = theta <= (tmax-tmin)

        return circmask*anglemask

    def degrade_psf_to_fake(self,smp_psfi,fake_psfi):
        smp_fwhm = self.get_fwhm_of_2d_psf(smp_psfi)
        sigma = .5
        degraded_smp_psfi_fwhm = 0.
        degraded_smp_psfi = copy(smp_psfi)
        while fake_psfi-degraded_smp_psfi_fwhm > .05:
            degraded_smp_psfi = scipy.ndimage.filters.gaussian_filter(degraded_smp_psfi, sigma)
            degraded_smp_psfi_fwhm = self.get_fwhm_of_2d_psf(degraded_smp_psfi)
            #print abs(fake_psfi-degraded_smp_psfi_fwhm)
            #raw_input()
        return degraded_smp_psfi

    def create_model(self,image,psf,scales):
        psfmax = max(psf.ravel())
        galaxy_peak_pix = np.argmax(image.ravel())
        galaxy_peak_val = np.max(image.ravel())

        smallmodel = copy(image.ravel())/2.
        #smallmodel[galaxy_peak_pix] = float(galaxy_peak_val)/(float(psfmax))
        model = np.append(smallmodel,scales)
        return model,galaxy_peak_pix

    def weighted_avg_and_std(self,values,weights):

        """
        Return the weighted average and standard deviation.
        values, weights -- Numpy ndarrays with the same shape.
        """
        average = np.average(values, weights=weights)
        variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
        return (average, variance**.5)

    def build_psfex(self, psffile,x,y,imfile,dogalsim=False,stop=False,psfexworked=True):
        if psfexworked:

            a = psfex.PSFEx(psffile)
            im = a.get_rec(y, x)[3:-4, 3:-4]
            im /= np.sum(im.ravel())

            return im, (round(x), round(y))
        else:
            stampsize=30
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
            psfout = np.zeros((stampsize, stampsize))
            for x, y, p in zip(ix, iy, psfval):
                psfout[y, x] = p

            return (psfout), (IMAGE_CENTERX, IMAGE_CENTERY)

    def build_psfex_old(self, psffile,x,y,imfile,dogalsim=False,stop=False):
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
        ### psf = os.popen("dump_psfex -inFile_psf %s -xpix %s -ypix %s -gridSize %s"%(psffile,x,y,
        ###                                                                           self.params.substamp)).read()
        #print "bash -c source dump_psfex.c -inFile_psf %s -xpix %s -ypix %s -gridSize %s"%(psffile,x,y,35)
        #print 'inside build psfex'
        #print 'xpix',x,'ypix',y
        xo = copy(x)
        yo = copy(y)
        psf = os.popen("dump_psfex -inFile_psf %s -xpix %s -ypix %s -gridSize %s"%(psffile,x,y,
                                                                                   30)).readlines()

        #ix, iy, psfval = np.genfromtxt(psffile, usecols = (1,2,5), skip_footer = 4)
        xin = copy(x)
        yin = copy(y)
        readdata,readheader = False,True
        ix,iy,psfval = [],[],[]
        IMAGE_CORNERX = 0
        IMAGE_CORNERY = 0
        for line in psf:
            #print line
            line = line.replace('\n','')
            if line.startswith('PSF:'):
                #linelist = filter(None,line.split(' '))
                linelist = line.split()
                ix += [int(linelist[1])]; iy += [int(linelist[2])]; psfval += [float(linelist[5])]

            elif line.startswith("IMAGE_CENTER"):
                linelist = line.split()
                IMAGE_CENTERX = float(linelist[1]); IMAGE_CENTERY = float(linelist[2])
            #elif line.startswith("IMAGE_CORNER"):
            #    linelist = line.split()
            #    IMAGE_CORNERX = float(linelist[1]); IMAGE_CORNERY = float(linelist[2])

        #IMAGE_CENTERX -= IMAGE_CORNERX; IMAGE_CENTERY -= IMAGE_CORNERY
        ix,iy,psfval = np.array(ix),np.array(iy),np.array(psfval)
        #psfout = np.zeros((2*self.params.fitrad + 1,2*self.params.fitrad + 1))
        psfout = np.zeros((self.params.substamp,self.params.substamp))
        for x,y,p in zip(ix,iy,psfval):
            #if x >= (35 - 2*self.params.fitrad -1)/2 and y >= (35 - 2*self.params.fitrad -1)/2 and x < (2*self.params.fitrad +1) and y < (2*self.params.fitrad + 1):
            #    psfout[y-(35 - 2*self.params.fitrad - 1)/2 ,x-(35 - 2*self.params.fitrad -1)/2 ] = p
            psfout[y,x] = p

        # print 'psfvalmax',np.max(psfval)
        # print 'psfshape',psfout.shape
        # print 'psfmax',np.max(psfout)
        # print 'psffile',psffile
        # print 'imfile',imfile
        # self.tmpwriter.savefits(psfout,'/pnfs/des/scratch/pysmp/test/aaapsf.fits')
        # imstamp = pf.getdata(imfile)
        # print 'imstamp shape',imstamp.shape
        # print xo,yo
        # xo = np.floor(xo)
        # yo = np.floor(yo)
        # imstamp = imstamp[yo-17:yo+17+1,xo-17-1:xo+17]
        # imstamp = imstamp/np.sum(imstamp)
        # print 'imstamp shape',imstamp.shape
        # self.tmpwriter.savefits(imstamp, '/pnfs/des/scratch/pysmp/test/aaaim.fits')
        # self.tmpwriter.savefits(imstamp-psfout, '/pnfs/des/scratch/pysmp/test/aaadms.fits')

        # imstamp = pf.getdata(imfile)
        # imstamp = imstamp[xo - 35/2 + 1:xo + 35/2, yo - 35/2 + 1:yo + 35/2]
        # imstamp = imstamp / np.sum(imstamp)
        #
        # self.tmpwriter.savefits(imstamp, '/pnfs/des/scratch/pysmp/test/aaaim2.fits')
        # self.tmpwriter.savefits(imstamp - psfout, '/pnfs/des/scratch/pysmp/test/aaadms2.fits')
        #if stop:
        #    sys.exit()
        if dogalsim:
            print imfile
            wcs = galsim.FitsWCS(imfile) #read in wcs
            psfstamp = galsim.ImageF(array=copy(psfout*0.)) #make emtpy stamp
            des_psfex = galsim.des.DES_PSFEx(psffile) #read in psf file
            psfcenter = galsim.PositionD(xin,yin) #set galsim position
            gs_psfex = des_psfex.getPSF(psfcenter) #get psf at position
            worldpsf = wcs.toWorld(gs_psfex,image_pos=psfcenter) #convert coords
            off = galsim.PositionD(-IMAGE_CENTERX+xin,-IMAGE_CENTERY+yin) #set subpixel offset
            gs_psfex.drawImage(image=psfstamp,offset=off) #draw psf to stamp at offset
            psfout = psfstamp.array

        return(psfout), (IMAGE_CENTERX, IMAGE_CENTERY)

def scene_check(p,x=None,y=None,fjac=None,params=None,err=None):
    """Scene modeling function, but with error 
    measurements and optionally saves stamps"""
    status = 0

    Nimage = len(x[:,0,0])
    substamp = float(params.substamp)
        
    model = np.zeros([Nimage,substamp,substamp])
    galaxy = p[0:substamp**2.].reshape(substamp,substamp)

    chi2 = np.zeros(Nimage)
    for i in range(Nimage):
        conv_prod = scipy.ndimage.convolve(galaxy,x[i,:,:])
        # model = scale + convolution + sky
        model[i,:,:] = p[substamp**2.+i]*x[i,:,:] + conv_prod + p[substamp**2+Nimage+i]

        xx = np.where(err < 10000.0)
        chi2[i]=np.sqrt(np.sum((model[xx]-y[xx])**2/err[xx]**2.)/float(len(xx)))
        if debug:
            import os
            import astropy.io.fits as pf
            if not os.path.exists('Stamps'):
                os.makedirs('Stamps')
            pf.writeto('Stamps/image{0}.fits'.format(i), y[i,:,:], clobber=True)
            pf.writeto('Stamps/model{0}.fits'.format(i), model[i,:,:], clobber=True)
            pf.writeto('Stamps/diff{0}.fits'.format(i), model[i,:,:] - y[i, :,:], clobber=True)
            pf.writeto('Stamps/chi2_{0}.fits'.format(i), (model[i,:,:]-y[i,:,:])**2/err[i,:,:]**2/float(len(model[i,:,:].flatten())), clobber=True)
            pf.writeto('Stamps/error{0}.fits'.format(i), err[i,:,:], clobber=True)
            pf.writeto('Stamps/psf{0}.fits'.format(i), x[i,:,:], clobber=True)
            pf.writeto('Stamps/sn{0}.fits'.format(i), p[substamp**2.+i]*x[i,:,:], clobber=True)
            pf.writeto('Stamps/gal{0}.fits'.format(i), conv_prod, clobber=True)
    return(chi2)



def scene(p,x=None,y=None,fjac=None,params=None,err=None):
    b = open('/global/u1/d/dbrout/PySMP/smpout_bigrun/didstart.txt','w')
    b.write('yes')
    b.close()
    """Scene modeling function given to mpfit"""
    status = 0

    Nimage = len(x[:,0,0])
    substamp = float(params.substamp)
        
    model = np.zeros([Nimage,substamp,substamp])
    galaxy = p[0:substamp**2.].reshape(substamp,substamp)
    conv_prod = np.zeros([Nimage,substamp,substamp])
    for i in range(Nimage):
        conv_prod[i] = scipy.ndimage.convolve(galaxy,x[i,:,:])
        # model = scale + convolution + sky
    model = (p[substamp**2.:substamp**2 +Nimage]*x.T + conv_prod.T + p[substamp**2+Nimage:]).T
    return(status, (y.reshape(Nimage*substamp*substamp)-model.reshape(Nimage*substamp*substamp))/err.reshape(Nimage*substamp*substamp))

def resid(param, psf, im, sigma, fitrad, sky, psfmag):
    model = psf * param / 10 ** (-0.4 * (psfmag - 25)) + sky
    residsig = (im - model) / sigma
    return np.array(residsig.ravel())

def starresid(param, psf, im, se, fitrad, sky, gain):
    weight = 1. / (se**2 + np.max([0, float(param)]) / gain + 1.)**.5
    model = psf * param  + sky
    residsig = (im - model) * weight * fitrad
    return np.array(residsig.ravel())


def simstamp(param, psf, im, sigma, fitrad, sky, psfmag):
    model = psf * param / 10 ** (-0.4 * (psfmag - 25)) + sky
    return model


if __name__ == "__main__":
    print 'ENTERING SMP'*100
    import sys,getopt
    # read in arguments and options
    try:
        if os.path.exists("default.config"):
            args = open("default.config", 'r').read().split()
        else:
            args = sys.argv[1:]
        
        #print args
        opt,arg = getopt.getopt(
            args,"hs:p:r:f:o:m:v:i:d:s",
            longopts=["help","snfile=","params=","rootdir=",
                      "filter=","nomask","nodiff","nozpt", "outfile=",
                      "mergeno=", "loadzpt","usefake","savezptstamps",
                      "debug","verbose","clearzpt",
                      "psf_model=","ismultiple",
                      "gal_model=","index=","diffimzpt","idlsky",
                      "dontgalfit","dontsnfit","dogalsimfit","dogalsimpixfit",
                      "fixgalzero","floatallepochs","dailyoff","snradecfit","dontglobalstar",
                      "snfilepath=","bigstarcatalog=",
                      "galaxyfoldername=",
                      "snfilelist=","files_split_by_filter","maskandnoise","stardumppsf",
                      "dosextractor","useweights","fermigrid","zptoutpath=",
                      "embarrasinglyParallelEnvVar=","fermigriddir=","worker",
                      "lcfilepath=","fermilog","isdonedir=","oldformat","continue","continuedir="])


        #print opt
        #print arg
    except getopt.GetoptError as err:
        print str(err)
        print "Error : incorrect option or missing argument."
        print __doc__
        sys.exit(1)

    try:
        args = sys.argv[1:]
        
        opt_command,arg = getopt.getopt(
            args,"hs:p:r:f:o:m:v:i:d:s",
            longopts=["help","snfile=","params=","rootdir=",
                      "filter=","nomask","nodiff","nozpt", "outfile=",
                      "mergeno=", "loadzpt","usefake","savezptstamps",
                      "debug","verbose","clearzpt",
                      "psf_model=","ismultiple",
                      "gal_model=","index=","diffimzpt","idlsky",
                      "dontgalfit","dontsnfit","dogalsimfit","dogalsimpixfit",
                      "fixgalzero","floatallepcohs","dailyoff","snradecfit","dontglobalstar",
                      "snfilepath=","bigstarcatalog=",
                      "galaxyfoldername=",
                      "snfilelist=","files_split_by_filter","maskandnoise","stardumppsf",
                      "dosextractor","useweights","fermigrid","zptoutpath=",
                      "embarrasinglyParallelEnvVar=","fermigriddir=","worker",
                      "lcfilepath=","fermilog","isdonedir","oldformat","continue","continuedir="])


        #print opt
        #print arg
    except getopt.GetoptError as err:
        print "No command line arguments"


    verbose,nodiff,debug,clear_zpt,psf_model,root_dir,mergeno,loadzpt,ismultiple,dogalfit,dosnfit,dogalsimfit,dogalsimpixfit = False,False,False,False,False,False,False,False,False,True,True,False,False
    oldformat = False
    savezptstamps = False
    fixgalzero,floatallepochs = False,False
    dailyoff = False
    usediffimzpt = False
    useidlsky = False
    snradecfit = False
    doglobalstar = True
    index = None
    gal_model = None
    snfile,param_file,outfile,filt = '','','','r'
    nomask,nozpt = 'none',False
    mergeno = 0
    snfilepath = None
    bigstarcatalog=None
    galaxyfoldername=''
    snfilelist = None
    files_split_by_filter = False
    useweights = False
    stardumppsf = False
    dosextractor=False
    fermigrid = False
    zptoutpath = './zpts/'
    isEmbarrasinglyParallel = False
    parallelvar = None
    fermigriddir = None
    worker = False
    lcfilepath=snfilepath
    fermilog = False
    isdonedir = None
    continu = False
    continudir = None

    dobigstarcat = True

    usefake = False

    print opt

    for o,a in opt:
        if o in ["-h","--help"]:
            print __doc__
            sys.exit(0)
        elif o in ["-s","--snfile"]:
            snfile = a
        elif o in ["-p","--params"]:
            param_file = a
        elif o in ["-r","--rootdir"]:
            root_dir = a
        elif o in ["-f","--filter"]:
            filt = a
        elif o in ["-v","--verbose"]:
            verbose = True
        elif o in ["-o","--outfile"]:
            outfile = a
            out_dir = a
        elif o in ["-m","--mergeno"]:
            mergeno = int(a)
        elif o in ["--loadzpt"]:
            loadzpt = True
        elif o == "--nomask":
            nomask = True
        elif o == "--usefake":
            usefake = True
        elif o == "--nodiff":
            nodiff = True
        elif o == "--nozpt":
            nozpt = True
        elif o == "--debug":
            debug = True
        elif o == "--ismultiple":
            ismultiple = True
        elif o == "--psf_model":
            psf_model = a.lower()
        elif o in ["-g","--gal_model"]:
            gal_model = a        
        elif o in ["-i","--index"]:
            index = a
        elif o in ["-d","--diffimzpt"]:
            usediffimzpt = True
        elif o in ["-s","--idlsky"]:
            useidlsky = True
        elif o in ["--dontgalfit"]:
            dogalfit = False
        elif o in ["--dontsnfit"]:
            dosnfit = False
        elif o in ["--dogalsimfit"]:
            dogalsimfit = True
        elif o in ["--dogalsimpixfit"]:
            dogalsimpixfit = True
        elif o in ["--fixgalzero"]:
            fixgalzero = True
        elif o in ["--floatallepochs"]:
            floatallepochs = True
        elif o in ["--dailyoff"]:
            dailyoff = True
        elif o in ["--snradecfit"]:
            snradecfit = True
        elif o in ["--dontglobalstar"]:
            doglobalstar = False
        elif o in ["--snfilepath"]:
            snfilepath = a
        elif o in ["--lcfilepath"]:
            lcfilepath = a
        elif o in["--bigstarcatalog"]:
            bigstarcatalog = a
        elif o in ["--galaxyfoldername"]:
            galaxyfoldername = a
        elif o in ["--snfilelist"]:
            snfilelist = a
        elif o == "--files_split_by_filter":
            files_split_by_filter = True
        elif o == "--maskandnoise":
            useweights = False
        elif o == "--stardumppsf":
            stardumppsf = True
        elif o == "--dosextractor":
            dosextractor = True
        elif o == "--useweights":
            useweights = True
        elif o == "--fermigrid":
            fermigrid = True
        elif o == "--zptoutpath":
             zptoutpath = a
        elif o == "--fermigriddir":
            fermigriddir = a
        elif o == "--isdonedir":
            isdonedir = a
        elif o == "--embarrasinglyParallelEnvVar":
            isEmbarrasinglyParallel = True
            parallelvar= a
        elif o == "--worker":
            worker = True
        elif o == "--savezptstamps":
            savezptstamps = True
        elif o == "--fermilog":
            fermilog = True
        elif o == "--oldformat":
            oldformat = True
        elif o == "--continue":
            continu = True
        elif o == "--continuedir":
            continudir = a
        else:
            print "Warning: option", o, "with argument", a, "is not recognized"


    #THESE WiLL OVERRIDE THE DEFAULTS
    for o,a in opt_command:
        if o in ["-h","--help"]:
            print __doc__
            sys.exit(0)
        elif o in ["-s","--snfile"]:
            snfile = a
        elif o in ["-p","--params"]:
            param_file = a
        elif o in ["-r","--rootdir"]:
            root_dir = a
        elif o in ["-f","--filter"]:
            filt = a
        elif o in ["-v","--verbose"]:
            verbose = True
        elif o in ["-o","--outfile"]:
            out_dir = a
            outfile = a
        elif o in ["-m","--mergeno"]:
            mergeno = int(a)
        elif o in ["--loadzpt"]:
            loadzpt = True
        elif o == "--nomask":
            nomask = True
        elif o == "--usefake":
            usefake = True
        elif o == "--nodiff":
            nodiff = True
        elif o == "--nozpt":
            nozpt = True
        elif o == "--debug":
            debug = True
        elif o == "--ismultiple":
            ismultiple = True
        elif o == "--psf_model":
            psf_model = a.lower()
        elif o in ["-g","--gal_model"]:
            gal_model = a        
        elif o in ["-i","--index"]:
            index = a
        elif o in ["-d","--diffimzpt"]:
            usediffimzpt = True
        elif o in ["-s","--idlsky"]:
            useidlsky = True
        elif o in ["--dontgalfit"]:
            dogalfit = False
        elif o in ["--dontsnfit"]:
            dosnfit = False
        elif o in ["--dogalsimfit"]:
            dogalsimfit = True
        elif o in ["--dogalsimpixfit"]:
            dogalsimpixfit = True
        elif o in ["--fixgalzero"]:
            fixgalzero = True
        elif o in ["--floatallepochs"]:
            floatallepochs = True
        elif o in ["--dailyoff"]:
            dailyoff = True
        elif o in ["--snradecfit"]:
            snradecfit = True
        elif o in ["--dontglobalstar"]:
            doglobalstar = False
        elif o in ["--snfilepath"]:
            snfilepath = a
        elif o in ["--lcfilepath"]:
            lcfilepath = a
        elif o in["--bigstarcatalog"]:
            bigstarcatalog = a
        elif o in ["--galaxyfoldername"]:
            galaxyfoldername = a
        elif o in ["--snfilelist"]:
            snfilelist = a
        elif o == "--files_split_by_filter":
            files_split_by_filter = True
        elif o == "--maskandnoise":
            useweights = False
        elif o == "--stardumppsf":
            stardumppsf = True
        elif o == "--dosextractor":
            dosextractor = True
        elif o == "--useweights":
            useweights = True
        elif o == "--fermigrid":
            fermigrid = True
        elif o == "--fermigriddir":
            fermigriddir = a
        elif o == "--zptoutpath":
            zptoutpath = a
        elif o == "--embarrasinglyParallelEnvVar":
            isEmbarrasinglyParallel = True
            parallelvar= a
        elif o == "--worker":
            worker = True
        elif o == "--savezptstamps":
            savezptstamps = True
        elif o == "--fermilog":
            fermilog = True
        elif o == "--isdonedir":
            isdonedir = a
        elif o == "--oldformat":
            oldformat = True
        elif o == "--continue":
            continu = True
        elif o == "--continuedir":
            continudir = a
        else:
            print "Warning: option", o, "with argument", a, "is not recognized"


    #worker = True
    if isEmbarrasinglyParallel:
        index = os.environ[parallelvar]

    if not os.path.exists(zptoutpath):
        if fermigrid & worker:
            if zptoutpath.split('/')[1] != 'pnfs':
                raise ValueError('--zptoutpath must be located at /pnfs/des/persistent/desdm/ for fermigrid running')
            os.system( 'ifdh mkdir '+zptoutpath)
        else:
            os.makedirs(zptoutpath)
    # if fermigrid:
    #     print '4200'
    #     param_file = os.path.join(fermigriddir, param_file)
    #     os.system('ifdh cp ' + param_file + ' .')
    #     param_file = param_file.split('/')[-1]
    #     print 'paramfile',param_file
    #     os.system('ls -ltr')
    if bigstarcatalog is None:
        dobigstarcat = False

    if galaxyfoldername is None:
        raise NameError("Must provide " +
                        "--galaxyfoldername=/location/to/previous/run in default.config \n Exiting now...")

    hostname = os.popen('hostname').read()
    if 'dsp053' in hostname:
        snfilelist = '/home/dscolnic/testdir/runfile_53.txt'
    if snfilelist is None:
        raise NameError("Must provide " +
                        "--snfilelist=/location/to/a/list/of/snfiles in default.config \n Exiting now...")

    if isdonedir is None:
        isdonedir = os.path.join(out_dir,'isdone')

    if 'dsp053' in hostname:
        isdonedir = '/home/dscolnic/testdir/isdone'

    if 'dsp057' in hostname:
        isdonedir = '/home/dscolnic/testdir/isdone'

    if not index is None:
        if index == 'all':
            for iii in np.arange(0,5000):

                a = open(snfilelist,'r')
                files = a.readlines()
                print 'files',files
                print 'index',index
                snfile = files[int(iii)].rstrip()
                a.close()
                if not snfilepath is None:
                    snfile = os.path.join(snfilepath, snfile.split('/')[-1])
                else:
                    lcfilepath = '.'
                print 'Index '+str(iii)

                print 'SN File '+snfile

                if not snfile or not param_file:
                    print("Error : snfile and params  must be provided")
                    print(__doc__)
                    sys.exit(1)

                if not mergeno:
                    mergeno = int(0)

                if not root_dir:
                    print("root_dir not specified. Assuming same directory as snfile...")
                    try:
                        root_dir = snfile.split('/')[:-1].join()
                    except:
                        root_dir = './'
                if not psf_model:
                    print("psf_model not specified. Assuming psfex...")
                    psf_model = 'psfex'


                snparams = get_snfile(snfile, root_dir, useweights)
                # if fermigrid:
                #     print 'ifdh 4258'
                #     os.system('ifdh cp '+param_file)
                print 'getting params'
                params = get_params(param_file)

                if nomask == 'none':
                    if params.mask_type.lower() == 'none':
                        nomask = True
                    else: nomask = False
                if nozpt == 'none':
                    if params.find_zpt.lower() == 'yes':
                        nozpt = True
                    else: nozpt = False
                if not filt:
                    print("Filt not defined.  Using all...")
                    filt = snparams.filters
                if not outfile:
                    print "Output file name not defined. Using /path/to/snfile/test.out ..."
                    try:
                        out_dir = snfile.split('/')[:-1].join()
                    except:
                        out_dir = './'
                    outfile = os.path.join(out_dir,'test.out')
                if not mergeno:
                    mergeno = 0
                try:
                    scenemodel = smp(snparams,params,root_dir,psf_model)
                    scenemodel.main(nodiff=nodiff,nozpt=nozpt,nomask=nomask,debug=debug,outfile=outfile,rootdir=root_dir,outdir=out_dir
                                 ,verbose=verbose,clear_zpt=True, mergeno=mergeno,usefake=usefake,snfile=snfile,
                                 gal_model=gal_model,stardumppsf=stardumppsf,dogalfit=dogalfit,dosnfit=dosnfit,
                                 dogalsimfit=dogalsimfit,dogalsimpixfit=dogalsimpixfit,dosnradecfit=snradecfit,
                                 usediffimzpt=usediffimzpt,useidlsky=useidlsky,fixgalzero=fixgalzero,floatallepochs=floatallepochs,
                                 dailyoff=dailyoff,doglobalstar=doglobalstar,bigstarcatalog=bigstarcatalog,dobigstarcat=dobigstarcat,
                                 galaxyfoldername=galaxyfoldername,isdonedir=isdonedir,
                                 useweights=useweights,dosextractor=dosextractor,fermigrid=fermigrid,zptoutpath=zptoutpath,
                                 fermigriddir=fermigriddir,worker=worker,lcfilepath=lcfilepath,savezptstamps=savezptstamps,
                                    fermilog=fermilog,oldformat=oldformat,continu=continu,continudir=continudir)
                    #scenemodel.afterfit(snparams,params,donesn=True)
                    print "SMP Finished!"
                except:
                    pass
            sys.exit()

        else:
            a = open(snfilelist,'r')
            files = a.readlines()

            snfile = files[int(index)].rstrip()

            a.close()
            if not snfilepath is None:
                osnfile = copy(snfile.split('/')[-1])
                snfile = os.path.join(snfilepath, snfile.split('/')[-1])

            print 'Index '+str(index)
            print 'SN File '+snfile
            print 'not stopping here because dont know what happens after'
            if fermigrid & worker:
                print 'ifdh 4366'
                os.system('ifdh cp --force=xrootd '+snfile+' .')
                snfile = osnfile

    '''if ismultiple:
        pbsint = os.environ['PBS_ARRAYID']
        a = open('snfiles_mjd.txt','r')
        files = a.readlines()
        snfile = files[int(pbsint)].rstrip()
        a.close()
        #b = open('snfiles_'+str(pbsint)+'.txt','w')
        #b.write(snfile)
        #b.close()
    '''
    if not snfile or not param_file:
        print("Error : snfile and params  must be provided")
        print(__doc__)
        sys.exit(1)

    if not mergeno:
        mergeno = int(0)

    if not root_dir:
        print("root_dir not specified. Assuming same directory as snfile...")
        try:
            root_dir = snfile.split('/')[:-1].join()
        except:
            root_dir = './'


    if files_split_by_filter:
        filt = snfile.split('_')[1].split('.')[0]
    print 'getting snparams'
    print snfile,root_dir
    try:
        snparams = get_snfile(snfile, root_dir, useweights)
    except:
        print sys.exc_info()
        if not os.path.exists(isdonedir):
            os.mkdir(isdonedir)
        os.system('touch ' + os.path.join(isdonedir, snfile.split('/')[-1].split('.')[0] + '.done'))
        sys.exit()
    print 'getting params'
    params = get_params(param_file)
    print 'done with params'
    if not params.psf_model:
        print("psf_model not specified. Assuming psfex...")
        psf_model = 'psfex'
    else:
        psf_model = params.psf_model
    if nomask == 'none':
        if params.mask_type.lower() == 'none':
            nomask = True
        else: nomask = False
    if nozpt == 'none':
        if params.find_zpt.lower() == 'yes':
            nozpt = True
        else: nozpt = False
    if not filt:
        print("Filt not defined.  Using all...")
        filt = snparams.filters
    if not outfile:
        print "Output file name not defined. Using /path/to/snfile/test.out ..."
        try:
            out_dir = snfile.split('/')[:-1].join()
        except:
            out_dir = './'
        outfile = os.path.join(out_dir,'test.out')
    if not mergeno:
        mergeno = 0
    print 'beginning smp'
    #sys.exit()
    scenemodel = smp(snparams,params,root_dir,psf_model)
    #print out_dir
    #print dobigstarcat
    #print bigstarcatalog
    #raw_input()
    if snparams.survey == 'PS1111':
        #print 'isps1'
        #raw_input()
        try:
            scenemodel.main(nodiff=nodiff, nozpt=nozpt, nomask=nomask, debug=debug, outfile=outfile, rootdir=root_dir,
                            outdir=out_dir,
                            verbose=verbose, clear_zpt=True, mergeno=mergeno, usefake=usefake, snfile=snfile,
                            gal_model=gal_model, stardumppsf=stardumppsf, dogalfit=dogalfit, dosnfit=dosnfit,
                            dogalsimfit=dogalsimfit, dogalsimpixfit=dogalsimpixfit, dosnradecfit=snradecfit,
                            usediffimzpt=usediffimzpt, useidlsky=useidlsky, fixgalzero=fixgalzero,
                            floatallepochs=floatallepochs,
                            dailyoff=dailyoff, doglobalstar=doglobalstar, bigstarcatalog=bigstarcatalog,
                            dobigstarcat=dobigstarcat,
                            galaxyfoldername=galaxyfoldername, isdonedir=isdonedir,
                            useweights=useweights, dosextractor=dosextractor, fermigrid=fermigrid,
                            zptoutpath=zptoutpath,
                            fermigriddir=fermigriddir, worker=worker, savezptstamps=savezptstamps,
                            fermilog=fermilog,oldformat=oldformat,continu=continu,continudir=continudir)
        except:
            print sys.exc_info()
            if not os.path.exists(isdonedir):
                os.mkdirs(isdonedir)
            os.system('touch '+os.path.join(isdonedir,snparams.snfile.split('/')[-1].split('.')[0] + '.done'))

    else:
        scenemodel.main(nodiff=nodiff,nozpt=nozpt,nomask=nomask,debug=debug,outfile=outfile,rootdir=root_dir,outdir=out_dir,
                     verbose=verbose,clear_zpt=True, mergeno=mergeno,usefake=usefake,snfile=snfile,
                     gal_model=gal_model,stardumppsf=stardumppsf,dogalfit=dogalfit,dosnfit=dosnfit,
                     dogalsimfit=dogalsimfit,dogalsimpixfit=dogalsimpixfit,dosnradecfit=snradecfit,
                     usediffimzpt=usediffimzpt,useidlsky=useidlsky,fixgalzero=fixgalzero,floatallepochs=floatallepochs,
                     dailyoff=dailyoff,doglobalstar=doglobalstar,bigstarcatalog=bigstarcatalog,dobigstarcat=dobigstarcat,
                     galaxyfoldername=galaxyfoldername,isdonedir=isdonedir,
                     useweights=useweights,dosextractor=dosextractor,fermigrid=fermigrid,zptoutpath=zptoutpath,
                     fermigriddir=fermigriddir,worker=worker,savezptstamps=savezptstamps,
                    fermilog=fermilog,oldformat=oldformat,continu=continu,continudir=continudir)
        #scenemodel.afterfit(snparams,params,donesn=True)
    print "SMP Finished!"
     
