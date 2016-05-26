#!/usr/bin/env python

stardeltasfolder = 'smp_y1y2_shallow_v3_40globalstars'
foldername = 'smp_y1y2_shallow_v3_58newsky_newzpt_exactpos_notgalsim'
galaxyfoldername = 'smp_y1y2_shallow_v3_40globalstars'
pixstart = None

"""
Scene modeling pipeline for DES and PanSTARRS.

Usage:

smp.py -s supernova_file -p parameter_file -s snfile -o outfile\ 
    --nomask --nodiff --nozpt -r root_dir -f filter --psf_model=psf_model

Default parameters are in parentheses.

-s/--supernova_file      : Filename with all the information about
                           the supernova. (Required)
-p/--params              : Parameter file (Required)
-r/--root_dir            : images root directory (/path/to/snfile)
-f/--filter              : observation filter, use 'all' for all filters 
                           ('all')
-o/--outfile             : output file (/path/to/snfile/test.out) 
-m/--mergeno             : create 2x2 merged pixels 'mergeno' times
-v                       : Verbose mode. Prints more information to terminal.
--nomask                 : set if no mask image exists (one will be 
                           created).
--nodiff                 : set if no difference image exists
--nozpt                  : set if zeropoints have not been measured
--debug                  : debug flag saves intermediate products
                           and prints additional information
--psf_model              : Name of psf model. Currently works for 
                           psfex and daophot. ('psfex')

"""
# TO DO: figure out how p gets passed into scene/mpfit...
import numpy as np
import exceptions
import os
#os.environ['SHELL'] = 'csh'
import sys
#os.system("alias funpack '/global/u1/d/dbrout/cfitsio/funpack'")
sys.path.append("/global/homes/d/dbrout/GalSim-1.3.0")
sys.path.append("/global/homes/d/dbrout/GalSim-1.3.0/lib")
import scipy.ndimage
import matplotlib as m
#import mcmc_emceen as mcmc
import mcmc as mcmc3
#import mcmctest as mcmc3
#import mcmcfloat as mcmcfloat
#import mcmc3galsim
import mcmcgalsim as mcmc3galsimpixshift
#import mcmcgtest as mcmc3galsimpixshift
m.use('Agg')
import matplotlib.pyplot as plt
import time
import pyfits as pf
import scipy.signal
import scipy.ndimage as nd
from copy import copy
import galsim
import galsim.des
import time
from astropy.io import fits
from scipy.interpolate import UnivariateSpline
import sigma_clip
import meanclip
from scipy.optimize import curve_fit
#import testsim
#import rdcol
#import getgalmodel
import scipy.interpolate as interpol
import cntrd,aper,getpsf,rdpsf
#from PythonPhot import cntrd,aper,getpsf,rdpsf
import pkfit_norecent_noise_smp
import addcoltoDESlightcurve as lc
import runsextractor

#from matplotlib.backends.backend_pdf import PdfPages

snkeywordlist = {'SURVEY':'string','SNID':'string','FILTERS':'string',
                 'PIXSIZE':'float','NXPIX':'float','NYPIX':'float',
                 'ZPFLUX':'float','RA':'string', 'FAKE_RA':'string','FAKE_DEC':'string',
                 'DECL':'string','PEAKMJD':'float','WEIGHT_BADPIXEL':'string',
                 'HOSTGAL_SB_FLUXCAL':[],
                 'STARCAT':'string', 'PSF_UNIT':'string', 'NOBS':'float'}
snvarnameslist = {'ID_OBS':'string','MJD':'float','BAND':'string',
                  'IMAGE_NAME_SEARCH':'string','IMAGE_NAME_WEIGHT':'string',
                  'FILE_NAME_PSF':'string','FAKE_TRUEMAG':'float','ZP':'float',
                  'FLUX':'float','FLUXERR':'float','PHOTFLAG':'string','SKYSIG':'float'}
paramkeywordlist = {'STAMPSIZE':'float','RADIUS1':'float',
                    'RADIUS2':'float','SUBSTAMP':'float',
                    'MAX_MASKNUM':'float','RDNOISE_NAME':'string',
                    'GAIN_NAME':'string','FWHM_MAX':'float',
                    'PSF_MAX':'float','WEIGHT_TYPE':'string',
                    'MASK_TYPE':'string','MJDPLUS':'float','MJDMINUS':'float',
                    'BUILD_PSF':'string','CNTRD_FWHM':'float','FITRAD':'float',
                    'FORCERADEC':'string','FRA':'float','FDEC':'float',
                    'FIND_ZPT':'string','PIXELATION_FACTOR':'float','SKYERR_RADIUS':'float',
                    'NEARBY_STARS_PIXEL_RAD':'float'
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
                #xvals = (bf+bs)/2.
                ww = [(x>bs)&(x<bf)]
                try:
                        medians[i] = np.median(y[ww])
                        mads[i] = 1.48*np.median(abs(y[ww]-medians[i]))*1/np.sqrt(len(y[ww]))
                        #xvals[i] = xv
                except IndexError:
                        medians[i] = np.nan
                        mads[i] = np.nan
                        #xvals[i] = np.nan
        return xvals,medians,mads
def parabola(x,a,b,c):
    y = a*x**2+b*x+c
    return y

class get_snfile:
    def __init__(self,snfile, rootdir):
        varnames = ''
        fin = open(snfile,'r')
        for line in fin:
            line = line.replace('\n','')
            if not line.startswith('#') and line.replace(' ',''):
                if not line.replace(' ','').startswith('OBS:') and \
                        not line.replace(' ','').startswith('VARNAMES:'):
                    key,val = line.split('#')[0].split(':')
                    key = key.replace(' ','')
                    if key.upper() == 'HOSTGAL_SB_FLUXCAL':
                        val = val.split()
                        self.__dict__[key.lower()] = val
                    elif key.upper() != 'WEIGHT_BADPIXEL' and (key.upper() != 'STARCAT' or not 'des' in snfile):             
                        val = val.split()[0]
                        val = val.replace(' ','')
                        self.__dict__[key.lower()] = val
                    elif key.lower() == 'starcat' and 'des' in snfile:
                        catfilter = val.split()[0]                        
                        if filt.lower() == catfilter.lower():
                            print val
                            self.__dict__["starcat"] = {catfilter.lower(): os.path.join(rootdir,val.split()[1])}
                        elif filt.lower() == 'all':
                            if "starcat" in self.__dict__:
                                self.__dict__["starcat"][val.split()[0]] = os.path.join(rootdir,val.split()[1])
                            else:
                                self.__dict__["starcat"] = {}
                                self.__dict__["starcat"][val.split()[0]] = os.path.join(rootdir,val.split()[1])
                    else:
                        try:
                            self.__dict__[key.lower()] = np.array(val.split()).astype('float')
                        except:
                            raise exceptions.RuntimeError("Error : WEIGHT_BADPIXEL cannot be parsed!")
                
                #elif line.replace(' ','').startswith('VARLIST:'):
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
                        #print var
                        #print val
                        #print snfile
                        #raw_input()
                        self.__dict__[var.lower()] = np.append(self.__dict__[var.lower()],val)

        catalog_exists = True
        for p in snkeywordlist.keys():
            #print p
            #raw_input()
            #print self.__dict__.keys()
            if not self.__dict__.has_key(p.lower()):
                if p.lower() != 'starcat':
                    raise exceptions.RuntimeError("Error : keyword %s doesn't exist in supernova file!!!"%p)
                else:
                    catalog_exists = False
            if snkeywordlist[p] == 'float':
                try:
                    self.__dict__[p.lower()] = float(self.__dict__[p.lower()])
                except:
                    raise exceptions.RuntimeError('Error : keyword %s should be set to a number!'%p)

        for p in snvarnameslist.keys():
            #print p
            #raw_input()
            if not self.__dict__.has_key(p.lower()):
                if p.lower() != 'starcat':
                    raise exceptions.RuntimeError("Error : field %s doesn't exist in supernova file!!!"%p)
                elif catalog_exists == False:
                    raise exceptions.RuntimeError("Error : field %s doesn't exist in supernova file!!!"%p)
            if snvarnameslist[p] == 'float':
                try:
                    self.__dict__[p.lower()] = self.__dict__[p.lower()].astype('float')
                except:
                    raise exceptions.RuntimeError('Error : keyword %s should be set to a number!'%p)
        #print self.__dict__.keys()
        #raw_input()

class get_params:
    def __init__(self,paramfile):
        fin = open(paramfile,'r')
        for line in fin:
            line = line.replace('\n','')
            if not line.startswith('#') and line.replace(' ',''):
                try:
                    key,val = line.split('#')[0].split(':')
                except:
                    raise exceptions.RuntimeError('Invalid format!  Should be key: value')
                key = key.replace(' ','')
                val = val.replace(' ','')
                self.__dict__[key.lower()] = val
        for p in paramkeywordlist.keys():
            #print paramfile
            #print p
            if not self.__dict__.has_key(p.lower()):
                raise exceptions.RuntimeError("Error : keyword %s doesn't exist in parameter file!!!"%p)
            if paramkeywordlist[p] == 'float':
                try:
                    self.__dict__[p.lower()] = float(self.__dict__[p.lower()])
                    #print float(self.__dict__[p.lower()])
                except:
                    raise exceptions.RuntimeError('Error : keyword %s should be set to a number!'%p)
            #raw_input()

class smp:
    def __init__(self,snparams,params,rootdir,psf_model):
        self.snparams = snparams
        self.params = params
        self.rootdir = rootdir
        self.psf_model = psf_model
        #print 'pdb'
        #pdb.set_trace()
        

    def main(self,nodiff=False,nozpt=False,
             nomask=False,outfile='',debug=False,
             verbose=False, clear_zpt=False,clear_checkstars=True,mergeno=0,
             mpfit_or_mcmc='mpfit',usefake=False,
             snfile='/test.dat',gal_model=None,stardumppsf=True,
             dogalfit=True,dosnfit=True,dogalsimfit=True, dogalsimpixfit=True,dosnradecfit=True,
             usediffimzpt=False,useidlsky=False,fixgalzero=True,floatallepochs=False,dailyoff=False,
             doglobalstar=True,exactpos=True):

        tstart = time.time()
        from txtobj import txtobj
        from astropy import wcs
        import astropy.io.fits as pyfits
        self.outfile = outfile
        self.checkstarfile = os.path.join(outfile,stardeltasfolder+'/SNe/'+snfile.split('/')[-1].split('.')[0] + '/'+filt+'/standardstarfits.txt')
        self.snfn = snfile.split('/')[-1].split('.')[0]
        if nozpt:
            self.zpt_fits = './zpts/zpt_plots.txt'
            self.big_zpt = './zpts/big_zpt'
            #self.checkstarfile = os.path.join(outfile,foldername+'/SNe/'+snfile.split('/')[-1].split('.')[0] + '/'+filt+'/standardstarfits.txt')
            print 'checkstarfile',self.checkstarfile
            #self.checkstarfile = self.params.checkstarfile
            if not os.path.exists('/'.join(self.checkstarfile.split('/')[:-1])):
                os.makedirs('/'.join(self.checkstarfile.split('/')[:-1]))
            if not os.path.exists('./zpts'):
                os.makedirs('./zpts/')
            a = open(self.zpt_fits,'w')
            a.write('ZPT FILE LOCATIONS\n')
            a.close()
            if clear_zpt:
                big = open(self.big_zpt+'.txt','w')
                big.write('Exposure Num\tRA\tDEC\tCat Zpt\tMPFIT Zpt\tMPFIT Zpt Err\tMCMC Zpt\tMCMC Zpt Err\tMCMC Model Errors Zpt\tMCMC Model Errors Zpt Err\tCat Mag\tMP Fit Mag\tMCMC Fit Mag\tMCMC Model Errors Fit Mag\tMCMC Analytical Simple\tMCMC Analytical Weighted\n')
                big.close()
            if clear_checkstars:
                big = open(self.checkstarfile,'w')
                big.write('Exposure Num\tMJD\tRA\tDEC\tCat Zpt\tMPFIT Zpt\tMPFIT Zpt Err\tFit Flux\tFit Flux Err\tCat Mag\n')
                big.close()
        self.verbose = verbose
        params,snparams = self.params,self.snparams
        print "FAKE TRUE MAGS"
        print snparams.fake_truemag
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

        self.bigcatalog = pf.open('/global/homes/d/dbrout/PySMP/SNscampCatalog/DES-SN_v2.cat')[2].data
        self.bigcatalogmags = self.bigcatalog['mag']
        self.bigcatalogras = self.bigcatalog['x_world']
        self.bigcatalogdecs = self.bigcatalog['y_world']
        #print self.bigcatalog.data
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

        if snparams.psf_model == 'psfex' and not snparams.__dict__.has_key('psf'):
            raise exceptions.RuntimeError('Error : PSF must be provided in supernova file!!!')
        if filt != 'all':
            snparams.nvalid = 0
          
            for b in snparams.band:
                if b in filt:
                    snparams.nvalid +=1
        else:
            snparams.nvalid = snparams.nobs

        smp_im = np.zeros([snparams.nvalid,params.substamp,params.substamp])
        smp_noise = np.zeros([snparams.nvalid,params.substamp,params.substamp])
        smp_psf = np.zeros([snparams.nvalid,params.substamp,params.substamp])
#        smp_bigim = np.zeros([snparams.nvalid,params.stampsize,params.stampsize])
#        smp_bignoise = np.zeros([snparams.nvalid,params.stampsize,params.stampsize])
#        smp_bigpsf = np.zeros([snparams.nvalid,params.stampsize,params.stampsize])
        

        smp_dict = {'scale':np.zeros(snparams.nvalid),
                    'scale_err':np.zeros(snparams.nvalid),
                    'mcmc_scale':np.zeros(snparams.nvalid),
                    'mcmc_scale_err':np.zeros(snparams.nvalid),
                    'image_scalefactor':np.zeros(snparams.nvalid),
                    'snx':np.zeros(snparams.nvalid),
                    'sny':np.zeros(snparams.nvalid),
                    'fwhm_arcsec':np.zeros(snparams.nvalid),
                    'sky':np.zeros(snparams.nvalid),
                    'skyerr':np.zeros(snparams.nvalid),
                    'flag':np.ones(snparams.nvalid),
                    'fitflag':np.ones(snparams.nvalid),
                    'psf':np.zeros(snparams.nvalid),
                    'psf_fwhm':np.zeros(snparams.nvalid),
                    'fakepsf':np.zeros(snparams.nvalid),
                    'zpt':np.zeros(snparams.nvalid),
                    'zpterr':np.zeros(snparams.nvalid),   
                    'skysig':np.zeros(snparams.nvalid),
                    'total_skyerr':np.zeros(snparams.nvalid),         
                    'mjd':np.zeros(snparams.nvalid),
                    'mjd_flag':np.zeros(snparams.nvalid),
                    'mjdoff':[],
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
                    'snra':np.zeros(snparams.nvalid),
                    'sndec':np.zeros(snparams.nvalid),
                    'notbrightflag':np.ones(snparams.nvalid)

                    }
        smp_scale = np.zeros(snparams.nvalid)
        smp_sky = np.zeros(snparams.nvalid)
        smp_flag = np.zeros(snparams.nvalid)
        for i in np.arange(snparams.nvalid):
            smp_dict['image_filename'][i] = 'na'

        snparams.cat_zpts = {}

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
        starcatras = []
        stardecs = []
        cntrs = 0
        cols = None

        filename = snparams.snfile.split('/')[-1].split('.')[0] +'_'+ filt
        outdir = os.path.join(outfile,stardeltasfolder+'/stardata/'+filt+'/')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        star_offset_file = os.path.join(outdir,filename+'band_starGlobalOffsets.npz')

        if not nozpt:
            try:
                staroffsets = np.load(star_offset_file)
            except:
                print 'Could not find star offset file. Calculating...'
                nozpt = True

        
        #############################################################################################################################
        ################################################# GET STAR GLOBAL OFFSETS ###################################################
        
        for imfile,noisefile,psffile,band,faketruemag, j in \
                zip(snparams.image_name_search,snparams.image_name_weight,snparams.file_name_psf,snparams.band,snparams.fake_truemag, range(len(snparams.band))):
            
            if not doglobalstar:
                continue
            if snparams.mjd[j] == 0:
                continue
            if not nozpt:
                continue
            skysig=np.nan
            nozpt = copy(orig_nozpt)

            self.ccdnum = imfile.split('/')[1].split('_')[1]
            self.field = imfile.split('/')[0].split('-')[1]

            if filt != 'all' and band not in filt:
                if verbose: print('filter %s not in filter list for image file %s'%(band,filt,imfile))
                #print 'filter %s,%s not in filter list for image file %s'%(band,filt,imfile)
                continue
            imfile,noisefile,psffile = os.path.join(self.rootdir,imfile),\
                os.path.join(self.rootdir,noisefile),os.path.join(self.rootdir,psffile)
            print imfile
            if not os.path.exists(imfile):
                if not os.path.exists(imfile+'.fz'):
                    print('Error : file %s does not exist'%imfile)
                    continue
                    print('Error : file %s does not exist'%imfile)
                    raise exceptions.RuntimeError('Error : file %s does not exist'%imfile)
                else:
                    os.system('/global/u1/d/dbrout/cfitsio/funpack %s.fz'%imfile)
            if not os.path.exists(noisefile):
                os.system('gunzip %s.gz'%noisefile)
                if not os.path.exists(noisefile):
                    os.system('/global/u1/d/dbrout/cfitsio/funpack %s.fz'%noisefile)
                    if not os.path.exists(noisefile):
                        raise exceptions.RuntimeError('Error : file %s does not exist'%noisefile)
            if not os.path.exists(psffile):
                if not os.path.exists(psffile+'.fz'):
                    raise exceptions.RuntimeError('Error : file %s does not exist'%psffile)
                else:
                    os.system('/global/u1/d/dbrout/cfitsio/funpack %s.fz'%psffile)

            if not nomask:
                maskfile = os.path.join(self.rootdir,snparams.image_name_search[j])

            fakeim = ''.join(imfile.split('.')[:-1])+'+fakeSN.fits'
            if not os.path.exists(fakeim):
                os.system('/global/u1/d/dbrout/cfitsio/funpack %s.fz'%fakeim)
                os.system('/global/u1/d/dbrout/cfitsio/gunzip %s.gz'%fakeim)
            if self.usefake:
                try:
                    im = pyfits.getdata(fakeim)
                except:
                    print fakeim+' is EMPTY, skipping star...'
                    continue
                hdr = pyfits.getheader(imfile)
            else:
                im = pyfits.getdata(imfile)
                hdr = pyfits.getheader(imfile)
            fakeim_hdr = pyfits.getheader(fakeim)
            snparams.cat_zpts[imfile] = fakeim_hdr['HIERARCH DOFAKE_ZP']
            snparams.platescale = hdr['PIXSCAL1']
            snparams.airmass = hdr['AIRMASS']

            noise = pyfits.getdata(noisefile)

            psf = pyfits.getdata(psffile)

            if params.weight_type.lower() == 'ivar':
                print 'ivar'
                raw_input()
                noise = np.sqrt(1/noise)
            elif params.weight_type.lower() != 'noise':
                raise exceptions.RuntimeError('Error : WEIGHT_TYPE value %s is not a valid option'%params.WEIGHT_TYPE)
            if nomask:
                mask = np.zeros(np.shape(noise))
                maskcols = np.where((noise < 0) |
                                    (np.isfinite(noise) == False))
                mask[maskcols] = 100.0
            else:
                mask = pyfits.getdata(maskfile)

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
                print wcsinfo
                ra1, dec1 = results[0]*radtodeg, results[1]*radtodeg
                results2 =  wcsinfo.tran([[snparams.nxpix-1], [snparams.nypix-1]])
                ra2, dec2 =results2[0]*radtodeg, results2[1]*radtodeg

            ra_high = np.max([ra1,ra2])
            ra_low = np.min([ra1,ra2])
            dec_high = np.max([dec1,dec2])
            dec_low = np.min([dec1,dec2])    

            if type(snparams.starcat) == np.array:
                if os.path.exists(snparams.starcat[j]):
                    starcat = txtobj(snparams.starcat[j],useloadtxt=True)
                    if not starcat.__dict__.has_key('mag'):
                        try:
                            starcat.mag = starcat.__dict__[band]
                            starcat.dmag = starcat.__dict__['d%s'%band]
                        except:
                            raise exceptions.RuntimeError('Error : catalog file %s has no mag column!!'%snparams.starcat[j])
                else: 
                    raise exceptions.RuntimeError('Error : catalog file %s does not exist!!'%snparams.starcat[j])
            elif type(snparams.starcat) == dict and 'des' in snfile:
                starcatfile = None
                starcatloc = '/'.join(imfile.split('/')[0:-1])+'/'
                #print 'starcatfile',starcatfile
                #print 'file should be here',os.listdir(starcatloc)
                for fl in os.listdir(starcatloc):
                    print fl
                    if 'STARCAT' in fl:
                        starcatfile = fl
                print starcatloc+starcatfile
                if os.path.exists(starcatloc+starcatfile):
                    starcat = txtobj(starcatloc+starcatfile,useloadtxt=True, des=True)
                    if not starcat.__dict__.has_key('mag_%s'%band):
                        try:
                            print starcat.__dict__
                            starcat.mag = starcat.__dict__[band]
                            starcat.dmag = starcat.__dict__['d%s'%band]
                        except:
                            raise exceptions.RuntimeError('Error : catalog file %s has no mag column!!'%snparams.starcat[band])
                else:
                    raise exceptions.RuntimeError('Error : catalog file %s does not exist!!'%snparams.starcat[band])
            else:
                if os.path.exists(snparams.starcat[j]):
                    starcat = txtobj(snparams.starcat[j],useloadtxt=True)
                    if not starcat.__dict__.has_key('mag'):
                        try:
                            starcat.mag = starcat.__dict__[band]
                            starcat.dmag = starcat.__dict__['d%s'%band]
                        except:
                            print snparams.starcat
                            raise exceptions.RuntimeError('Error : catalog file %s has no mag column!!'%snparams.starcat[j])

                else: 
                    raise exceptions.RuntimeError('Error : catalog file %s does not exist!!'%snparams.starcat[j])


            if nozpt:
                self.rdnoise = hdr[params.rdnoise_name]
                self.gain = hdr[params.gain_name]
                
                cols = np.where((starcat.ra > ra_low) & 
                                (starcat.ra < ra_high) & 
                                (starcat.dec > dec_low) & 
                                (starcat.dec < dec_high))[0]


                if not len(cols):
                    raise exceptions.RuntimeError("Error : No stars in image!!")



                if wcsworked:
                    coords = zip(*w.wcs_world2pix(np.array(zip(starcat.ra,starcat.dec)),0))
                    x_star,y_star = [],[]
                    for xval,yval in zip(*coords):
                        x_star += [xval]
                        y_star += [yval]

                else:
                    coords = wcsinfo.tran([starcat.ra/radtodeg,starcat.dec/radtodeg],False)
                    #print coords
                   #print coords[0]
                    #print len(coords)
                    #print len(starcat.ra[cols])
                    #print len(coords[0])
                    #print coords.shape
                    #print 'cooooords'
                    #raw_input()
                #print starcat.ra[cols][starcat.ra[cols] , starcat.dec[cols]
                #print coords
                #print 'stop'
                #raw_input()

                x_star,y_star = [],[]
                for xval,yval in zip(*coords):
                    x_star += [xval]
                    y_star += [yval]
                
                #print x_star
                #raw_input()
                #print y_star
                #raw_input()

                #print len(coords)
                #print coords[0][16],coords[1][16]
                #print x_star[16],y_star[16]
                #print 'stop'
                #raw_input()

                x_star1,y_star1 = np.array(x_star),np.array(y_star)
                x_star,y_star = cntrd.cntrd(im,x_star1,y_star1,params.cntrd_fwhm)
                newra,newdec = zip(*w.wcs_pix2world(np.array(zip(x_star,y_star)),0))
                for rrr in starcat.objid:
                    starids.append(rrr)
                for rrr,zzz in zip(newra,newdec):
                    starras.append(rrr)
                    stardecs.append(zzz)
                for rrr in starcat.ra:
                    starcatras.append(rrr)
                cntrs += 1

        
        if nozpt:
            starids = np.array(starids)
            starras = np.array(starras)
            stardecs = np.array(stardecs)
            np.savez(star_offset_file,starras=starras,stardecs=stardecs,starids=starids)

        staroffsets = np.load(star_offset_file)
        starras = staroffsets['starras']
        stardecs = staroffsets['stardecs']
        starids = staroffsets['starids']
        #print star_offset_file
        #raw_input()
        starglobalids = []
        starglobalras = []
        starglobaldecs = []
        #starcatras = []
        for ide in np.unique(np.array(starids)):
            ww = (starids == ide)
            starglobalids.append(ide)
            starglobalras.append(np.median(starras[ww]))
            starglobaldecs.append(np.median(stardecs[ww]))
            #starcatras.append(np.median(starcatras[ww]))

        starglobalids = np.array(starglobalids)
        starglobalras = np.array(starglobalras)
        starglobaldecs = np.array(starglobaldecs)
        starcatras = np.array(starcatras)
        
        scampra,scampdec = self.getProperCatRaDec(starglobalras,starglobaldecs)
        offsetra = np.array(starglobalras) - np.array(scampra)
        offsetdec = np.array(starglobaldecs) - np.array(scampdec)
        

        print starglobalras
        #print starras
        #print starcatras - starglobalras
        #print 'checking global offsets'
        #raw_input()


        #############################################################################################################################
        #############################################################################################################################

        cccc = 0
        for imfile,noisefile,psffile,band,faketruemag, j in \
                zip(snparams.image_name_search,snparams.image_name_weight,snparams.file_name_psf,snparams.band,snparams.fake_truemag, range(len(snparams.band))):
            if snparams.mjd[j] == 0:
                continue
            #cccc += 1
            #if cccc > 8:
            #    continue
            #print imfile
            #raw_input()
            skysig=np.nan
            nozpt = copy(orig_nozpt)

            self.ccdnum = imfile.split('/')[1].split('_')[1]
            self.field = imfile.split('/')[0].split('-')[1]
            self.rickfakestarfile = 'data/fixmagCoords_SN-'+self.field+'.dat'
            if filt != 'all' and band not in filt:
                if verbose: print('filter %s not in filter list for image file %s'%(band,filt,imfile))
                #print 'filter %s,%s not in filter list for image file %s'%(band,filt,imfile)
                continue
            imfile,noisefile,psffile = os.path.join(self.rootdir,imfile),\
                os.path.join(self.rootdir,noisefile),os.path.join(self.rootdir,psffile)
            print imfile
            if not os.path.exists(imfile):
                if not os.path.exists(imfile+'.fz'):
                    print('Error : file %s does not exist'%imfile)
                    continue
                    print('Error : file %s does not exist'%imfile)
                    raise exceptions.RuntimeError('Error : file %s does not exist'%imfile)
                else:
                    os.system('/global/u1/d/dbrout/cfitsio/funpack %s.fz'%imfile)
            if not os.path.exists(noisefile):
                os.system('gunzip %s.gz'%noisefile)
                if not os.path.exists(noisefile):
                    os.system('/global/u1/d/dbrout/cfitsio/funpack %s.fz'%noisefile)
                    if not os.path.exists(noisefile):
                        raise exceptions.RuntimeError('Error : file %s does not exist'%noisefile)
            if not os.path.exists(psffile):
                if not os.path.exists(psffile+'.fz'):
                    raise exceptions.RuntimeError('Error : file %s does not exist'%psffile)
                else:
                    os.system('/global/u1/d/dbrout/cfitsio/funpack %s.fz'%psffile)

            if not nomask:
                maskfile = os.path.join(self.rootdir,snparams.image_name_search[j])
                #if not os.path.exists(maskfile):
                #    os.system('gunzip %s.gz'%maskfile)
                #    if not os.path.exists(maskfile):
                #        raise exceptions.RuntimeError('Error : file %s does not exist'%maskfile)
 
            # read in the files
            fakeim = ''.join(imfile.split('.')[:-1])+'+fakeSN.fits'
            if not os.path.exists(fakeim):
                os.system('/global/u1/d/dbrout/cfitsio/funpack %s.fz'%fakeim)
                os.system('/global/u1/d/dbrout/cfitsio/gunzip %s.gz'%fakeim)
            if self.usefake:
                try:
                    im = pyfits.getdata(fakeim)
                    #hdulist = pyfits.open(fakeim)
                except:
                    print fakeim+' is EMPTY, skipping star...'
                    #print 'stopped'
                    #raw_input()
                    continue
                hdr = pyfits.getheader(imfile)
            else:
                im = pyfits.getdata(imfile)
                #hdulist = pyfits.open(imfile)
                hdr = pyfits.getheader(imfile)
            fakeim_hdr = pyfits.getheader(fakeim)
            snparams.cat_zpts[imfile] = fakeim_hdr['HIERARCH DOFAKE_ZP']
            snparams.platescale = hdr['PIXSCAL1']
            snparams.airmass = hdr['AIRMASS']
            #snparams.snr = hdr['SNR']
            #print hdulist[1].data.shape
            #raw_input()
            noise = pyfits.getdata(noisefile)
            #noise[noise < 1e-5] = 0.
            print noise.shape
            print np.max(noise)
            print np.min(noise)
            print np.mean(noise)
            #raw_input()
            #raw_input()
            psf = pyfits.getdata(psffile)

            if params.weight_type.lower() == 'ivar':
                print 'ivar'
                raw_input()
                noise = np.sqrt(1/noise)
            elif params.weight_type.lower() != 'noise':
                raise exceptions.RuntimeError('Error : WEIGHT_TYPE value %s is not a valid option'%params.WEIGHT_TYPE)
            if nomask:
                mask = np.zeros(np.shape(noise))
                maskcols = np.where((noise < 0) |
                                    (np.isfinite(noise) == False))
                mask[maskcols] = 100.0
            else:
                mask = pyfits.getdata(maskfile)

            print noise.shape
            print np.max(noise)
            print np.min(noise)
            print np.mean(noise)
            #raw_input()

            #wcs = astWCS.WCS(imfile)
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

            #ra1,dec1 = wcs.pix2wcs(0,0)
            #ra2,dec2 = wcs.pix2wcs(snparams.nxpix-1,
            #                       snparams.nypix-1)
            
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
                print wcsinfo
                #print 'wcsinfo'
                #raw_input()
                #print results
                #print results[0]
                #print results[1]
                ra1, dec1 = results[0]*radtodeg, results[1]*radtodeg
                results2 =  wcsinfo.tran([[snparams.nxpix-1], [snparams.nypix-1]])
                ra2, dec2 =results2[0]*radtodeg, results2[1]*radtodeg

            ra_high = np.max([ra1,ra2])
            ra_low = np.min([ra1,ra2])
            dec_high = np.max([dec1,dec2])
            dec_low = np.min([dec1,dec2])

            #print ra_high*radtodeg,ra_low*radtodeg
            #print snparams.ra
            #print 'lowhitest'
            try:
                if self.exactpos:
                    snparams.RA = float(snparams.fake_ra)
                    snparams.DECL = float(snparams.fake_dec)
                else:
                    snparams.RA = float(snparams.ra)
                    snparams.DECL = float(snparams.decl)
            except:
                try:
                    if self.exactpos:
                        snparams.RA = astCoords.hms2decimal(snparams.fake_ra,':')
                        snparams.DECL = astCoords.dms2decimal(snparams.fake_dec,':')
                    else:
                        snparams.RA = astCoords.hms2decimal(snparams.ra,':')
                        snparams.DECL = astCoords.dms2decimal(snparams.decl,':')

                except:
                    raise exceptions.RuntimeError('Error : RA/Dec format unrecognized!!')


            #xsn,ysn = wcs.wcs2pix(snparams.RA,snparams.DECL)
            #if wcsworked:

            if params.forceradec.lower() == 'true':
                print float(params.fra),float(params.fdec)
                print 'fake ra dec'
                #raw_input()
                xsn,ysn = zip(*w.wcs_world2pix(np.array([[float(params.fra),float(params.fdec)]]), 0))
            else:
                xsn,ysn = zip(*w.wcs_world2pix(np.array([[snparams.RA,snparams.DECL]]), 0))

            xsn = xsn[0]
            ysn = ysn[0]
            testoffset = False
            offsetx = .6
            offsety = .3
            if testoffset:
                xsn += offsetx
                ysn += offsety

            #else:
            #    print snparams.RA
            #    print snparams.DECL
            #    print 'GOGOGOGOG'
            #    rr = wcsinfo.tran([[snparams.RA],[snparams.DECL]],False)
            #    print rr
            #    xsn = rr[0]
            #    ysn = rr[1]
            #    print xsn
            #    print ysn
            #    print 'xsnysn'
            #    raw_input()



            #xsn-=200
            #ysn-=200
            #print snparams.starcat[j]
            #raw_input()

            print xsn
            print ysn
            #print snparams.npix
            #raw_input()
            if xsn < 0 or ysn < 0 or xsn > snparams.nxpix-1 or ysn > snparams.nypix-1:
                #raise exceptions.RuntimeError("Error : SN Coordinates %s,%s are not within image"%(snparams.ra,snparams.decl))
                print "Error : SN Coordinates %s,%s are not within image"%(snparams.ra,snparams.decl)
                badflag = 1

            if type(snparams.starcat) == np.array:
                if os.path.exists(snparams.starcat[j]):
                    starcat = txtobj(snparams.starcat[j],useloadtxt=True)
                    if not starcat.__dict__.has_key('mag'):
                        try:
                            starcat.mag = starcat.__dict__[band]
                            starcat.dmag = starcat.__dict__['d%s'%band]
                        except:
                            raise exceptions.RuntimeError('Error : catalog file %s has no mag column!!'%snparams.starcat[j])
                else: 
                    raise exceptions.RuntimeError('Error : catalog file %s does not exist!!'%snparams.starcat[j])
            elif type(snparams.starcat) == dict and 'des' in snfile:
                starcatfile = None
                starcatloc = '/'.join(imfile.split('/')[0:-1])+'/'
                #print 'starcatfile',starcatfile
                #print 'file should be here',os.listdir(starcatloc)
                for fl in os.listdir(starcatloc):
                    print fl
                    if 'STARCAT' in fl:
                        starcatfile = fl
                print starcatloc+starcatfile
                if os.path.exists(starcatloc+starcatfile):
                    starcat = txtobj(starcatloc+starcatfile,useloadtxt=True, des=True)
                    if not starcat.__dict__.has_key('mag_%s'%band):
                        try:
                            print starcat.__dict__
                            starcat.mag = starcat.__dict__[band]
                            starcat.dmag = starcat.__dict__['d%s'%band]
                        except:
                            raise exceptions.RuntimeError('Error : catalog file %s has no mag column!!'%snparams.starcat[band])
                else:
                    raise exceptions.RuntimeError('Error : catalog file %s does not exist!!'%snparams.starcat[band])
            else:
                if os.path.exists(snparams.starcat[j]):
                    starcat = txtobj(snparams.starcat[j],useloadtxt=True)
                    if not starcat.__dict__.has_key('mag'):
                        try:
                            starcat.mag = starcat.__dict__[band]
                            starcat.dmag = starcat.__dict__['d%s'%band]
                        except:
                            print snparams.starcat
                            raise exceptions.RuntimeError('Error : catalog file %s has no mag column!!'%snparams.starcat[j])

                else: 
                    raise exceptions.RuntimeError('Error : catalog file %s does not exist!!'%snparams.starcat[j])
                    
            #print snparams.psf_model.lower()
            #raw_input()
            if snparams.psf_model.lower() == 'daophot':
                if params.build_psf == 'yes':
                    cols = np.where((starcat.ra > ra_low) & 
                                    (starcat.ra < ra_high) & 
                                    (starcat.dec > dec_low) & 
                                    (starcat.dec < dec_high))[0]
                    if not len(cols):
                        raise exceptions.RuntimeError("Error : No stars in image!!")
                    
                    mag_star = starcat.mag[cols]
                    #x_star,y_star = wcs.wcs2pix(starcat.ra[cols],starcat.dec[cols])
                    x_star,y_star = zip(*w.wcs_world2pix(np.array(zip(starcat.ra[cols],starcat.dec[cols])),0))
                    if not dontcentroid:
                        x_star,y_star = cntrd.cntrd(im,x_star,y_star,params.cntrd_fwhm)
                        newra,newdec = zip(*w.wcs_pix2world(np.array(zip(xstar_,y_star)),0))
                    #print starcat.ra[cols] - newra
                    #raw_input()
                    mag,magerr,flux,fluxerr,sky,skyerr,badflag,outstr = \
                        aper.aper(im,x_star,y_star,apr = params.fitrad)

                    if badflag == 1:
                        print 'aper1 badflag'
                        #raw_input()

                    self.rdnoise = hdr[params.rdnoise_name]
                    self.gain = hdr[params.gain_name]
                    if not os.path.exists(psffile) or params.clobber_psf == 'yes':
                        gauss,psf,magzpt = getpsf.getpsf(im,x_star,y_star,mag,sky,
                                                         hdr[params.rdnoise_name],hdr[params.gain_name],
                                                         range(len(x_star)),params.fitrad,
                                                         psffile)
                        hpsf = pyfits.getheader(psffile)
                        #self.gauss = gauss
                    else:
                        print('PSF file exists.  Not clobbering...')
                        hpsf = pyfits.getheader(psffile)
                        magzpt = hpsf['PSFMAG']

                        #self.gauss = [hpsf['GAUSS1'],hpsf['GAUSS2'],hpsf['GAUSS3'],hpsf['GAUSS4'],hpsf['GAUSS5']]
                elif nozpt:
                    self.rdnoise = hdr[params.rdnoise_name]
                    self.gain = hdr[params.gain_name] #1

                    cols = np.where((starcat.ra > ra_low) & 
                                    (starcat.ra < ra_high) & 
                                    (starcat.dec > dec_low) & 
                                    (starcat.dec < dec_high))[0]

                    if not len(cols):
                        raise exceptions.RuntimeError("Error : No stars in image!!")
                    
                    mag_star = starcat.mag[cols]
                    #coords = wcs.wcs2pix(starcat.ra[cols],starcat.dec[cols])
                    coords = zip(*w.wcs_world2pix(np.array(zip(starcat.ra[cols],starcat.dec[cols])),0))
                    x_star,y_star = [],[]
                    for c in coords:
                        x_star += [c[0]]
                        y_star += [c[1]]
                    x_star,y_star = np.array(x_star),np.array(y_star)
                    if not doncentroid:
                        x_star,y_star = cntrd.cntrd(im,x_star,y_star,params.cntrd_fwhm)
                        newra,newdec = zip(*w.wcs_pix2world(np.array(zip(xstar_,y_star)),0))
                        
                    #print starcat.ra[cols] - newra
                    #print 'hehehehehehe'
                    #raw_input()
                    mag,magerr,flux,fluxerr,sky,skyerr,badflag,outstr = \
                        aper.aper(im,x_star,y_star,apr = params.fitrad)
                    if badflag == 1:
                        print 'aper2 badflag'
                        #raw_input()
                    hpsf = pyfits.getheader(psffile)
                    magzpt = hpsf['PSFMAG']

                    #self.gauss = [hpsf['GAUSS1'],hpsf['GAUSS2'],hpsf['GAUSS3'],hpsf['GAUSS4'],hpsf['GAUSS5']]
                else:

                    hpsf = pyfits.getheader(psffile)
                    magzpt = hpsf['PSFMAG']
                    #print hpsf
                    #print magzpt
                    #raw_input()
                    ##self.gauss = [hpsf['GAUSS1'],hpsf['GAUSS2'],hpsf['GAUSS3'],hpsf['GAUSS4'],hpsf['GAUSS5']]
                    self.rdnoise = hdr[params.rdnoise_name]
                    self.gain = hdr[params.gain_name]


                #fwhm = 2.355*self.gauss[3]

            # begin taking PSF stamps

            if snparams.psf_model.lower() == 'psfex':
                #hpsf = pyfits.getheader(psffile)
                #print hpsf
                #print hpsf['PSF_FWHM']
                #raw_input()
                hdulist = fits.open(psffile)
                hdulist.info()
                print hdulist[1].header['PSF_FWHM']
                psf_fwhm = hdulist[1].header['PSF_FWHM']
                self.psf, self.psfcenter= self.build_psfex(psffile,xsn,ysn,imfile)
                self.psf = self.psf/np.sum(self.psf)
                print 'snpsfcenter',self.psfcenter
                #fitpsf = self.get_fwhm_of_2d_psf(self.psf)
                #self.psf *= self.sector_mask(self.psf.shape,(self.psf.shape[0]/2.,self.psf.shape[1]/2.),2.*fitpsf/self.snparams.platescale)
            elif snparams.psf_model.lower() == 'daophot':
                self.psf = rdpsf.rdpsf(psffile)[0]/10.**(0.4*(25.-magzpt))
            else:
                raise exceptions.RuntimeError("Error : PSF_MODEL not recognized!")

            self.rdnoise = hdr[params.rdnoise_name]
            self.gain = hdr[params.gain_name]

            if not nozpt:
                try:
                    if doglobalstar:
                        zpt_file = imfile.split('.')[-2] + '_'+str(filt)+'band_dillonzptinfo_globalstar.npz'
                        #print zpt_file
                        #raw_input()
                    else:
                        zpt_file = imfile.split('.')[-2] + '_'+str(filt)+'band_dillonzptinfo.npz'
                    zptdata = np.load(zpt_file) #load previous zpt information
                    #zpt = zptdata['mcmc_me_zpt'] #set zpt to mcmc zpt and continue
                    #print zpt
                    #zpt_err = zptdata['mcmc_me_zpt_std']
                    zpt = zptdata['mpfit_zpt']
                    zpterr = zptdata['mpfit_zpt_std']
                    mjdoff = zptdata['mjdoff']
                    mjdslopeinteroff = zptdata['mjdslopeinteroff']
                except:
                    print('Warning : IMAGE_ZPT field does not exist!  Calculating')
                    nozpt = True
            if nozpt:
                self.rdnoise = hdr[params.rdnoise_name]
                self.gain =  hdr[params.gain_name]

                cols = np.where((starcat.ra > ra_low) & 
                                (starcat.ra < ra_high) & 
                                (starcat.dec > dec_low) & 
                                (starcat.dec < dec_high))[0]

                #checkstars = cols[:params.numcheckstars]
                #usecols = cols[params.numheckstars:]

                if not len(cols):
                    raise exceptions.RuntimeError("Error : No stars in image!!")
                try:
                    if band.lower() == 'g':
                        mag_star = starcat.mag_g[cols]
                    elif band.lower() == 'r':
                        mag_star = starcat.mag_r[cols]
                    elif band.lower() == 'i':
                        mag_star = starcat.mag_i[cols]
                    elif band.lower() == 'z':
                        mag_star = starcat.mag_z[cols]
                    else:
                        raise Exception("Throwing all instances where mag_%band fails to mag. Should not appear to user.")
                except:
                    mag_star = starcat.mag[cols]
                #coords = wcs.wcs2pix(starcat.ra[cols],starcat.dec[cols])



                if wcsworked:
                    x_starold,y_starold = [],[]
                    x_star,y_star = [],[]
                    coords = zip(*w.wcs_world2pix(np.array(zip(starcat.ra[cols],starcat.dec[cols])),0))

                    for xval,yval in zip(*coords):
                        x_starold += [xval]
                        y_starold += [yval]
                    if doglobalstar:
                        for ide in starcat.objid[cols]:
                            #print starglobalras
                            tra = starglobalras[starglobalids == ide]
                            tdec = starglobaldecs[starglobalids == ide]
                            
                            coords = zip(*w.wcs_world2pix(np.array(zip(tra,tdec)),0))
                            for xval,yval in zip(*coords):
                                x_star += [xval]
                                y_star += [yval]
                    else:
                        x_star = x_starold
                        y_star = y_starold
                    #print np.array(x_starold) - np.array(x_star)
                    #print 'checking new star positions'
                    #raw_input()
                else:
                    coords = wcsinfo.tran([starcat.ra[cols]/radtodeg,starcat.dec[cols]/radtodeg],False)
                    #print coords
                    #print coords[0]
                    #print len(coords)
                    #print len(starcat.ra[cols])
                    #print len(coords[0])
                    #print coords.shape
                    #print 'cooooords'
                    #raw_input()
                #print starcat.ra[cols][starcat.ra[cols] , starcat.dec[cols]
                #print coords
                #print 'stop'
                #raw_input()

                #x_star,y_star = [],[]
                #for xval,yval in zip(*coords):
                #    x_star += [xval]
                #    y_star += [yval]
                
                #print x_star
                #raw_input()
                #print y_star
                #raw_input()

                #print len(coords)
                #print coords[0][16],coords[1][16]
                #print x_star[16],y_star[16]
                #print 'stop'
                #raw_input()

                x_star1,y_star1 = np.array(x_star),np.array(y_star)
                mag,magerr,flux,fluxerr,sky,skyerr,badflag,outstr = \
                    aper.aper(im,x_star1,y_star1,apr = params.fitrad)
                #I REMOVED CENTROIDING BECAUSE WE NOW FIND A GLOBAL RA AND DEC FOR THE STAR SIMILARLY TO THE SN
                newx_star,newy_star = cntrd.cntrd(im,x_star1,y_star1,params.cntrd_fwhm)

                if wcsworked:
                    newra,newdec = zip(*w.wcs_pix2world(np.array(zip(newx_star,newy_star)),0))
                else:
                    ccc = wcsinfo.tran([x_star,y_star])
                    newra,newdec = ccc[0]*radtodeg,ccc[1]*radtodeg

                
                #print x_star,y_star
                #print 'yoyoyo'
                #raw_input()
                
                catra,catdec = self.getProperCatRaDec(starcat.ra[cols],starcat.dec[cols])
                #print starcat.ra[cols]-catra
                #print 'stopped'
                #raw_input()
                deltara = catra - newra
                deltadec = catdec - newdec
                deltamjd = copy(deltara)*0. + snparams.mjd[j]

                jjj = abs(deltara) < 1
                deltaram = deltara[jjj]
                print np.array(newra), deltara, jjj
                ram = np.array(newra)[jjj]
                mm, s, iii = meanclip.meanclip( deltaram, clipsig = 3., maxiter = 8, returnSubs=True)
                bra = np.polyfit(ram[iii], deltaram[iii], 0)
                mra,cra = np.polyfit(ram[iii], deltaram[iii], 1)

                jjj = abs(deltadec) < 1
                deltadecm = deltadec[jjj]
                decm = np.array(newdec)[jjj]
                mm, s, iii = meanclip.meanclip( deltadecm, clipsig = 3., maxiter = 8, returnSubs=True)
                bdec = np.polyfit(decm[iii], deltadecm[iii], 0)
                mdec,cdec = np.polyfit(decm[iii], deltadecm[iii], 1)

            
                mjdoff = [bra[0],bdec[0]]
                mjdslopeinteroff = [[mra,cra],[mdec,cdec]]
                
                '''badflagd = 0
                if dailyoff:
                    if bra[0] > .5:
                        badflagd = 1
                    if bdec[0] > .5:
                        badflagd = 1
                    print bra[0],bdec[0]
                    print xsn,ysn
                    xsn,ysn = zip(*w.wcs_world2pix(np.array([[snparams.RA-bra[0],snparams.DECL-bdec[0]]]), 0))
                    print xsn[0],ysn[0]
                    opsf = copy(self.psf)
                    opsfcenter = copy(self.psfcenter)
                    self.psf, self.psfcenter= self.build_psfex(psffile,xsn[0],ysn[0],imfile)
                    self.psf = self.psf/np.sum(self.psf)
                    print opsf - self.psf
                    print opsfcenter-self.psfcenter
                    raw_input()
                '''
                #print starcat.ra[cols]
                #print newra
                #print deltadec
                #print ' deltaaaaaaa'
                #raw_input()
                self.deltaras.extend(deltara)
                self.deltadecs.extend(deltadec)
                self.deltamjds.extend(deltamjd)
                self.ras.extend(catra)
                self.decs.extend(catdec)
                self.x_stars.extend(x_star)
                self.y_stars.extend(y_star)
                #print snparams.airmass
                #print 'airmasssss'
                #raw_input()
                self.airmasses.extend(starcat.ra[cols]*0. + round(snparams.airmass,2))
                #print self.airmasses
                #raw_input()
                #print len(self.deltaras)
                #raw_input()
                
                #mag,magerr,flux,fluxerr,sky,skyerr,badflag,outstr = \
                #    aper.aper(im,x_star,y_star,apr = params.fitrad)
                
                #if badflagd == 1:
                #    badflag = 1
                #if badflag == 1:
                #    print 'aper3 badflag'
                    #raw_input()
                #print 'ra,dec'
                #print starcat.ra[cols][22],starcat.dec[cols][22]
                #print x_star1[22],y_star1[22]
                #print x_star[22],y_star1[22]
                #raw_input()
                self.psf = self.psf/np.sum(self.psf)
                #print snparams.mjd[j]
                #raw_input()
                skipactualzeropoint = False
                if not skipactualzeropoint:
                    zpt,zpterr,zpt_file = self.getzpt(x_star,y_star,starcat.ra[cols], starcat.dec[cols],starcat,mag,sky,skyerr,snparams.mjd[j],
                                         badflag,mag_star,im,noise,mask,psffile,imfile,snparams,params.substamp,mjdoff,mjdslopeinteroff,
                                         psf=self.psf)    
                else:
                    if doglobalstar:
                        zpt_file = imfile.split('.')[-2] + '_'+str(filt)+'band_dillonzptinfo_globalstar.npz'
                    else:
                        zpt_file = imfile.split('.')[-2] + '_'+str(filt)+'band_dillonzptinfo.npz'

                    zptdata = np.load(zpt_file) #load previous zpt information                                                                                                                       
                    zpt = zptdata['mpfit_zpt']
                    zpterr = zptdata['mpfit_zpt_std']
                    mjdoff = zptdata['mjdoff']
                    mjdslopeinteroff = zptdata['mjdslopeinteroff']
                dotestoff = False
                if dotestoff:
                    self.teststarpos(self.rickfakestarfile,w,zpt,sky,skyerr,im,noise,mask,psffile,imfile,snparams,params.substamp,snparams.zp[j],psf=self.psf)
                fakestarflux,fakestarfluxerr,fakestarzpt = self.fake20magstars(self.rickfakestarfile,w,zpt,sky,skyerr,im,noise,mask,psffile,imfile,snparams,params.substamp,psf=self.psf)
                self.fakestarfluxes.extend(fakestarflux)
                self.fakestarfluxerrs.extend(fakestarfluxerr)
                self.fakestarzpts.extend(fakestarzpt)
            #idlmjdnearest = min(idlmjd, key=lambda x:abs(x-float(snparams.mjd[j])))
            #zpt = idlzpt[idlmjd == idlmjdnearest]

            if not ('firstzpt' in locals()): firstzpt = 31. ####firstzpt = zpt
            if self.usediffimzpt:
                scalefactor = 10**(-.4*(snparams.zp[j]-firstzpt))
            else:
                if zpt != 0.0 and np.min(self.psf) > -10000:
                    scalefactor = 10.**(-0.4*(zpt-firstzpt))
                if zpt == 0.:
                    badflag = 1
                    scalefactor = 0.
            print scalefactor
            im *= scalefactor
            im[np.where(mask != 0)] =-999999.0
            #im[np.where(mask != 0)] =np.median(im.ravel())
            '''if snparams.fake_truemag[j] < 90.:
                print "FAKEMAG"
                print snparams.fake_truemag[j]
                print badflag
                print fwhm_arcsec < params.fwhm_max 
                print np.min(im[ysn-2:ysn+3,xsn-2:xsn+3]) != np.max(im[ysn-2:ysn+3,xsn-2:xsn+3])
                print len(np.where(mask[ysn-25:ysn+26,xsn-25:xsn+26] != 0)[0]) < params.max_masknum
                print np.max(psf_stamp[params.substamp/2+1-3:params.substamp/2+1+4,params.substamp/2+1-3:params.substamp/2+1+4]) == np.max(psf_stamp[:,:])
                raw_input()
            '''
            print 'heeeeee'
            print xsn
            print ysn
            print snparams.nxpix-25
            print snparams.nypix-25
            print scalefactor

            #mjdoff = [mjdoff[0],bdec[0]]
            #mjdslopeinteroff = [[mra,cra],[mdec,cdec]]
            badflagd = 0
            if dailyoff:
                if mjdoff[0] > .5:
                    badflagd = 1
                if mjdoff[1] > .5:
                    badflagd = 1
                print mjdoff[0],mjdoff[1]
                print xsn,ysn
                xsn,ysn = zip(*w.wcs_world2pix(np.array([[snparams.RA+mjdoff[0],snparams.DECL+mjdoff[1]]]), 0))
                print xsn[0],ysn[0]
                xsn = xsn[0]
                ysn = ysn[0]
                opsf = copy(self.psf)
                opsfcenter = copy(self.psfcenter)
                self.psf, self.psfcenter= self.build_psfex(psffile,xsn,ysn,imfile)
                self.psf = self.psf/np.sum(self.psf)
                #print np.sum(opsf - self.psf)
                #print opsfcenter,self.psfcenter
                #raw_input()



            if xsn > 25 and ysn > 25 and xsn < snparams.nxpix-25 and ysn < snparams.nypix-25 and np.isfinite(scalefactor):
                index += 1
                radius1 = 5
                radius2 = 8 
                fff = float(snparams.psf[j])
                skyrad=[radius1*fff,radius2*fff]

                magsn,magerrsn,fluxsn,fluxerrsn,skysn,skyerrsn,badflag,outstr = aper.aper(im,xsn,ysn,apr = params.fitrad)#,skyrad=skyrad)
                print skyerrsn,np.sqrt(skysn/self.gain)
                mygain = ((1/skyerrsn**2)*skysn)
                #print mygain,hdr['GAINA'],hdr['GAINB']
                #raw_input()
                if badflagd == 1:
                    badflag = 1
                if badflag == 1:
                    print 'aper4 badflag'
                    #raw_input()
                #skysn = idlsky[idlmjd == idlmjdnearest]

                if np.sum(mask[ysn-params.fitrad:ysn+params.fitrad+1,xsn-params.fitrad:xsn+params.fitrad+1]) != 0:
                    badflag = 1
                    print 'mask badflag'
                    #raw_input()
                if skysn < -1e5:
                    badflag = 1
                    print 'skysn badflag'
                    #raw_input()
                if not badflag:
                    print noise.shape
                    print np.max(noise)
                    print np.min(noise)
                    print np.mean(noise)
                    #raw_input()

                    #SKY SIMGA
                    #save_fits_image(im[xsn-200:xsn+200,ysn-200:ysn+200],'./testskysig/'+str(snparams.mjd[i])+'.fits')
                    print xsn,ysn
                    print im.shape
                    stampsize = 256
                    print 'hahhahah'
                    #raw_input()
                    if ysn-stampsize < 0:
                        ylow = 0
                    else:
                        ylow = ysn-stampsize
                    if ysn+stampsize>im.shape[0]:
                        yhi = im.shape[0]-1
                    else:
                        yhi = ysn+stampsize
                    if xsn-stampsize < 0:
                        xlow = 0
                    else:
                        xlow = xsn-stampsize
                    if xsn+stampsize>im.shape[1]:
                        xhi = im.shape[1]-1
                    else:
                        xhi = xsn+stampsize
                    print 'yoyoyo'
                    #raw_input()
                    print ylow,yhi
                    print xlow,xhi,xsn
                    print np.median(im)
                    #raw_input()sigma_clip.meanclip( skyvals, clipsig = 3, maxiter = 8 ) 
                    mean,st,vals = sigma_clip.meanclip(im[ylow:yhi,xlow:xhi],clipsig = 4, maxiter = 8)
                    skysig=1.48*np.median(abs(vals-np.median(vals)))
                    mysky = np.median(vals)
                    mygain = (np.sqrt(mysky)/(skysig))**2
                    mygainsn =  (np.sqrt(skysn)/(skyerrsn))**2
                    print mygain,mygainsn,hdr['GAINA'],hdr['GAINB']
                    #import runsextractor
                    #print im
                    sexsky,sexrms = runsextractor.getsky_and_skyerr(imfile,xlow,xhi,ylow,yhi)
                    print mysky,skysig,skysn,skyerrsn
                    print sexsky,sexrms,snparams.skysig[j]
                    #print 'hahahahahahahaha'
                    #raw_input()
                    # for sss in im[ylow:yhi,xlow:xhi].ravel():
                    #     print sss
                    # print snparams.mjd[j]
                    # raw_input()
                    skyvals = im[ylow:yhi,xlow:xhi].ravel()
                    #m, s, clipped_skyvals = sigma_clip.meanclip( skyvals, clipsig = 3, maxiter = 8 )
                    #Freedman Diaconis binwidth
                    #IQR = np.subtract(*np.percentile(clipped_skyvals, [75, 25]))
                    #bw = (max(clipped_skyvals)-min(clipped_skyvals)) / 2 * IQR / len(clipped_skyvals)**(1./3.)
                    '''bw = 3.
                    hist, bin_edges = np.histogram(clipped_skyvals,bins=np.arange(min(clipped_skyvals),max(clipped_skyvals),bw))
                    bin_means = (bin_edges[:-1]+bin_edges[1:])/2.
                    print hist.shape
                    print 'bin_means'
                    print bin_means[100:250]
                    print hist[100:250]
                    ww = (bin_means > skysn - 4*skysig) & (bin_means < skysn +4*skysig)
                    bin_means = bin_means[ww]
                    hist = hist[ww]
                    guess = [1,skysn,skysig]
                    xpopt, xpcov = curve_fit(gaussian, bin_means,hist,p0=guess)


                    gauss_sky = xpopt[1]
                    gauss_skysig = xpopt[2]

                    print xpopt
                    xs = np.arange(min(clipped_skyvals),max(clipped_skyvals),.1)

                    thissky = np.median(im[ylow:yhi,xlow:xhi])
                    plt.clf()
                    plt.hist(clipped_skyvals,bins=np.arange(min(clipped_skyvals),max(clipped_skyvals),bw))
                    plt.plot(bin_means,hist,color='green')
                    plt.plot(xs,gaussian(xs, xpopt[0], xpopt[1], xpopt[2]),color='black')
                    plt.axvline(xpopt[1],color='magenta')
                    plt.axvline(thissky,color='yellow')
                    plt.axvline(skysn,color='red')
                    plt.xlim(skysn-3*skysig,skysn+3*skysig)
                    #plt.xlim(gauss_sky-3*gauss_skysig,gauss_sky+3*gauss_skysig)
                    plt.xlabel('Pixel Value')
                    plt.savefig('backgroundtest.png')
                    print skysig
                    '''



                    #raw_input()
                    pk = pkfit_norecent_noise_smp.pkfit_class(im,self.psf,self.psfcenter,self.rdnoise,self.gain,noise,mask)
                    #pk = pkfit_norecent_noise_smp.pkfit_class(im,self.gauss,self.psf,self.rdnoise,self.gain,noise,mask)
                    try:
                        errmag,chi,niter,scale,iylo,iyhi,ixlo,ixhi,image_stamp,noise_stamp,mask_stamp,psf_stamp = pk.pkfit_norecent_noise_smp(1,xsn,ysn,skysn,skyerrsn,params.fitrad,returnStamps=True,stampsize=params.substamp)
                    except ValueError:
                        raise ValueError('SN too close to edge of CCD!')

                    print noise_stamp.shape
                    print np.max(noise_stamp)
                    print np.min(noise_stamp)
                    print np.mean(noise_stamp)
                    #raw_input()
                    #mcmc_scale, mcmc_scale_err = pk.pkfit_norecent_noise_smp(1,xsn,ysn,skysn,skyerrsn,self.params.fitrad,
                    ##                                                        mpfit_or_mcmc='mcmc',
                    #                                                        counts_guess=scale,model_errors=False,
                    #                                                        maxiter=10000,gewekenum=100,useskyerr=True)

                    msk = copy(image_stamp)
                    msk[msk!=0.] = 1
                    model = np.append(np.zeros(len(image_stamp.ravel())),scale)
                    newsub = int(image_stamp.shape[0])
                    stdev = np.zeros(len(model))
                    stdev[-1] = np.sqrt(model[-1])


                    '''result  = mcmc.metropolis_hastings(
                                                model=model,psfs=psf_stamp,
                                                data=image_stamp,weights=1./noise_stamp,
                                                substamp=newsub, stdev=stdev, sky=[skysn],
                                                model_errors=False,mask=msk,
                                                Nimage=1,maxiter=10000,mjd=[float(snparams.mjd[j])],gewekenum=500,
                                                skyerr=skyerrsn,useskyerr=True,shiftpsf=False)
                    #try:
                    model, uncertainty, history, sims, xhistory, yhistory = result.get_params()
                    mcmc_scale = model[-1]
                    mcmc_scale_err = uncertainty[-1]

                    print len(model)
                    '''
                    print "mag sn pkfit"
                    print 31 - 2.5*np.log10(model[-1])


                if snparams.psf_model.lower() == 'psfex':
                    fwhm = float(snparams.psf[j])
                if snparams.psf_unit.lower() == 'arcsec':
                    fwhm_arcsec = fwhm
                elif snparams.psf_unit.lower().startswith('sigma-pix') or snparams.psf_unit.lower().startswith('pix'):
                    print snparams.psf_model.lower()
                    fwhm_arcsec = fwhm*snparams.platescale
                    print imfile
                    print fwhm
                    print snparams.platescale
                    #raw_input()
                else:
                    raise exceptions.RuntimeError('Error : FWHM units not recognized!!')

                '''if snparams.fake_truemag[j] < 90.:
                    print "FAKEMAG2"
                    print snparams.fake_truemag[j]
                    print badflag
                    print fwhm_arcsec < params.fwhm_max 
                    print np.min(im[ysn-2:ysn+3,xsn-2:xsn+3]) != np.max(im[ysn-2:ysn+3,xsn-2:xsn+3])
                    print len(np.where(mask[ysn-25:ysn+26,xsn-25:xsn+26] != 0)[0]) < params.max_masknum
                    print np.max(psf_stamp[params.substamp/2+1-3:params.substamp/2+1+4,params.substamp/2+1-3:params.substamp/2+1+4]) == np.max(psf_stamp[:,:])
                    raw_input()
                '''
                
                if not badflag:
                    if not np.isfinite(skysig):
                        print 'infinite skysig'
                        badflag = 1
                    if skysig < 1:
                        print 'skysig less than one'
                        badflag = 1
                print badflag

                try:
                    fname = self.checkstarfile.split('.')[0]+'_deltaradec.npz'
                    self.deltastarsfile = fname
                    df = np.load(self.deltastarsfile)
                except:
                    print 'MJD Skipped, could not load deltastarfile'
                    badflag = 1

                badflags.append(badflag)
                if not badflag:
                    if fwhm_arcsec < params.fwhm_max:
                        if np.min(im[ysn-2:ysn+3,xsn-2:xsn+3]) != np.max(im[ysn-2:ysn+3,xsn-2:xsn+3]):
                            if len(np.where(mask[ysn-25:ysn+26,xsn-25:xsn+26] != 0)[0]) < params.max_masknum:
                                if np.max(psf_stamp[params.substamp/2+1-3:params.substamp/2+1+4,params.substamp/2+1-3:params.substamp/2+1+4]) == np.max(psf_stamp[:,:]):
                                    
                                    noise_stamp[noise_stamp > 0.] = 1

                                    noise_stamp[noise_stamp <= 0.] = 0

                                    smp_im[i,:,:] = image_stamp
                                    # ccr = -1
                                    # for sss in image_stamp.ravel():
                                    #     ccr += 1
                                    #     if msk.ravel()[ccr] == 1:
                                    #         print sss
                                    # print snparams.mjd[j]
                                    # raw_input()
                                    smp_noise[i,:,:] = noise_stamp*1/(skysig**2)
                                    print float(snparams.mjd[j])
                                    print skysig
                                    print 1/(skysig**2)
                                    #raw_input()
                                    #if float(snparams.mjd[j]) == 56536.243:
                                    #    raw_input() 
                                    smp_psf[i,:,:] = psf_stamp/np.sum(psf_stamp)

                                    c = 20
                                    psa = self.snparams.platescale

                                    rmat = np.zeros(smp_psf[i,:,:].shape)
                                    for a in range(0,40):
                                        for b in range(0,40):
                                            r = np.sqrt((a-20)**2+(b-20)**2)
                                            rmat[a,b] = r
                                    #print rmat[20,:]

                                    #Following formalism of 
                                    #http://das.sdss2.org/ge/sample/sdsssn/SNANA-PUBLIC/doc/snana_manual.pdf
                                    #page 8
                                    integral = 2*np.pi*np.trapz((smp_psf[i,:,:]**2*rmat).ravel())
                                    total_skyerr = np.sqrt(integral**(-1)*skysn)


                                    smp_dict['total_skyerr'][i] = total_skyerr

                                    smp_dict['scale'][i] = scale
                                    #smp_dict['mcmc_scale'][i] = mcmc_scale
                                    smp_dict['scale_err'][i] = errmag
                                    #smp_dict['mcmc_scale_err'][i] = mcmc_scale_err
                                    #smp_dict['sky'][i] = skysn
                                    #smp_dict['sky'][i] = mysky
                                    #smp_dict['skyerr'][i] = skyerrsn
                                    smp_dict['sky'][i] = sexsky
                                    smp_dict['skyerr'][i] = sexrms
                                    smp_dict['flag'][i] = 0
                                    #CHECK FOR DIFFIM FLAGS
                                    if (int(snparams.photflag[j]) & 1016) > 0:
                                        smp_dict['flag'][i] = 2
                                        smp_dict['scale'][i] = np.nan
                                        smp_dict['scale_err'][i] = np.nan
                                    smp_dict['zpt'][i] = zpt
                                    smp_dict['zpterr'][i] = zpterr
                                    smp_dict['mjd'][i] = float(snparams.mjd[j])
                                    smp_dict['mjdoff'].append( mjdoff )
                                    smp_dict['mjdslopeinteroff'].append(mjdslopeinteroff)
                                    smp_dict['image_scalefactor'][i] = scalefactor
                                    smp_dict['snx'][i] = xsn
                                    smp_dict['sny'][i] = ysn
                                    smp_dict['skysig'][i] = skysig
                                    smp_dict['imwcs'].append(w)
                                    msk = copy(image_stamp)
                                    msk[msk!=0.] = 1
                                    smp_dict['mask'].append(msk)
                                    smp_dict['fwhm_arcsec'][i] = fwhm_arcsec
                                    smp_dict['image_filename'][i] = imfile
                                    smp_dict['zpt_file'][i] = zpt_file
                                    smp_dict['psf_filename'][i] = psffile
                                    smp_dict['psf_fwhm'][i] = psf_fwhm
                                    smp_dict['fakepsf'][i] = snparams.psf[j]
                                    smp_dict['weight_filename'][i] = noisefile
                                    smp_dict['fakemag'][i] = snparams.fake_truemag[j]
                                    smp_dict['fakezpt'][i] = snparams.zp[j]
                                    smp_dict['diffim_flux'][i] = snparams.flux[j]
                                    smp_dict['diffim_fluxerr'][i] = snparams.fluxerr[j]
                                    smp_dict['id_obs'][i] = snparams.id_obs[j]
                                    fs = snparams.flux
                                    brightlimit = fs[np.argsort(fs)][::-1][:15]
                                    print brightlimit
                                    brightlimit = brightlimit[-1]
                                    print brightlimit
                                    print snparams.flux[j]
                                    #raw_input()
                                    if snparams.flux[j] > brightlimit:
                                        smp_dict['notbrightflag'][i] = 0
                                    else:
                                        smp_dict['notbrightflag'][i] = 1
                                    #print 'sb flux', snparams.hostgal_sb_fluxcal
                                    #raw_input()
                                    #smp_dict['hostgal_sbmag'][i] = 1


                                    #START HERE TOMORROW
                                    if not nozpt:
                                        fname = self.checkstarfile.split('.')[0]+'_deltaradec.npz'
                                        self.deltastarsfile = fname
                                        df = np.load(self.deltastarsfile)
                                        self.usedeltaras = np.array(df['deltaras'])
                                        self.usedeltadecs = np.array(df['deltadecs'])
                                        self.usemjds = np.array(df['mjds'])
                                        self.useras = np.array(df['ras'])
                                        self.usedecs = np.array(df['decs'])
                                        self.usexstar = np.array(df['x_star'])
                                        self.useystar = np.array(df['y_star'])
                                    else:    
                                        self.usedeltaras = np.array(copy(self.deltaras))
                                        self.usedeltadecs = np.array(copy(self.deltadecs))
                                        self.usemjds = np.array(copy(self.deltamjds))
                                        self.useras = np.array(copy(self.ras))
                                        self.usedecs = np.array(copy(self.decs))
                                        self.usexstar = np.array(copy(self.x_stars))
                                        self.useystar = np.array(copy(self.y_stars))
            
                                    srad = params.nearby_stars_pixel_rad
                                    rad = ((self.usexstar-xsn)**2+(self.useystar-ysn)**2)**.5

                                    nearbystars_onthisCCDepoch_indices = [(rad < srad) & (self.usemjds == float(snparams.mjd[j]))]
                                    nearby_xstar = self.usexstar[nearbystars_onthisCCDepoch_indices]
                                    nearby_ystar = self.useystar[nearbystars_onthisCCDepoch_indices]
                                    print xsn
                                    print nearby_xstar
                                    print self.usedeltaras[nearbystars_onthisCCDepoch_indices]
                                    print np.mean(self.usedeltaras[nearbystars_onthisCCDepoch_indices])
                                    print 'hshshshshshshshs'
                                    #raw_input()

                                    #NOW FIND NEARbY STARS AND GRAB OFFSETS AND APPLY TO SN (MAJE A NEW SMP_DICT ELEMENT)
                                        

                                    if filt == 'u':
                                        smp_dict['hostgal_mag'][i] = snparams.fake_hostmag_u
                                        smp_dict['hostgal_sbmag'][i] = -99
                                    if filt == 'g':
                                        smp_dict['hostgal_mag'][i] = snparams.fake_hostmag_g
                                        smp_dict['hostgal_sbmag'][i] = 27.5 - 2.5*np.log10(float(snparams.hostgal_sb_fluxcal[0]))
                                    if filt == 'r':
                                        smp_dict['hostgal_mag'][i] = snparams.fake_hostmag_r
                                        smp_dict['hostgal_sbmag'][i] = 27.5 - 2.5*np.log10(float(snparams.hostgal_sb_fluxcal[1]))
                                    if filt == 'i':
                                        smp_dict['hostgal_mag'][i] = snparams.fake_hostmag_i
                                        smp_dict['hostgal_sbmag'][i] = 27.5 - 2.5*np.log10(float(snparams.hostgal_sb_fluxcal[2]))
                                    if filt == 'z':
                                        smp_dict['hostgal_mag'][i] = snparams.fake_hostmag_z
                                        smp_dict['hostgal_sbmag'][i] = 27.5 - 2.5*np.log10(float(snparams.hostgal_sb_fluxcal[3]))
                                    if filt== 'y':
                                        smp_dict['hostgal_mag'][i] = snparams.fake_hostmag_y
                                        smp_dict['hostgal_sbmag'][i] = -99

                                    if smp_dict['mjd'][i] < snparams.peakmjd - params.mjdminus or \
                                        smp_dict['mjd'][i] > snparams.peakmjd + params.mjdplus:
                                        smp_dict['mjd_flag'][i] = 1
                                    if self.dogalfit:
                                        if smp_dict['mjd'][i] > snparams.peakmjd + params.mjdplus:
                                            smp_dict['fitflag'][i] = 0
                                        if smp_dict['mjd'][i] < snparams.peakmjd - params.mjdminus:
                                            smp_dict['fitflag'][i] = 0


                                    i += 1
        if mergeno == 0:
            zeroArray = np.zeros(smp_noise.shape)
            largeArray = zeroArray + 1E10
            #largeArray = copy(zeroArray)
            #smp_noise = np.fmin(smp_noise,largeArray)
            #smp_psfWeight = np.fmin(smp_psf,largeArray) 
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
        # Now all the images are in the arrays
        # Begin the fitting
        #badnoisecols = np.where(smp_noise <= 1)
        badnoisecols = np.where(smp_noise < 1e-5)
        #smp_noise[badnoisecols] = 1e10
        smp_noise[badnoisecols] = 0.
        badpsfcols = np.where(smp_psf < 0)
        #smp_noise[badpsfcols] = 1e10
        smp_noise[badpsfcols] = 0.0
        smp_psf[badpsfcols] = 0.0

#        badnoisecols = np.where(smp_bignoise <= 1)
#        smp_bignoise[badnoisecols] = 1e10
#        badpsfcols = np.where(smp_bigpsf < 0)
#        smp_bignoise[badpsfcols] = 1e10
#        smp_bigpsf[badpsfcols] = 0.0
        # data can't be sky subtracted with this cut in place
        infinitecols = np.where((smp_im == 0) | (np.isfinite(smp_im) == 0) | (np.isfinite(smp_noise) == 0))
        #smp_noise[infinitecols] = 1e10
        smp_noise[infinitecols] = 0.0
        #smp_noise[smp_noise > 1e9] = 0.0
        smp_im[infinitecols] = 0
        mpparams = np.concatenate((np.zeros(float(params.substamp)**2.),smp_dict['scale'],smp_dict['sky']))

        mpdict = [{'value':'','step':0,
                  'relstep':0,'fixed':0, 'xtol': 1E-15} for i in range(len(mpparams))]
        # provide an initial guess - CHECK
        #First Guess
        #maxcol = np.where(smp_im[0,:,:].reshape(params.substamp**2.) == np.max(smp_im[0,:,:]))[0][0]
        #mpparams[maxcol+1] = np.max(smp_im[0,:,:])/np.max(smp_psf[0,:,:])
        #End First Guess
        #badpsfweightcols = np.where(smp_psf == 0)
        #smp_psf_weight = np.copy(smp_psf)
        #smp_psf_weight[badpsfweightcols] = 1E10
        #mpparams[:params.substamp**2] = (smp_im[0,:,:]/smp_psf_weight[0,:,:]).flatten()
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
        #mpdict[1012]['value'] = 10**((31-19.033)/2.5)
        #mpdict[1012]['fixed'] = 1
        for col in range(int(params.substamp)**2+len(smp_dict['scale'])):
            #mpdict[col]['step']=np.max(smp_dict['scale'])
            mpdict[col]['step']=np.sqrt(np.max(smp_dict['scale']))
        #for i in range(len(mpparams)):
        #    mpdict[i]['xtol'] = (np.fmax(0.1, np.sqrt(mpdict[i]['value'])/10.0))  
        #Setting parameter values for all galaxy pixels with at least one valid psf and image pixel
        #mpdict[:]['value'] = mpparams[:]
        #Fixing parameter values for all galaxy pixels with no valid psf or galaxy pixel
        #mpdict[mpparams != mpparams]['value'] = 0.0
        #mpdict[mpparams != mpparams]['fixed'] = 1
        #mpdict[mpparams >1E307]['value'] = 0.0
        #mpdict[mpparams >1E307]['fixed'] = 1
        #Fixing parameter values for all epochs that were flagged
        for col in np.where((smp_dict['mjd_flag'] == 1) | (smp_dict['flag'] == 1))[0]+int(params.substamp)**2:
            print 'flagged '+str(col)
            mpdict[col]['fixed'] = 1
            mpdict[col]['value'] = 0
        #Setting parameter values for all good epochs

        #Setting step values for all parameters
        #Temporarily setting to zero to tell mpfit to calculate automatically.
        #mpdict[range(int(params.substamp)**2+len(smp_dict['scale']))]['step']=0#np.max(smp_dict['scale'])
        
        #Setting other arguments to scene and scene_check for mpfit
        #mpargs = {'x':smp_psf,'y':smp_im,'err':smp_noise,'params':params}

        #Setting final iteration tolerance of mpfit to sqrt(value)/10.0
        #mpdict[range(len(mpparams))]['xtol'] = (np.fmax(0.1, np.sqrt(mpdict[range(len(mpparams))]['value'])/10.0))  


        if verbose: print('Creating Initial Scene Model')
        #first_result = mpfit(scene,parinfo=mpdict,functkw=mpargs, debug = True, quiet=False)
        #print smp_psf.shape
        #print mpdict[:]['value']
        #raw_input()
        
        #make pixelized galaxy model from data outside mjd range (no supernova)
        #ww = np.where((smp_dict['mjd_flag'] == 1))[-1][1]
        

        arg = -1
        mn = 999999999
        for sky in smp_dict['sky']:
            arg += 1
            if smp_dict['flag'][arg] == 0:
                if smp_dict['mjd_flag'][arg] == 1:
                    if sky < mn:
                        if sky > 0.:
                            if np.any(smp_im[arg,:,:]):
                                mn = sky
                                usearg = arg
                        else:
                            print 'herehhere neggggggg'
                            smp_dict['flag'][arg] = 1

        #Make sure the psf is not zero
        for i in np.arange(len(smp_dict['sky'])):
            if np.max(smp_psf[i,:,:]) == np.min(smp_psf[i,:,:]):
                print 'hererererere psffsfsffsfsf',smp_dict['mjd'][i]
                smp_dict['flag'][i] = 1
                smp_dict['mjd'][i]


        #print badflags
        #raw_input()

        try:
            ww = usearg
        except:
            ww = 0
        #print np.where((smp_dict['mjd_flag'] == 0))
        print ww

        model, modelpeak = self.create_model(smp_im[ww,:,:]-smp_dict['sky'][ww],smp_psf[ww,:,:],smp_dict['scale'])

        #stdev = np.sqrt(np.append(smp_im[ww,:,:]-smp_dict['sky'][ww],smp_dict['scale']))
        
        #model = model*0.0
        stdev = np.sqrt(copy(model))
        #stdev[modelpeak] = np.sqrt(model[modelpeak])
        newsub = int(smp_im[ww].shape[0])

        modelvec = np.zeros(len(smp_dict['scale']))
        modelstd = np.zeros(len(smp_dict['scale']))

        for i,scale in enumerate(smp_dict['scale']):
            if i in np.where((smp_dict['mjd_flag'] == 1) | (smp_dict['flag'] == 1))[0]:
                model[newsub**2+i] = 0
                #stdev[newsub**2+i] = 0.
                modelvec[i] = 0.
                #modelstd[i] = 0.
            else:
                if scale < 25.:
                    #model[newsub**2+i] = scale/4.
                    model[newsub**2+i] = 0.
                    modelvec[i] = 0.
                    #stdev[newsub**2+i] = 0.
                else:
                    modelvec[i] = scale
                #modelstd[i] = np.sqrt(scale)

        #SET GALAXY MODEL HERE!!!! TO IMAGE - SKY
        galmodel = smp_im[ww,:,:]-smp_dict['sky'][ww]

        pkyerr = -2.5*np.log10(smp_dict['mcmc_scale']) + 2.5*np.log10(smp_dict['mcmc_scale'] + smp_dict['mcmc_scale_err'])

        





        outfolder = os.path.join(outfile,foldername)
        out = os.path.join(outfile,foldername+'/SNe/'+snparams.snfile.split('/')[-1].split('.')[0] + '/'+filt+'/')
        outimages = os.path.join(out,'image_stamps/')

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

        save_fits_image(model[:params.substamp**2].reshape(params.substamp,params.substamp),os.path.join(outimages,'initialmodel.fits'))

        imodel = copy(model)
        #self.gal_model = '/global/scratch2/sd/dbrout/smp_y1y2_pix/SNe/'+snparams.snfile.split('/')[-1].split('.')[0]+'/'+filt+'/image_stamps/finalmodel.fits'

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

        self.scaled_diffim_flux = scaled_diffim_flux
        self.scaled_diffim_fluxerr = scaled_diffim_fluxerr

        filename = snparams.snfile.split('/')[-1].split('.')[0] +'_'+ filt
        outdir = os.path.join(outfile,foldername+'/np_data/'+filt+'/')
        galaxyoutdir = os.path.join(outfile,galaxyfoldername+'/np_data/'+filt+'/')
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        maxiter = 1200
        print os.path.join(outdir,filename+'_mcmc_input.npz')
        #np.savez(os.path.join(outdir,filename+'_smpDict.npz'),**smp_dict)
        print outimages
        print filename
        filename = snparams.snfile.split('/')[-1].split('.')[0] +'_'+ filt
        lightcurves = os.path.join(outfile,foldername+'/lightcurves/'+filt+'/')  
        if not os.path.exists(lightcurves):
            os.makedirs(lightcurves)  


        for i,scale in enumerate(smp_dict['scale']):
            if i in np.where((smp_dict['flag'] == 1))[0]:
                pass
            else:
                pass
                #smp_psf[i] = self.degrade_psf_to_fake(copy(smp_psf[i]),copy(smp_dict['fakepsf'][i]))
                
    
        #print 'MJD','\t','\t','FIT_PSF','\t','FAKE_PSF'

        for i,scale in enumerate(smp_dict['scale']):
            if i in np.where((smp_dict['flag'] == 1))[0]:
                pass
            else:
                #fitpsf = self.get_fwhm_of_2d_psf(smp_psf[i])
                #smp_psf[i,:,:] *= self.sector_mask(smp_psf[i,:,:].shape,(smp_psf[i,:,:].shape[0]/2.,smp_psf[i,:,:].shape[1]/2.),2.*fitpsf/self.snparams.platescale)
                
                fakepsf = smp_dict['fakepsf'][i]
                #print smp_dict['mjd'][i],'\t',round(fitpsf,2),'\t','\t',fakepsf
                #print 'fit-fake psf', fitpsf-fakepsf
                #raw_input()

        print smp_dict['image_filename'][-1]        
        print 'MJD','\t','BAND','\t','FIT_ZPT','\t','FAKE_ZPT','\t','PSF','\t','SKY','\t','Skyerr','\t','Skysig','\t','IMAGE_FILENAME',''
        psfs = []
        for i,scale in enumerate(smp_dict['scale']):
            if i in np.where((smp_dict['flag'] == 1))[0]:
                psfs.append(-9.)
                print smp_dict['mjd'][i],'\t','FLAGGED'
            else:
                fitzpt = smp_dict['zpt'][i]
                #CuT PSF OFF AT 2*FWHM
                fakezpt = smp_dict['fakezpt'][i]
                psfs.append(round(self.get_fwhm_of_2d_psf(smp_psf[i]),2))
                print smp_dict['mjd'][i],'\t',filt,round(fitzpt,2),'\t','\t',round(fakezpt,2),'\t',round(self.get_fwhm_of_2d_psf(smp_psf[i]),2),round(smp_dict['sky'][i],2),round(smp_dict['skyerr'][i],2),round(smp_dict['skysig'][i],2),smp_dict['image_filename'][i]
                #if abs(fitzpt - fakezpt) > .025:
                #    smp_dict['fitflag'][i] = 1
                #if smp_dict['mjd'][i] < snparams.peakmjd +100:
                #    smp_dict['flag'][i] = 1
                if smp_dict['skyerr'][i] < .001:
                    smp_dict['flag'][i] = 1
                if smp_dict['skyerr'][i] > 1000.:
                    smp_dict['flag'][i] = 1


        zptnpz = os.path.join(outdir,filename+'_imagezpts.npz')

        np.savez(zptnpz
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
        #save_fits_image(galmodel,'./test/initalgalmodel.fits')
        
        fname = self.checkstarfile.split('.')[0]+'_deltaradec.npz'
        self.deltastarsfile = fname
        
        if nozpt:
            np.savez(self.deltastarsfile,deltaras=self.deltaras,deltadecs=self.deltadecs,mjds=self.deltamjds,ras=self.ras,decs=self.decs,airmasses=self.airmasses,x_star=self.x_stars,y_star=self.y_stars)
            print fname,'SAVED'
        else:
            dsf  = np.load(self.deltastarsfile)
            self.deltaras = dsf['deltaras']
            self.deltadecs = dsf['deltadecs']
            self.deltamjds = dsf['mjds']
            self.ras = dsf['ras']
            self.decs = dsf['decs']
            self.airmasses = dsf['airmasses']
        fname = self.checkstarfile.split('.')[0]+'_20magfakes.npz'
        self.moneyfile = fname
        if nozpt:
            np.savez(self.moneyfile,flux=self.fakestarfluxes,fluxerr=self.fakestarfluxerrs,zpt=self.fakestarzpts)
            print fname,'SAVED'

        #self.plotcheckstars()
        #self.plotallfake20staroffsets()

        print 'skyerr',smp_dict['skyerr']
        print 'flag',smp_dict['flag']
        print 'notbrightflag',smp_dict['notbrightflag']
        #raw_input()
        print 'mjdflag',smp_dict['mjd_flag']
        print 'fitflag',smp_dict['fitflag']
        print snparams.peakmjd
        print smp_dict['mjd']
        print os.path.join(outdir,filename+'_mcmc_input.npz')

        print 'mjdoff',smp_dict['mjdoff']
        print 'mjdslopeinteroff',smp_dict['mjdslopeinteroff']

        print os.path.join(outdir,filename+'_mcmc_input.npz')

        np.savez( os.path.join(outdir,filename+'_mcmc_input.npz'), 
                galmodel = galmodel
                , modelvec = modelvec*0.
                , galstd = np.sqrt(galmodel)*2.
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
                lcout = lightcurves+filename,
                model_data_index = ww,
                mjdoff = smp_dict['mjdoff'],
                mjdslopeinteroff = smp_dict['mjdslopeinteroff'],
                star_offset_file = star_offset_file,
                zpt_files = smp_dict['zpt_file'],
                starglobalids = starglobalids,
                globalraoffsets = offsetra,
                globaldecoffsets = offsetdec
                )
        
        np.savez(os.path.join(outdir,filename+'_smpDict.npz'),**smp_dict)
        #raw_input()

        if self.dogalfit:
            aaa = mcmc3.metropolis_hastings( 
                    galmodel = galmodel
                    , modelvec = modelvec*0.
                    , galstd = np.sqrt(galmodel)/3
                    , modelstd = modelvec*0.
                    , data = smp_im
                    , psfs = smp_psf
                    , weights = smp_noise
                    , substamp = params.substamp
                    , Nimage = len(smp_dict['sky'])
                    , maxiter = 200000
                    , mask = None
                    , sky=smp_dict['sky']
                    , mjd=smp_dict['mjd']
                    , gewekenum=9999999
                    , skyerr=smp_dict['skyerr']
                    , useskyerr = True
                    , flags = smp_dict['flag']
                    , fitflags = smp_dict['fitflag']
                    , psf_shift_std = .0
                    , shiftpsf = False
                    , fileappend = ''
                    , stop = False
                    , skyerr_radius = 16.
                    , outpath = outimages
                    , compressionfactor = 100
                    , fix_gal_model = None
                    , pixelate_model = 1.
                    , burnin = .75
                    , lcout = lightcurves+filename
                    , chainsnpz = os.path.join(outdir,filename+'_nosn.npz')
                    )



            modelvec, modelvec_uncertainty, galmodel_params, galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, sims, xhistory,yhistory,accepted_history,pix_stamp,chisqhist,redchisqhist  = aaa.get_params()

            print 'TOTAL Galfit SMP TIME ',time.time()-tstart

            print os.path.join(outdir,filename+'_nosn.npz')



        self.outdir = outdir
        self.galaxyoutdir = galaxyoutdir
        self.filename = filename
        self.foldername = foldername
        self.filt = filt
        self.smp_dict = smp_dict
        self.smp_im = smp_im
        self.smp_psf = smp_psf
        self.smp_noise = smp_noise

        if self.dosnradecfit:
            if not self.dogalfit:
                chains = np.load(os.path.join(galaxyoutdir,filename+'_nosn.npz'))
                galmodel_params = chains['galmodel_params']
                galmodel_uncertainty = chains['galmodel_uncertainty']
            galmodel = galmodel_params
            galstd = np.sqrt(abs(galmodel))/5.
            tstart = time.time()
            modelvec = scaled_diffim_flux
            modelstd = scaled_diffim_fluxerr/5.
            if not self.floatallepochs:
                modelvec[smp_dict['mjd_flag'] == 1] = 0
                modelstd[smp_dict['mjd_flag'] == 1] = 0

            fixgal = True
            print 'fitting SN RA DEC'
            aaa = mcmc3.metropolis_hastings(
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
                , useskyerr = False
                , usesimerr = True
                , flags = smp_dict['flag']+smp_dict['mjd_flag']
                , fitflags = smp_dict['fitflag']*0.
                , psf_shift_std = .001
                , xoff = 0.
                , yoff = 0.
                , shiftpsf = True
                , fileappend = ''
                , stop = False
                , skyerr_radius = 16.
                , outpath = outimages
                , compressionfactor = 100
                , fix_gal_model = fixgal
                , pixelate_model = 1.
                , burnin = .75
                , lcout = lightcurves+filename
                , chainsnpz = os.path.join(outdir,filename+'_withSngetRADEC.npz')
                , platescale = .27
                , mjdoff = smp_dict['mjdoff']
                , fitradec = True
                )
            modelvec, modelvec_uncertainty, galmodel_params, galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, sims, xhistory,yhistory,accepted_history,pix_stamp,chisqhist,redchisqhist  = aaa.get_params()

            xoff = np.mean(xhistory[int(3*len(xhistory)/4.):])/.27
            yoff = np.mean(yhistory[int(3*len(yhistory)/4.):])/.27
            print xoff,yoff
            print 'xoff,yoff'
            #raw_input()
        if self.dosnfit:
            if not self.dogalfit:
                chains = np.load(os.path.join(galaxyoutdir,filename+'_nosn.npz'))
                galmodel_params = chains['galmodel_params']
                galmodel_uncertainty = chains['galmodel_uncertainty']
            if not self.dosnradecfit:
                try:
                    print os.path.join(outdir,filename+'_withSngetRADEC.npz')
                    chains = np.load(os.path.join(outdir,filename+'_withSngetRADEC.npz'))
                    xhistory = chains['xhistory']
                    yhistory = chains['yhistory']
                    decoff = chains['decoff']
                    xoff = np.mean(xhistory[int(3*len(xhistory)/4.):])/.27
                    yoff = np.mean(yhistory[int(3*len(yhistory)/4.):])/.27
                    modelvec = chains['modelvec']
                    modelstd = chains['modelvec_uncertainty']/5.
                except:
                    print 'could not find ra dec file, setting to zero...'
                    #sys.exit()
                    xoff = 0.
                    yoff = 0.
                    modelvec = scaled_diffim_flux
                    modelstd = scaled_diffim_fluxerr/7.
            galmodel = galmodel_params
            galstd = np.sqrt(abs(galmodel))/10.
            tstart = time.time()
            #if not self.dosnradecfit:
            #    modelvec = scaled_diffim_flux
            #    modelstd = scaled_diffim_fluxerr/5.
            if not self.floatallepochs:
                modelvec[smp_dict['mjd_flag'] == 1] = 0
                modelstd[smp_dict['mjd_flag'] == 1] = 0

            print 'snfit modelvec',modelvec
            print 'snfit modelstd',modelstd
            print 'galmodelshape', galmodel.shape
            print 'xoff,yoff',xoff,yoff
            #print 'decoff',decoff
            #if self.fixgalzero:
            #    galmodel = galmodel*0.
            #    galstd = galstd*0.
            #    fixgal = True
            #else:
            #    fixgal = False
            fixgal = False
            #print smp_dict['sky']
            #print np.sqrt(smp_dict['sky']/4.)
            #print smp_dict['skyerr']
            #print 'skycheck',self.gain
            
            aaa = mcmc3.metropolis_hastings( 
                    galmodel = galmodel
                    , modelvec = modelvec
                    , galstd = galstd
                    , modelstd = modelstd
                    , data = smp_im
                    , psfs = smp_psf
                    , weights = smp_noise
                    , substamp = params.substamp
                    , Nimage = len(smp_dict['sky'])
                    , maxiter = 300000
                    , mask = None
                    , sky=smp_dict['sky']
                    , mjd=smp_dict['mjd']
                    , gewekenum=9999999
                    , skyerr=smp_dict['skyerr']
                    , useskyerr = False
                    , usesimerr = True
                    , flags = smp_dict['flag']
                    , fitflags = smp_dict['fitflag']*0.
                    , psf_shift_std = None#.0001
                    , xoff = xoff
                    , yoff = yoff
                    , shiftpsf = False
                    , fileappend = ''
                    , stop = False
                    , skyerr_radius = 16.
                    , outpath = outimages
                    , compressionfactor = 100
                    , fix_gal_model = fixgal
                    , pixelate_model = 1.
                    , burnin = .75
                    , lcout = lightcurves+filename
                    , chainsnpz = os.path.join(outdir,filename+'_withSn.npz')
                    , mjdoff = smp_dict['mjdoff']
                    )
            print 'modelvec before',modelvec
            modelveco = copy(modelvec)
            
            modelvec, modelvec_uncertainty, galmodel_params, galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, sims, xhistory,yhistory,accepted_history,pix_stamp,chisqhist,redchisqhist  = aaa.get_params()
            print 'modelvec after',modelvec
            print 'modelvec after-before',modelvec-modelveco
            #print modelvec
            #raw_input()
            print modelvec_nphistory[:,-5]
            #raw_input()
            print 'TOTAL SMP SN TIME ',time.time()-tstart

            #np.savez(os.path.join(outdir,filename+'_withSn.npz'),modelvec=modelvec, modelvec_uncertainty=modelvec_uncertainty, galmodel_params=galmodel_params, galmodel_uncertainty=galmodel_uncertainty, modelvec_nphistory=modelvec_nphistory, galmodel_nphistory=galmodel_nphistory, sims=sims,data=smp_im,accepted_history=accepted_history,chisqhist=chisqhist,redchisqhist=redchisqhist)
            print os.path.join(outdir,filename+'_withSn.npz')
            print smp_dict['image_filename']
            print smp_dict['image_filename'].shape
            #raw_input()
            #raw_input()

        if self.dogalsimfit:
            #fixmodelvec = self.afterfit(self.snparams,self.params,donesn=False)
            #print fixmodelvec.shape
            #print modelvec.shape
            #raw_input()
            if not self.dogalfit:
                chains = np.load(os.path.join(galaxyoutdir,filename+'_nosn.npz'))
                galmodel_params = chains['galmodel_params']
                galmodel_uncertainty = chains['galmodel_uncertainty']
                
            galmodel = galmodel_params
            #galstd = galmodel_uncertainty
            galstd = np.sqrt(abs(galmodel))/5.
            #print galmodel[10:20,10:20]
            #print galstd[10:20,10:20]
            #print modelvec/2.
            #raw_input()
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
            modelvec[smp_dict['mjd_flag'] == 1] = 0
            modelstd[smp_dict['mjd_flag'] == 1] = 0

            print modelstd
            print 'galmodelshape', galmodel.shape
            aaa = mcmc3galsim.metropolis_hastings( 
                    galmodel = galmodel
                    , modelvec = fixmodelvec
                    , galstd = galstd
                    , modelstd = modelstd
                    , data = smp_im
                    , psfs = smp_psf
                    , weights = smp_noise
                    , substamp = params.substamp
                    , Nimage = len(smp_dict['sky'])
                    , maxiter = 1000
                    , mask = None
                    , sky=smp_dict['sky']
                    , mjd=smp_dict['mjd']
                    , gewekenum=9999999
                    , skyerr=smp_dict['skyerr']
                    , useskyerr = True
                    , flags = smp_dict['flag']
                    , fitflags = smp_dict['fitflag']*0.
                    , psf_shift_std = .0005
                    , shiftpsf = True
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
                    , burnin = 100
                    , model_pixel_scale = .27
                    , lcout = lightcurves+filename
                    , chainsnpz = os.path.join(outdir,filename+'_withSnAndGalsim.npz') 
                    )

            print 'modelvec before',modelvec
            modelveco = copy(modelvec)
            
            modelvec, modelvec_uncertainty, galmodel_params, galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, sims, xhistory,yhistory,accepted_history,pix_stamp,chisqhist  = aaa.get_params()
            print 'modelvec after',modelvec
            print 'modelvec after-before',modelvec-modelveco
            #print modelvec
            #raw_input()
            print modelvec_nphistory[:,-5]
            #raw_input()
            print 'TOTAL SMP SN TIME ',time.time()-tstart

            np.savez(os.path.join(outdir,filename+'_withSnAndGalsim.npz'),modelvec=modelvec, modelvec_uncertainty=modelvec_uncertainty, galmodel_params=galmodel_params, galmodel_uncertainty=galmodel_uncertainty, modelvec_nphistory=modelvec_nphistory, galmodel_nphistory=galmodel_nphistory, sims=sims,data=smp_im,accepted_history=accepted_history,chisqhist=chisqhist)
            print os.path.join(outdir,filename+'_withSnAndGalsim.npz')
        
        if self.dogalsimpixfit:    
            if not self.dogalfit:
                chains = np.load(os.path.join(galaxyoutdir,filename+'_nosn.npz'))
                galmodel_params = chains['galmodel_params']
                galmodel_uncertainty = chains['galmodel_uncertainty']
            if not pixstart == None:
                usedir = os.path.join(outfile,pixstart+'/np_data/'+filt+'/')
                chains = np.load(os.path.join(usedir,filename+'_nosn.npz'))
                modelvec = chains['modelvec']
                modelstd = scaled_diffim_fluxerr/5.
                galmodel = chains['galmodel_params']
            else:
                modelvec = scaled_diffim_flux
                modelstd = scaled_diffim_fluxerr/5.
                galmodel = galmodel_params
            #galstd = galmodel_uncertainty
            galstd = np.sqrt(abs(galmodel))/10.
            #print galmodel[10:20,10:20]
            #print galstd[10:20,10:20]
            #print modelvec/2.
            #raw_input()
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
            extraflag[-10:] = 1

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
                    , psf_shift_std = None#.00008
                    , shiftpsf = False
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
                    , model_pixel_scale = .6
                    , lcout = lightcurves+filename
                    , chainsnpz = os.path.join(outdir,filename+'_withSnAndGalsimPix.npz')
                    , platescale = .27
                    , snraoff = 0.
                    , sndecoff = 0.
                    )

            modelveco = copy(modelvec)
            
            modelvec, modelvec_uncertainty, galmodel_params, galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, sims, xhistory,yhistory,accepted_history,pix_stamp,chisqhist,rahistory,dechistory  = aaa.get_params()
            print 'modelvec after',modelvec
            print  'modelvec before',modelveco
            print 'modelvec after-before',modelvec-modelveco
            #print modelvec
            #raw_input()
            #print modelvec_nphistory[:,-5]
            #raw_input()
            print 'TOTAL SMP SN TIME ',time.time()-tstart

            #np.savez(os.path.join(outdir,filename+'_withSnAndGalsimPix.npz'),modelvec=modelvec, modelvec_uncertainty=modelvec_uncertainty, galmodel_params=galmodel_params, galmodel_uncertainty=galmodel_uncertainty, modelvec_nphistory=modelvec_nphistory, galmodel_nphistory=galmodel_nphistory, sims=sims,data=smp_im,accepted_history=accepted_history,chisqhist=chisqhist,snrahistory=rahistory,sndechistory=dechistory)
            print os.path.join(outdir,filename+'_withSnAndGalsimPix.npz')

        self.outdir = outdir
        self.filename = filename
        self.foldername = foldername
        self.filt = filt
        self.smp_dict = smp_dict
        self.smp_im = smp_im
        self.smp_psf = smp_psf
        self.smp_noise = smp_noise
        sys.exit()
        return

    def closest_node(self,ra,dec):
        tra = self.bigcatalogras*0. + ra
        tdec = self.bigcatalogdecs*0. + dec

        dist_2 = ((self.bigcatalogras-tra)**2+(self.bigcatalogdecs-tdec)**2)**.5
        return np.argmin(dist_2)


    def getProperCatRaDec(self,ra,dec):
        properra = np.zeros(len(ra))
        properdec = np.zeros(len(dec))
        #self.bigcatalog
        for i in np.arange(0,len(ra)):
            j = self.closest_node(ra[i],dec[i])
            properra[i] = self.bigcatalogras[j]
            properdec[i] = self.bigcatalogdecs[j]
        return properra,properdec
        


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
            plt.savefig(fname)
            print fname
        
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
        plt.savefig(fname)
        print fname

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
        plt.savefig(fname)
        print fname
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
        plt.savefig(fname)
        print fname
        
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
        plt.savefig(fname)
        print fname

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
        plt.savefig(fname)
        print fname

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
        plt.savefig(fname)
        print fname

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
                SNscale,SNscale_std,chisq,dms = self.getfluxsmp(im,psf,sky,weight,fitrad,galmodel_wosn,mjd,None)
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

        np.savez(os.path.join(self.outdir,snparams.snfile.split('/')[-1].split('.')[0]+'_'+self.filt+'_finalresults.npz'),
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
            np.savez(os.path.join(self.outdir,snparams.snfile.split('/')[-1].split('.')[0]+'_'+self.filt+'_nosnresults.npz'),
                fixedgal_scale = np.array(final_mcmc_fixfluxes),
                fixedgal_std = np.array(final_mcmc_fixstd)
            )
        print 'SMP3 Completed Successfully'
        try:
            print np.array(final_mcmc_fixfluxes)[:,0]
            rr = np.array(final_mcmc_fixfluxes)[:,0]
        except:
            rr = np.array(final_mcmc_fixfluxes)
        return rr



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

        lightcurves = os.path.join(outfile,foldername+'/lightcurves/'+filt+'/') 
        plt.savefig(lightcurves+filename+'_chisqlike.pdf')
        
        print lightcurves+filename+'_chisqlike.pdf'
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
    
        plt.savefig(lightcurves+filename+'_lightcurve.pdf')
        print lightcurves+filename+'_lightcurve.pdf'
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
        np.savez( fileout
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
        np.savez( fileout
                ,history = modelhistory
                ,substamp = substamp)

        #fileout = os.path.join(outdir,filename+'_galfluxhistory.npz')
        #np.savez( fileout
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

    def getfluxsmp(self,im,psf,sky,weight,fitrad,gal,mjd,guess_scale):

        chisqvec = []
        fluxvec = []
        
        galconv = scipy.signal.fftconvolve(gal,psf,mode='same')

        radius = 12
        substamp = galconv.shape[0]
        #Make a mask with radius
        fitrad = np.zeros([substamp,substamp])
        for x in np.arange(substamp):   
            for y in np.arange(substamp):
                if np.sqrt((substamp/2. - x)**2 + (substamp/2. - y)**2) < radius:
                    fitrad[int(x),int(y)] = 1.


        if guess_scale is None:
            for i in np.arange(-10000,200000,5):
                sim = galconv + sky + i*psf
                chisqvec.append(np.sum((im-sim)**2*weight*fitrad))
                fluxvec.append(i)
        else: 
            for i in np.arange(guess_scale-2000,guess_scale+2000,1):
                sim = galconv + sky + i*psf
                chisqvec.append(np.sum((im-sim)**2*weight*fitrad))
                fluxvec.append(i)

        ii = fitrad.ravel()
        i = ii[ii != 0]
        
        ndof = len(i) + 1

        fluxvec = np.array(fluxvec)
        chisqvec = np.array(chisqvec)
        hh = chisqvec*0 + min(chisqvec)
        mchisq = min(chisqvec)
        idx = np.isclose(chisqvec, hh, atol=1.)
        # print len(fluxvec[idx])
        #
        # plt.clf()
        # plt.plot(fluxvec,np.array(chisqvec),color='black')
        # plt.axhline(min(chisqvec)+1.0,color='red')
        # plt.axvline(fluxvec[idx][0])
        # plt.axvline(fluxvec[idx][-1])
        # plt.xlabel('Flux')
        # plt.ylabel('Chisq')
        # plt.ylim(min(chisqvec),min(chisqvec)+10.)
        # plt.xlim(-1000,10000)
        # plt.savefig('./chisqmin/testflux_'+str(mjd)+'.png')
        # print './chisqmin/testflux_'+str(mjd)+'.png'
        #raw_input()


        #if mjd > 
        
        #raw_input()
        sim = galconv + sky + fluxvec[chisqvec == min(chisqvec)]*psf
        sum_data_minus_sim = np.sum(im-sim)
        return fluxvec[chisqvec == min(chisqvec)], fluxvec[chisqvec == min(chisqvec)] - fluxvec[idx][0], mchisq/ndof, sum_data_minus_sim

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
            aper.aper(im,x_star1,y_star1,apr = params.fitrad)
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
                            psf, psfcenter = self.build_psfex(psffile,x,y,imfile)
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
                        flux_star[i,j] = copy(scale)
                        flux_star_std[i,j] = copy(cscale_std)

                    except:
                        continue


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
                        psf, psfcenter = self.build_psfex(psffile,x,y,imfile)
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



    def getzpt(self,xstar,ystar,ras, decs,starcat,mags,sky,skyerr,thismjd,
                badflag,mag_cat,im,noise,mask,psffile,imfile,snparams,substamp,
                mjdoff,mjdslopeinteroff,psf='',
                mpfit_or_mcmc='mpfit',cat_zpt=-999):
        """Measure the zeropoints for the images"""
        import pkfit_norecent_noise_smp
        #from PythonPhot import iterstat
        import astropy.io.fits as pyfits
        #from PythonPhot import pkfit_norecent_noise
        counter = 0

        flux_star = np.array([-999.]*len(xstar))        
        flux_star_std = np.array([-999.]*len(xstar))
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

        radius = 10.
        fitrad = np.zeros([substamp,substamp])
        for x in np.arange(substamp):   
            for y in np.arange(substamp):
                if np.sqrt((substamp/2. - x)**2 + (substamp/2. - y)**2) < radius:
                    fitrad[int(x),int(y)] = 1.
        print 'xstarrrrrrr',len(xstar)
        for x,y,m,s,se,mc,i in zip(xstar,ystar,mags,sky,skyerr,mag_cat,range(len(xstar))):
            #print 'xstar',xstar
            #raw_input()
            #y -= 2.
            if i < float(params.numcheckstars):
                isnotcheckstars[i] = 0

            if mc > 21:
                continue
            if x > 51 and y > 51 and x < self.snparams.nxpix-51 and y < self.snparams.nypix-51:
                if self.stardumppsf:
                    if self.snparams.psf_model.lower() == 'psfex':
                        psf, psfcenter = self.build_psfex(psffile,x,y,imfile)
                        print psf.shape
                    elif psf == '':
                        raise exceptions.RuntimeError("Error : PSF array is required!")
                else:
                    psf, psfcenter = self.psf, self.psfcenter
                    print 'psfcenter',psfcenter

                
                counter += 1
                pk = pkfit_norecent_noise_smp.pkfit_class(im,psf/np.sum(psf),psfcenter,self.rdnoise,self.gain,noise,mask)
                #Run for MPFIT

                try:
                    errmag, chi, niter, scale, iylo, iyhi, ixlo, ixhi, image_stamp, noise_stamp, mask_stamp, psf_stamp = \
                        pk.pkfit_norecent_noise_smp(1, x, y, s, se, params.fitrad, returnStamps=True,
                                                    stampsize=params.substamp)
                    noise_stamp[noise_stamp > 0.] = 1
                    noise_stamp[noise_stamp <= 0.] = 0
                    sexsky, sexrms = runsextractor.getsky_and_skyerr(imfile, ixlo, ixhi, iylo, iyhi)
                    # noise_stamp = noise_stamp*1/(se**2)
                    noise_stamp = noise_stamp * 1 / (sexrms ** 2)
                    gal = np.zeros(image_stamp.shape)
                    mjd = 000.
                    #oldcscale, cscale_std, chisq, dms = self.getfluxsmp(image_stamp, psf_stamp, s, noise_stamp, fitrad, gal,
                    #                                                    mjd, scale)
                    cscale, cscale_std, chisq, dms = self.getfluxsmp(image_stamp, psf_stamp, sexsky, noise_stamp, fitrad,
                                                                     gal, mjd, scale)
                    #print 'checking!!!', cscale, oldcscale
                    # print 'DIFFFFFF',scale,cscale
                    scale = cscale
                    #raw_input()
                except:
                    print 'skipped star...'
                    continue
                #print self.params.fitrad
                #resid = imstamp - psf*scale - s
                #print 
                print mc
                print 2.5*np.log10(scale)  
                print errmag 
                print 'aaa'
                '''
                if abs(mc +2.5*np.log10(scale) - 30.61) > 0.1:

                    #if i == 22:
                    #    save_fits_image(imstamp,imfile.split('.')[-2] + '_'+str(filt)+'band_starfit'+str(i)+'_im.fits')
                    #    print 'hhh'
                    #    raw_input()
                    #if i == 23:
                    #    save_fits_image(imstamp,imfile.split('.')[-2] + '_'+str(filt)+'band_starfit'+str(i)+'_im.fits')
                    #    print 'hhh'
                    #    raw_input()
                    plt.clf()
                    plt.imshow(resid)
                    plt.colorbar()
                    plt.title('Residual')
                    plt.savefig(imfile.split('.')[-2] + '_'+str(filt)+'band_starfit_'+str(i)+'_resid_br.png')
                    print imfile.split('.')[-2] + '_'+str(filt)+'band_starfit_'+str(i)+'_resid_br.png'
                    plt.clf()
                    plt.imshow(imstamp)
                    plt.colorbar()
                    plt.title('Image')
                    plt.savefig(imfile.split('.')[-2] + '_'+str(filt)+'band_starfit_'+str(i)+'_im_br.png')
                    print imfile.split('.')[-2] + '_'+str(filt)+'band_starfit_'+str(i)+'_im_b.png'
                    plt.clf()

                    plt.imshow(psf*scale+s)
                    plt.colorbar()
                    plt.title('Scale*PSF + sky')
                    plt.savefig(imfile.split('.')[-2] + '_'+str(filt)+'band_starfit_'+str(i)+'_psf_br.png')
                    print imfile.split('.')[-2] + '_'+str(filt)+'band_starfit_'+str(i)+'_psf_b.png'
                    print 'stopped'
                    #raw_input()
                '''
                #print psf.shape
                #print psfcenter
                #print x,y
                #print scale
                #raw_input()
                #raw_input()
                flux_star[i] = scale #write file mag,magerr,pkfitmag,pkfitmagerr and makeplots
                flux_star_std[i] = cscale_std
                #if abs(m + 2.5*np.log10(scale) - 30.61) > .05:
                #    raw_input()

                #print x
                #print y
                print scale
                #raw_input()

                #THIS IS THE MCMC... UNCOMMENT TO RUN
                #show = False
                #gain = 1.0
                '''if scale < 60000.:
                    # MCMC without Model Errors
                    #val, std = pk.pkfit_norecent_noise_smp(1,x,y,s,se,self.params.fitrad,mpfit_or_mcmc='mcmc',counts_guess=scale,show=show,gain=gain)
                    #flux_star_mcmc[i] = val
                    #flux_star_std_mcmc[i] = std
                    #mcmc_mag_std[i] = abs(2.5*np.log10(val)-2.5*np.log10(val+std))
                    
                    #MCMC With Model Errors
                    valb, std = pk.pkfit_norecent_noise_smp(1,x,y,s,se,self.params.fitrad,mpfit_or_mcmc='mcmc',counts_guess=scale,show=show,gain=gain,model_errors=True)
                    flux_star_mcmc_modelerrors[i] = valb
                    flux_star_std_mcmc_modelerrors[i] = std
                    mcmc_me_mag_std[i] = abs(2.5*np.log10(valb)-2.5*np.log10(valb+std))
                    
                    # Analytical simple scale=sum(pix)/sum(psf)
                    #valsimple = pk.pkfit_norecent_noise_smp(1,x,y,s,se,self.params.fitrad,mpfit_or_mcmc='mcmc',analytical='simple',counts_guess=scale,show=show,gain=gain,model_errors=True)
                    #flux_star_mcmc_me_simple[i] = valsimple

                    #valweighted = pk.pkfit_norecent_noise_smp(1,x,y,s,se,self.params.fitrad,mpfit_or_mcmc='mcmc',analytical='weighted',counts_guess=scale,show=show,gain=gain,model_errors=True)
                    #flux_star_mcmc_me_weighted[i] = valweighted
                else:
                    #flux_star_mcmc[i] = 0.0
                    flux_star_mcmc_modelerrors[i] = 0.0
                    #flux_star_std_mcmc[i] = 0.0
                    flux_star_std_mcmc_modelerrors[i] = 0.0
                    #flux_star_mcmc_me_simple[i] = 0.0
                    #flux_star_mcmc_me_simple[i] = 0.0
                    #flux_star_std_mcmc_me_simple[i] = 0.0
                    #flux_star_std_mcmc_me_simple[i] = 0.0
                    #flux_star_mcmc_me_weighted[i] = 0.0
                    #flux_star_mcmc_me_weighted[i] = 0.0
                    #flux_star_std_mcmc_me_weighted[i] = 0.0
                    #flux_star_std_mcmc_me_weighted[i] = 0.0
                
                
                '''
                ##Run for MCMC
                #errmag_mcmc,chi_mcmc,niter_mcmc,scale_mcmc = pk.pkfit_norecent_noise_smp(1,x,y,s,se,self.params.fitrad,mpfit_or_mcmc='mcmc')
                #flux_star_mcmc[i] = scale_mcmc

        badflag = badflag.reshape(np.shape(badflag)[0])
        
        #check for only good fits MPFIT        
        goodstarcols = np.where((mag_cat != 0) & 
                                (mag_cat < 21.5) &
                                (flux_star != 1) & 
                                (flux_star < 1e7) &
                                #(flux_star_mcmc < 1e7) &
                                #(flux_star_mcmc != 0) &
                                #(flux_star_mcmc_modelerrors != 0) &
                                #(flux_star_mcmc_modelerrors < 1e7) &
                                #(flux_star_std_mcmc > 1.0) &
                                #(flux_star_std_mcmc_modelerrors > 1.0) &
                                (np.isfinite(mag_cat)) &
                                (np.isfinite(flux_star)) &
                                (flux_star > 0) &
                                (badflag == 0) &
                                (isnotcheckstars == 1))[0]


        checkstarcols = np.where((mag_cat != 0) &
                                (mag_cat < 21.5) &
                                (flux_star != 1) &
                                (flux_star < 1e7) &
                                (np.isfinite(mag_cat)) &
                                (np.isfinite(flux_star)) &
                                (flux_star > 0) &
                                (badflag == 0) &
                                (isnotcheckstars == 0))[0]

        #NEED TO MAKE A PLOT HERE!
        if len(goodstarcols) > 10:
            md,std,num = self.iterstat(mag_cat[goodstarcols]+2.5*np.log10(flux_star[goodstarcols]),
                                       startMedian=True,sigmaclip=1.5,iter=10)
            
            print 'zpt',md
            print 'std',std

            dstd = 1.48*np.median(abs(mag_cat[goodstarcols]+2.5*np.log10(flux_star[goodstarcols])- np.ones(len(flux_star[goodstarcols]))*md))/np.sqrt(len(flux_star[goodstarcols]))
            std = float(std)/float(num**.5)
            print 'reduced std', std
            print 'dan std',dstd
            #raw_input()
            #mcmc_md, mcmc_std = self.weighted_avg_and_std(mag_cat[goodstarcols]+2.5*np.log10(flux_star_mcmc[goodstarcols]),1.0/(mcmc_mag_std[goodstarcols])**2)
            mcmc_md = -999.
            mcmc_std = -999.
            #mcmc_md,mcmc_std = iterstat.iterstat(mag_cat[goodstarcols]+2.5*np.log10(flux_star_mcmc[goodstarcols]),
            #                           startMedian=True,sigmaclip=3.0,iter=10)
            
            mcmc_me_md,mcmc_me_std = self.weighted_avg_and_std(mag_cat[goodstarcols]+2.5*np.log10(flux_star_mcmc_modelerrors[goodstarcols]),1.0/(mcmc_me_mag_std[goodstarcols])**2)

            #mcmc_me_md,mcmc_me_std = iterstat.iterstat(mag_cat[goodstarcols]+2.5*np.log10(flux_star_mcmc_modelerrors[goodstarcols]),
            #                           startMedian=True,sigmaclip=3.0,iter=10)
            zpt_plots_out = mag_compare_out = imfile.split('.')[-2] + '_zptPlots'
            exposure_num = imfile.split('/')[-1].split('_')[1]
            #print cat_zpt
            #raw_input()
            #self.make_zpt_plots(zpt_plots_out,goodstarcols,mag_cat,flux_star,flux_star_mcmc_modelerrors,mcmc_me_mag_std,md,mcmc_me_md,starcat,cat_zpt=float(cat_zpt))
            if nozpt:
                if os.path.isfile(self.big_zpt+'.txt'):
                    b = open(self.big_zpt+'.txt','a')
                else:
                    b = open(self.big_zpt+'.txt','w')
                    b.write('Exposure Num\tRA\tDEC\tCat Zpt\tMPFIT Zpt\tMPFIT Zpt Err\tMCMC Zpt\tMCMC Zpt Err\tMCMC Model Errors Zpt\tMCMC Model Errors Zpt Err\tCat Mag\tMP Fit Mag\tMCMC Fit Mag\tMCMC Model Errors Fit Mag\tMCMC Analytical Simple\tMCMC Analytical Weighted\n')
                for i in goodstarcols:
                    b.write(str(exposure_num)+'\t'+str(ras[i])+'\t'+str(decs[i])+'\t'+str(cat_zpt)+'\t'+str(md)+'\t'+str(std)\
                        +'\t'+str(mcmc_md)+'\t'+str(mcmc_std)+'\t'+str(mcmc_me_md)+'\t'+str(mcmc_me_std)+'\t'+str(mag_cat[i])\
                        +'\t'+str(-2.5*np.log10(flux_star[i]))+'\t'+str(-2.5*np.log10(flux_star_mcmc[i]))\
                        +'\t'+str(-2.5*np.log10(flux_star_mcmc_modelerrors[i]))\
                        +'\t'+str(-2.5*np.log10(flux_star_mcmc_me_simple[i]))
                        +'\t'+str(-2.5*np.log10(flux_star_mcmc_me_weighted[i]))
                        +'\n')
                b.close()

                if os.path.isfile(self.big_zpt+'.txt'):
                    b = open(self.checkstarfile,'a')
                else:
                    b = open(self.checkstarfile,'w')
                    b.write('Exposure Num\tMJD\tRA\tDEC\txstar\tystar\tCat Zpt\tMPFIT Zpt\tMPFIT Zpt Err\tFit Flux\tFit Flux Err\tCat Mag\n')

                for i in checkstarcols:
                    b.write(str(exposure_num)+'\t'+str(thismjd)+'\t'+str(ras[i])+'\t'+str(decs[i])+'\t'+str(xstar[i])+'\t'+str(ystar[i])+'\t'+str(cat_zpt)+'\t'+str(md)+'\t'+str(std)+'\t'+str(flux_star[i])+'\t'+str(flux_star_std[i])+'\t'+str(mag_cat[i])+'\n')

            '''i = rdcol.read('/global/homes/d/dbrout/smppro/idlzpttest.txt',1,2,' ')
            istarmags = np.array(i['starmags'])
            istarcats = np.array(i['catmags'])

            istarmags = istarmags[istarcats < 21.5]
            istarcats = istarcats[istarcats < 21.5]

            print istarmags
            print istarcats
            '''

            #r = open('resids.txt','a')


            hh = mag_cat[goodstarcols]+2.5*np.log10(flux_star[goodstarcols]) - np.ones(len(flux_star[goodstarcols]))*md
            hh = hh[abs(hh < .25)]

            plt.clf()
            #plt.hist([mag_cat[goodstarcols]+2.5*np.log10(flux_star[goodstarcols]) - np.ones(len(flux_star[goodstarcols]))*md,istarmags-istarcats+30.6198],bins=np.arange(-.25,.25,.04),label=['python','idl'])
            plt.hist(mag_cat[goodstarcols]+2.5*np.log10(flux_star[goodstarcols]) - np.ones(len(flux_star[goodstarcols]))*md,bins=np.arange(-.25,.25,.04),label='mean: '+str(np.mean(hh))+' std: '+str(np.std(hh)))
            plt.xlabel('cat mag + 2.5log10(flux) - zeropoint')
            plt.ylabel('counts')
            plt.xlim(-.25,.25)
            #plt.legend()
            plt.savefig(imfile.split('.')[-2] + '_'+str(filt)+'band_starfitresids1s.png')
            print imfile.split('.')[-2] + '_'+str(filt)+'band_starfitresids1s.png'
            #r.write(imfile.split('.')[-2] + '_'+str(filt)+'band_starfitresids1s.png\n')
            #r.close()

            plt.clf()
            #plt.hist([mag_cat[goodstarcols]+2.5*np.log10(flux_star[goodstarcols]) - np.ones(len(flux_star[goodstarcols]))*md,istarmags-istarcats+30.6198],bins=np.arange(-.25,.25,.04),label=['python','idl'])
            print 'scatter'
            plt.scatter(mag_cat[goodstarcols], -2.5*np.log10(flux_star[goodstarcols]))
            print 'plot'
            plt.plot([min(mag_cat[goodstarcols]),max(mag_cat[goodstarcols])],[min(mag_cat[goodstarcols]),max(mag_cat[goodstarcols])]-md,color='black')
            plt.xlabel('cat mag')
            plt.ylabel('-2.5log10(flux)')
            #plt.legend()
            print 'saving'
            print mag_cat[goodstarcols].shape
            plt.savefig(imfile.split('.')[-2] + '_'+str(filt)+'band_starfit_zptplot.png')
            print imfile.split('.')[-2] + '_'+str(filt)+'band_starfit_zptplot.png'
            #raw_input()

            '''print 'mean python', np.mean(hh)
            print 'mean idl ',np.mean(istarmags-istarcats+30.6198)
            print 'std python', np.std(hh)
            print 'std idl ',np.std(istarmags-istarcats+30.6198)
            '''
            #raw_input()
            #Writing mags out to file .zpt in same location as image
            print 'saving npz'
            if doglobalstar:
                mag_compare_out = imfile.split('.')[-2] + '_'+str(filt)+'band_dillonzptinfo_globalstar.npz'
            else:
                mag_compare_out = imfile.split('.')[-2] + '_'+str(filt)+'band_dillonzptinfo.npz'
            np.savez( mag_compare_out
                ,ra = 1
                )
            print 'ttttt'
            np.savez( mag_compare_out
                ,ra = ras[goodstarcols]
                ,dec = decs[goodstarcols]
                ,cat_mag = mag_cat[goodstarcols]
                ,mpfit_mag = -2.5*np.log10(flux_star[goodstarcols])
                ,mcmc_me_fit_mag = -2.5*np.log10(flux_star_mcmc_modelerrors[goodstarcols])
                ,mcmc_me_fit_mag_std = mcmc_me_mag_std[goodstarcols]
                ,mpfit_zpt = md
                ,mpfit_zpt_std = std
                ,mcmc_me_zpt = mcmc_me_md
                ,mcmc_me_zpt_std = mcmc_me_std
                ,cat_zpt = cat_zpt
                ,mjd = mjd
                ,mjdoff=mjdoff
                ,mjdslopeinteroff=mjdslopeinteroff
                )
            
            print 'done saving npz'
            print mag_compare_out
            #raw_input()
            #print 'cat,fit'
            #print mag_cat,mpfit_mag

            '''print 'mag cat'
            print mag_cat[goodstarcols]
            print 'mpfit mag'
            print -2.5*np.log10(flux_star[goodstarcols])
            print mcmc_me_md
            print cat_zpt
            print mag_compare_out
            print 'zeropoint', md
            '''
            print md
            print 'hereeeeeeeeeee'
            #raw_input()
        else:
            print len(goodstarcols)
            print len(checkstarcols)
            print isnotcheckstars
            print params.numcheckstars
            raise exceptions.RuntimeError('Error : not enough good stars to compute zeropoint!!!')

        if self.verbose:
            print('measured ZPT: %.3f +/- %.3f'%(md,std))

        return(md,std,mag_compare_out)

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

    def build_psfex(self, psffile,x,y,imfile,dogalsim=False):
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
        psf = os.popen("dump_psfex -inFile_psf %s -xpix %s -ypix %s -gridSize %s"%(psffile,x,y,
                                                                                   35)).readlines()
        #ix, iy, psfval = np.genfromtxt(psffile, usecols = (1,2,5), skip_footer = 4)
        xin = copy(x)
        yin = copy(y)
        readdata,readheader = False,True
        ix,iy,psfval = [],[],[]
        IMAGE_CORNERX = 0
        IMAGE_CORNERY = 0
        for line in psf:
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
        psfout = np.zeros((2*self.params.fitrad + 1,2*self.params.fitrad + 1))
        for x,y,p in zip(ix,iy,psfval):
            if x >= (35 - 2*self.params.fitrad -1)/2 and y >= (35 - 2*self.params.fitrad -1)/2 and x < (2*self.params.fitrad +1) and y < (2*self.params.fitrad + 1):
                psfout[y-(35 - 2*self.params.fitrad - 1)/2,x-(35 - 2*self.params.fitrad -1)/2] = p 
            #psfout[y,x] = p


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

if __name__ == "__main__":

    import sys,getopt
    # read in arguments and options
    try:
        if os.path.exists("./defaults/default.config"):
            args = open("./defaults/default.config", 'r').read().split()
        else:
            args = sys.argv[1:]
        
        print args
        opt,arg = getopt.getopt(
            args,"hs:p:r:f:o:m:v:i:d:s",
            longopts=["help","snfile=","params=","rootdir=",
                      "filter=","nomask","nodiff","nozpt", "outfile=",
                      "mergeno=", "loadzpt",
                      "debug","verbose","clearzpt",
                      "psf_model=","ismultiple",
                      "gal_model=","index=","diffimzpt","idlsky",
                      "dontgalfit","dontsnfit","dontgalsimfit","dontgalsimpixfit",
                      "fixgalzero","floatallepochs","dailyoff","snradecfit","dontglobalstar"])


        print opt
        print arg
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
                      "mergeno=", "loadzpt",
                      "debug","verbose","clearzpt",
                      "psf_model=","ismultiple",
                      "gal_model=","index=","diffimzpt","idlsky",
                      "dontgalfit","dontsnfit","dontgalsimfit","dontgalsimpixfit",
                      "fixgalzero","floatallepcohs","dailyoff","snradecfit","dontglobalstar"])


        print opt
        print arg
    except getopt.GetoptError as err:
        print "No command line arguments"


    verbose,nodiff,debug,clear_zpt,psf_model,root_dir,mergeno,loadzpt,ismultiple,dogalfit,dosnfit,dogalsimfit,dogalsimpixfit = False,False,False,False,False,False,False,False,False,True,True,True,True
    fixgalzero,floatallepochs = False,False
    dailyoff = False
    usediffimzpt = False
    useidlsky = False
    snradecfit = False
    doglobalstar = True
    index = None
    gal_model = None
    snfile,param_file,outfile,filt = '','','',''
    nomask,nozpt = 'none',False
    mergeno = 0 


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
        elif o in ["-m","--mergeno"]:
            mergeno = int(a)
        elif o in ["--loadzpt"]:
            loadzpt = True
        elif o == "--nomask":
            nomask = True
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
        elif o in ["--dontgalsimfit"]:
            dogalsimfit = False        
        elif o in ["--dontgalsimpixfit"]:
            dogalsimpixfit = False
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
            outfile = a
        elif o in ["-m","--mergeno"]:
            mergeno = int(a)
        elif o in ["--loadzpt"]:
            loadzpt = True
        elif o == "--nomask":
            nomask = True
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
        elif o in ["--dontgalsimfit"]:
            dogalsimfit = False
        elif o in ["--dontgalsimpixfit"]:
            dogalsimpixfit = False
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
        else:
            print "Warning: option", o, "with argument", a, "is not recognized"
        #elif o == "--clearzpt":
        #    clear_zpt = True

    if not index is None:
        if index == 'all':
            for iii in np.arange(16,50):

                a = open('./data/snfiles.txt','r')
                files = a.readlines()
                print 'files',files
                print 'index',index
                snfile = files[int(iii)].rstrip()
                a.close()
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

                #print snfile
                #raw_input()
                snparams = get_snfile(snfile, root_dir)
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
                    scenemodel.main(nodiff=nodiff,nozpt=nozpt,nomask=nomask,debug=debug,outfile=outfile
                                 ,verbose=verbose,clear_zpt=True, mergeno=mergeno,usefake=True,snfile=snfile,
                                 gal_model=gal_model,stardumppsf=True,dogalfit=dogalfit,dosnfit=dosnfit,
                                 dogalsimfit=dogalsimfit,dogalsimpixfit=dogalsimpixfit,dosnradecfit=snradecfit,
                                 usediffimzpt=usediffimzpt,useidlsky=useidlsky,fixgalzero=fixgalzero,floatallepochs=floatallepochs,
                                 dailyoff=dailyoff,doglobalstar=doglobalstar)
                    #scenemodel.afterfit(snparams,params,donesn=True)
                    print "SMP Finished!"
                except:
                    pass
            sys.exit()

        else:
            a = open('./data/snfiles.txt','r')
            files = a.readlines()
            print 'files',files
            print 'index',index
            print len(files)
            snfile = files[int(index)].rstrip()
            a.close()
            print 'Index '+str(index)
            print 'SN File '+snfile

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
    if not psf_model:
        print("psf_model not specified. Assuming psfex...")
        psf_model = 'psfex'

    #print snfile
    #raw_input()
    snparams = get_snfile(snfile, root_dir)
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

    scenemodel = smp(snparams,params,root_dir,psf_model)
    scenemodel.main(nodiff=nodiff,nozpt=nozpt,nomask=nomask,debug=debug,outfile=outfile
                     ,verbose=verbose,clear_zpt=True, mergeno=mergeno,usefake=True,snfile=snfile,
                     gal_model=gal_model,stardumppsf=True,dogalfit=dogalfit,dosnfit=dosnfit,
                     dogalsimfit=dogalsimfit,dogalsimpixfit=dogalsimpixfit,dosnradecfit=snradecfit,
                     usediffimzpt=usediffimzpt,useidlsky=useidlsky,fixgalzero=fixgalzero,floatallepochs=floatallepochs,
                     dailyoff=dailyoff,doglobalstar=doglobalstar)
    scenemodel.afterfit(snparams,params,donesn=True)
    print "SMP Finished!"
     
