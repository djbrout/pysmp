#!/usr/bin/env python
# Dillon Brout 3/10/2015
# dbrout@physics.upenn.edu

"""
Usage:
import mcmc
a = mcmc.metropolis_hastings( model, data, psfs, weights, substamp, , Nimage )
a.run_d_mc()

1D arrays (all of same size)
model                 : contains all model parameters

2D Stamps (all of same size) 
data                  : data stamps (1 for each epoch)
psfs                  : psf stamps (1 for each epoch)
weights               : uncertainty stamps (1 for each epoch)

Integers
substamp              : size of one edge of a stamp
Nimage                : Number of epochs


To do list: 
only test for convergence with supernova/star parameters
check autocorr with galsim_iter.py and figure out when to stop mcmc
figure out how to calculate mean and uncertainty properly
calculate covariance

Geweke is slow for lots of iter

GEWEKE COVARIANCE

"""




import numpy as np
import scipy.ndimage
import scipy.signal
import scipy.ndimage as nd
#from . import Matplot, flib
#from .utils import autocorr, autocov
from copy import copy
import pdb
from numpy import corrcoef, sum, log, arange
from numpy.random import rand
#from pylab import pcolor, show, colorbar, xticks, yticks
#import pylab as plt
import time
import pyfits as pf
import os
import math
import galsim
import matplotlib as m
m.use('Agg')
import gc

import matplotlib.pyplot as plt
#import pyfftw
import dilltools as dt
import multiprocessing
from multiprocessing import Process, Queue, current_process, cpu_count
ncpu = cpu_count()

#from pympler.tracker import SummaryTracker

class metropolis_hastings():

    def __init__(self
                , galmodel = None
                , modelvec = None
                , galstd = None
                , modelstd = None
                , data = None
                , psfs = None
                , weights = None
                , substamp = 0
                , Nimage = 1
                , maxiter = 100000
                , gain = 1.0
                , model_errors = False
                , readnoise = 5.
                , analytical = 'No'
                , mask = None
                , fix = None
                , sky=None
                , mjd=None
                , gewekenum=1000
                , skyerr=None
                , useskyerr = False
                , flags = None
                , fitflags = None
                , psf_shift_std = .0005
                , shiftpsf = False
                , fileappend = ''
                , stop = False
                , skyerr_radius = 16.
                , outpath = './'
                , compressionfactor = 1
                , fix_gal_model = False
                , pixelate_model = None
                , imagefiles = None
                , psffiles = None
                , weightfiles = None
                , snra = None
                , sndec = None
                , skysig = None
                , model_pixel_scale = .27
                , lcout = None
                , chainsnpz = None
                , platescale = .27
                , snraoff = 0.
                , sndecoff = 0.
                , fitradius = 15
                , isfermigrid = False
                , psfcenterx = None
                 ,psfcentery = None
                ):


        self.galmodel = galmodel
        self.galmodel = self.galmodel
        self.modelvec = modelvec
        self.galstd = galstd
        self.modelstd = modelstd
        self.isfermigrid=isfermigrid
        self.fermigrid = self.isfermigrid
        #self.galdeltas = copy(self.galstd)
        self.modeldeltas = copy(self.modelstd)

        #self.deltas = copy(self.stdev) #this vec will change for each iter
        self.substamp = substamp
        self.Nimage = Nimage
        self.maxiter = maxiter
        self.gain = gain
        self.model_errors = model_errors
        self.readnoise = readnoise
        self.sky = sky
        self.mjd = mjd
        self.flags = flags
        self.fitflags = fitflags
        self.imagefiles = imagefiles
        self.psffiles = psffiles
        self.weightfiles = weightfiles
        self.snra = snra
        self.sndec = sndec
        self.skysig = skysig
        self.model_pixel_scale = model_pixel_scale

        self.gewekenum = gewekenum
        self.fix_gal_model = fix_gal_model
        #self.skyerr = skyerr
        self.psf_shift_std = psf_shift_std
        self.current_x_offset = 0.
        self.current_y_offset = 0.
        self.compressioncounter = 0
        self.shiftpsf = shiftpsf
        self.stop = stop
        self.outpath = outpath
        self.compressionfactor = compressionfactor
        self.pixelate_model = pixelate_model

        self.lcout = lcout
        self.chainsnpz = chainsnpz

        self.useskyerr = useskyerr
        self.psfcenterx = psfcentery
        self.psfcentery = psfcenterx

        self.inmask = mask

        self.didtimeout = False
        #self.comboerr = True

        self.pixelation_factor = model_pixel_scale/platescale
        self.pixelation_factor == 1
        if not self.pixelation_factor == 1:
            self.galaxy_model = self.pixelate(galmodel,self.pixelation_factor)
        else:
            self.galaxy_model = galmodel
        #newfitrad = self.pixelate(fitradius,self.pixelation_factor)
        #self.galaxy_model = np.zeros((fitradius*2+1,fitradius*2+1)) + np.mean(galmodel)
        self.kicked_galaxy_model = copy(self.galaxy_model)
        self.kicked_galmodel = copy(self.galaxy_model)
        #print 'gmshape',self.galaxy_model.shape
        #raw_input()
        #self.galaxy_model = np.ones()
        #self.model_radius = int(np.floor(self.galaxy_model.shape[0]/2.))

        if self.isfermigrid:
            #print 'we have correct tempwriter'
            #raw_input()
            self.tmpwriter = dt.tmpwriter(tmp_subscript='snfit_', useifdh=True)
        else:
            self.tmpwriter = dt.tmpwriter(tmp_subscript=self.chainsnpz.split('/')[-1].split('.')[0])

        self.data = data
        #for d in range(Nimage):
        #    self.data[d,:,:] = self.data[d,:,:].T
        #self.psfs = psfs
        #self.original_psfs = copy(psfs)
        self.weights = weights
        self.fiducial_coord = galsim.CelestialCoord(self.snra * galsim.degrees, self.sndec * galsim.degrees)
        self.psfs = []
        #self.weights = []
        self.imagestamps = []
        #self.imagestampsformodel = []
        self.simstamps = []
        self.snobjs = []
        self.snoffsets = []
        self.snraoff = snraoff
        self.sndecoff = sndecoff
        self.fitradius = fitradius
        #radius = self.galmodel.shape[0]/2.


        self.baseim = galsim.fits.read(self.imagefiles[( self.flags == 0 ) & (self.modelstd == 0.)][0])
        print self.imagefiles[( self.flags == 0 ) & (self.modelstd == 0.)][0]
        print self.baseim
        #raw_input()

        for i in np.arange(self.Nimage):
            #print self.imagefiles[i]

            if self.fermigrid :
                imfile = self.imagefiles[i]
                psffile = self.psffiles[i]
                noisefile = self.weightfiles[i]
                ifdhls = os.popen('ifdh lss '+imfile).read()
                if (len(ifdhls) > 1):
                    if (int(ifdhls.split()[-1]) > 0) :
                        print 'Copying over',imfile
                        #####os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp '+imfile+' .').read()
                        #imfilel = copy(imfilel)
                        self.imagefiles[i] = imfile.split('/')[-1]
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
                        #####os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp '+noisefile+' .').read()
                        self.weightfiles[i] = noisefile.split('/')[-1]
                        weightsfile = noisefile
                        #print 'ifdh cp ' + psffile + ' .'
                        #####os.popen('IFDH_CP_MAXRETRIES=1; ifdh cp ' + psffile + ' .').read()
                        self.psffiles[i] = psffile.split('/')[-1]
                else:
                    self.imagefiles[i] = 'na'
            if self.imagefiles[i] == 'na':
                self.psfs.append(np.nan)
                #self.imagestamps.append(np.nan)
                self.simstamps.append(np.nan)
                self.snobjs.append(np.nan)
                self.snoffsets.append(np.nan)
                self.flags[i] = 1
                #self.snras.append(np.nan)
                #self.sndecs.append(np.nan)
            else:



                #addd astrometric fit for each night just for the galaxy model


                full_data_image = galsim.fits.read(self.imagefiles[i])
                #full_weights_image = galsim.fits.read(self.weightfiles[i])
                stamp_center = full_data_image.wcs.posToImage(self.fiducial_coord)
                cx = int(round(stamp_center.x))
                cy = int(round(stamp_center.y))
                bigim_pix_center = galsim.PositionD(cx,cy)
                des_psfex = galsim.des.DES_PSFEx(self.psffiles[i])
                thispsf = des_psfex.getPSF(stamp_center)

                im = self.baseim[galsim.BoundsI( self.psfcenterx[i] - substamp / 2.,self.psfcenterx[i] + substamp / 2. -1,
                                                          self.psfcentery[i] - substamp / 2.,self.psfcentery[i] + substamp / 2. -1)]

                #self.data[i,:,:] = im.array
                #im = full_data_image[ galsim.BoundsI( cx-self.fitradius,cx+self.fitradius-1,
                #                                      cy-self.fitradius,cy+self.fitradius-1 ) ]
                
                self.psfs.append(im.wcs.toWorld(thispsf,image_pos=stamp_center))
                #self.imagestamps.append(im)
                #print np.median(im.array),self.sky[i], np.median(self.data[i])
                # if self.galaxy_model.shape[0] % 2 == 0:
                #     modelim = full_data_image[galsim.BoundsI( cx-self.model_radius,cx+self.model_radius-1,
                #                                               cy-self.model_radius,cy+self.model_radius-1 )]
                # else:
                if self.flags[i]  == 0:
                    if self.modelstd[i] == 0:
                        self.modelim = full_data_image[galsim.BoundsI( self.psfcenterx[i] - substamp / 2.,self.psfcenterx[i] + substamp / 2. - 1,
                                                          self.psfcentery[i] - substamp / 2.,self.psfcentery[i] + substamp / 2. -1 )]*0.0

                #[galsim.BoundsI( cx-self.fitradius,cx+self.fitradius-1,
                #                                                cy-self.fitradius,cy+self.fitradius-1 ) ]*0.
                #self.imagestampsformodel.append(modelim)

                #print self.modelvec
                #print i
                #raw_input()
                # self.simstamps.append(full_data_image[ galsim.BoundsI( cx-fitradius,cx+fitradius-1,
                #                                                        cy-fitradius,cy+fitradius-1 ) ] * 0.0)
                #
                self.simstamps.append(full_data_image[ galsim.BoundsI( self.psfcenterx[i] - substamp / 2.,self.psfcenterx[i] + substamp / 2. -1,
                                                          self.psfcentery[i] - substamp / 2.,self.psfcentery[i] + substamp / 2. -1) ] * 0.0)


                self.snobjs.append(galsim.Gaussian(sigma = 1.e-8, flux = self.modelvec[i]))
                self.snoffsets.append(im.wcs.toWorld(im.trueCenter()).project(self.fiducial_coord))

                if self.fermigrid:
                    # print 'cleaning up copied files'
                    os.popen('rm '+self.imagefiles[i]).read()
                    os.popen('rm '+self.psffiles[i]).read()
                    os.popen('rm '+self.weightfiles[i]).read()

                #self.snras.append(0.)
                #self.sndecs.append(0.)

        self.model_pixel_scale_galsim = self.model_pixel_scale * galsim.arcsec
        #self.model_wcs = galsim.PixelScale(self.model_pixel_scale_galsim/galsim.arcsec)
        self.model_wcs = im.wcs
        #print galsim.GSParams().__dict__
        #raw_input()
        self.big_fft_params = galsim.GSParams(maximum_fft_size=2024000,folding_threshold=1.e-1,maxk_threshold=1.e-1)
        self.psfparams = galsim.GSParams(maximum_fft_size=20500)

        self.kicked_snraoff = copy(self.snraoff)
        self.kicked_sndecoff = copy(self.sndecoff)

        #self.xshift = copy(self.snras)*0.
        #self.yshift = copy(self.snras)*0.

        self.z_scores_say_keep_going = True

        self.sims = np.zeros([Nimage,substamp,substamp])

        #tempgalmodel = copy(self.galaxy_model)*0.

        print 'hhhhhh'
        
        self.skyerr = np.zeros([Nimage,substamp,substamp]) 
        self.mask = np.zeros([substamp,substamp]) 
        self.skyerr = self.skyerr + 99999999.
        for i in np.arange(Nimage):
            for x in np.arange(substamp):
                for y in np.arange(substamp):
                    if np.sqrt((substamp/2. - x)**2 + (substamp/2. - y)**2) < skyerr_radius:
                        self.skyerr[i,int(x),int(y)] = skyerr[i]
                        #tempgalmodel[int(x),int(y)] = copy(self.galaxy_model[int(x),int(y)])
                        self.mask[int(x),int(y)] = 1.



        #self.galaxy_model = copy(tempgalmodel)

        self.platescale = platescale
        self.model_substamp = self.galaxy_model.shape[0]
        #self.galstd = self.pixelate(self.galstd,model_pixel_scale/platescale)
        self.galdeltas = np.sqrt(self.galaxy_model)/2.#copy(self.galstd)
        self.run_d_mc()

    #@profile
    def run_d_mc( self ):
        print 'running'
        self.lastchisq = 9999999999.9
        self.chisq = []
        self.chisq.append(self.lastchisq)
        self.galhistory = []
        self.modelvechistory = []
        self.xhistory = []
        self.yhistory = []
        self.snrahistory = []
        self.sndechistory = []
        self.accepted_history = 0
        self.accepted_int = 0
        self.t1 = time.time()
        self.counter = 0
        #plt.imshow(self.data)
        #plt.show()
        #self.t2 = time.time()
        #tracker = SummaryTracker()

        while self.z_scores_say_keep_going:
            #self.t2 = time.time()
            self.counter += 1
            #print self.counter
            self.accepted_int += 1
            self.mcmc_func()
            #print self.counter
            #if (self.counter % 10) == 0:#every 100 iterations
            #    print np.array(self.csv) / len(self.mask[self.mask>0.].ravel())
            #    collected = gc.collect()
            #    #print "Garbage collector: collected %d objects." % (collected)
            #    #tracker.print_diff()


            # if self.counter == 20000:
            #     mn, st, num = dt.iterstat(np.array(self.csv)[np.array(self.csv) > 0.] / len(self.mask[self.mask>0.].ravel()),
            #                               startMedian=True, sigmaclip=3, iter=3)
            #
            #     print np.array(self.csv)/ len(self.mask[self.mask>0.].ravel()) - mn
            #     print st
            #     #raw_input()
            #     self.flags[np.array(self.csv)/ len(self.mask[self.mask>0.].ravel())-mn > 3*st + 5 ] = 1
            #     self.modelvec[np.array(self.csv)/ len(self.mask[self.mask>0.].ravel())-mn > 3*st + 5]=0.
            #     self.modelstd[np.array(self.csv)/ len(self.mask[self.mask>0.].ravel())-mn > 3*st + 5]=0.
            #
            #     self.modelvec[np.array(self.kicked_modelvec) < -50000] = 0.
            #     self.modelstd[np.array(self.kicked_modelvec) < -50000] = 0.
            #     self.flags[np.array(self.kicked_modelvec) < -50000] = 1.

            #Check Geweke Convergence Diagnostic every 5000 iterations
            if (self.counter % self.gewekenum) == self.gewekenum-1: 
                self.check_geweke()
                self.last_geweke = self.counter

            print 'psf position', self.kicked_snraoff, self.kicked_sndecoff,round(self.thischisq/len(self.mask[self.mask>0.].ravel())/len(self.flags[self.flags==0]),3)
            if (self.counter % 100) ==0:
                self.t2 = time.time()
                print 'Total Time: ' + str( self.t2 - self.t1 )
                print 'Num Iterations: ' + str( self.counter )
                print 'Accepted Percentage: ' + str( self.accepted_history )
                print 'Seconds per iteration: '+str(float(( self.t2 - self.t1 )/self.counter))

                print 'Chi Square: '+str(round(self.thischisq/len(self.mask[self.mask>0.].ravel())/len(self.flags[self.flags==0]),3))
                #print 'Chisqvec',self.chisqvec
                if (self.counter % 5000) == 0:
                    self.plotchains()
                    self.plotstamps()
                #self.savechains()

                # np.savez('test.npz',galmodel=self.modelim,
                #          wcs=self.model_wcs,
                #          snoffset=self.snoffsets[0],
                #          psf=self.psfs[0],
                #          gssimstamp=self.simstamps[0],
                #          sky=self.sky[0])
                #print 'saved test.npz'

                #raw_input()
            if self.counter > self.maxiter:
                self.z_scores_say_keep_going = False#GETOUT
                self.didtimeout = True
            #plt.imshow(self.data[20,self.substamp/2.-14.:self.substamp/2.+14.,self.substamp/2.-14.:self.substamp/2.+14.])
            #plt.show()

        self.summarize_run()
        self.model_params()

        self.t2 = time.time()

    def summarize_run( self ):
        self.t2 = time.time()
        print 'Total Time: ' + str( self.t2 - self.t1 )
        print 'Num Iterations: ' + str( self.counter )
        print 'Accepted Percentage: ' + str( self.accepted_history )
        print 'Seconds per iteration: '+str(float(( self.t2 - self.t1 )/self.counter))
        #np.savez(self.results_npz, pixel_history = self.pixel_history
        #                        , simulated_stamps = self.simulated_images
        #                        , data_stamps = self.real_data_stamps_trimmed
        #                        , sn_flux_history  = self.sn_flux_history
        #                        )
        self.plotstamps()

    #@profile
    def mcmc_func( self ):

        #t1 = time.time()
        #print 'adjusting'
        self.adjust_model()
        #t2 = time.time()

        #if self.shiftpsf:
        #    self.float_sn_pos()

        # Contains the convolution
        #print 'interpolating'

        ###new_gal_model = galsim.InterpolatedImage(self.modelim + self.kicked_galaxy_model,k_interpolant='linear')
        ###gs_model = galsim.Image(ncol=self.modelim.array.shape[1], nrow=self.modelim.array.shape[0], wcs=self.model_wcs)
        ###new_gal_model.drawImage(image=gs_model, method='no_pixel')

        self.gs_model_interp = galsim.InterpolatedImage(self.modelim + self.kicked_galaxy_model, x_interpolant='lanczos3',k_interpolant='linear',
                                                   calculate_stepk=False, calculate_maxk=False, gsparams=self.psfparams)

        #self.mapkernel()
        #self.kernel()
        #self.mapkernel()
        #print 'simming'
        # q = Queue()
        #
        # jobs = []
        # for i in range(len(self.sky)):
        #     if self.flags[i] == 0:
        #         p = multiprocessing.Process(target=self.mapkernel, args=(q, i,self.flags[i],self.fitflags[i],
        #                                                                  self.kicked_modelvec[i], self.snoffsets[i],
        #                                                                  self.psfs[i], self.simstamps[i], self.sky[i],))
        #         jobs.append(p)
        #         p.start()
        #         sim, ind = q.get()
        #         self.sims[ind,:,:] = sim
        #
        # for j in jobs:
        #     j.join()
        #     #print '%s.exitcode = %s' % (j.name, j.exitcode)


        # task_queue = Queue()
        # for i in range(len(self.sky)):
        #     if self.flags[i] == 0:
        #         task_queue.put(((i,self.flags[i],self.fitflags[i],
        #                       self.kicked_modelvec[i], self.snoffsets[i],
        #                       self.psfs[i], self.simstamps[i], self.sky[i],),))


        # for i in range(len(self.sky)):
        #     if self.flags[i] == 0:
        #         p = Process(target=f, args=(q,))
        #         p.start()
        #         q.get()  # prints "[42, None, 'hello']"
        #         p.join()
        #
        # done_queue = Queue()
        # for k in range(ncpu):
        #     Process(target=self.poolkernel, args=(task_queue, done_queue)).start()
        #
        # for i in range(len(self.sky)):
        #     if self.flags[i] == 0:
        #         psim, proc = done_queue.get()
        #         self.sims[proc,:,:] = psim
        #
        # for k in range(len(self.sky)):
        #     if self.flags[k] == 0:
        #         task_queue.put('STOP')
        #

        # pool = multiprocessing.Pool(processes=32)
        # pool.map(self.mapkernel, (range(len(self.sky)),self.flags,self.fitflags, self.kicked_modelvec, self.snoffsets,
        #          self.psfs, self.simstamps, self.sky,))

        self.sims = map(self.mapkernel, self.flags,self.fitflags, self.kicked_modelvec, self.snoffsets, self.psfs, self.simstamps, self.sky)

        #t3 = time.time()
        #print 'kernel',t3-t2

        #Calculate Chisq over all epochs
        #t4 = time.time()
        #self.thischisq = self.chisq_sim_and_real()
        #print 'chisqing'
        self.csv = map(self.mapchis, self.sims, self.data, self.flags, self.fitflags, self.skyerr, self.sky, self.gain, self.weights,self.inmask)
        #chsqs = self.csv

        self.thischisq = np.sum(self.csv)
        #t5 = time.time()

        #print self.thischisq
        #print 'chisq',t5-t4
        #raw_input()
        #t4 = time.time()

        #decide whether to accept new values
        #print self.lastchisq,self.thischisq
        #raw_input()
        #t4 = time.time()
        accept_bool = self.accept(self.lastchisq,self.thischisq)
        if self.counter == 1:
            accept_bool = False
        #print self.thischisq/len(self.mask[self.mask>0.].ravel())/len(self.flags[self.flags==0]),self.lastchisq/len(self.mask[self.mask>0.].ravel())/len(self.flags[self.flags==0]), accept_bool

        #t5 = time.time()
        #print 'accept',t5-t4
        '''
        print 'kernel'
        print t3 - t2
        print 'chisq'
        print t4 - t3
        print 'Accept bool'
        print t5 - t4
        raw_input()
        '''
        #t6 = time.time()
        if accept_bool:
            #print 'accepted'
            self.lastchisq = self.thischisq
            self.accepted_history = ( self.accepted_history * self.accepted_int + 1.0 ) / ( self.accepted_int + 1 )
            self.copy_adjusted_image_to_model()
            #self.copy_shifted_psf()
            self.update_history()
            self.chisq.append( self.thischisq )
        else:
            self.accepted_history = ( self.accepted_history * self.accepted_int ) / ( self.accepted_int + 1 )
            self.update_unaccepted_history()
            self.chisq.append(self.lastchisq)
        #t7 = time.time()
        #print 'update',t7-t6
        #print 'galval',self.galhistory[-1][10,10]
        #raw_input()
        #t6 = time.time()
        #print 'adjust model '+str(t2-t1)
        #print 'kernel '+str(t3-t2)
        #print 'chisq '+str(t4-t3)
        #print 'accept bool '+str(t5-t4)
        #print 'history update '+str(t6-t5)
        #raw_input()

    def adjust_model( self ):
        
        for i in np.arange( self.galdeltas.shape[0]  ):
            for j in np.arange( self.galdeltas.shape[1] ):
                if self.galstd[i,j] > 0.:
                    #print self.galstd.shape
                    #print self.galdeltas.shape
                    self.galdeltas[i,j] = np.random.normal( scale= self.galstd[ i, j ] )
                    '''if i == 10:
                        if j == 10:
                            print self.galstd[ i, j ]
                            print self.galdeltas[i,j]
                            raw_input()
                    '''
                #except:
                #    self.deltas[ i ] = np.random.normal( scale= self.stdev[ i ] )
                else:
                    self.galdeltas[ i, j ] = 0.0

        #print 'modelstd', self.modelstd
        #raw_input()
        for i in np.arange(len(self.modelstd)):

            if self.modelstd[i] > 0.:
                self.modeldeltas[i] = np.random.normal(scale=self.modelstd[i])
            else:
                self.modeldeltas[i] = 0.

        self.kicked_galmodel = self.galaxy_model + self.galdeltas
        self.kicked_modelvec = self.modelvec + self.modeldeltas
        
        if not self.psf_shift_std is None:
            #for i in np.arange(self.Nimage):
            if self.psf_shift_std > 0:
                self.xshift = np.random.normal(scale=self.psf_shift_std)
                self.yshift = np.random.normal(scale=self.psf_shift_std)
            else:
                self.xshift = 0.
                self.yshift = 0.
            self.kicked_snraoff = self.snraoff + self.xshift
            self.kicked_sndecoff = self.sndecoff + self.yshift
                
                #fiducial_coord = galsim.CelestialCoord(self.kicked_snras[i] * galsim.degrees, self.kicked_sndecs[i] * galsim.degrees)
                #im = self.imagestamps[i]
                #try:
                #    self.snoffsets[i] = im.wcs.toWorld(im.trueCenter()).project(fiducial_coord)
                #except:
                #    self.snoffsets[i] = np.nan

        #print self.modeldeltas
        #raw_input()
        #self.kicked_model = self.model + self.deltas

        #if not self.pixelate_model is None:
        #    self.kicked_galaxy_model = self.unpixelate(self.kicked_galmodel,self.pixelate_model,self.substamp)
        #else:
        self.kicked_galaxy_model = self.kicked_galmodel
        
        #self.kicked_galaxy_model = self.kicked_model[ 0 : self.substamp**2. ].reshape( self.substamp, self.substamp )
        return

    def float_sn_pos( self ):
        self.x_pix_offset = np.random.normal( scale= self.psf_shift_std )
        self.y_pix_offset = np.random.normal( scale= self.psf_shift_std ) 
        self.shiftPSF(x_offset=self.x_pix_offset,y_offset=self.y_pix_offset)

    def poolkernel(self,input,output):
        for (args) in iter(input.get):
            print 'argsshapeeeeeeeee',args.shape
            index, flags, fitflags, kicked_modelvec, snoffsets, psfs, simstamps, sky = args
            #= args

            sims = simstamps
            if flags == 0:
                if fitflags == 0.:
                    sn = galsim.Gaussian(sigma=1.e-8, flux=kicked_modelvec, gsparams=self.psfparams)
                    sn = sn.shift(snoffsets)  # arcsec (relative to galaxy center)
                    if not self.psf_shift_std is None:
                        sn = sn.shift(self.kicked_snraoff, self.kicked_sndecoff)
                    total_model = self.gs_model_interp + sn
                    conv = galsim.Convolve(total_model, psfs, gsparams=self.psfparams)

                    conv.drawImage(image=simstamps,
                                   method='no_pixel')  # ,offset=offset)#Draw my model to the stamp at new wcs

                    sims = simstamps.array + sky

            output.put((sims, index))

    #@profile
    def mapkernel(self, flags, fitflags, kicked_modelvec ,snoffsets, psfs, simstamps, sky ):
        #q, index,
        #self.psfparams = galsim.GSParams(maximum_fft_size=2024000,kvalue_accuracy=1.e-3,folding_threshold=1.e-1,maxk_threshold=1.e-1)

        #print sky
        sims = simstamps
        if flags == 0:
            if fitflags == 0.:
                #sn = galsim.Gaussian(sigma=1.e-8, flux=kicked_modelvec )
                #sn = sn.shift(snoffsets)  # arcsec (relative to galaxy center)
                #if not self.psf_shift_std is None:
                #    sn = sn.shift(self.kicked_snraoff, self.kicked_sndecoff)


                total_model = self.gs_model_interp


                conv1 = galsim.Convolve(total_model, psfs.shift(snoffsets), gsparams=self.psfparams)
                #conv2 = psfs.withFlux(kicked_modelvec).shift(snoffsets)

                conv1.drawImage(image=simstamps, method='no_pixel')  # ,offset=offset)#Draw my model to the stamp at new wcs
                #conv2.drawImage(image=simstamps, method='no_pixel',add_to_image=True)  # ,offset=offset)#Draw my model to the stamp at new wcs

                sims = simstamps.array + sky
        #self.sims[index,:,:] = sims

        #output.put((sims, index))

        #q.put((sims,index))
        return sims


    def kernel( self,  ):
        #t1 = time.time()


        new_gal_model = galsim.InterpolatedImage(self.modelim + self.kicked_galaxy_model,k_interpolant='linear')
        gs_model = galsim.Image(ncol=self.modelim.array.shape[1], nrow=self.modelim.array.shape[0], wcs=self.model_wcs)
        new_gal_model.drawImage(image=gs_model,method='no_pixel')

        gs_model_interp = galsim.InterpolatedImage(image=gs_model, x_interpolant='lanczos3',k_interpolant='linear',
                                                               calculate_stepk=False, calculate_maxk=False)



        #t2 = time.time()
        #print 'setting up galmodel',t2-t1
        a = []
        b = []
        #totshiftime = 0.
        #totconvtime = 0.
        #totdrawtime = 0.
        #t3 = time.time()
        for epoch in np.arange( self.Nimage ):

            #if self.fix_gal_model:
            #    star_conv = self.kicked_modelvec[ epoch ] * self.kicked_psfs[ epoch,:,:]
            #    self.sims[ epoch,:,:] =  (star_conv + self.gal_conv[epoch] + self.sky[epoch])*self.mask
            #else:
            if self.flags[epoch] == 0:
                if self.fitflags[epoch] == 0.:
                    #print 'hhhshshshsshshshshshshshshs',self.kicked_modelvec[epoch]
                    #raw_input()
                    sn = galsim.Gaussian(sigma = 1.e-8, flux = self.kicked_modelvec[epoch])
                    sn = sn.shift(self.snoffsets[epoch])  # arcsec (relative to galaxy center)
                    #t6 = time.time()
                    if not self.psf_shift_std is None:
                        sn = sn.shift(self.kicked_snraoff, self.kicked_sndecoff)
                    #t7 = time.time()
                    #totshiftime += t7-t6

                    total_model = gs_model_interp + sn

                    #print 'convolving'
                    #t4 = time.time()
                    conv = galsim.Convolve(total_model,self.psfs[epoch],gsparams=self.psfparams)

                    #print 'drawing'
                    #t5 = time.time()
                    conv.drawImage(image=self.simstamps[epoch],method='no_pixel')#,offset=offset)#Draw my model to the stamp at new wcs
                    #t6 = time.time()
                    self.sims[epoch,:,:] = self.simstamps[epoch].array + self.sky[epoch]
                    #t7 = time.time()
                    #totshiftime += t7-t6 
                    #totdrawtime += t6-t5
                    #totconvtime += t5-t4
                    #galaxy_conv = scipy.signal.fftconvolve(self.kicked_galaxy_model, self.centered_psfs[ epoch,:,:],mode='same')
                    
                    #oldgal = scipy.signal.fftconvolve(self.galaxy_model, self.centered_psfs[ epoch,:,:],mode='same')
                    #star_conv = self.kicked_modelvec[ epoch ] * self.kicked_psfs[ epoch,:,:]
                    #self.sims[ epoch,:,:] =  (star_conv + galaxy_conv + self.sky[epoch])*self.mask
        #t4 = time.time()
        #print 'convolving and drawing',t4-t3
        #print 'tot shift time',totshiftime
        #print 'tot conv time ',totconvtime
        #print 'tot draw time',totdrawtime


    def mapchis(self, sims, data, flags, fitflags, skyerr, sky, gain,weights,inmask):
        chisq = 0
        if flags == 0:
            if fitflags == 0:
                #self.readnoise = self.sky * 0. + 1
                self.gain = self.sky * 0 + 1.
                chisq += np.sum(((sims - data) ** 2 * self.mask.T * inmask.T / (
                skyerr**2 + ((sims - sky) ** 2) ** .5 / gain +
                    (self.readnoise / gain) ** 2)).ravel())
        return chisq

    def chisq_sim_and_real( self, model_errors = False ):
        chisq = np.float64(0.0)
        dms = np.float64(0.0)
        chisq_pixels = []
        self.chisqarr = np.zeros(self.sims.shape)
        self.chisqvec = np.zeros(len(self.sims))
        #print self.skyerr,self.skyerr**2
        #if self.Nimage == 1:
        #    if self.model_errors:
        #        chisq += np.sum( ( (self.sims[ 0, :,:] - self.data[ 0, :,:])**2 / (self.sims[ 0,:,:]/self.gain + (self.readnoise/self.gain)**2) ).ravel() )
        #    else:
        #        if self.useskyerr:
        #            chisq += np.sum( ( (self.sims[ 0, :,:] - self.data[ 0, :,:])**2 / self.skyerr**2).ravel() )
        #        else:
        #            chisq += np.sum( ( (self.sims[ 0, :,:] - self.data[ 0, :,:])**2 * (self.weights[ 0,:,:])**2).ravel() )
        #else:
        for epoch in np.arange( self.Nimage ):
            #print epoch
            #print chisq
            if self.flags[epoch] == 0:
                if self.fitflags[epoch] == 0:
                    self.model_errors = True
                    self.readnoise = self.sky*0. + 1
                    self.gain = self.sky*0+1.
                    if self.model_errors:
                        #chisq += np.sum( ( (self.sims[ epoch, :,:] - self.data[ epoch, : ,: ])**2 / (self.sims[ epoch,:,:]/self.gain + (self.readnoise/self.gain)**2) ).ravel() )
                        chisq += np.sum( ( (self.sims[ epoch, :,:] - self.data[ epoch, : ,: ])**2 *self.mask / (self.skyerr[epoch,:,:]**2 + ((self.sims[ epoch,:,:]-self.sky[epoch])**2)**.5/self.gain[epoch] + (self.readnoise[epoch]/self.gain[epoch])**2) ).ravel() )
                        #print tchisq
                        #if tchisq > 0:
                    else:
                        if self.useskyerr:
                            a = np.sum( ( (self.sims[ epoch, :,:] - self.data[ epoch, :,: ])**2 / self.skyerr[epoch]**2 * self.mask).ravel() )
                            chisq_pixels.extend(( (self.sims[ epoch, :,:] - self.data[ epoch, :,: ])**2 / self.skyerr[epoch]**2 * self.mask).ravel() )
                            self.chisqarr[epoch,:,:] = (self.sims[ epoch, :,:] - self.data[ epoch, :,: ])**2 / self.skyerr[epoch]**2 * self.mask
                            self.chisqvec[epoch] = a
                            #print a/len(self.sims[epoch,:,:].ravel())
                            if np.isnan(a):
                                chisq += 0.
                            else:
                                chisq += a
                            #chisq += np.float64(np.sum( ( (self.sims[ epoch, :,:] - self.data[ epoch, :,:])**2 / self.skyerr[epoch]**2).ravel() ))
                        else:
                            chisq += np.sum( ( (self.sims[ epoch, :,:] - self.data[ epoch, :,: ])**2 * (self.weights[ epoch ] ) * self.mask).ravel() )
                            #chisq_pixels.extend(( (self.sims[ epoch, :,:] - self.data[ epoch, :,: ])**2 * (self.weights[ epoch ] ) * self.mask).ravel())
                            dms +=  np.sum( self.data[ epoch,:,: ] - self.sims[ epoch, :,:])
        
        chisq_pixels = np.array(chisq_pixels)
        #print chisq_pixels.shape
        #chisq_pixels = chisq_pixels[~np.nan(chisq_pixels)]
        #chisq_pixels = chisq_pixels[~np.inf(chisq_pixels)]
        # plt.clf()
        # plt.hist(chisq_pixels,bins=np.arange(0,25000,1000))
        # chisqmeans = np.mean(chisqarr,axis=0)
        # print chisqmeans.shape
        # plt.clf()
        # plt.imshow(chisqmeans)
        # plt.colorbar()
        # plt.savefig('chisqs.png')
        # print 'saved chisq hist'
        # plt.clf()
        # plt.imshow(np.mean(self.data,axis=0))
        # plt.colorbar()
        # plt.savefig('data.png')
        # plt.clf()
        # plt.imshow(np.mean(self.sims[0],axis=0))
        # plt.colorbar()
        # plt.savefig('sim.png')
        # raw_input()
        #print 'chisq', chisq/len(self.mask[self.mask>0.].ravel())
        #print 'dms',dms
        #raw_input()
        if self.stop:
            print 'chisq map here'
            for epoch in np.arange(self.Nimage):
                save_fits_image((self.sims[ epoch, :,:] - self.data[ epoch, :,:])**2 / self.skyerr[0]**2 * self.mask,'./test/'+str(epoch)+'chisq.fits')
                save_fits_image(self.sims[ epoch, :,:]*self.mask,'./test/'+str(epoch)+'sim.fits')
                save_fits_image(self.data[ epoch, :,:],'./test/'+str(epoch)+'data.fits')
                #save_fits_image(self.psfs[epoch, :, :],'./test/'+str(epoch)+'psf.fits')
                save_fits_image(self.skyerr[epoch,:,:],'./test/'+str(epoch)+'skyerr.fits')
            save_fits_image(self.galaxy_model,'./test/galmodel.fits')
            
                #save_fits_image(self.weights[ epoch,:,:],'./weights.fits')
                #save_fits_image(self.model[:self.substamp**2].reshape(self.substamp,self.substamp),'./model.fits')
                #save_fits_image((self.sims[ epoch,:,:]/self.gain + self.readnoise/self.gain**2),'./modelerrors.fits')
            
            raw_input()
        #print 'Chisquare: '+str(chisq)
        #print 'Flux ' +str(self.model[-1])
        #raw_input()

        return chisq

    def accept( self, last_chisq, this_chisq ):
        alpha = np.exp( (last_chisq - this_chisq)/2. )
        #print ' alpha',alpha
        return_bool = False
        if alpha >= 1:
            return_bool = True
        else:
            if np.random.rand() < alpha:
                return_bool = True
        return return_bool


    def copy_adjusted_image_to_model( self ):
        self.modelvec = copy( self.kicked_modelvec )
        self.galaxy_model = copy( self.kicked_galmodel )
        if not self.psf_shift_std is None:
            self.snraoff = copy( self.kicked_snraoff )
            self.sndecoff = copy( self.kicked_sndecoff )

        return

    #def copy_shifted_psf( self ):
    #    self.psfs = copy(self.kicked_psfs)

    def update_history( self ):
        self.compressioncounter += 1
        if self.compressioncounter % self.compressionfactor == 0:
            #print 'len gal history', len(self.galhistory)
            self.galhistory.append( self.kicked_galmodel )
            self.modelvechistory.append(self.kicked_modelvec)
            
            if not self.psf_shift_std is None:
                self.snrahistory.append(self.kicked_snraoff)
                self.sndechistory.append(self.kicked_sndecoff)

            # if self.shiftpsf:
            #     self.current_x_offset += self.x_pix_offset
            #     self.current_y_offset += self.y_pix_offset
            #     self.xhistory.append(self.current_x_offset)
            #     self.yhistory.append(self.current_y_offset)
        return

    def update_unaccepted_history( self ):
        self.compressioncounter += 1
        if self.compressioncounter % self.compressionfactor == 0:
            self.galhistory.append( self.galaxy_model )
            self.modelvechistory.append( self.modelvec )

            if not self.psf_shift_std is None:
                self.snrahistory.append(self.snraoff)
                self.sndechistory.append(self.sndecoff)

            # if self.shiftpsf:
            #     self.xhistory.append(self.current_x_offset)
            #     self.yhistory.append(self.current_y_offset)
        return

    def model_params( self ):
        self.make_history()
        burn_in = int(self.modelvec_nphistory.shape[0]*.5)
        self.modelvec_params = copy(self.modelvec)
        self.galmodel_params = copy(self.galaxy_model)        
        self.modelvec_uncertainty = copy(self.modelvec)
        self.galmodel_uncertainty = copy(self.galaxy_model)
        #self.model_params = copy( self.model )
        #self.model_uncertainty = copy( self.model )
        for i in np.arange( len( self.modelvec ) ):
            self.modelvec_params[ i ] = np.median( self.modelvec_nphistory[ burn_in : , i ] )
            self.modelvec_uncertainty[ i ] = np.std( self.modelvec_nphistory[ burn_in : , i ] )            
        for i in np.arange(self.galaxy_model.shape[0]):
            for j in np.arange(self.galaxy_model.shape[1]):
                self.galmodel_params[ i, j ] = np.median( self.galmodel_nphistory[ burn_in : , i, j ] )
                self.galmodel_uncertainty[ i, j ] = np.std( self.galmodel_nphistory[ burn_in : , i, j ] )


    def autocorr( self, x ):
        result = np.correlate( x, x, mode='full' )
        return result[ result.size / 2 : ]

    def plotstamps(self):
        from matplotlib.backends.backend_pdf import PdfPages
        if self.isfermigrid:
            pdf_pages = PdfPages('galsimstamps.pdf')
        else:
            pdf_pages = PdfPages(self.lcout + '_galsimstamps.pdf')
        fig = plt.figure(figsize=(25, 10))
        for i in range(self.Nimage):
            #self.x_pix_offset = -0.49250885143
            #self.y_pix_offset = .627071191203

            self.model_params()
            #print self.x_pix_offset,self.y_pix_offset
            #raw_input('sssss')
            #map(self.mapshiftPSF, np.arange(self.Nimage))
            #self.sims = map(self.mapkernel, self.modelvec_params, self.kicked_psfs, self.centered_psfs, self.sky,
            #            self.flags, self.fitflags, self.sims, self.gal_conv)
            #self.kernel()
            wmask = copy(self.weights[i,:,:])
            wmask[wmask > 0] = 1
            v = ((self.sims[i] - self.data[i,:,:]) ** 2 *self.mask*self.inmask[i] / (self.skyerr[i]**2 + (self.sims[i] - self.sky[i]) / self.gain[i] + self.readnoise/self.gain[i])).ravel()  # hardcoded gain, hardcoded readnoise
            # v = np.real(v)
            chisq = np.sum(v[(v > 0.) & (v < 99999999.)])
            tchi = chisq/len(self.mask[self.mask>0.].ravel())
            #tchi = np.sum((self.data[i,:,:] - self.sims[i]) ** 2 / self.skyerr[i]**2 * self.mask)/len(self.mask[self.mask>0.].ravel())
            if not tchi > -1.:
                continue
            if self.flags[i] == 1:
                continue
            #fig = plt.figure(figsize=(20, 10))
            plt.clf()
            axgm = plt.subplot(151)
            axim = plt.subplot(152)
            axpsf = plt.subplot(153)
            axdiff = plt.subplot(154)
            axchi = plt.subplot(155)
            for ax, title in zip([axgm, axim, axpsf, axdiff, axchi], ['pgalmodel','image MJD '+str(round(self.mjd[i])), 'model', 'resid', 'chisq: '+str(round(tchi,2))]):
                ax.set_title(title)
            axs = axgm.imshow(self.galaxy_model * self.mask,cmap='gray',interpolation='nearest')
            cbar = fig.colorbar(axs, ax=axgm)
            axs = axim.imshow(self.data[i,:,:] * self.mask, cmap='gray', interpolation='nearest',vmin=np.min(self.sky[i]-self.sky[i]/3.),vmax=np.max(self.data[i,:,:]))
            cbar = fig.colorbar(axs, ax=axim)
            axs = axpsf.imshow(self.sims[i] * self.mask, cmap='gray', interpolation='nearest',vmin=np.min(self.sky[i]-self.sky[i]/3.),vmax=np.max(self.data[i,:,:]))
            cbar = fig.colorbar(axs, ax=axpsf)
            resid = (self.data[i,:,:] - self.sims[i])*self.mask
            md = np.median(resid[resid!=0.].ravel())
            std = np.std(resid[resid!=0.].ravel())
            axs = axdiff.imshow((self.data[i,:,:] - self.sims[i]) * self.mask, cmap='gray', interpolation='nearest',vmin=-4*std,vmax=4*std)
            cbar = fig.colorbar(axs, ax=axdiff)
            axs = axchi.imshow((self.data[i,:,:] - self.sims[i]) ** 2 / self.skyerr[i]**2 * self.mask, cmap='gray', interpolation='nearest', vmin=0, vmax=6.)
            cbar = fig.colorbar(axs, ax=axchi)
            # plt.imshow((subim-scaledpsf)/imhdr['SKYSIG'],cmap='gray',interpolation='nearest')
            # plt.colorbar()
            plt.title(title)
            pdf_pages.savefig(fig)
        pdf_pages.close()
        plt.close()
        gc.collect()
        if self.isfermigrid:
            print os.popen('ifdh rm ' + self.lcout + '_galsimstamps.pdf').read()
            print os.popen('ifdh cp galsimstamps.pdf '+self.lcout+'_galsimstamps.pdf').read()
        print 'Saved', self.lcout + '_galsimstamps.pdf'

    def get_params( self ):
        #save_fits_image(self.data[ 0, :,:],'./data.fits')
        #if self.didtimeout:

        #filelist = [ f for f in os.listdir("./out")]
        #for f in filelist:
        #    os.remove('./out/'+f)
        #print len(self.imagestamps)
        #print len(self.sims)
        #print self.Nimage
        #print 'hhhhhhh'
        #raw_input()

        datastamps = []
        simstamps = []
        galmodelstamps = []
        weightstamps = []
        psfstamps = []
        chisqstamps = []

        for i in np.arange(self.Nimage): 
            #try:
            if True:
                self.tmpwriter.savefits(self.data[i],os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_pixfluxdata.fits'))
                datastamps.append(os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_pixfluxdata.fits'))
                self.tmpwriter.savefits(self.sims[i],os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_pixfluxsim.fits'))
                simstamps.append(os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_pixfluxsim.fits'))
                self.tmpwriter.savefits(self.data[i]-self.sky[i],os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_pixfluxdataminussky.fits'))
                self.tmpwriter.savefits(self.galaxy_model,os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_pixgalmodel.fits'))
                galmodelstamps.append(os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_pixgalmodel.fits'))
                self.tmpwriter.savefits(self.weights[i,:,:],os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_pixfluxnoise.fits'))
                weightstamps.append(os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_pixfluxnoise.fits'))
                self.tmpwriter.savefits(self.data[i]-self.sims[i],os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_pixdataminussim.fits'))

                #save_fits_image(self.original_psfs[i,:,:],os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_psf.fits'))
                self.tmpwriter.savefits(self.skyerr[i,:,:],os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_pixskyerr.fits'))
                psfstamps.append('nan')
                epoch = int(i)
                chisqst = ( (self.sims[ epoch] - self.data[ epoch])**2 *self.mask / (self.skyerr[epoch]**2 + ((self.sims[ epoch]-self.sky[epoch])**2)**.5/self.gain[epoch] + (self.readnoise/self.gain[epoch])**2) )
                self.tmpwriter.savefits(chisqst,os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_chisq.fits'))
                chisqstamps.append(os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_chisq.fits'))
                # a = open(os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_pixskyval.txt'),'w')
                # a.write(str(self.sky[i]))
                # a.close()
                #print os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_pixdataminussim.fits')
            #except:
            #    pass
            #print 'didnt save image'
            ##print self.sims.shape
            #return self.model_params,self.model_uncertainty,self.nphistory, self.sims, np.asarray(self.xhistory),np.asarray(self.yhistory)
            #return np.zeros(len(self.model_params))+1e8,np.zeros(len(self.model_params))+1e9,self.nphistory
        #save_fits_image(self.data[0,:,:],'./out/MDJ'+str(self.mjd)+'data.fits')
        #return self.modelvec_params, self.modelvec_uncertainty, self.galmodel_params, self.galmodel_uncertainty, self.modelvec_nphistory, self.galmodel_nphistory, self.sims,np.asarray(self.xhistory),np.asarray(self.yhistory),self.accepted_history,self.model_substamp,self.chisq,self.ra_nphistory,self.dec_nphistory # size: self.history[num_iter,len(self.model_params)]

        stamps = [datastamps, simstamps, galmodelstamps, weightstamps, psfstamps, chisqstamps]
        chsqs = np.asarray(self.csv) / np.float64(len(self.mask[self.mask > 0.].ravel()))
        return self.modelvec_params, self.modelvec_uncertainty, self.galmodel_params, self.galmodel_uncertainty, self.modelvec_nphistory, self.galmodel_nphistory, self.sims,np.asarray(self.xhistory),np.asarray(self.yhistory),self.accepted_history, self.model_substamp, self.chisq, self.chisq, stamps, chsqs  # size: self.history[num_iter,len(self.model_params)]

    def get_params_analytical_weighted( self ):
        burn_in = int(self.nphistory.shape[0]*.5)
        model_params = copy( self.model )
        model_uncertainty = copy( self.model )
        
        for i in np.arange( len( self.model ) ):
            model_params[ i ] = np.mean( self.nphistory[ burn_in : , i ] )
            model_uncertainty[ i ] = np.std( self.nphistory[ burn_in : , i ] )

        sim = self.model_params[:self.substamp**2] + self.psfs[0,:,:].ravel()*model_params[self.substamp**2]

        sum_numer = np.sum(sim.ravel()*self.psfs[0,:,:].ravel()*self.weights[0,:,:].ravel())
        sum_denom = np.sum(self.psfs[0,:,:].ravel()*self.psfs[0,:,:].ravel()*self.weights[0,:,:].ravel())

        scale = sum_numer/sum_denom

        #compute an image of modle params and then compute sum.
        
        return scale
        
    def get_params_analytical_simple( self ):
        burn_in = int(self.nphistory.shape[0]*.5)
        model_params = copy( self.model )
        model_uncertainty = copy( self.model )
        
        for i in np.arange( len( self.model ) ):
            model_params[ i ] = np.mean( self.nphistory[ burn_in : , i ] )
            model_uncertainty[ i ] = np.std( self.nphistory[ burn_in : , i ] )

        sim = self.model_params[:self.substamp**2] + self.psfs[0,:,:].ravel()*model_params[self.substamp**2]

        sum_numer = np.sum(sim.ravel())
        sum_denom = np.sum(self.psfs[0,:,:].ravel())

        scale = sum_numer/sum_denom

        return scale

    def make_history( self ):
        num_iter = len( self.galhistory )
        self.galmodel_nphistory = np.zeros( (num_iter , self.galaxy_model.shape[0], self.galaxy_model.shape[1]))
        self.modelvec_nphistory = np.zeros( (num_iter , len(self.modelvec)))
        self.ra_nphistory = np.array(self.snrahistory)
        self.dec_nphistory = np.array(self.sndechistory)

        for i in np.arange( num_iter ):
            self.galmodel_nphistory[ i , : , : ] = self.galhistory[ i ]
            self.modelvec_nphistory[ i, : ] = self.modelvechistory[ i ]

    def savefig(self, fname):
        if self.isfermigrid:
            plt.savefig('tmp.png')
            os.popen('ifdh rm ' + fname).read()
            print os.popen('ifdh cp tmp.png '+fname).read()
            os.popen('rm tmp.png')
        else:
            plt.savefig(fname)

        print 'saved',fname

    def plotchains( self ):
        self.model_params()
        numepochs = self.modelvec_nphistory.shape[1]
        #print self.modelvec_nphistory.shape
        #raw_input()
        plt.clf()
        fig = plt.figure(1, figsize=(10,7))
        for e in np.arange(numepochs):
            plt.plot(np.arange(0,len(self.modelvec_nphistory[:,e])*self.compressionfactor,self.compressionfactor),self.modelvec_nphistory[::1,e])
            plt.xlabel('Step')
            plt.ylabel('SN Flux')
        self.savefig(str(self.lcout)+'_SNchainsg.png')
        print str(self.lcout)+'_SNchainsg.png'

        plt.clf()
        plt.close(1)

        fig = plt.figure(1,figsize=(10,7))
        #for e in np.arange(numepochs):
        plt.plot(np.arange(0,len(self.ra_nphistory)*self.compressionfactor,self.compressionfactor),self.ra_nphistory,label='RA Offset')
        plt.plot(np.arange(0,len(self.dec_nphistory)*self.compressionfactor,self.compressionfactor),self.dec_nphistory,label='DEC Offset')
        plt.legend()
        plt.xlabel('Step')
        plt.ylabel('SN Shift (arcsec)')
        self.savefig(str(self.lcout)+'_SNradechistoryg.png')
        print str(self.lcout)+'_SNradechistoryg.png'

        plt.clf()
        plt.close(1)
        '''
        fig = plt.figure(figsize=(10,7))
        for e in np.arange(numepochs):
            plt.plot(np.arange(0,len(self.dec_nphistory[:,e])*self.compressionfactor,self.compressionfactor),self.dec_nphistory[::1,e])
            plt.xlabel('Step')
            plt.ylabel('SN Shift')
        plt.savefig(str(self.lcout)+'_SNdechistory.png')
        print str(self.lcout)+'_SNdechistory.png'
        '''

    def savechains( self ):
        modelvec, modelvec_uncertainty, galmodel_params, galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, sims, xhistory,yhistory,accepted_history,model_substamp,chisqhist,self.ra_nphistory,self.dec_nphistory  = self.get_params()
        np.savez(self.chainsnpz,modelvec=modelvec, modelvec_uncertainty=modelvec_uncertainty, galmodel_params=galmodel_params, galmodel_uncertainty=galmodel_uncertainty, modelvec_nphistory=modelvec_nphistory, galmodel_nphistory=galmodel_nphistory, sims=sims,data=self.data,accepted_history=accepted_history,chisqhist=chisqhist,snrahistory=self.ra_nphistory,sndechistory=self.dec_nphistory)
    #DIAGNOSTICS
    # def check_geweke( self, zscore_mean_crit=1, zscore_std_crit=1.0 ):
    #     #print 'making history'
    #     self.make_history()
    #     #print 'geweke'
    #     zscores = self.geweke( self.nphistory[:, self.pix_stamp**2 : ] )
    #     #print 'done'
    #     #If abs(mean) of zscores is less than .5 and if stdev lt 1.0 then stop and calculate values and cov
    #     means = np.mean(zscores[1,:,:], axis=0)
    #     print means
    #     stdevs = np.std(zscores[1,:,:], axis=0)
    #     print stdevs
    #     alltrue = True
    #     for mean in means:
    #         if alltrue:
    #             if (abs(mean) > zscore_mean_crit) or (math.isnan(mean)):
    #                 alltrue = False
    #     if alltrue:
    #         for std in stdevs:
    #             if alltrue:
    #                 if (std > zscore_std_crit) or (math.isnan(std)):
    #                     alltrue = False
    #     if alltrue:
    #         self.z_scores_say_keep_going = False
    #         print 'Zscores computed and convergence criteria has been met'
    #     else:
    #         print 'Zscores computed and convergence criteria have not been met, mcmc will continue...'

    #     return

    # def geweke( self, x_in, first = .1, last = .5, intervals = 20, maxlag = 20):
    #     """Return z-scores for convergence diagnostics.
    #     Compare the mean of the first percent of series with the mean of the last percent of
    #     series. x is divided into a number of segments for which this difference is
    #     computed. If the series is converged, this score should oscillate between
    #     -1 and 1.
    #     Parameters
    #     ----------
    #     x : array-like, size x[num_params,num_iter]
    #       The trace of some stochastic parameter.
    #     first : float
    #       The fraction of series at the beginning of the trace.
    #     last : float
    #       The fraction of series at the end to be compared with the section
    #       at the beginning.
    #     intervals : int
    #       The number of segments.
    #     maxlag : int
    #       Maximum autocorrelation lag for estimation of spectral variance
    #     Returns
    #     -------

    #     """
    
    #     # Filter out invalid intervals
    #     if first + last >= 1:
    #         raise ValueError(
    #             "Invalid intervals for Geweke convergence analysis",
    #             (first, last))

    #     #if its 1d make it 2d so all code can be the same
    #     ndim = np.ndim(x_in)
    #     if ndim == 1:
    #         x = np.array(x_in.shape[0],1)
    #         x[:,0] = x_in
    #     else:
    #         x = x_in
    #     starts = np.linspace(0, int(x[:,0].shape[0]*(1.-last)), intervals).astype(int)


    #     # Initialize list of z-scores
    #     zscores = [None] * intervals 
    #     zscores = np.zeros((2,len(starts),x.shape[1]))


    #     # Loop over start indices
    #     #print len(starts)
    #     for i,s in enumerate(starts):

    #         # Size of remaining array
    #         x_trunc = x[s:,:]
    #         #print x_trunc.shape
    #         n = x_trunc.shape[0]

    #         # Calculate slices
    #         first_slice = x_trunc[ :int(first * n),:]
    #         last_slice = x_trunc[ int(last * n):,:]

    #         z = (first_slice.mean(axis=0) - last_slice.mean(axis=0))
            
    #         #spectral density
    #         z /= np.sqrt(np.fft.rfft(first_slice,axis=0)[0]/first_slice.shape[0] +
    #                  np.fft.rfft(last_slice,axis=0)[0]/last_slice.shape[0])
            
    #         #print zscores.shape
    #         #print x.shape[0]
    #         zscores[0,i,:] = np.ones(x.shape[1])*x.shape[0] - n
    #         #print z.shape
    #         zscores[1,i,:] = z
    #         #print zscores[1,:,:]

    #     #print zscores[1,:,:]
    #     #raw_input()
    #     return zscores


    def shiftPSF(self,y_offset=0.0,x_offset=0.0):

        #psf_shape = self.psfs[0,:,:].shape
        #xvals = np.arange(psf_shape[0])
        #yvals = np.arange(psf_shape[1])

        for epoch in np.arange(self.Nimage):
            #self.kicked_psfs[epoch,:,:] = self.psfs[epoch,:,:]
            #spline = scipy.interpolate.RectBivariateSpline(xvals, yvals, self.psfs[epoch,:,:])
            int_spline = np.zeros(self.psf_shape)

            #For some reason things are flipped
            x_off = y_offset
            y_off = x_offset
            #################################

            #Interpolate spline at offset
            for x,val in enumerate(self.xvals):
                #use_spline[x] = spline.ev(xvals*0 + x,yvals*0 + y)
                int_spline[x] = self.splines[epoch].ev(self.xvals*0 + x + x_off,self.yvals+y_off)

            self.kicked_psfs[epoch,:,:] = int_spline
        return

    def plot_covar( self, data ):


        # generating some uncorrelated data
        data = rand(10,100) # each row of represents a variable

        # creating correlation between the variables
        # variable 2 is correlated with all the other variables
        data[2,:] = sum(data,0)
        # variable 4 is correlated with variable 8
        data[4,:] = log(data[8,:])*0.5

        # plotting the correlation matrix
        R = corrcoef(data)
        pcolor(R)
        colorbar()
        yticks(arange(0.5,10.5),range(0,10))
        xticks(arange(0.5,10.5),range(0,10))
        show()

    def pad(self,a,new_stampsize):
        xwidth = a.shape[0]
        zeros_needed = new_stampsize - xwidth
        hz = np.floor(zeros_needed/2.)

        out = np.zeros((new_stampsize,new_stampsize))

        out[hz:-(hz+1),hz:-(hz+1)] = a[:,:]
            
        return out

    def pixelate(self,matrix,pixelation_factor):
        zmatrix = nd.interpolation.zoom(matrix, 1./float(pixelation_factor))
        return zmatrix
    
    def unpixelate(self,matrix,pixelation_factor,substamp):
        bigmat = nd.interpolation.zoom(matrix, float(pixelation_factor))
        if bigmat.shape[0] == substamp:
            outmat = bigmat
        else:
            outmat = self.pad(bigmat,substamp)
        return outmat

    '''def pixelate_galmodel(self,matrix,pixelation_factor,substamp):
        zmatrix = nd.interpolation.zoom(matrix, 1./float(pixelation_factor))
        bigmat = nd.interpolation.zoom(zmatrix, float(pixelation_factor))
        if matrix.shape[0] == substamp:
            outmat = bigmat
        else:
            outmat = self.pad(bigmat,substamp)

        return outmat
    '''
class CustomFFTConvolution(object):

    def __init__(self, A, B, threads=1):

        shape = (np.array(A.shape) + np.array(B.shape))-1
        #shape = np.array(A.shape)
        if np.iscomplexobj(A) and np.iscomplexobj(B):
            self.fft_A_obj = pyfftw.builders.fftn(
                    A, s=shape, threads=threads)
            self.fft_B_obj = pyfftw.builders.fftn(
                    B, s=shape, threads=threads)
            self.ifft_obj = pyfftw.builders.ifftn(
                    self.fft_A_obj.get_output_array(), s=shape,
                    threads=threads)

        else:
            self.fft_A_obj = pyfftw.builders.rfftn(
                    A, s=shape, threads=threads)
            self.fft_B_obj = pyfftw.builders.rfftn(
                    B, s=shape, threads=threads)
            self.ifft_obj = pyfftw.builders.irfftn(
                    self.fft_A_obj.get_output_array(), s=shape,
                    threads=threads)

    def __call__(self, A, B):

        fft_padded_A = self.fft_A_obj(A)
        fft_padded_B = self.fft_B_obj(B)

        return self.ifft_obj(fft_padded_A * fft_padded_B)


def save_fits_image(image,filename):
    hdu = pf.PrimaryHDU(image)
    if os.path.exists(filename):
        os.remove(filename)
    hdu.writeto(filename)
    return

if __name__ == "__main__":
    '''
    #TEST DATA
    # 4 by for image with 4 supernova epochs initalized to 1
    Nepochs = 4
    substamp = 5
    model = np.array([250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,0,119900,160000,200000])
    stdev = np.array([20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,0.,25.,25.,25.])
    

    data = np.zeros((4,5,5))
    a = [250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250]
    ina = np.asarray(a)
    x = ina.reshape(5,5)
    data[0,:,:] = x
    a = [250,250,250,250,250,250,250,250,250,250,250,250,3250,250,250,250,250,250,250,250,250,250,250,250,250]
    ina = np.asarray(a)
    x = ina.reshape(5,5)
    data[1,:,:] = x
    a = [250,250,250,250,250,250,250,250,250,250,250,250,4250,250,250,250,250,250,250,250,250,250,250,250,250]
    ina = np.asarray(a)
    x = ina.reshape(5,5)
    data[2,:,:] = x
    a = [250,250,250,250,250,250,250,250,250,250,250,250,5250,250,250,250,250,250,250,250,250,250,250,250,250]
    ina = np.asarray(a)
    x = ina.reshape(5,5)
    data[3,:,:] = x

    psfs = np.ones((4,5,5))/1000.

    #psf = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]).reshape(5,5)
    #for epoch in np.arange(Nepochs):
    #    psfs[epoch,:,:] = psf

    weights = 1/(np.ones((4,5,5))+10)
    #weight = 1./(np.ones(25).reshape(5,5)+4.)
    #for epoch in np.arange(Nepochs):
    #    weights[epoch,:,:] = weight


    a = metropolis_hastings( model = model
        , stdev = stdev
        , data = data
        , psfs = psfs
        , weights = weights
        , substamp = substamp
        , Nimage = Nepochs
        )

    model, uncertainty, history = a.get_params()

    print 'FINAL MODEL'
    print model
    print 'MODEL Uncertainty'
    print uncertainty
    '''


    f = np.load('/scratch2/scratchdirs/dbrout/smp_y1y2_shallow62/np_data/r/des_fake_00224387_r_mcmc_input.npz')
    t1 = time.time()
    a = metropolis_hastings(galmodel = f['galmodel']
                , modelvec = f['modelvec']
                , galstd = f['galstd']
                , modelstd = f['modelstd']
                , data = f['data']
                , psfs = f['psfs']
                , weights = f['weights']
                , substamp = f['substamp']
                , Nimage = f['Nimage']
                , maxiter = 100.
                , mask = None
                , sky= f['sky']
                , mjd= f['mjd']
                , gewekenum= f['gewekenum']
                , skyerr= f['skyerr']
                , useskyerr = True
                , flags = f['flags']
                , psf_shift_std = .0005
                , shiftpsf = False
                , fileappend = ''
                , stop = False
                , skyerr_radius = 16.
                , outpath = f['outpath']
                , compressionfactor = 1
                , fix_gal_model = False
                , pixelate_model = 1.
                )
    t2 = time.time()
    print 'seconds per iter',(t2-t1)/100.

