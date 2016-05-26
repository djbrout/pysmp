import os
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(threshold=np.nan)
import sys
import pyfits as pf

filt = 'r'
zpt = 27.5
nominal_zpt = 31.
isdeep = False
#smpzpt = 31.
dopix = False
loadstardeltas = False
doglobalsn = False
dooff = False
dostardeltas = False
dogalsim = True
usechisq = False
#outfolder = '/global/cscratch1/sd/dbrout/smp_y1y2_shallow79_galsim11'
outfolder = '/global/cscratch1/sd/dbrout/v3/smp_y1y2_shallow_v3_57testfixedgalsimfloatpos'
pixoutfolder = '/global/cscratch1/sd/dbrout/v3/smp_y1y2_shallow_v3_40globalstars'
nothingfolder = '/global/cscratch1/sd/dbrout/v3/smp_y1y2_shallow_v3_35'
nodailyofffolder = '/global/cscratch1/sd/dbrout/v3/smp_y1y2_shallow_v3_40globalstars'
dov3=False
files = os.listdir(os.path.join(outfolder,'np_data/'+filt+'/'))
#pixfiles =  os.listdir(os.path.join(pixoutfolder,'np_data/'+filt+'/'))
#nfiles = os.listdir(os.path.join(nothingfolder,'np_data/'+filt+'/'))
#dofiles = os.listdir(os.path.join(nodailyofffolder,'np_data/'+filt+'/'))
#print files
#raw_input()
	

snnames = []
for fle in files:
	if fle.split('_')[-1] == 'input.npz':
		snnames.append('_'.join(fle.split('_')[0:4]))
#print sns
#raw_input()

bigcat = pf.open('/global/homes/d/dbrout/PySMP/SNscampCatalog/DES-SN_v2.cat')[2].data
bigcatras = bigcat['x_world']
#print bigcatras
#raw_input()
bigcatdecs = bigcat['y_world']
bigcatids = bigcat['source_number'] 
finalresultsfiles = []
mcmc_inputs = []
withsns = []
withsngalsims = []
withsngalsimspix = []
nosns = []
justnames = []
withsnsRADEC = []
for sn in snnames:
	print sn.split('/')[-1]
	finalresultsfiles.append(sn+'_plotresults.npz')
	mcmc_inputs.append(sn+'_mcmc_input.npz')
	withsns.append(sn+'_withSn.npz')
	withsngalsims.append(sn+'_withSnAndGalsim.npz')
	nosns.append(sn+'_nosn.npz')
	withsngalsimspix.append(sn+'_withSnAndGalsimPix.npz')
	justnames.append(sn.split('/')[-1])
	withsnsRADEC.append(sn+'_withSngetRADEC.npz')

print 'Number of SNe ',len(snnames)
print 'press enter to continue'
raw_input()

def madstd(x):
	madstd = 1.48*np.median(abs(x-np.median(x)))
	return madstd

def bindata(x,y,bins):
	medians = np.zeros(len(bins)-1)
	mads = np.zeros(len(bins)-1)
	for i in np.arange(len(bins)-1):
		bs = bins[i]
		bf = bins[i+1]
		ww = [(x>bs)&(x<bf)]
		yhere = y[ww]
		ss = [abs(yhere) < 3*np.std(yhere)]
		try:
			medians[i] = np.median(yhere[ss])
			mads[i] = 1.48*np.median(abs(yhere[ss]-medians[i]))*1/np.sqrt(len(yhere[ss]))
		except IndexError:
			medians[i] = np.nan
			mads[i] = np.nan
	xvals = (bins[1:] + bins[:-1])/2.
	return xvals,medians,mads

def bindataweighted(x,y,ye,bins):
        medians = np.zeros(len(bins)-1)
        mads = np.zeros(len(bins)-1)
        for i in np.arange(len(bins)-1):
                bs = bins[i]
                bf = bins[i+1]
                ww = [(x>bs)&(x<bf)]
		yhere = y[ww]
		yehere = ye[ww]
		ss = [abs(yhere) < 3*np.std(yhere)]
                try:
                        medians[i] = np.average(yhere[ss],weights=yehere[ss])
                        mads[i] = 1.48*np.median(abs(yhere[ss]-medians[i]))*1/np.sqrt(len(yhere[ss]))
                except:
                        medians[i] = np.nan
                        mads[i] = np.nan
        xvals = (bins[1:] + bins[:-1])/2.
        return xvals,medians,mads


def bindata2(x,y,nbins):
	n, b = np.histogram(x, bins=nbins)
	sy, b = np.histogram(x, bins=nbins, weights=y)
	sy2, b = np.histogram(x, bins=nbins, weights=y*y)
	n = np.float64(n)
	n[n==0] = np.nan
	mean = sy / n
	std = np.sqrt(sy2/n - mean*mean)
	xvals = (b[1:] + b[:-1])/2
	return xvals,mean,std

def numfivesigma(x,y,nbins):
	n, b = np.histogram(x, bins=nbins)
	sy, b = np.histogram(x, bins=nbins, weights=y)
	sy2, b = np.histogram(x, bins=nbins, weights=y*y)
	n = np.float64(n)
	n[n==0] = np.nan
	mean = sy / n
	std = np.sqrt(sy2/n - mean*mean)
	xvals = (b[1:] + b[:-1])/2

	num = 0
	tot = 0
	for i in np.arange(len(b)-1):
		xstart = b[i]
		xfin = b[i+1]
		ww = [(x>xstart) & (x<xfin)]
		thisbit = y[ww]
		outliers = thisbit[(abs(thisbit - mean[i]) > 5*std[i])]
		num += len(outliers)
		tot += len(thisbit)
	try:
		p =  float(num)/float(tot)
	except:
		p = 0
	return num, p
def closest_node(ra,dec,catalogras,catalogdecs):
        tra = catalogras*0. + ra
        tdec = catalogdecs*0. + dec
	
        dist_2 = ((catalogras-tra)**2+(catalogdecs-tdec)**2)**.5
	return np.argmin(dist_2)
def getProperCatRaDec(ra,dec,catalogras,catalogdecs):
        properra = np.zeros(len(ra))
        properdec = np.zeros(len(dec))
        for i in np.arange(0,len(ra)):
            j = closest_node(ra[i],dec[i],catalogras,catalogdecs)
            properra[i] = catalogras[j]
            properdec[i] = catalogdecs[j]
        return properra,properdec


from scipy.stats import gaussian_kde
def kde_scipy(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scipy"""
    # Note that scipy weights its bandwidth by the covariance of the
    # input data.  To make the results comparable to the other methods,
    # we divide the bandwidth by the sample standard deviation here.
    kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
    return kde.evaluate(x_grid)

#os.system('rm '+os.path.join(outfolder,'lightcurves/r/*.pdf'))
#npzs = []
#inputs = []
#for f in files:
#	if '.npz' in f:#
#		npzs.append(f)


#print npzs
#plt.hold(False)
mjds = []
chis = []
chisqs = []
fakechis = []
fakechisqs = []
mcmc_fixmags = []
mcmc_fixmag_errs = []
mcmc_fixfluxes = []
mcmc_fixfluxerrs = []
mcmc_floatmags = []
mcmc_floatmag_errs = []
dofluxes = []
dofluxerrs = []
domags = []
domagerrs = []
nfluxes = []
nfluxerrs = []
nmags = []
nmagerrs = []
mcmc_gfloatmags = []
mcmc_gfloatmag_errs = []
mcmc_gpixfloatmags = []
mcmc_gpixfloatmag_errs = []

snraoff = []
sndecoff = []

mcmc_floatfluxes = []
mcmc_floatfluxerrs = []
mcmc_gfloatfluxes = []
mcmc_gfloatfluxerrs = []
mcmc_gpixfloatfluxes = []
mcmc_gpixfloatfluxerrs = []
diffim_mags = []
diffim_mag_errs = []
diffim_fluxes = []
diffim_fluxerrs = []

input_diffim_fluxerrs = []
input_fakeflux = []
input_fakemag = []
input_flags = []
fakemags = []
peakmjds = []
fake_fluxes = []
fake_zpts = []
zpt_errs = []
psfs = []
flags = []
skyerrs = []
tses = []
skysigs= []
hostgal_bandmags = []
fakemags_plus = []
fakemags_minus = []
fake_fluxes_plus = []
fake_fluxes_minus = []
fake_flux_errs = []
zpt_offsets = []
acceptedrate = []
mjdflags = []
sumchisqs = []
snnamess = []
snraoffs = []
sndecoffs = []
globalraoffsets = []
globaldecoffsets = []
globalsnra = []
globalsndec = []
toomuch = -1
i = -1
totscounter = 0
#print snnames
for sn in np.arange(len(snnames)):
	toomuch += 1
	if totscounter > 100000:
		continue
	i+=1
	#print i
	# print os.path.join(outfolder,'np_data/'+filt+'/'+finalresultsfiles[sn])
	# print os.path.join(outfolder,'np_data/'+filt+'/'+finalresultsfiles[sn])
	# finalresults = np.load(os.path.join(outfolder,'np_data/'+filt+'/'+finalresultsfiles[sn]))
	# mcmc_input = np.load(os.path.join(outfolder,'np_data/'+filt+'/'+finalresultsfiles[sn]))
	
	#try:
	#	print os.path.join(outfolder,'np_data/'+filt+'/'+withsns[sn])
	#finalresults = np.load(os.path.join(outfolder,'np_data/'+filt+'/'+finalresultsfiles[sn]))
	skip = False
	try:
		mcmc_input = np.load(os.path.join(outfolder,'np_data/'+filt+'/'+mcmc_inputs[sn]))
	except:
		continue
	print 'mjdoff',np.max(abs(mcmc_input['mjdoff'])),mcmc_inputs[sn]
	#raw_input()
	try:	
		if dogalsim:
			print 'galsimmm'
			withsn = np.load(os.path.join(outfolder,'np_data/'+filt+'/'+withsngalsimspix[sn]))
			print withsn.keys()
			#print withsn['modelvec_nphistory'].shape
			#raw_input()
		else:
			withsn = np.load(os.path.join(outfolder,'np_data/'+filt+'/'+withsns[sn]))
	except:
		print 'Could not find ',os.path.join(outfolder,'np_data/'+filt+'/'+withsns[sn])
		skip = True
		#raw_input()
                #withsngalsim = np.load(os.path.join(outfolder,'np_data/'+filt+'/'+withsngalsims[sn]))
	if dopix:
		try:
			withsngalsimpix = np.load(os.path.join(pixoutfolder,'np_data/'+filt+'/'+withsngalsimspix[sn]))
		except:
		        #print 'Skipped '+os.path.join(outfolder,'np_data/'+filt+'/'+finalresultsfiles[sn])+' Could not find npz file'
			print 'Could not find',os.path.join(pixoutfolder,'np_data/'+filt+'/'+withsngalsimspix[sn])
			skip = True
	if dov3:
		try:
			#nf = np.load(os.path.join(nothingfolder,'np_data/'+filt+'/'+withsns[sn]))
			dof = np.load(os.path.join(nodailyofffolder,'np_data/'+filt+'/'+withsns[sn]))
		except:
			print 'Could not find '+os.path.join(nodailyofffolder,'np_data/'+filt+'/'+withsns[sn])
                        skip = True

	if skip:
		continue
	#print finalresults.keys()
	#print withsngalsimpix.keys()
	if doglobalsn:
		try:
			snradec = np.load(os.path.join(outfolder,'np_data/'+filt+'/'+withsnsRADEC[sn]))
			raoff = snradec['raoff'][0]
			decoff = snradec['decoff'][0]
		except:
			print 'Could not find '+os.path.join(outfolder,'np_data/'+filt+'/'+withsnsRADEC[sn])
			#raw_input()
			skip = True
	if skip:
		continue
	if dopix:
		if withsngalsimpix['accepted_history'] > .7:
			print snnames[sn]+' skipped too many steps accepted:', withsngalsimpix['accepted_history']
			raw_input()
			continue
		if withsngalsimpix['accepted_history'] < .15:
			raw_input()
			print snnames[sn]+' skipped too many steps skipped:',withsngalsimpix['accepted_history']
			continue
	#print withsn.keys()
	#print withsn['chisqvec'].shape
	#print withsn['chisqhist'].shape
	#print withsn['redchisqhist'].shape
	#print withsn['modelvec_nphistory']
	mm = np.argmax(withsn['modelvec'])
	wsn_val = []
	wsnh = withsn['modelvec_nphistory'][:,mm]
	if len(wsnh) < 300:
		print snnames[sn]+' Not enough steps for convergence '+str(len(wsnh))+'... Skipping...'
		raw_input()
		continue
	
	#print 'heyyyyyy'
	#raw_input()
	nval = []
	doval = []
	'''
	if dostardeltas:
		print outfolder+'/stardata/'+filt+'/'+justnames[sn]+'band_starGlobalOffsets.npz'
		staroffsets = np.load(outfolder+'/stardata/'+filt+'/'+justnames[sn]+'band_starGlobalOffsets.npz')
		starras = staroffsets['starras']
		stardecs = staroffsets['stardecs']
		starids = staroffsets['starids']
		starglobalids = []
		starglobalras = []
		starglobaldecs = []

		starcatras = []
		starcatdecs= []
		
		for ide in np.unique(np.array(starids)):
			ww = (starids == ide)
			starglobalids.append(ide)
			starglobalras.append(np.median(starras[ww]))
			starglobaldecs.append(np.median(stardecs[ww]))
			nra,ndec = getProperCatRaDec([starglobalras[-1]],[starglobaldecs[-1]],bigcatras,bigcatdecs)
		#print so.keys()
		scampra,scampdec = getProperCatRaDec(starglobalras,starglobaldecs,bigcatras,bigcatdecs)
	offsetra = np.array(starglobalras) - np.array(scampra)
	offsetdec = np.array(starglobaldecs) - np.array(scampdec)
	globalraoffsets.append(offsetra)
	globaldecoffsets.append(offsetdec)
	'''
	#raw_input()
	'''dh = dof['modelvec_nphistory'][:,mm]
	nh = nf['modelvec_nphistory'][:,mm]
	for i in range(0,len(wsnh),10):
		wsn_val.append(np.mean(wsnh[i:]))
		nval.append(np.mean(dh[i:]))
		doval.append(np.mean(nh[i:]))
	plt.clf()
	plt.plot(wsn_val,color='blue')
        #plt.plot(nval)
        #plt.plot(doval)
	plt.axhline(withsn['modelvec'][mm],color='red')
	plt.axhline(np.median(withsn['modelvec_nphistory'][:,mm]),color='blue')
	plt.savefig('fluxhistory.png')
	print 'fluxhistory.png'
	print withsn['modelvec'].shape
	print withsn['modelvec_nphistory'].shape
	if 31.-2.5*np.log10(withsn['modelvec'][mm]) < 21.:
		raw_input()
	'''
	#print 'all',np.mean(withsn['redchisqhist'][1000:])
	#print 'doff',np.mean(dof['redchisqhist'][1000:])
	#print 'soff',np.mean(nf['redchisqhist'][1000:])

	#raw_input()

	#print snnames[sn]
	#print withsn['modelvec']
	#raw_input()

	#ipt = np.load(os.path.join(outfolder,'np_data/'+filt+'/'+mcmc_inputs[sn]))

	mcmc_fixflux = []
	mcmc_fixfluxerr = []
	mcmc_floatflux = []
	mcmc_floatfluxerr = []
	mcmc_floatchi = []
	mcmc_gfloatflux = []
	mcmc_gfloatfluxerr = []
	mcmc_gpixfloatflux = []
	mcmc_gpixfloatfluxerr = []
	doflux = []
	dofluxerr =[]
	nflux = []
	nfluxerr =[]

	snraoff = []
	sndecoff = []

	zpt = []
	fakezpt = []
	chisq = []
	chi = []
	diffim_flux = []
	diffim_fluxerr = []
	mjd = []
	fakemag = []
	hostmags = []
	flag = []
	skyerr = []
	skysig = []
	tse = []
	sns = []
	psf = []
	try:
		a = mcmc_input['skysig']
	except:
		print 'could not find skysig'
		continue
	if usechisq:
		if withsn['redchisqhist'][-1] > 1.05:
			print 'REDUCED CHI SQARED IS POOR'
			continue
	#print withsngalsim.keys()
	#raw_input()
	#try:
	for i in np.arange(len(mcmc_input['mjd'])):
		mp = withsn['modelvec_nphistory'][:,i]
		#print (np.mean(mp[-1*int(len(mp)/2):])-np.mean(mp[-1*int(len(mp)/4):]))/np.mean(mp[-1*int(len(mp)/2):]),np.mean(mp),len(mp)
		#raw_input()
		if not float(madstd(withsn['modelvec_nphistory'][:,i])) != 0.:
			mcmc_floatflux.append(np.nan)
			mcmc_floatfluxerr.append(np.nan)
		        #elif (np.mean(mp[-1*int(len(mp)/2):])-np.mean(mp[-1*int(len(mp)/4):]))/np.mean(mp[-1*int(len(mp)/2):]) :
		elif abs(np.mean(mp[-1*int(len(mp)/2):])-np.mean(mp[-1*int(len(mp)/4):]))/np.mean(mp[-1*int(len(mp)/2):]) < 500.:
			mcmc_floatflux.append(float(np.median(withsn['modelvec_nphistory'][-1*int(len(mp)/4):,i])))
			mcmc_floatfluxerr.append(float(madstd(withsn['modelvec_nphistory'][-1*int(len(mp)/4):,i])))
		else:
			mcmc_floatflux.append(np.nan)
                        mcmc_floatfluxerr.append(np.nan)


		#mcmc_gfloatflux.append(float(withsngalsim['galsim_scales'][i]))
		#mcmc_gfloatfluxerr.append(float(withsngalsim['galsim_stds'][i]))
		#mcmc_gfloatflux.append(float(withsngalsim['modelvec'][i]))
		#mcmc_gfloatfluxerr.append(float(withsngalsim['modelvec_uncertainty'][i]))
		if dopix:
			mcmc_gpixfloatflux.append(float(withsngalsimpix['modelvec'][i]))
			mcmc_gpixfloatfluxerr.append(float(withsngalsimpix['modelvec_uncertainty'][i]))
			
		if dov3:
			doflux.append(float(np.median(dof['modelvec_nphistory'][:,i])))
                        dofluxerr.append(float(madstd(dof['modelvec_nphistory'][:,i])))

			#doflux.append(dof['modelvec'][i])
			#dofluxerr.append(dof['modelvec_uncertainty'][i])

			#nflux.append(nf['modelvec'][i])
			#nfluxerr.append(nf['modelvec_uncertainty'][i])

			#nflux.append(float(np.median(nf['modelvec_nphistory'][:,i])))
	                #nfluxerr.append(float(madstd(nf['modelvec_nphistory'][:,i])))
			nflux.append(0)
			nfluxerr.append(0)
		zpt.append(float(mcmc_input['zpt'][i]))
		if doglobalsn:
			snraoff.append(float(withsn['raoff']))
			sndecoff.append(float(withsn['decoff']))
		fakezpt.append(float(mcmc_input['fakezpt'][i]))
		# chisq.append(float(finalresults['chisq'][i]))
		# chi.append(float(finalresults['dms'][i]))
		mjd.append(float(mcmc_input['mjd'][i]))
		psf.append(float(mcmc_input['fwhm'][i]))
		flag.append(float(mcmc_input['flags'][i]))
		# diffim_flux.append(float(finalresults['diffim_flux'][i]))
		# diffim_fluxerr.append(float(finalresults['diffim_flux'][i]))
		fakemag.append(float(mcmc_input['fakemag'][i]))
		skyerr.append(float(mcmc_input['skyerr'][i]))
		skysig.append(float(mcmc_input['skysig'][i]))
		tse.append(float(mcmc_input['total_skyerr'][i]))
		diffim_flux.append(float(mcmc_input['diffim_flux'][i]))
		diffim_fluxerr.append(float(mcmc_input['diffim_fluxerr'][i]))
		hostmags.append(float(mcmc_input['hostgal_sbmag'][i]))
		#print mcmc_input['hostgal_sbmag']
		#raw_input()
		sns.append(snnames[sn])
	#except:
	#	print 'gogogogogogogog'
	#	continue

	#print mcmc_gpixfloatflux
	#print mcmc_floatflux
	#raw_input()	

	#thissn = []
	#for jj in mcmc_flux:
	#	thissn.append(sn)
	if dopix:
		acceptedrate.append(withsngalsimpix['accepted_history'])
	else:
		acceptedrate.append(1)

	mjd = np.asarray(mjd)
	mcmc_fixflux = np.asarray(mcmc_fixflux)
	mcmc_fixfluxerr = np.asarray(mcmc_fixfluxerr)
	mcmc_floatflux = np.asarray(mcmc_floatflux)
	mcmc_floatfluxerr = np.asarray(mcmc_floatfluxerr)
	mcmc_gfloatflux = np.asarray(mcmc_gfloatflux)
	mcmc_gfloatfluxerr = np.asarray(mcmc_gfloatfluxerr)
	mcmc_gpixfloatflux = np.asarray(mcmc_gpixfloatflux)
	mcmc_gpixfloatfluxerr = np.asarray(mcmc_gpixfloatfluxerr)
	if dov3:
		nflux = np.asarray(nflux)
		nfluxerr = np.asarray(nfluxerr)
		doflux = np.asarray(doflux)
		dofluxerr = np.asarray(dofluxerr)
	diffim_flux = np.asarray(diffim_flux)
	diffim_fluxerr = np.asarray(diffim_fluxerr)
	fakemag = np.asarray(fakemag)
	hostmags = np.asarray(hostmags)
	#thissn = np.asarray(thissn,dtype='string')

	mjds.extend(mjd)
	psfs.extend(psf)
	flags.extend(flag)
	#chis.extend(chi)
	#chisqs.extend(chisq)
	#sns.extend(thissn)
	#mjdflags.extend(mjdflag)
	fakemags.extend(fakemag)
	skyerrs.extend(skyerr)
	skysigs.extend(skysig)
	tses.extend(tse)
	snnamess.extend(sns)

	mcmc_fixfluxes.extend(mcmc_fixflux)
	mcmc_fixfluxerrs.extend(mcmc_fixfluxerr)
	mcmc_floatfluxes.extend(mcmc_floatflux)
	mcmc_floatfluxerrs.extend(mcmc_floatfluxerr)
	mcmc_gfloatfluxes.extend(mcmc_gfloatflux)
	mcmc_gfloatfluxerrs.extend(mcmc_gfloatfluxerr)
	mcmc_gpixfloatfluxes.extend(mcmc_gpixfloatflux)
	mcmc_gpixfloatfluxerrs.extend(mcmc_gpixfloatfluxerr)
	hostgal_bandmags.extend(hostmags)
			       
	snraoffs.extend(snraoff)
	sndecoffs.extend(sndecoff)
	
	if dov3:
		nfluxes.extend(nflux)
		nfluxerrs.extend(nfluxerr)
		dofluxes.extend(doflux)
		dofluxerrs.extend(dofluxerr)

	# if isdeep:
	# 	prevday = 0.
	# 	currentvec = []
	# 	indexes = []
	# 	sumt = 0
	# 	sume = 0
	# 	smpf = []
	# 	smpfe = []
	# 	index = -1
	# 	started = False
	# 	#print mjd
	# 	for day in mjd:
	# 		#print day
	# 		index += 1
	# 		if day - prevday < 1.:
	# 			sumt += diffim_flux[index]
	# 			sume += diffim_fluxerr[index]**2
	# 			indexes.append(index)
	# 			if not started:
	# 				sumt += diffim_flux[index-1]
	# 				sume += diffim_fluxerr[index-1]**2
	# 				indexes.append(index - 1)
	# 				started = True
	# 		else:
	# 			'''try:
	# 				print smpf
	# 				print np.average(np.asarray(smpf),weights=1./np.asarray(smpfe))
	# 				raw_input()
	# 			except:
	# 				continue'''
	# 			if len(indexes) > 0.:
	# 				flux = np.average(np.asarray(mcmc_floatflux)[indexes],weights=1./np.asarray(mcmc_floatfluxerr)[indexes])
	# 				fluxerr = np.sqrt(np.sum(np.asarray(mcmc_floatfluxerr)[indexes]**2))
	# 				#print flux
	# 				#print fluxerr
	# 			for i in indexes:

	# 				#print 'before',mjd[i],diffim_flux[i]
	# 				diffim_flux[i] = sumt
	# 				diffim_fluxerr[i] = np.sqrt(sume)

	# 				mcmc_floatflux[i] = flux
	# 				mcmc_floatfluxerr[i] = fluxerr

	# 				#diffim_flux[i+1] = sumt
	# 				#diffim_fluxerr[i+1] = np.sqrt(sume)
	# 				#print 'after',mjd[i],diffim_flux[i]
	# 			#print mjd[indexes]
	# 			indexes = []
	# 			sume = 0
	# 			sumt = 0
	# 			smpf = []
	# 			smpfe = []
	# 			started = False
	# 			#raw_input()
	# 		#print indexes
	# 		#raw_input()
	# 		prevday = day




	# diffim_mag = 31.-2.5*np.log10(diffim_flux)
	# diffim_mag_err = -2.5*np.log10(diffim_flux)+2.5*np.log10(diffim_flux+diffim_fluxerr)

	if dov3:
		nmag = 31.-2.5*np.log10(nflux)
		nmagerr = -2.5*np.log10(nflux)+2.5*np.log10(nflux+nfluxerr)
		domag = 31.-2.5*np.log10(doflux)
		domagerr= -2.5*np.log10(doflux)+2.5*np.log10(doflux+dofluxerr)
		nmags.extend(nmag)
		nmagerrs.extend(nmagerr)
		domags.extend(domag)
		domagerrs.extend(domagerr)
	mcmc_fixmag = 31.-2.5*np.log10(mcmc_fixflux)
	mcmc_fixmag_err = -2.5*np.log10(mcmc_fixflux)+2.5*np.log10(mcmc_fixflux+mcmc_fixfluxerr)
	mcmc_fixmags.extend(mcmc_fixmag)
	mcmc_fixmag_errs.extend(mcmc_fixmag_err)	
	mcmc_floatmag = 31.-2.5*np.log10(mcmc_floatflux)
	mcmc_floatmag_err = -2.5*np.log10(mcmc_floatflux)+2.5*np.log10(mcmc_floatflux+mcmc_floatfluxerr)
	mcmc_floatmags.extend(mcmc_floatmag)
	mcmc_floatmag_errs.extend(mcmc_floatmag_err)

	mcmc_gfloatmag = 31.-2.5*np.log10(mcmc_gfloatflux)
	mcmc_gfloatmag_err = -2.5*np.log10(mcmc_gfloatflux)+2.5*np.log10(mcmc_gfloatflux+mcmc_gfloatfluxerr)
	mcmc_gfloatmags.extend(mcmc_gfloatmag)
	mcmc_gfloatmag_errs.extend(mcmc_gfloatmag_err)


	mcmc_gpixfloatmag = 31.-2.5*np.log10(mcmc_gpixfloatflux)
	mcmc_gpixfloatmag_err = -2.5*np.log10(mcmc_gpixfloatflux)+2.5*np.log10(mcmc_gpixfloatflux+mcmc_gpixfloatfluxerr)
	mcmc_gpixfloatmags.extend(mcmc_gpixfloatmag)
	mcmc_gpixfloatmag_errs.extend(mcmc_gpixfloatmag_err)

	#print np.mean(mcmc_floatmag[np.isfinite(mcmc_floatmag) & np.isfinite(domag)]-domag[np.isfinite(mcmc_floatmag) & np.isfinite(domag)])
	#raw_input('jjj')

	fake_flux = 10.**(.4*(31.-fakemag))
	
	#SCALE FAKE FLUX BY DIFFERENCE BETWEEN RICK (PLOPPED) AND SMP ZEROPOINTS

	fakezpt = np.asarray(fakezpt)
	zpt = np.asarray(zpt)
	scalefactor = 10.**(0.4*(fakezpt-zpt))
	fake_flux_plus = fake_flux*10.**(0.4*(fakezpt-zpt))
	fake_flux_minus = fake_flux*10.**(-0.4*(fakezpt-zpt))

	fakemag_minus = fakemag - (fakezpt-zpt)
	fakemag_plus = fakemag + (fakezpt-zpt)
	#mcmc_fixmag_plus = mcmc_fixmag + (fakezpt - zpt)
	#mcmc_fixmag_minus = mcmc_fixmag - (fakezpt - zpt)
	#mcmc_floatmag_plus = mcmc_floatmag + (fakezpt - zpt)
	#mcmc_floatmag_minus = mcmc_floatmag - (fakezpt - zpt)

	fakemags_plus.extend(fakemag_plus)
	fakemags_minus.extend(fakemag_minus)
	fake_fluxes_plus.extend(fake_flux_plus)
	fake_fluxes_minus.extend(fake_flux_minus)
	zpt_offsets.extend(fakezpt-zpt)


	diffim_fluxes.extend(diffim_flux)
	diffim_fluxerrs.extend(diffim_fluxerr)
	# diffim_mags.extend(diffim_mag)
	# diffim_mag_errs.extend(diffim_mag_err)
	fake_fluxes.extend(fake_flux)
	if doglobalsn:
		#print np.median(np.array(snradec['xhistory']))
		globalsnra.extend(fakemag*0. + raoff)
		globalsndec.extend(fakemag*0 + decoff)
	if dostardeltas:
	        print outfolder+'/stardata/'+filt+'/'+justnames[sn]+'band_starGlobalOffsets.npz'
                staroffsets = np.load(outfolder+'/stardata/'+filt+'/'+justnames[sn]+'band_starGlobalOffsets.npz')
                starras = staroffsets['starras']
                stardecs = staroffsets['stardecs']
                starids = staroffsets['starids']
                starglobalids = []
                starglobalras = []
                starglobaldecs = []

		starcatras = []
                starcatdecs= []
		
		for ide in np.unique(np.array(starids)):
                        ww = (starids == ide)
                        starglobalids.append(ide)
                        starglobalras.append(np.median(starras[ww]))
                        starglobaldecs.append(np.median(stardecs[ww]))
                        #nra,ndec = getProperCatRaDec([starglobalras[-1]],[starglobaldecs[-1]],bigcatras,bigcatdecs)
                scampra,scampdec = getProperCatRaDec(starglobalras,starglobaldecs,bigcatras,bigcatdecs)
		offsetra = np.array(starglobalras) - np.array(scampra)
		offsetdec = np.array(starglobaldecs) - np.array(scampdec)
		globalraoffsets.append(offsetra)
		globaldecoffsets.append(offsetdec)
	#print bigcatras[0:10]
	h = np.zeros(len(fakemag))
	#peakmjd = ipt['peakmjd']
	#cs = np.array(finalresults['chisq'])
	#print cs
	#cs = cs[cs<9999999999999]
	#cs = cs[cs>-1]
	#print np.sum(cs)
	#sumchisqs.extend(cs)
	totscounter += 1
	print totscounter

#if dov3:
#	nmags = 31.-2.5*np.log10(nfluxes)
#	nmagerrs =  -2.5*np.log10(nfluxes)+2.5*np.log10(nfluxes+nfluxerrs)
#	domags = 31.-2.5*np.log10(dofluxes)
#	domagerrs = -2.5*np.log10(dofluxes)+2.5*np.log10(dofluxes+dofluxerrs)

print 'PASSED PERCENTAGE:',float(totscounter)/float(len(snnames))
print 'total sne potential',len(snnames)
raw_input()
#START PLOTTING
plt.clf()

out = os.path.join(outfolder,'plots/'+filt+'/')
if not os.path.exists(out):
	os.makedirs(out)

mcmc_fixmags = np.asarray(mcmc_fixmags)
mcmc_fixmag_errs = np.asarray(mcmc_fixmag_errs)
mcmc_floatmags = np.asarray(mcmc_floatmags)
mcmc_floatmag_errs = np.asarray(mcmc_floatmag_errs)

mcmc_gfloatmags = np.asarray(mcmc_gfloatmags)
mcmc_gfloatmag_errs = np.asarray(mcmc_gfloatmag_errs)

mcmc_gpixfloatmags = np.asarray(mcmc_gpixfloatmags)
mcmc_gpixfloatmag_errs = np.asarray(mcmc_gpixfloatmag_errs)

if dov3:
	domags = np.asarray(domags)
	domagerrs = np.asarray(domagerrs)
	nmags = np.asarray(nmags)
	nmagerrs = np.asarray(nmagerrs)
#print np.mean(mcmc_floatmags[np.isfinite(mcmc_floatmags) & np.isfinite(domags)]-domags[np.isfinite(mcmc_floatmags) & np.isfinite(domags)]) 
#raw_input('lll')
# diffim_mags = np.asarray(diffim_mags)
# diffim_mag_errs = np.asarray(diffim_mag_errs)
fakemags = np.asarray(fakemags)
skyerrs = np.asarray(skyerrs)
skysigs = np.asarray(skysigs)
tses = np.asarray(tses)	
#print snnamess
snnamess = np.asarray(snnamess)
#print snnamess
#raw_input()
#mcmc_fluxes = ((smpzpt -mcmc_mags)/2.5)**10.
#diffim_fluxes = ((smpzpt - diffim_mags)/2.5)**10.
#fakefluxes = ((fakezpt - fakemags)/2.5)**10.
hostgal_bandmags = np.asarray(hostgal_bandmags)
mcmc_fixfluxes = np.asarray(mcmc_fixfluxes)
mcmc_fixfluxerrs = np.asarray(mcmc_fixfluxerrs)
mcmc_floatfluxes = np.asarray(mcmc_floatfluxes)
mcmc_floatfluxerrs = np.asarray(mcmc_floatfluxerrs)

mcmc_gfloatfluxes = np.asarray(mcmc_gfloatfluxes)
mcmc_gfloatfluxerrs = np.asarray(mcmc_gfloatfluxerrs)

mcmc_gpixfloatfluxes = np.asarray(mcmc_gpixfloatfluxes)
mcmc_gpixfloatfluxerrs = np.asarray(mcmc_gpixfloatfluxerrs)

diffim_fluxes = np.asarray(diffim_fluxes)
diffim_fluxerrs = np.asarray(diffim_fluxerrs)

diffim_mags = 31.-2.5*np.log10(diffim_fluxes)
#diffim_mag_errs = -2.5*np.log10(diffim_fluxes)+2.5*np.log10(diffim_fluxes)

fake_fluxes = np.asarray(fake_fluxes)
zpt_errs = np.asarray(zpt_errs)
mjds = np.asarray(mjds)
psfs = np.asarray(psfs)
#chis = np.asarray(chis)
#chisqs = np.asarray(chisqs)
#mjdflags = np.asarray(mjdflags)
#fakechis = np.asarray(fakechis)
#fakechisqs = np.asarray(fakechisqs)
#sns = np.asarray(sns,dtype='string')
#psfs = np.asarray(psfs)
fakemags_plus = np.asarray(fakemags_plus)
fakemags_minus = np.asarray(fakemags_minus)
fake_fluxes_plus = np.asarray(fake_fluxes_plus)
fake_fluxes_minus = np.asarray(fake_fluxes_minus)
flags = np.asarray(flags)
zpt_offsets = np.asarray(zpt_offsets)
snraoffs = np.asarray(snraoffs)
sndecoffs = np.asarray(sndecoffs)

ww = (skysigs < 1.) & (flags == 0)
#print snnamess[ww], mjds[ww], skysigs[ww]
#raw_input()

os.system('rm '+out+'*.pdf')
os.system('rm '+out+'*.png')

'''
print 'here1'
plt.hold(True)

i=i+5
plt.figure(i+1)
a = (diffim_mags-fakemags)
b = (mcmc_mags-fakemags)
ax,ay,aystd = bindata(fakemags,a,np.arange(18,30,1))
anumout,apout = numfivesigma(fakemags,a,np.arange(18,30,1))
bx,by,bystd = bindata(fakemags,b,np.arange(18,30,1))
bnumout,bpout = numfivesigma(fakemags,b,np.arange(18,30,1))
plt.scatter(fakemags,a,alpha=.4,color='navy')
plt.scatter(fakemags,b,alpha=.4,color='indianred')
plt.errorbar(ax,ay,aystd,color='navy',label='Diffim #Outliers '+str(anumout)+', '+str(round(apout*100.,3))+'$\%$',fmt='o')
plt.errorbar(bx,by,bystd,color='indianred',label='SMP #Outliers '+str(bnumout)+', '+str(round(bpout*100.,3))+'$\%$',fmt='o')
plt.ylabel('(Fit Mag- Fake Mag)')
plt.xlabel('Fake Mag')
plt.title(filt+' Band')
plt.xlim(19,24)
plt.ylim(-.25,.25)
plt.plot([19,24],[0,0],color='black')
plt.legend(fontsize = 'small')
plt.savefig(out+'fit_mag_resid_'+filt+'.pdf')
print 'here2'
'''
'''
# (Recoved Flux (scene) - Fake Flux)/(Fake Flux) versus Fake Mag
plt.figure(i+2)
a = (mcmc_mags-fakemags)/fakemags
plt.scatter(fakemags,a,alpha=.3)
plt.ylabel('(SMP Mag - Fake Mag)/Fake Mag')
plt.xlabel('Fake Mag')
plt.xlim(18,30)
plt.title(filt+' Band')
plt.savefig(out+'smp_mag_resid_'+filt+'.pdf')
'''
'''
#(Recoved Flux (scene) -Recoved Flux (diff) )/(Fake Flux)) versus Fake Mag
plt.figure(i+3)
a = (mcmc_fluxes-diffim_fluxes)
ax,ay,aystd = bindata(fakemags,a,np.arange(18,30,1))
anumout,apout = numfivesigma(fakemags,a,np.arange(18,30,1))
plt.scatter(fakemags,a,alpha=.2,color='black')
plt.errorbar(ax,ay,aystd,color='blue',fmt='o',label='#Outliers '+str(anumout)+', '+str(round(apout*100.,3))+'$\%$')
plt.ylabel('(SMP Flux - Diffim Flux)')
plt.xlabel('Fake Mag')
plt.xlim(19,24)
#plt.ylim(-.2,.2)
plt.legend(fontsize = 'small')
plt.plot([19,24],[0,0],color='black')
plt.title(filt+' Band')
plt.savefig(out+'smpdiffim_flux_resid_'+filt+'.pdf')
print 'here3'
'''
'''
# (Recoved Flux (diff) - Fake Flux)/(Fake Flux)) versus Fake Mag                                                                                                    
plt.figure(i+4)
a = (diffim_fluxes-fake_fluxes)/diffim_fluxerrs
b = (mcmc_fluxes-fake_fluxes)/mcmc_fluxerrs
ax,ay,aystd = bindata(fakemags,a,np.arange(18,30,1))
anumout,apout = numfivesigma(fakemags,a,np.arange(18,30,1))
bx,by,bystd = bindata(fakemags,b,np.arange(18,30,1))
bnumout,bpout = numfivesigma(fakemags,b,np.arange(18,30,1))
plt.scatter(hostgal_bandmags,a,alpha=.4,color='navy')
plt.scatter(hostgal_bandmags,b,alpha=.4,color='indianred')
plt.errorbar(ax,ay,aystd,color='navy',label='Diffim #Outliers '+str(anumout)+', '+str(round(apout*100.,3))+'$\%$',fmt='o')
plt.errorbar(bx,by,bystd,color='indianred',label='SMP #Outliers '+str(bnumout)+', '+str(round(bpout*100.,3))+'$\%$',fmt='o')
plt.ylabel('(Fit Flux - Fake Flux)/Fit Flux Err')
plt.xlabel('Fake Mag')
plt.xlim(19,24)
plt.ylim(-3.,3.)
plt.title(filt+' Band')
plt.legend(fontsize = 'small')
plt.plot([19,24],[0,0],color='black')
plt.savefig(out+'mag_resid_vs_host_'+filt+'.pdf')
'''
'''
#(Recoved Flux (scene) - Fake Flux)/(Fake Flux)) versus Fake Mag                                                                                                   
plt.figure(i+5)
a = (mcmc_mags-fakemags)/fakemags
plt.scatter(hostgal_bandmags,a,alpha=.3)
plt.ylabel('SMP Mag - Fake Mag)/Fake Mag)')
plt.xlabel('Host Mag')
plt.xlim(18,30)
plt.title(filt+' Band')
plt.savefig(out+'smp_mag_resid_vs_host_'+filt+'.pdf')
'''

#(Recoved Flux (scene) -Recoved Flux (diff) )/(Fake Flux)) versus Fake Mag                                                      
'''
plt.figure(i+6)
a = np.asarray(mcmc_mags-diffim_mags)
ax,ay,aystd = bindata(hostgal_bandmags,a,np.arange(18,30,1))
anumout,apout = numfivesigma(hostgal_bandmags,a,np.arange(18,30,1))
plt.scatter(hostgal_bandmags,a,alpha=.2,color='black')
plt.errorbar(ax,ay,aystd,color='blue',fmt='o',label='#Outliers '+str(anumout)+', '+str(round(apout*100.,3))+'$\%$')
plt.ylabel('(SMP Mag - Diffim Mag)')
plt.xlabel('Host Mag')
plt.xlim(19,24)
plt.legend(fontsize = 'small')
plt.title(filt+' Band')
plt.plot([19,24],[0,0],color='black')
plt.savefig(out+'smpdiffim_mag_resid_vs_host_'+filt+'.pdf')
'''
'''
plt.hold(True)
plt.figure(i+7)
plt.hist(mcmc_fluxes[(mcmc_fluxes > 0.) & (mcmc_fluxes < 9999999. )],100,alpha=.5,label='SMP',color='green')
plt.xlabel('Flux')
plt.hist(diffim_fluxes[(diffim_fluxes > 0.) & (diffim_fluxes < 9999999. )],100,alpha=.5,label='Diffim',color='blue')
plt.ylabel('counts')
plt.legend()
plt.title(filt+' Band')
plt.savefig(out+'flux_hist_'+filt+'.pdf')
'''

#mcmc_std = np.std(mcmc_mags)
#mcmc_mean = np.mean(mcmc_mags)
#diffim_std = np.std(diffim_mags)
#diffim_mean = np.mean(diffim_mags)

#mcmc_out_num = len(mcmc_mags[abs(mcmc_mags-mcmc_mean) > 5.*mcmc_std])
#mcmc_tot = len(mcmc_fluxes)
#mcmc_out_pct = float(mcmc_out_num)/float(mcmc_tot)
#diffim_out_num = len(diffim_mags[abs(diffim_mags-diffim_mean) > 5.*diffim_std])
#diffim_tot = len(diffim_fluxes)
#diffim_out_pct = float(diffim_out_num)/float(diffim_tot)
#print 'here4'

fig = plt.figure(199,figsize=(8,8))
plt.savefig(out+'dump.pdf')
plt.clf()


'''
fig = plt.figure(198,figsize=(8,8))
plt.scatter(np.log10(input_fakeflux),(
	-input_fakeflux)/input_diffim_fluxerrs,alpha=.15,color='red',label='Diffim')
ax,ay,aystd = bindata(np.log10(input_fakeflux),(
	-input_fakeflux)/input_diffim_fluxerrs,np.arange(min(fake_fluxes_plus),max(fake_fluxes_plus),.25))
plt.errorbar(ax,ay,aystd,color='red',fmt='o')
plt.xlabel('Log 10 Fake Flux')
plt.ylim(-5,5)
plt.legend()
plt.plot([min(np.log10(input_fakeflux)),max(np.log10(input_fakeflux))],[0,0],color='black')
plt.ylabel('(Fit Flux - Fake Flux Adjusted IDEAL) / Fit Error')
plt.title(filt+' Band')
plt.savefig(out+'FluxResidStd_TEST_'+filt+'.png')

print out+'FluxResidStd_TEST_'+filt+'.png'
raw_input()
'''
#plt.hold(True)
'''
fig = plt.figure(i+8)
bins = np.arange(18,24,.5)
plt.hist([fakemags,mcmc_mags,diffim_mags],label=['Fake','SMP #Outliers '+str(mcmc_out_num)+', '+str(round(mcmc_out_pct*100.,3))+'$\%$',
		'Diffim #Outliers '+str(diffim_out_num)+', '+str(round(diffim_out_pct*100.,3))+'$\%$'],color=['magenta','red','blue'] ,histtype='bar',bins=bins)
#plt.hist(mcmc_mags[(mcmc_mags > 0.) & (mcmc_mags < 999999.)],20,alpha=1,label='SMP',color='red', histtype='bar', stacked=True,bins=bins)
#plt.hist(diffim_mags[(diffim_mags > 0.) & (diffim_mags < 99999.)],20,alpha=1,label='Diffim',color='green', histtype='bar', stacked=True,bins=bins)
plt.xlabel('Mag')
plt.ylabel('counts')
plt.legend(fontsize = 'small')
plt.title(filt+' Band')
plt.savefig(out+'mag_hist_'+filt+'.png')
'''
print 'here5'


rates = []
for l in acceptedrate:
	if not np.isnan(l):
		rates.append(float(l))

fig = plt.figure(200)
plt.clf()
plt.hist(rates,label='acceptedrate',color='blue' ,histtype='bar',bins=20)
plt.xlabel('Acceptence Rate')
plt.ylabel('Counts')
plt.title(filt+' Band')
plt.xlim(0,1)
plt.savefig(out+'acceptance_hist_'+filt+'.pdf')

'''
fig = plt.figure(201,figsize=(8,8))
plt.scatter(fakemags,(mcmc_fluxes-fake_fluxes_plus)/mcmc_fluxerrs,alpha=.15,color='red',label='SMP')
ax,ay,aystd = bindata(fakemags,(mcmc_fluxes-fake_fluxes_plus)/mcmc_fluxerrs,np.arange(min(fake_fluxes_plus),max(fake_fluxes_plus),.25))
plt.errorbar(ax,ay,aystd,color='red',fmt='o')
plt.scatter(fakemags,(diffim_fluxes-fake_fluxes)/diffim_fluxerrs,alpha=.15,color='blue',label='Diffim')
ax,ay,aystd = bindata(fakemags,(diffim_fluxes-fake_fluxes)/diffim_fluxerrs,np.arange(min(fake_fluxes_plus),max(fake_fluxes_plus),.25))
plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
plt.xlabel('Fake Mag')
plt.ylim(-5,5)
plt.legend()
plt.xlim(18,30)
plt.plot([min(fakemags),max(fakemags)],[0,0],color='black')
plt.ylabel('(Fit Flux - Fake Flux Adjusted IDEAL) / Fit Error')
plt.title(filt+' Band')
plt.savefig(out+'FluxResidStd_vs_true_adjusted_'+filt+'.png')
'''


'''fig = plt.figure(201,figsize=(8,8))
plt.scatter(fakemags,(mcmc_fluxes-fake_fluxes_plus),alpha=.15,color='red',label='SMP')
ax,ay,aystd = bindata(fakemags,(mcmc_fluxes-fake_fluxes_plus),np.arange(min(fake_fluxes_plus),max(fake_fluxes_plus),.25))
plt.errorbar(ax,ay,aystd,color='red',fmt='o')
plt.scatter(fakemags,(diffim_fluxes-fake_fluxes),alpha=.15,color='blue',label='Diffim')
ax,ay,aystd = bindata(fakemags,(diffim_fluxes-fake_fluxes),np.arange(min(fake_fluxes_plus),max(fake_fluxes_plus),.25))
plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
plt.xlabel('Fake Mag')
plt.ylim(-5,5)
plt.legend()
plt.xlim(18,30)
plt.plot([min(fakemags),max(fakemags)],[0,0],color='black')
plt.ylabel('(Fit Flux - Fake Flux Adjusted IDEAL)')
plt.title(filt+' Band')
plt.savefig(out+'FluxResid_vs_true_adjusted_'+filt+'.png')
'''

# fig = plt.figure(30)
# plt.scatter(np.sqrt(tses[fakemags < 24]),mcmc_floatfluxerrs[fakemags<24],alpha=.2)
# plt.plot([0,1400],[0,1400])
# plt.ylabel('mcmc float flux errors')
# plt.xlabel('total sky sigma')
# plt.savefig(out+'tse_vs_fluxerr'+filt+'.png')
# print out+'tse_vs_fluxerr'+filt+'.png'
# fig = plt.figure(29)
# plt.scatter(np.sqrt(tses[fakemags < 24]),mcmc_fixfluxerrs[fakemags<24],alpha=.2)
# plt.plot([0,1400],[0,1400])
# plt.ylabel('mcmc fix flux errors')
# plt.xlabel('total sky sigma')
# plt.savefig(out+'tse_vs_fixfluxerr'+filt+'.png')
# print out+'tse_vs_fixfluxerr'+filt+'.png'

# fig = plt.figure(31)
# plt.scatter(skyerrs[fakemags < 24],mcmc_floatfluxerrs[fakemags<24],alpha=.2)
# plt.plot([10,60],[10,60])
# plt.xlim(10,60)
# plt.ylim(10,450)
# plt.ylabel('mcmc float flux errors')
# plt.xlabel('skyerrs')
# plt.savefig(out+'skyerr_vs_fluxerr'+filt+'.png')
# print out+'skyerr_vs_fluxerr'+filt+'.png'

# fig, ax = plt.subplots(4,figsize=(8,12.5))

# mcmc_fixstdes = (mcmc_fixfluxes-fake_fluxes)/mcmc_fixfluxerrs
# mcmc_floatstdes = (mcmc_floatfluxes-fake_fluxes)/mcmc_floatfluxerrs
# mcmc_gfloatstdes = (mcmc_gfloatfluxes-fake_fluxes)/mcmc_gfloatfluxerrs
# mcmc_gpixfloatstdes = (mcmc_gpixfloatfluxes-fake_fluxes)/mcmc_gpixfloatfluxerrs

# ww1 = (fakemags > 50) & (flags == 0) & (psfs < 1.5) & (mcmc_floatfluxes != 0) & (~np.isnan(mcmc_floatstdes)) & (~np.isinf(mcmc_floatstdes))
# ww2 = (fakemags > 50) & (flags == 0) & (psfs > 1.5) & (psfs < 2.) & (mcmc_floatfluxes != 0) & (~np.isnan(mcmc_floatstdes)) & (~np.isinf(mcmc_floatstdes))
# ww3 = (fakemags > 50) & (flags == 0) & (psfs > 2.) & (psfs < 2.5) & (mcmc_floatfluxes != 0) & (~np.isnan(mcmc_floatstdes)) & (~np.isinf(mcmc_floatstdes))
# ww4 = (fakemags > 50) & (flags == 0) & (psfs > 2.5) & (mcmc_floatfluxes != 0) & (~np.isnan(mcmc_floatstdes)) & (~np.isinf(mcmc_floatstdes))

# axs = ax.ravel()
# axs[0].hist(mcmc_floatstdes[ww1],bins=np.arange(-5,5,.25),color='green')
# axs[1].hist(mcmc_floatstdes[ww2],bins=np.arange(-5,5,.25),color='green')
# axs[2].hist(mcmc_floatstdes[ww3],bins=np.arange(-5,5,.25),color='green')
# axs[3].hist(mcmc_floatstdes[ww4],bins=np.arange(-5,5,.25),color='green')

# plt.xlabel('(Fit Flux - Fake Flux)/Fit Error')
# plt.ylabel('Normalized Counts')
# axs[0].set_title('fwhm<1.5')
# axs[1].set_title('1.5<fwhm<2')
# axs[2].set_title('2<fwhm<2.5')
# axs[3].set_title('2.5<fwhm<3')
# plt.savefig(out+'fwhm_stdHist_'+filt+'.png')


'''
fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,figsize=(8,12.5))
from scipy.stats.kde import gaussian_kde

mcmc_fixstdeso = (mcmc_fixfluxes-fake_fluxes)/mcmc_fixfluxerrs
mcmc_floatstdeso = (mcmc_floatfluxes-fake_fluxes)/mcmc_floatfluxerrs
#mcmc_gfloatstdeso = (mcmc_gfloatfluxes-fake_fluxes)/mcmc_gfloatfluxerrs
if dopix:
	mcmc_gpixfloatstdeso = (mcmc_gpixfloatfluxes-fake_fluxes)/mcmc_gpixfloatfluxerrs
else:
	mcmc_gpixfloatstdeso = mcmc_fixstdeso*0.-999
# mcmc_fixstdes = (mcmc_fixfluxes-fake_fluxes)/tses
# mcmc_floatstdes = (mcmc_floatfluxes-fake_fluxes)/tses
# mcmc_gfloatstdes = (mcmc_gfloatfluxes-fake_fluxes)/tses
# mcmc_gpixfloatstdes = (mcmc_gpixfloatfluxes-fake_fluxes)/tses

#for i in hostgal_bandmags:
#	print i
#print 'sopp'
#raw_input()
hostgal_bandmags = np.array(hostgal_bandmags)
hostgal_bandmags_fix = hostgal_bandmags[(fakemags < 50) & (flags == 0)]
hostgal_bandmags_float = hostgal_bandmags[(fakemags < 50) & (flags == 0) & (mcmc_floatfluxes != 0)]
if dopix:
	hostgal_bandmags_pix = hostgal_bandmags[(fakemags < 50) & (flags == 0) & (mcmc_gpixfloatfluxes != 0)]
else:
	hostgal_bandmags_pix = hostgal_bandmags_fix*0.-999

mcmc_fixstdes = mcmc_fixstdeso[(fakemags < 50) & (flags == 0)]

mcmc_floatstdes = mcmc_floatstdeso[(fakemags < 50) & (flags == 0) & (mcmc_floatfluxes != 0)]

#mcmc_gfloatstdes = mcmc_gfloatstdeso[(fakemags < 50) & (flags == 0) & (mcmc_gfloatfluxes != 0)]
#print mcmc_gfloatstdes

if dopix:
	mcmc_gpixfloatstdes = mcmc_gpixfloatstdeso[(fakemags < 50) & (flags == 0) & (mcmc_gpixfloatfluxes != 0)]
else:
	mcmc_gpixfloatstdes = mcmc_fixstdes*0.-999

# diffim_stdes = (diffim_fluxes - fake_fluxes)/diffim_fluxerrs
hostgal_bandmags_fix = hostgal_bandmags_fix[~np.isnan(mcmc_fixstdes)]
mcmc_fixstdes = mcmc_fixstdes[~np.isnan(mcmc_fixstdes)]
hostgal_bandmags_fix = hostgal_bandmags_fix[~np.isinf(mcmc_fixstdes)]
mcmc_fixstdes = mcmc_fixstdes[~np.isinf(mcmc_fixstdes)]

hostgal_bandmags_float = hostgal_bandmags_float[~np.isnan(mcmc_floatstdes)]
mcmc_floatstdes = mcmc_floatstdes[~np.isnan(mcmc_floatstdes)]
hostgal_bandmags_float = hostgal_bandmags_float[~np.isinf(mcmc_floatstdes)]
mcmc_floatstdes = mcmc_floatstdes[~np.isinf(mcmc_floatstdes)]



#mcmc_gfloatstdes = mcmc_gfloatstdes[~np.isnan(mcmc_gfloatstdes)]
#mcmc_gfloatstdes = mcmc_gfloatstdes[~np.isinf(mcmc_gfloatstdes)]

mcmc_gpixfloatstdes = mcmc_gpixfloatstdes[~np.isnan(mcmc_gpixfloatstdes)]
mcmc_gpixfloatstdes = mcmc_gpixfloatstdes[~np.isinf(mcmc_gpixfloatstdes)]





# diffim_stdes = diffim_stdes[~np.isnan(diffim_stdes)]
# diffim_stdes = diffim_stdes[~np.isinf(diffim_stdes)]

#mcmc_fixpdf = gaussian_kde(mcmc_fixstdes,)
#mcmc_floatpdf = gaussian_kde(mcmc_floatstdes,)
#mcmc_gfloatpdf = gaussian_kde(mcmc_gfloatstdes,)
#mcmc_gpixfloatpdf = gaussian_kde(mcmc_gpixfloatstdes,)



#diffim_pdf = gaussian_kde(diffim_stdes,)
x_grid = np.linspace(-5,5,num=200)

# plt.plot(x_grid,mcmc_fixpdf(x_grid),'r',label='SMP Brute Force',color='red')
# plt.plot(x_grid,mcmc_floatpdf(x_grid),'r',label='SMP Float SN',color='green')
# plt.plot(x_grid,mcmc_gfloatpdf(x_grid),'r',label='SMP Float SN + .27 Resampling',color='blue')
# plt.plot(x_grid,mcmc_gpixfloatpdf(x_grid),'r',label='SMP Float SN + .5 Resampling',color='gold')
#ax1.hist(mcmc_fixstdes,bins=np.arange(-5,5,.25),color='red',label='SMP Brute Force',normed=True)
ax2.hist(mcmc_floatstdes,bins=np.arange(-5,5,.25),color='green',label='SMP Float SN',normed=True)
#ax3.hist(mcmc_gfloatstdes,bins=np.arange(-5,5,.25),color='blue',label='SMP Float SN + .27 Resampling',normed=True)
#if dopix:
#	ax4.hist(mcmc_gpixfloatstdes,bins=np.arange(-5,5,.25),color='gold',label='SMP Galsim .5" Fit RA+DEC',normed=True)

# plt.plot(x_grid,diffim_pdf(x_grid),'r',label='Diffim',color='blue')

import matplotlib.mlab as mlab
import math
mean = 0
variance = 1
sigma = math.sqrt(variance)
x = np.arange(-5,5,.1)
ax1.plot(x,mlab.normpdf(x,mean,sigma),color='black',label='Gaussian Normal')
ax2.plot(x,mlab.normpdf(x,mean,sigma),color='black')
ax3.plot(x,mlab.normpdf(x,mean,sigma),color='black')
ax4.plot(x,mlab.normpdf(x,mean,sigma),color='black')



plt.xlabel('(Fit Flux - Fake Flux)/Fit Error')
plt.ylabel('Normalized Counts')
ax1.set_title('SMP Brute Force')
ax2.set_title('SMP Float SN')
ax3.set_title('SMP Float SN + .27 Resampling')
ax4.set_title('SMP Float SN + .5 Resampling')
plt.savefig(out+'stdHist_onlyFlux'+filt+'.png')

'''
'''
fig, axs = plt.subplots(4,2,figsize=(10,20))
axs = axs.ravel()
axs[0].hist(mcmc_fixstdes[(hostgal_bandmags_fix > 16) & (hostgal_bandmags_fix < 20)],bins=np.arange(-5.175,5,.25),color='red',label='SMP Brute Force',normed=True)
axs[1].hist(mcmc_floatstdes[(hostgal_bandmags_float > 16) & (hostgal_bandmags_float < 20)],bins=np.arange(-5.175,5,.25),color='green',label='SMP Float SN',normed=True)
axs[2].hist(mcmc_fixstdes[(hostgal_bandmags_fix > 20) & (hostgal_bandmags_fix < 22)],bins=np.arange(-5.175,5,.25),color='red',label='SMP Brute Force',normed=True)
axs[3].hist(mcmc_floatstdes[(hostgal_bandmags_float > 20) & (hostgal_bandmags_float < 22)],bins=np.arange(-5.175,5,.25),color='green',label='SMP Float SN',normed=True)
axs[4].hist(mcmc_fixstdes[(hostgal_bandmags_fix > 22) & (hostgal_bandmags_fix < 24)],bins=np.arange(-5.175,5,.25),color='red',label='SMP Brute Force',normed=True)
axs[5].hist(mcmc_floatstdes[(hostgal_bandmags_float > 22) & (hostgal_bandmags_float < 24)],bins=np.arange(-5.175,5,.25),color='green',label='SMP Float SN',normed=True)
axs[6].hist(mcmc_fixstdes[(hostgal_bandmags_fix > 24) & (hostgal_bandmags_fix < 200)],bins=np.arange(-5.175,5,.25),color='red',label='SMP Brute Force',normed=True)
axs[7].hist(mcmc_floatstdes[(hostgal_bandmags_float > 24) & (hostgal_bandmags_float < 200)],bins=np.arange(-5.175,5,.25),color='green',label='SMP Float SN',normed=True)
plt.legend(loc=4)
axs[0].set_xlabel('stdev 16<hostmag<20')
axs[1].set_xlabel('stdev 16<hostmag<20')
axs[2].set_xlabel('stdev 20<hostmag<22')
axs[3].set_xlabel('stdev 20<hostmag<22')
axs[4].set_xlabel('stdev 22<hostmag<24')
axs[5].set_xlabel('stdev 22<hostmag<24')
axs[6].set_xlabel('stdev 24<hostmag<26')
axs[6].legend(loc=4)
axs[7].set_xlabel('stdev 24<hostmag<26')
axs[0].set_ylabel('Normalized Counts')
axs[1].set_ylabel('Normalized Counts')
axs[2].set_ylabel('Normalized Counts')
axs[3].set_ylabel('Normalized Counts')
axs[4].set_ylabel('Normalized Counts')
axs[5].set_ylabel('Normalized Counts')
axs[6].set_ylabel('Normalized Counts')
axs[7].set_ylabel('Normalized Counts')

#ax3.hist(mcmc_gfloatstdes,bins=np.arange(-5,5,.25),color='blue',label='SMP Float SN + .27 Resampling',normed=True)
#ax4.hist(mcmc_gpixfloatstdes,bins=np.arange(-5,5,.25),color='gold',label='SMP Float SN + .5 Resampling')

# plt.plot(x_grid,diffim_pdf(x_grid),'r',label='Diffim',color='blue')

import matplotlib.mlab as mlab
import math
mean = 0
variance = 1
sigma = math.sqrt(variance)
x = np.arange(-5,5,.1)
axs[0].plot(x,mlab.normpdf(x,mean,sigma),color='black')
axs[1].plot(x,mlab.normpdf(x,mean,sigma),color='black')
axs[2].plot(x,mlab.normpdf(x,mean,sigma),color='black')
axs[3].plot(x,mlab.normpdf(x,mean,sigma),color='black')
axs[4].plot(x,mlab.normpdf(x,mean,sigma),color='black')
axs[5].plot(x,mlab.normpdf(x,mean,sigma),color='black')
axs[6].plot(x,mlab.normpdf(x,mean,sigma),color='black')
axs[7].plot(x,mlab.normpdf(x,mean,sigma),color='black')

axs[0].set_xlim(-4,4)
axs[1].set_xlim(-4,4)
axs[2].set_xlim(-4,4)
axs[3].set_xlim(-4,4)
axs[4].set_xlim(-4,4)
axs[5].set_xlim(-4,4)
axs[6].set_xlim(-4,4)
axs[7].set_xlim(-4,4)
#plt.xlabel('(Fit Flux - Fake Flux)/Fit Error')
#plt.ylabel('Normalized Counts')
# ax1.set_title('SMP Brute Force')
# ax2.set_title('SMP Float SN')
# ax3.set_title('SMP Float SN + .27 Resampling')
# ax4.set_title('SMP Float SN + .5 Resampling')
plt.savefig(out+'stdHist_vsHostgalMag_'+filt+'.png')
'''
np.savez(outfolder+'/np_data/v4_'+filt+'.npz'
	 , flux=mcmc_floatfluxes
	 , fluxerr = mcmc_floatfluxerrs
	 #, bfflux = mcmc_fixfluxes
	 #, bffluxerr = mcmc_fixfluxes
	 #, pixflux = mcmc_gpixfloatfluxes
	 #, pixfluxerr = mcmc_gpixfloatfluxerrs
	 , fakemags = fakemags_minus
	 , fakefluxes = fake_fluxes
	 , mjd = mjds
	 , sns = snnamess
	 , hostgal_bandmags = hostgal_bandmags
	 , diffim_fluxes = diffim_fluxes
	 , diffim_fluxerrs = diffim_fluxerrs
	 #, stds=mcmc_floatstdes
	 , globalraoffsets = globalraoffsets
	 , globaldecoffsets = globaldecoffsets
	 , globalsnra = globalsnra
	 , globalsndec = globalsndec
	 , snraoff = snraoffs
	 , sndecoff = sndecoffs
         #, nmags = nmags
	 #, nmagerrs = nmagerrs
	 #, domags = domags
	 #, domagerrs = domagerrs
	 #, nfluxes = nfluxes
	 #, nfluxerrs = nfluxerrs
	 #, dofluxes = dofluxes
	 #, dofluxerrs = dofluxerrs
	 )

fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,figsize=(8,12.5))
from scipy.stats.kde import gaussian_kde

mcmc_fixresid = (mcmc_fixfluxes-fake_fluxes)
mcmc_floatresid = (mcmc_floatfluxes-fake_fluxes)
#mcmc_gfloatresid = (mcmc_gfloatfluxes-fake_fluxes)
#mcmc_gpixfloatresid = (mcmc_gpixfloatfluxes-fake_fluxes)

# mcmc_fixstdes = (mcmc_fixfluxes-fake_fluxes)/tses
# mcmc_floatstdes = (mcmc_floatfluxes-fake_fluxes)/tses
# mcmc_gfloatstdes = (mcmc_gfloatfluxes-fake_fluxes)/tses
# mcmc_gpixfloatstdes = (mcmc_gpixfloatfluxes-fake_fluxes)/tses

#print snnamess[ (mcmc_floatstdeso > 3) & (mcmc_floatresid < 1000) & (flags == 0) & (fakemags > 50) & (mcmc_floatfluxes != 0) ][-3], mjds[ (mcmc_floatstdeso > 3) & (mcmc_floatresid < 1000) & (flags == 0) & (fakemags > 50) & (mcmc_floatfluxes != 0) ][-3], mcmc_floatresid[(mcmc_floatstdeso > 3) & (mcmc_floatresid < 1000) & (flags == 0) & (fakemags > 50) & (mcmc_floatfluxes != 0)][-3]
#print mcmc_floatstdeso[ (mcmc_floatstdeso > 3) & (mcmc_floatresid < 1000) & (flags == 0) & (fakemags > 50) & (mcmc_floatfluxes != 0) ][-3]
#print mcmc_floatfluxes[ (mcmc_floatstdeso > 3) & (mcmc_floatresid < 1000) & (flags == 0) & (fakemags > 50) & (mcmc_floatfluxes != 0) ][-3]
#print mcmc_floatfluxerrs[ (mcmc_floatstdeso > 3) & (mcmc_floatresid < 1000) & (flags == 0) & (fakemags > 50) & (mcmc_floatfluxes != 0) ][-3]
#raw_input()

#hostgal_bandmags = np.array(hostgal_bandmags)
hostgal_bandmags_fix = hostgal_bandmags[(fakemags > -50) & (flags == 0)]
hostgal_bandmags_float = hostgal_bandmags[(fakemags > -50) & (flags == 0) & (mcmc_floatfluxes != 0)]

mcmc_fixresid = mcmc_fixresid[(fakemags > -50) & (flags == 0)]

mcmc_floatresid = mcmc_floatresid[(fakemags > -50) & (flags == 0) & (mcmc_floatfluxes != 0)]

#mcmc_gfloatresid = mcmc_gfloatresid[(fakemags > 50) & (flags == 0) & (mcmc_gfloatfluxes != 0)]

#mcmc_gpixfloatresid = mcmc_gpixfloatresid[(fakemags > 50) & (flags == 0) & (mcmc_gpixfloatfluxes != 0)]





# diffim_stdes = (diffim_fluxes - fake_fluxes)/diffim_fluxerrs
print mcmc_floatresid
try:
	hostgal_bandmags_fix = hostgal_bandmags_fix[~np.isnan(mcmc_fixstdes)]
	mcmc_fixstdes = mcmc_fixresid[~np.isnan(mcmc_fixstdes)]
except:
	pass
try:
	hostgal_bandmags_fix = hostgal_bandmags_fix[~np.isinf(mcmc_fixstdes)]
	mcmc_fixstdes = mcmc_fixresid[~np.isinf(mcmc_fixstdes)]
except:
	pass
try:
	hostgal_bandmags_float = hostgal_bandmags_float[~np.isnan(mcmc_floatstdes)]
	mcmc_floatstdes = mcmc_floatresid[~np.isnan(mcmc_floatstdes)]
except:
	pass
try:
	hostgal_bandmags_float = hostgal_bandmags_float[~np.isinf(mcmc_floatstdes)]
	mcmc_floatstdes = mcmc_floatresid[~np.isinf(mcmc_floatstdes)]
except:
	pass
#mcmc_gfloatstdes = mcmc_gfloatresid[~np.isnan(mcmc_gfloatstdes)]
#mcmc_gfloatstdes = mcmc_gfloatresid[~np.isinf(mcmc_gfloatstdes)]

#mcmc_gpixfloatresid = mcmc_gpixfloatresid[~np.isnan(mcmc_gpixfloatresid)]
#mcmc_gpixfloatresid = mcmc_gpixfloatresid[~np.isinf(mcmc_gpixfloatresid)]


ax1.hist(mcmc_fixresid,bins=np.arange(-1000,1000,50),color='red',label='SMP Brute Force')
ax2.hist(mcmc_floatresid,bins=np.arange(-1000,1000,50),color='green',label='SMP Float SN')
#ax3.hist(mcmc_gfloatresid,bins=np.arange(-1000,1000,50),color='blue',label='SMP Float SN + .27 Resampling')
#ax4.hist(mcmc_gpixfloatresid,bins=np.arange(-1000,1000,50),color='gold',label='SMP Float SN + .5 Resampling')


plt.xlabel('(Fit Flux - Fake Flux)')
plt.ylabel('Normalized Counts')
ax1.set_title('SMP Brute Force')
ax2.set_title('SMP Float SN')
ax3.set_title('SMP Float SN + .27 Resampling')
ax4.set_title('SMP Float SN + .5 Resampling')
plt.savefig(out+'residHist_'+filt+'.png')





np.savez('reallyv2float.npz',flux=mcmc_floatfluxes
					  , fluxerr = mcmc_floatfluxerrs
					  , fakemags = fakemags
					  , mjd = mjds
					  , sns = snnamess
					  , hostgal_bandmags = hostgal_bandmags
					  , resid = mcmc_floatresid
					  , stds=mcmc_floatstdes)














'''
fig, axs = plt.subplots(4,2,figsize=(10,20))
axs = axs.ravel()
#axs[0].hist(mcmc_fixstdes[(hostgal_bandmags_fix > 16) & (hostgal_bandmags_fix < 20)],bins=np.arange(-1025,1000,50),color='red',label='SMP Brute Force',normed=True)
axs[1].hist(mcmc_floatstdes[(hostgal_bandmags_float > 16) & (hostgal_bandmags_float < 20)],bins=np.arange(-1025,1000,50),color='green',label='SMP Float SN',normed=True)
#axs[2].hist(mcmc_fixstdes[(hostgal_bandmags_fix > 20) & (hostgal_bandmags_fix < 22)],bins=np.arange(-1025,1000,50),color='red',label='SMP Brute Force',normed=True)
axs[3].hist(mcmc_floatstdes[(hostgal_bandmags_float > 20) & (hostgal_bandmags_float < 22)],bins=np.arange(-1025,1000,50),color='green',label='SMP Float SN',normed=True)
#axs[4].hist(mcmc_fixstdes[(hostgal_bandmags_fix > 22) & (hostgal_bandmags_fix < 24)],bins=np.arange(-1025,1000,50),color='red',label='SMP Brute Force',normed=True)
axs[5].hist(mcmc_floatstdes[(hostgal_bandmags_float > 22) & (hostgal_bandmags_float < 24)],bins=np.arange(-1025,1000,50),color='green',label='SMP Float SN',normed=True)
#axs[6].hist(mcmc_fixstdes[(hostgal_bandmags_fix > 24) & (hostgal_bandmags_fix < 200)],bins=np.arange(-1025,1000,50),color='red',label='SMP Brute Force',normed=True)
axs[7].hist(mcmc_floatstdes[(hostgal_bandmags_float > 24) & (hostgal_bandmags_float < 200)],bins=np.arange(-1025,1000,50),color='green',label='SMP Float SN',normed=True)
plt.legend(loc=4)
axs[0].set_xlabel('residuals 16<hostmag<20')
axs[1].set_xlabel('residuals 16<hostmag<20')
axs[2].set_xlabel('residuals 20<hostmag<22')
axs[3].set_xlabel('residuals 20<hostmag<22')
axs[4].set_xlabel('residuals 22<hostmag<24')
axs[5].set_xlabel('residuals 22<hostmag<24')
axs[6].set_xlabel('residuals 24<hostmag<26')
axs[6].legend(loc=4)
axs[7].set_xlabel('residuals 24<hostmag<26')
axs[0].set_ylabel('Normalized Counts')
axs[1].set_ylabel('Normalized Counts')
axs[2].set_ylabel('Normalized Counts')
axs[3].set_ylabel('Normalized Counts')
axs[4].set_ylabel('Normalized Counts')
axs[5].set_ylabel('Normalized Counts')
axs[6].set_ylabel('Normalized Counts')
axs[7].set_ylabel('Normalized Counts')

mean = 0
variance = 200**2
sigma = math.sqrt(variance)
x = np.arange(-1000,1000,.1)
axs[0].plot(x,mlab.normpdf(x,mean,sigma),color='black')
axs[1].plot(x,mlab.normpdf(x,mean,sigma),color='black')
axs[2].plot(x,mlab.normpdf(x,mean,sigma),color='black')
axs[3].plot(x,mlab.normpdf(x,mean,sigma),color='black')
axs[4].plot(x,mlab.normpdf(x,mean,sigma),color='black')
axs[5].plot(x,mlab.normpdf(x,mean,sigma),color='black')
axs[6].plot(x,mlab.normpdf(x,mean,sigma),color='black')
axs[7].plot(x,mlab.normpdf(x,mean,sigma),color='black')

plt.savefig(out+'residHist_vsHostMag_'+filt+'.png')
'''








hostgal_bandmags_fix = hostgal_bandmags[(fakemags > -50) & (flags == 0)]
hostgal_bandmags_float = hostgal_bandmags[(fakemags > -50) & (flags == 0) & (mcmc_floatfluxes != 0)]

mcmc_fixfe = mcmc_fixfluxerrs[(fakemags > -50) & (flags == 0)]

mcmc_floatfe = mcmc_floatfluxerrs[(fakemags > -50) & (flags == 0) & (mcmc_floatfluxes != 0)]


hostgal_bandmags_fix = hostgal_bandmags_fix[~np.isnan(mcmc_fixfe)]
mcmc_fixfe = mcmc_fixfe[~np.isnan(mcmc_fixfe)]
hostgal_bandmags_fix = hostgal_bandmags_fix[~np.isinf(mcmc_fixfe)]
mcmc_fixfe = mcmc_fixfe[~np.isinf(mcmc_fixfe)]

hostgal_bandmags_float = hostgal_bandmags_float[~np.isnan(mcmc_floatfe)]
mcmc_floatfe = mcmc_floatfe[~np.isnan(mcmc_floatfe)]
hostgal_bandmags_float = hostgal_bandmags_float[~np.isinf(mcmc_floatfe)]
mcmc_floatfe = mcmc_floatfe[~np.isinf(mcmc_floatfe)]

fig, axs = plt.subplots(2,figsize=(15,25))
axs = axs.ravel()
axs[0].hist([mcmc_fixfe[(hostgal_bandmags_fix > 16) & (hostgal_bandmags_fix < 20)],
			mcmc_fixfe[(hostgal_bandmags_fix > 20) & (hostgal_bandmags_fix < 22)],
			mcmc_fixfe[(hostgal_bandmags_fix > 22) & (hostgal_bandmags_fix < 24)],
			mcmc_fixfe[(hostgal_bandmags_fix > 24) & (hostgal_bandmags_fix < 26)]]
			,bins=np.arange(0,500,100),
			label=['SMP Brute Force Err 16<hm<20','SMP Brute Force Err 20<hm<22',
			'SMP Brute Force Err 22<hm<24','SMP Brute Force Err 24<hm<26'],normed=True)
axs[1].hist([mcmc_floatfe[(hostgal_bandmags_float > 16) & (hostgal_bandmags_float < 20)],
			mcmc_floatfe[(hostgal_bandmags_float > 20) & (hostgal_bandmags_float < 22)],
			mcmc_floatfe[(hostgal_bandmags_float > 22) & (hostgal_bandmags_float < 24)],
			mcmc_floatfe[(hostgal_bandmags_float > 24) & (hostgal_bandmags_float < 26)]]
			,bins=np.arange(0,500,100),
			label=['SMP Float Err 16<hm<20','SMP Float Err 20<hm<22',
			'SMP Float Err 22<hm<24','SMP Float Err 24<hm<26'],normed=True)

axs[0].legend()
axs[0].set_xlabel('Uncertainty')
axs[0].set_ylabel('Normalized Counts')
axs[1].legend()
axs[1].set_xlabel('Uncertainty')
axs[1].set_ylabel('Normalized Counts')
plt.savefig(out+'err_vsHostMag_'+filt+'.png')







plt.clf()
plt.scatter(diffim_fluxerrs[fakemags<24],mcmc_floatfluxerrs[fakemags<24])
plt.plot([100,600],[100,600])
plt.ylim(100,600)
plt.xlim(100,600)
plt.xlabel('diffim flux errros')
plt.ylabel('mcmc_fixfluxerr')
plt.savefig(out+'fluxerr_comparo.png')


print 'here8'

plt.clf()
a = np.array(mcmc_fixfluxerrs)
print a
b = np.array(mcmc_floatfluxerrs)
#	c = np.array(mcmc_gfloatfluxerrs)
#d = np.array(mcmc_gpixfloatfluxerrs)

a = a[fakemags>50]
b = b[fakemags>50]
#c = c[fakemags>50]
#d = d[fakemags>50]

a = a[~np.isnan(a)]
a = a[~np.isinf(a)]
b = b[~np.isnan(b)]
b = b[~np.isinf(b)]
#c = c[~np.isnan(c)]
#c = c[~np.isinf(c)]
#d = d[~np.isnan(d)]
#d = d[~np.isinf(d)]



plt.hist([a,b],bins=np.arange(0,1000,100),label=['fix errors','float'],normed=True)
plt.legend()
plt.xlabel('Error in Counts')
plt.ylabel('#')
plt.savefig(out+'errHist_'+filt+'.png')
print out+'errHist_'+filt+'.png'

'''

aa = diffim_fluxerrs>0.
bb = mcmc_floatfluxerrs>0.

diffim_fluxerrs = diffim_fluxerrs[aa & bb]
mcmc_floatfluxerrs = mcmc_floatfluxerrs[aa & bb]
diffim_fluxes = diffim_fluxes[aa & bb]


plt.clf()
plt.scatter(diffim_fluxes,diffim_fluxerrs-mcmc_floatfluxerrs,alpha=.2)
plt.plot([min(diffim_fluxes),max(diffim_fluxes)],[0,0],color='black')
plt.xlabel('Diffim Flux')
plt.ylabel('Diffim Flux Error - SMP Flux Error')
plt.savefig(out+'uncertaintyScatter_'+filt+'.png')
#mcmc_floatfluxerrs = mcmc_floatfluxerrs[mcmc_floatfluxerrs>0.]
#mcmc_fixfluxerrs = mcmc_fixfluxerrs[mcmc_fixfluxerrs>0.]
#diffim_fluxerrs = diffim_fluxerrs[diffim_fluxerrs>0.]

plt.clf()
plt.hist([mcmc_floatfluxerrs,mcmc_fixfluxerrs,diffim_fluxerrs],bins=np.arange(0,600,25),label=['SMP Float Uncertainties','SMP Fix Uncertainties','Diffim Uncertainties'],normed=True)
plt.xlabel('Uncertainty (Flux)')
plt.ylabel('Count')
plt.legend()
plt.savefig(out+'uncertaintyHist_'+filt+'.png')
print out+'uncertaintyHist_'+filt+'.png'
'''

# sumchisqs = np.array(sumchisqs)
# print sumchisqs.shape
# sumchisqs = sumchisqs[~np.nan(sumchisqs)]
# sumchisqs = sumchisqs[~np.inf(sumchisqs)]
weights = np.ones_like(sumchisqs)/float(len(sumchisqs))

# plt.clf()
# plt.hist(sumchisqs,bins=np.arange(0,5,.1),label='SMP Chisqs',weights=weights)
# plt.plot(x,mlab.normpdf(x,mean,sigma),color='black',label='Gaussian Normal')
# plt.xlabel('Total Chisq For Each SNe Epoch')
# plt.ylabel('Count')
# plt.xlim(0,5)
# plt.legend()
# plt.savefig(out+'chisqHist_'+filt+'.png')
# print out+'chisqHist_'+filt+'.png'

# aa = (mcmc_fixfluxes != 0)
# bb = (mcmc_floatfluxes != 0)
# cc = (diffim_fluxes != 0)

# mcmc_fixfluxes = mcmc_fixfluxes[aa & bb & cc]
# mcmc_floatfluxes = mcmc_floatfluxes[aa & bb & cc]
# diffim_fluxes = diffim_fluxes[aa & bb & cc]

# plt.clf()
# plt.scatter(diffim_fluxes,diffim_fluxes-mcmc_floatfluxes,alpha=.2)
# plt.plot([min(mcmc_floatfluxes),max(mcmc_floatfluxes)],[0,0],color='black')
# plt.xlabel('SMP float galaxy float SN flux')
# plt.ylabel('Diffim flux - SMP Flux')
# plt.xlim(-5000,40000)
# plt.ylim(-1000,1000)
# plt.savefig(out+'diffimVsFloat_'+filt+'.png')
# print out+'diffimVsFloat_'+filt+'.png'
'''
fig = plt.figure(i+10)
bins = np.arange(0,5,.15)
plt.hist(diffim_mag_errs/mcmc_mag_errs,color=['red'] ,histtype='bar',bins=bins)
plt.xlabel('diffim_mag_errs/mcmc_mag_errs')
plt.ylabel('counts')
plt.xlim(0,5)
plt.title(filt+' Band')
plt.savefig(out+'errmag_comparo_hist_'+filt+'.pdf')
print 'here8'

fig = plt.figure(i+11)
plt.scatter(mcmc_fluxerrs,diffim_fluxerrs ,alpha=.3,color='black')
plt.xlabel('SMP Err (Counts)')
plt.ylabel('Diffim Err (Counts)')
plt.xlim(0.,50.)
plt.ylim(0.,50)
plt.plot([0,50],[0,50],color='black')
plt.title(filt+' Band')
plt.savefig(out+'errflux_scatter_'+filt+'.pdf')
print 'here9'

fig = plt.figure(i+12)
plt.scatter(fakemags,mcmc_mag_errs-diffim_mag_errs ,alpha=.15,color='black')
ax,ay,aystd = bindata(fakemags,mcmc_mag_errs-diffim_mag_errs,np.arange(18,30,1))
plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
plt.xlabel('Fake Mag')
plt.ylabel('SMP Err - Diffim Err')
#plt.xlim(-.,1.)
#plt.ylim(0.,1.)
plt.xlim(19,24)
plt.plot([19,24],[0,0],color='black')
plt.title(filt+' Band')
plt.savefig(out+'errmag_scatter2_'+filt+'.pdf')
print 'here10'



#fig = plt.figure(i+15)
#plt.scatter(fakemags,(mcmc_mags-diffim_mags)/mcmc_mag_errs,alpha=.15,color='black')
#ax,ay,aystd = bindata(fakemags,(mcmc_mags-diffim_mags)/mcmc_mag_errs,np.arange(18,30,1))
#plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
#plt.xlabel('Fake Mag')
#plt.ylabel('(SMP Mag - Diffim Mag)/SMP Mag Error')
#plt.xlim(-.,1.)
#plt.ylim(-5,5.)
#plt.xlim(19,24)
#plt.plot([19,24],[0,0],color='black')
#plt.title(filt+' Band')
#plt.savefig(out+'SMP_minus_Diffim_over_err_'+filt+'.pdf')
print 'here11'
'''
'''
fig = plt.figure(i+16)
bins = np.arange(-10,10,.3)
plt.hist((mcmc_mags-diffim_mags)/mcmc_mag_errs,color='red' ,histtype='bar',bins=bins,normed=True)
import matplotlib.mlab as mlab
import math
mean = 0
variance = 1
sigma = math.sqrt(variance)
x = np.arange(-10,10,.01)
plt.plot(x,mlab.normpdf(x,mean,sigma),color='black')
plt.xlabel('diffim_mag_errs/mcmc_mag_errs')
plt.ylabel('Normalized Counts')
#plt.xlim(0,5)
plt.title(filt+' Band')
plt.savefig(out+'SMP_minus_Diffim_over_err_hist_'+filt+'.pdf')
print 'here12'
'''
#make histograms of -2.5*alog10((Recoved Flux - Fake Flux)/(Fake Flux)) 
#for fake mag of (<19 mag),(19<fake<21mag),(21<fake<23mag),(23<fake).  
#Make a note of the number of 5sigma outliers for each histogram.

#print mcmc_fluxes.shape
#print diffim_fluxes.shape
#print mcmc_fluxerrs.shape
'''
fig = plt.figure(i+20)
plt.scatter(mjds,(mcmc_fluxes-diffim_fluxes)/mcmc_fluxerrs,alpha=.15,color='black')
ax,ay,aystd = bindata(mjds,(mcmc_fluxes-diffim_fluxes)/mcmc_fluxerrs,np.arange(min(mjds),max(mjds),1))
plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
plt.xlabel('MJD')
plt.ylabel('(SMP Flux - Diffim Flux)/SMP Flux Error')
#plt.xlim(-.,1.)
plt.ylim(-3.,3.)
#plt.xlim(19,24)
plt.plot([min(mjds),max(mjds)],[0,0],color='black')
plt.title(filt+' Band')
plt.savefig(out+'residstd_vs_mjd_'+filt+'.pdf')
print 'here13'
'''
'''
fig = plt.figure(i+21)
ww = fake_fluxes_plus > 20.
plt.scatter(psfs[ww],(mcmc_fluxes[ww]-fake_fluxes_plus[ww])/mcmc_fluxerrs[ww],alpha=.15,color='red',label='SMP')
ax,ay,aystd = bindata(psfs[ww],(mcmc_fluxes[ww]-fake_fluxes_plus[ww])/mcmc_fluxerrs[ww],np.arange(min(psfs),max(psfs),.1))
plt.errorbar(ax,ay,aystd,color='red',fmt='o')
plt.scatter(psfs[ww],(diffim_fluxes[ww]-fake_fluxes[ww])/diffim_fluxerrs[ww],alpha=.15,color='blue',label='Diffim')
ax,ay,aystd = bindata(psfs[ww],(diffim_fluxes[ww]-fake_fluxes[ww])/diffim_fluxerrs[ww],np.arange(min(psfs),max(psfs),.1))
plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
plt.xlabel('PSF')
plt.ylabel('(Fit Flux - Fake Flux Adjusted)/Fit Flux Error')
#plt.xlim(-.,1.)
plt.ylim(-20.,20.)
#plt.xlim(19,24)
plt.plot([min(psfs),max(psfs)],[0,0],color='black')
plt.title(filt+' Band')
plt.legend()
plt.savefig(out+'ResidStd_vs_psffwhm_'+filt+'.pdf')
print 'here14'
'''
'''fig = plt.figure(i+21)
plt.scatter(psfs,(mcmc_fluxes-fake_fluxes_plus)/mcmc_fluxerrs,alpha=.15,color='red',label='SMP')
ax,ay,aystd = bindata(psfs,(mcmc_fluxes-fake_fluxes_plus)/mcmc_fluxerrs,np.arange(min(psfs),max(psfs),.1))
plt.errorbar(ax,ay,aystd,color='red',fmt='o')
plt.scatter(psfs,(diffim_fluxes-fake_fluxes)/diffim_fluxerrs,alpha=.15,color='blue',label='Diffim')
ax,ay,aystd = bindata(psfs,(diffim_fluxes-fake_fluxes)/diffim_fluxerrs,np.arange(min(psfs),max(psfs),.1))
plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
plt.xlabel('PSF')
plt.ylabel('(Fit Flux - Fake Flux Adjusted)/Fit Flux Error')
plt.xlim(-.,1.)                                                                                                                                                                                         
plt.ylim(-20.,20.)
plt.plot([min(psfs),max(psfs)],[0,0],color='black')
plt.title(filt+' Band')
plt.legend()
plt.savefig(out+'ResidStd_vs_psffwhm_'+filt+'.pdf')
'''
'''
fig = plt.figure(i+22)
plt.scatter(fakemags,(mcmc_fluxes-fake_fluxes),alpha=.15,color='red',label='SMP')
print 'herea'
ax,ay,aystd = bindata(fakemags,(mcmc_fluxes-fake_fluxes),np.arange(min(fake_fluxes),max(fake_fluxes),.25))
print 'hereb'
plt.errorbar(ax,ay,aystd,color='red',fmt='o')
print 'herec'
plt.scatter(fakemags,(diffim_fluxes-fake_fluxes),alpha=.15,color='blue',label='Diffim')
print 'hered'
ax,ay,aystd = bindata(fakemags,(diffim_fluxes-fake_fluxes),np.arange(min(fake_fluxes),max(fake_fluxes),.25))
print 'heree'
plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
plt.xlabel('Fake Mag')
plt.ylim(-200,200)
plt.legend()
plt.xlim(20,30)
plt.plot([min(fakemags),max(fakemags)],[0,0],color='black')
plt.ylabel('Fit Flux - Fake Flux')
#plt.plot([min(psfs),max(psfs)],[0,0],color='black')
plt.title(filt+' Band')
plt.savefig(out+'FluxResid_vs_true_'+filt+'.png')
print 'here15'
'''
'''
fig = plt.figure(i+26)
plt.scatter(fakemags,(mcmc_fluxes-fake_fluxes_plus),alpha=.15,color='red',label='SMP')
ax,ay,aystd = bindata(fakemags,(mcmc_fluxes-fake_fluxes_plus),np.arange(min(fake_fluxes_plus),max(fake_fluxes_plus),.25))
plt.errorbar(ax,ay,aystd,color='red',fmt='o')
plt.scatter(fakemags,(diffim_fluxes-fake_fluxes_plus),alpha=.15,color='blue',label='Diffim')
ax,ay,aystd = bindata(fakemags,(diffim_fluxes-fake_fluxes_plus),np.arange(min(fake_fluxes_plus),max(fake_fluxes_plus),.25))
plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
plt.xlabel('Fake Mag')
#plt.ylim(-200,200)
plt.legend()
plt.xlim(20,30)
plt.plot([min(fakemags),max(fakemags)],[0,0],color='black')
plt.ylabel('Fit Flux - Fake Flux Adjusted (plus)')
plt.title(filt+' Band')
plt.savefig(out+'FluxResid_vs_true_adjustedplus_'+filt+'.png')
print 'here16'
'''


'''
fig = plt.figure(203)
plt.scatter(fakemags,(mcmc_fluxes-fake_fluxes_plus),alpha=.15,color='red',label='SMP')
ax,ay,aystd = bindata(fakemags,(mcmc_fluxes-fake_fluxes_plus),np.arange(min(fake_fluxes_plus),max(fake_fluxes_plus),.25))
plt.errorbar(ax,ay,aystd,color='red',fmt='o')
plt.scatter(fakemags,(diffim_fluxes-fake_fluxes),alpha=.15,color='blue',label='Diffim')
ax,ay,aystd = bindata(fakemags,(diffim_fluxes-fake_fluxes),np.arange(min(fake_fluxes_plus),max(fake_fluxes_plus),.25))
plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
plt.xlabel('Fake Mag')
plt.ylim(-200,200)
plt.legend()
plt.xlim(20,30)
plt.plot([min(fakemags),max(fakemags)],[0,0],color='black')
plt.ylabel('Fit Flux - Fake Flux Adjusted')
plt.title(filt+' Band')
plt.savefig(out+'FluxResid_vs_true_adjusted_'+filt+'.png')
print 'here17'
'''
'''
fig = plt.figure(204)
plt.scatter(fakemags,(mcmc_mags-fakemags_plus),alpha=.15,color='red',label='SMP')
ax,ay,aystd = bindata(fakemags,(mcmc_mags-fakemags_plus),np.arange(min(fake_fluxes_plus),max(fake_fluxes_plus),.25))
plt.errorbar(ax,ay,aystd,color='red',fmt='o')
plt.scatter(fakemags,(diffim_mags-fakemags),alpha=.15,color='blue',label='Diffim')
ax,ay,aystd = bindata(fakemags,(diffim_mags-fakemags),np.arange(min(fake_fluxes_plus),max(fake_fluxes_plus),.25))
plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
plt.xlabel('Fake Mag')
plt.ylim(-1,1)
plt.legend()
plt.xlim(20,30)
plt.plot([min(fakemags),max(fakemags)],[0,0],color='black')
plt.ylabel('Fit Mag - Fake Mag Adjusted')
plt.title(filt+' Band')
plt.savefig(out+'MagResid_vs_true_adjusted_'+filt+'.png')
print 'here18'

fakeflux = 10**(.4*(31.-fakemags))
mcmcflux = 10**(.4*(31.-mcmc_mags))
diffimflux = 10**(.4*(31.-diffim_mags))
'''
'''
fig = plt.figure(206)
plt.scatter(fakemags,(mcmcflux-fakeflux),alpha=.15,color='red',label='SMP')
ax,ay,aystd = bindata(fakemags,(mcmcflux-fakeflux),np.arange(min(fakemags),max(fakemags),.25))
plt.errorbar(ax,ay,aystd,color='red',fmt='o')
plt.scatter(fakemags,(diffimflux-fakeflux),alpha=.15,color='blue',label='Diffim')
ax,ay,aystd = bindata(fakemags,(diffimflux-fakeflux),np.arange(min(fakemags),max(fakemags),.25))
plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
plt.xlabel('Fake Mag')
plt.ylim(-1000,1000)
#plt.xlim(0,5000)
plt.legend()
#plt.xlim(20,30)
#plt.plot([0,20000],[0,0],color='black')
plt.plot([min(fakemags),max(fakemags)],[0,0],color='black')
plt.ylabel('Fit Flux - Fake Flux')
plt.title(filt+' Band')
plt.savefig(out+'INPUTFluxResid_vs_true_'+filt+'.png')

print 'here18'
print out+'INPUTFluxResid_vs_true_'+filt+'.png'
'''

#print chis[abs(chis)>5000]
#print sns[abs(chis)>5000]
#print 'eerererere'
#raw_input()

# fig = plt.figure(372)
# plt.scatter(chis[chis<20000],mcmc_fluxes[chis<20000]-fake_fluxes[chis<20000],color='black',alpha=.4)
# plt.plot([min(chis),max(chis)],[0,0],color='black')
# plt.xlim(-1000,1000)
# plt.ylim(-1000,1000)
# plt.xlabel('Chi')
# plt.ylabel('Flux Residual')
# plt.savefig(out+'FluxResid_vs_Chi_'+filt+'.png')

# fig = plt.figure(373)
# plt.scatter(chisqs,chis,color='black',alpha=.4)
# plt.plot([min(chisqs),max(chisqs)],[0,0],color='black')
# plt.xlabel('Chisq')
# plt.ylabel('Chi')
# plt.savefig(out+'Chi_vs_Chisq_'+filt+'.png')

# fig = plt.figure(374)
# plt.scatter(chisqs,mcmc_fluxes-fake_fluxes,color='black',alpha=.4)
# plt.plot([min(chisqs),max(chisqs)],[0,0],color='black')
# plt.xlabel('Chisq')
# plt.ylabel('Flux Residual')
# plt.savefig(out+'FluxResid_vs_Chisq_'+filt+'.png')

# fig = plt.figure(375)
# plt.hist([diffim_fluxes[mjdflags==0]-fake_fluxes[mjdflags==0],mcmc_fluxes[mjdflags==0]-fake_fluxes[mjdflags==0]],bins=range(-4000,4000,100),label=['Diffim','MCMC'])
# plt.savefig(out+'DiffimFluxResid_hist_'+filt+'.png')

# fig = plt.figure(376)
# #plt.hist([fakechis[mjdflags==0],chis[mjdflags==0]],bins=range(-2000,2000,200),normed=True,label=['Fake Inserted Total = '+str(round(sum(fakechis[abs(fakechis)<2000]),1)),'SMP, Total = '+str(round(sum(chis[abs(chis)<2000]),1))])
# plt.hist(chis[mjdflags==0],bins=range(-2000,2000,200),normed=True,label='SMP, Total = '+str(round(sum(chis[abs(chis)<2000]),1)))
# plt.xlabel('Sum(Data-Sim)')
# plt.legend(fontsize = 'small')
# plt.savefig(out+'Chi_hist_'+filt+'.png')

# fig = plt.figure(377)
# #plt.hist([fakechisqs[(fakechisqs<20)&(mjdflags==0)],chisqs[(chisqs<20)&(mjdflags==0)],chisqs[mjdflags==1]],bins=np.arange(0,4,.2),normed=True,label=['Fake Inserted, Total = '+str(round(sum(fakechisqs[fakechisqs<20]),1)),'SMP, Total = '+str(round(sum(chisqs[chisqs<20]),1)),'SMP No Supernova'])
# plt.hist([chisqs[(chisqs<20)&(mjdflags==0)],chisqs[mjdflags==1]],bins=np.arange(0,4,.2),normed=True,label=['Total = '+str(round(sum(fakechisqs[fakechisqs<20]),1))+'SMP', 'Total = '+str(round(sum(chisqs[chisqs<20]),1))+'SMP No Supernova'])
# plt.xlabel('Chisq')
# plt.legend(fontsize = 'small')
# plt.savefig(out+'Chisq_hist_'+filt+'.png')
'''
for i in mjdflags:
	print i

print chis.shape
print mjdflags.shape
print chis[mjdflags==1.].shape
print chisqs[mjdflags==1.].shape
'''

# fig = plt.figure(378)
# #plt.hist(chis[mjdflags==1],bins=range(-2000,2000,200))
# plt.xlabel('Sum(Data-Sim)')
# plt.legend(fontsize = 'small')
# plt.savefig(out+'Chi_hist_noflux'+filt+'.png')

# fig = plt.figure(379)
# #plt.hist(chisqs[mjdflags==1],bins=20)
# plt.xlabel('Chisq')
# plt.legend(fontsize = 'small')
# plt.savefig(out+'Chisq_hist_noflux'+filt+'.png')


#fakeflux = 10.**(.4*(31.-fakemags))
#mcmc_fluxa = 10.**(.4*(31.-mcmc_mags))
#diffim_flux = 10.**(.4*(31.-diffim_mags))

# fig = plt.figure(205)

# plt.scatter(fakemags,(mcmc_fixfluxes-fake_fluxes),alpha=.15,color='green',label='SMP Fix Gal')
# ax,ay,aystd = bindata(fakemags,(mcmc_fixfluxes-fake_fluxes),np.arange(min(fake_fluxes),max(fake_fluxes),.25))
# plt.errorbar(ax,ay,aystd,color='green',fmt='o')

# plt.scatter(fakemags,(mcmc_floatfluxes-fake_fluxes),alpha=.15,color='red',label='SMP Float Gal')
# ax,ay,aystd = bindata(fakemags,(mcmc_floatfluxes-fake_fluxes),np.arange(min(fake_fluxes),max(fake_fluxes),.25))
# plt.errorbar(ax,ay,aystd,color='red',fmt='o')

# plt.scatter(fakemags,(diffim_fluxes-fake_fluxes),alpha=.15,color='blue',label='Diffim')
# ax,ay,aystd = bindata(fakemags,(diffim_fluxes-fake_fluxes),np.arange(min(fake_fluxes),max(fake_fluxes),.25))
# plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
# plt.xlabel('Fake Mag')
# plt.ylim(-1000,1000)
# plt.plot([min(fakemags),max(fakemags)],[0,0],color='black',label='Total Sim SNe '+str(totscounter))
# plt.legend()
# plt.xlim(20,28)
# plt.ylabel('Fit Flux - Fake Flux')
# plt.title(filt+' Band')
# plt.savefig(out+'FluxResid_vs_TrueMag_'+filt+'.png')

'''
fig = plt.figure(209)
plt.scatter(fake_fluxes,(mcmc_fluxes-fake_fluxes),alpha=.15,color='green',label='SMP')
ax,ay,aystd = bindata(fake_fluxes,(mcmc_fluxes-fake_fluxes),np.arange(min(fake_fluxes)-250,max(fake_fluxes),250))
plt.errorbar(ax,ay,aystd,color='green',fmt='o')
plt.scatter(fake_fluxes,(diffim_fluxes-fake_fluxes),alpha=.15,color='blue',label='Diffim')
ax,ay,aystd = bindata(fake_fluxes,(diffim_fluxes-fake_fluxes),np.arange(min(fake_fluxes)-250,max(fake_fluxes),250))
plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
plt.xlabel('Fake Flux')
plt.ylim(-1000,1000)
plt.legend()
plt.xlim(-250,5000)
plt.plot([min(fake_fluxes),max(fake_fluxes)],[0,0],color='black')
plt.ylabel('Fit Flux - Fake Flux')
plt.title(filt+' Band')
plt.savefig(out+'FluxResid_vs_TrueFlux_'+filt+'.png')
'''

fig = plt.figure(210,figsize=(10,10))

# plt.scatter(fakemags_minus,(mcmc_fixmags-fakemags_minus),alpha=.07,color='green')
# ax,ay,aystd = bindata(fakemags_minus,(mcmc_fixmags-fakemags_minus),np.arange(min(fakemags),max(fakemags),.25))
# plt.errorbar(ax,ay,aystd,color='green',fmt='o',label='SMP Fix Gal')
'''
plt.scatter(fakemags_minus,(mcmc_fixmags-fakemags_minus),alpha=.1,color='red')
ax,ay,aystd = bindataweighted(fakemags_minus,(mcmc_fixmags-fakemags_minus),np.arange(min(fakemags_minus),max(fakemags_minus),.25))
plt.errorbar(ax,ay,aystd,markersize=10,color='red',fmt='o',label='SMP Brute Force')
'''
a = 1./mcmc_floatmag_errs
print 'weights'
ww = np.isfinite(a) 

if dov3:
	a = 1./domagerrs
	print 'weights'
	jj = np.isfinite(a)
	a = 1./nmagerrs
	kk = np.isfinite(a)
	
#raw_input()
plt.scatter(fakemags_minus[ww],(mcmc_floatmags[ww]-fakemags_minus[ww]),alpha=.1,color='green',label='DES Astrometry')
ax,ay,aystd = bindata(fakemags_minus[ww],(mcmc_floatmags[ww]-fakemags_minus[ww]),np.arange(min(fakemags_minus),max(fakemags_minus),.25))
plt.errorbar(ax,ay,aystd,markersize=10,color='green',fmt='o')

#plt.scatter(fakemags_minus,(mcmc_gfloatmags-fakemags_minus),alpha=.1,color='blue')
#ax,ay,aystd = bindata(fakemags_minus,(mcmc_gfloatmags-fakemags_minus),np.arange(min(fakemags_minus),max(fakemags_minus),.25))
#plt.errorbar(ax,ay,aystd,markersize=10,color='blue',fmt='o',label='SMP Float SN + .27 Resampling')
if dopix:
	plt.scatter(fakemags,(mcmc_gpixfloatmags-fakemags_plus),alpha=.1,color='blue')
	ax,ay,aystd = bindata(fakemags,(mcmc_gpixfloatmags-fakemags_plus),np.arange(min(fakemags_minus),max(fakemags_minus),.25))
	plt.errorbar(ax,ay,aystd,markersize=10,color='blue',fmt='o',label='SMP .5" Galsim Fit RA+DEC')

if dov3:
	plt.scatter(fakemags_minus[jj],(domags[jj]-fakemags_minus[jj]),alpha=.1,color='blue')
        ax,ay,aystd = bindata(fakemags_minus[jj],(domags[jj]-fakemags_minus[jj]),np.arange(min(fakemags_minus),max(fakemags_minus),.25))
        plt.errorbar(ax,ay,aystd,markersize=10,color='blue',fmt='o',label='Float Global Stars and Global SN')
        #plt.scatter(fakemags_minus[kk],(nmags[kk]-fakemags_minus[kk]),alpha=.1,color='blue')
        #ax,ay,aystd = bindata(fakemags_minus[kk],(nmags[kk]-fakemags_minus[kk]),np.arange(min(fakemags_minus),max(fakemags_minus),.25))
        #plt.errorbar(ax,ay,aystd,markersize=10,color='red',fmt='o',label='DES Astrometry')
#diffdata = np.load('diffimresults.npz')
#fakemags = diffdata['fakemags']
#diffim_mags = diffdata['diffimmags']
#plt.scatter(fakemags,(diffim_mags-fakemags),alpha=.1,color='blue')
#ax,ay,aystd = bindata(fakemags,(diffim_mags-fakemags),np.arange(min(fakemags),max(fakemags),.25))
#plt.errorbar(ax,ay,aystd,color='blue',fmt='o',label='Diffim')

plt.xlabel('Fake Mag')
plt.ylim(-.05,.05)
plt.plot([min(fakemags),max(fakemags)],[0,0],color='black',label='Total Sim SNe '+str(totscounter))
plt.legend()
plt.xlim(20.,25)
plt.ylabel('Fit Mag - Fake Mag')
plt.title(filt+' Band')
plt.savefig(out+'MagResid_vs_TrueMag_withoffset'+filt+'.png')
print out+'MagResid_vs_TrueMag_withoffset'+filt+'.png'





plt.clf()
fig = plt.figure(figsize=(10,7))
plt.plot([20,29],[1,1],linestyle='--',color='black')
plt.xlabel('sb mag (mag/asec^2)')
plt.ylabel('RMS(fitflux-fakeflux / fluxerr)')
diffimbins = []
smpbins = []
floatbins = []
brightsmpbins = []
for i in np.arange(20.5,29.5,.1):
	ww = [(hostgal_bandmags < i+.25) & (hostgal_bandmags > i - .25) & (fakemags <22.5)]
        xx = [(hostgal_bandmags < i+.25) & (hostgal_bandmags > i - .25) & (fakemags >22.5)]

	diffimstd = (diffim_fluxes[ww] - fake_fluxes[ww]) / diffim_fluxerrs[ww]
	smpstd = (mcmc_gpixfloatfluxes[ww] - fake_fluxes[ww]) / mcmc_gpixfloatfluxerrs[ww]
	diffimstd = diffimstd[~np.isinf(diffimstd)]
	smpstd = smpstd[~np.isinf(smpstd)]
	print 'bright',i,len(diffimstd),len(smpstd)
	diffim_rmspt = np.sqrt(np.nanmean(diffimstd**2))
	smp_rmspt = np.sqrt(np.nanmean(smpstd**2))

	diffimbins.append(diffim_rmspt)
	brightsmpbins.append(smp_rmspt)
	#if dopix:
	#	smpstd = (mcmc_gpixfloatfluxes[xx] - fake_fluxes[xx]) / mcmc_gpixfloatfluxerrs[xx]
        floatstd = (mcmc_floatfluxes[xx] - fake_fluxes[xx]) / mcmc_floatfluxerrs[xx]
	#if dopix:
	#	smpstd = smpstd[~np.isinf(smpstd)]
	#	smp_rmspt = np.sqrt(np.nanmean(smpstd**2))
	#	smpbins.append(smp_rmspt)

        floatstd = floatstd[~np.isinf(floatstd)]
        float_rmspt = np.sqrt(np.nanmean(floatstd**2))
	floatbins.append(float_rmspt)

	#print 'dim',i,len(diffimstd),len(smpstd)

#bbb = np.convolve(brightsmpbins, np.ones((4,))/4, mode='same')
#if dopix:
#	sss = np.convolve(smpbins, np.ones((6,))/6, mode='same')
fff = np.convolve(floatbins,np.ones((6,))/6, mode='same')
plt.plot(np.arange(21.5,29.5,.1),fff[10:],color='blue',linestyle='--',label='Fake Mag > 22.5')
#if dopix:
#	plt.plot(np.arange(21.5,29.5,.1),sss[10:],color='red',linestyle='--',label='Resampling Fake Mag > 22.5')
plt.ylim(0,8)
plt.xlim(20,29)
plt.legend()
plt.savefig(out+'_sbbins.png')
print out+'_sbbins.png'
raw_input()











plt.clf()
fig = plt.figure(211,figsize=(10,10))
plt.scatter(mcmc_fixmags,mcmc_floatmags,alpha=.2,label='SMP Float SN',color='green')
#plt.scatter(mcmc_fixmags,mcmc_gfloatmags,alpha=.2,label='SMP Float SN + .27 Resampling',color='blue')
#plt.scatter(mcmc_fixmags,mcmc_gpixfloatmags,alpha=.2,label='SMP Float SN + .5 Resampling',color='gold')
plt.plot([19,25],[19,25],color='black')

plt.xlabel('Brute Force SN Mag')
plt.ylabel('Float SN Mag')
plt.legend()
plt.ylim(19,25)
plt.xlim(19,25)
plt.title(filt+' Band')
plt.savefig(out+'Float_vs_FixMag_'+filt+'.png')
print out+'Float_vs_FixMag_'+filt+'.png'

print mcmc_floatmags-diffim_mags
plt.clf()
fig = plt.figure(234,figsize=(10,10))
aaa = mcmc_floatmags-diffim_mags
aaa = aaa[~np.isinf(aaa)]
aaa = aaa[~np.isnan(aaa)]
plt.hist(aaa,bins=np.arange(-1,1,.05))
plt.title(filt+' Band')
plt.savefig(out+'Diffim-Float_hist_'+filt+'.png')
print out+'Diffim-Float_hist_'+filt+'.png'
plt.clf()
fig = plt.figure(232,figsize=(10,10))
plt.scatter(fakemags,mcmc_floatmags-diffim_mags,alpha=.2,label='SMP Float SN',color='green')
ax,ay,aystd = bindata(fakemags,mcmc_floatmags-diffim_mags,np.arange(min(fakemags),max(fakemags),.25))
plt.errorbar(ax,ay,aystd,markersize=10,color='green',fmt='o')
#plt.scatter(mcmc_fixmags,mcmc_gfloatmags,alpha=.2,label='SMP Float SN + .27 Resampling',color='blue')
#plt.scatter(mcmc_fixmags,mcmc_gpixfloatmags,alpha=.2,label='SMP Float SN + .5 Resampling',color='gold')
plt.plot([19,25],[0,0],color='black')

plt.xlabel('Fake Mag Mag')
plt.ylabel('Float SN Mag - Diffim Mag')
plt.legend()
plt.ylim(-.2,.2)
plt.xlim(19,25)
plt.title(filt+' Band')
plt.savefig(out+'Diffim-Float_vs_FakeMag_'+filt+'.png')
print out+'Diffim-Float_vs_FakeMag_'+filt+'.png'


# fig = plt.figure(210)
# plt.scatter(fake_fluxes,(mcmc_fluxes-fake_fluxes),alpha=.15,color='green',label='SMP')
# ax,ay,aystd = bindata(fake_fluxes,(mcmc_fluxes-fake_fluxes),np.arange(min(fake_fluxes),max(fake_fluxes),250))
# plt.errorbar(ax,ay,aystd,color='green',fmt='o')
# plt.scatter(fake_fluxes,(diffim_fluxes-fake_fluxes),alpha=.15,color='blue',label='Diffim')
# ax,ay,aystd = bindata(fake_fluxes,(diffim_fluxes-fake_fluxes),np.arange(min(fake_fluxes),max(fake_fluxes),250))
# plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
# plt.xlabel('Fake Flux')
# plt.ylim(-1000,1000)
# plt.legend()
# plt.xlim(0,5000)
# plt.plot([min(fake_fluxes),max(fake_fluxes)],[0,0],color='black')
# plt.ylabel('Fit Flux - Fake Flux')
# plt.title(filt+' Band')
# plt.savefig(out+'FluxResid_vs_TrueFlux_'+filt+'.png')

print 'here19'
'''
fig = plt.figure(i+23)
plt.scatter(fakemags,(mcmc_fluxes-fake_fluxes)/mcmc_fluxerrs,alpha=.15,color='red',label='SMP')
ax,ay,aystd = bindata(fakemags,(mcmc_fluxes-fake_fluxes)/mcmc_fluxerrs,np.arange(min(fake_fluxes),max(fake_fluxes),.250))
plt.errorbar(ax,ay,aystd,color='red',fmt='o')
plt.scatter(fakemags,(diffim_fluxes-fake_fluxes)/diffim_fluxerrs,alpha=.15,color='blue',label='Diffim')
ax,ay,aystd = bindata(fakemags,(diffim_fluxes-fake_fluxes)/diffim_fluxerrs,np.arange(min(fake_fluxes),max(fake_fluxes),.250))
plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
plt.xlabel('Fake Mag')
plt.xlim(19,24)
plt.legend()
plt.plot([min(fakemags),max(fakemags)],[0,0],color='black')
plt.ylim(-20,20)
plt.ylabel('(Fit Flux - Fake Flux)/(Fit Flux Err)')
#plt.plot([min(psfs),max(psfs)],[0,0],color='black')
plt.title(filt+' Band')
plt.savefig(out+'FluxResidStd_vs_true_'+filt+'.png')
print 'here20'
'''
'''
fig = plt.figure(i+24)
plt.scatter(zpt_offsets,(fake_fluxes_plus - fake_fluxes_minus),alpha=.15,color='red',label='SMP')
#ax,ay,aystd = bindata(fake_fluxes,(mcmc_fluxes-fake_fluxes)/mcmc_fluxerrs,np.arange(min(fake_fluxes),max(fake_fluxes),50))
#plt.errorbar(ax,ay,aystd,color='red',fmt='o')
plt.xlabel('Zpt Offsets')
#plt.legend()
plt.plot([min(zpt_offsets),max(zpt_offsets)],[0,0])
plt.ylim(-20,20)
plt.xlim(min(zpt_offsets),max(zpt_offsets))
plt.ylabel('(Fake Flux Plus - Fake Flux Minus)')
#plt.plot([min(psfs),max(psfs)],[0,0],color='black')
plt.title(filt+' Band')
plt.savefig(out+'CheckZptOffset_vs_true_'+filt+'.png')

print 'here21'
'''
'''
fig = plt.figure(i+25)
plt.scatter(fake_fluxes,mcmc_fluxes,alpha=.15,color='red',label='SMP')
#ax,ay,aystd = bindata(fakemags,(mcmc_fluxes-fake_fluxes),np.arange(min(fake_fluxes),max(fake_fluxes),.25))
#plt.errorbar(ax,ay,aystd,color='red',fmt='o')
plt.scatter(fake_fluxes,diffim_fluxes,alpha=.15,color='blue',label='Diffim')
#ax,ay,aystd = bindata(fakemags,(diffim_fluxes-fake_fluxes),np.arange(min(fake_fluxes),max(fake_fluxes),.25))
#plt.errorbar(ax,ay,aystd,color='blue',fmt='o')
plt.xlabel('Fake Flux')
plt.ylim(-100,1500)
plt.legend()
plt.xlim(-100,1500)
plt.plot([-100,1500],[-100,1500],color='black')
plt.ylabel('Fit Flux')
#plt.plot([min(psfs),max(psfs)],[0,0],color='black')                                                                                                                                                        
plt.title(filt+' Band')
plt.savefig(out+'Fit_vs_true_'+filt+'.png')
'''
print 'Done'
sys.exit()

a = (diffim_mags-fakemags)
b = (mcmc_mags-fakemags)

bin_edges = [0.,19.,21.,23.,99.]
bins = np.arange(-.5,.5,.05)
f, axt = plt.subplots(2,2)
axs = axt.ravel()

if len(a[fakemags < 19.]) > 0:
	aa = a[fakemags < 19.]
	anum = len(aa[aa-np.mean(aa)>5.*np.std(aa)])
	atot = len(aa)
	ap = float(anum)/float(atot)
	bb = b[fakemags < 19.]
	bnum = len(bb[bb-np.mean(bb)>5.*np.std(bb)])
	btot = len(bb)
	bp = float(bnum)/float(btot)
	axs[0].hist([a[fakemags < 19.],b[fakemags < 19.]],bins=bins,color=['blue','red'],label=['Diffim #Outliers '+str(anum)+
				', '+str(round(ap*100.,3))+'$\%$'
				,'SMP #Outliers '+str(bnum)+', '+str(round(bp*100.,3))+'$\%$'])
#else:
#	axs[0].hist([[0],[0]],color=['blue','green'],label=['Diffim','SMP'])
axs[0].set_title('Fake Mag < 19')
axs[0].set_xlabel('(Fit Mag - Fake Mag)',size=10)
axs[0].legend(fontsize = 'xx-small')


aa = a[(fakemags < 21.) & (fakemags > 19.)]
anum = len(aa[aa-np.mean(aa)>5.*np.std(aa)])
atot = len(aa)
ap = float(anum)/float(atot)
bb = b[(fakemags < 21.) & (fakemags > 19.)]
bnum = len(bb[bb-np.mean(bb)>5.*np.std(bb)])
btot = len(bb)
bp = float(bnum)/float(btot)

axs[1].hist([a[(fakemags < 21.) & (fakemags > 19.)],b[(fakemags < 21.) & (fakemags > 19.)]],bins=bins,color=['blue','red'],label=['Diffim #Outliers '+str(anum)+
				', '+str(round(ap*100.,3))+'$\%$'
				,'SMP #Outliers '+str(bnum)+', '+str(round(bp*100.,3))+'$\%$'])
axs[1].set_title('19 < Fake Mag < 21')
axs[1].set_xlabel('(Fit Mag - Fake Mag)',size=10)
axs[1].legend(fontsize = 'xx-small')


aa = a[(fakemags < 23.) & (fakemags > 21.)]
anum = len(aa[aa-np.mean(aa)>5.*np.std(aa)])
atot = len(aa)
ap = float(anum)/float(atot)
bb = b[(fakemags < 23.) & (fakemags > 21.)]
bnum = len(bb[bb-np.mean(bb)>5.*np.std(bb)])
btot = len(bb)
bp = float(bnum)/float(btot)


axs[2].hist([a[(fakemags < 23.) & (fakemags > 21.)],b[(fakemags < 23.) & (fakemags > 21.)]],bins=bins,color=['blue','red'],label=['Diffim #Outliers '+str(anum)+
				', '+str(round(ap*100.,3))+'$\%$'
				,'SMP #Outliers '+str(bnum)+', '+str(round(bp*100.,3))+'$\%$'])
axs[2].set_title('21 < Fake Mag < 23')
axs[2].set_xlabel('(Fit Mag - Fake Mag)',size=10)
axs[2].legend(fontsize = 'xx-small')


print 'here23'

if len(a[fakemags > 23.]) > 0:
	aa = a[fakemags > 23.]
	anum = len(aa[aa-np.mean(aa)>5.*np.std(aa)])
	atot = len(aa)
	ap = float(anum)/float(atot)
	bb = b[fakemags > 23.]
	bnum = len(bb[bb-np.mean(bb)>5.*np.std(bb)])
	btot = len(bb)
	bp = float(bnum)/float(btot)
	axs[3].hist([a[fakemags > 23.],b[fakemags > 23.]],bins=bins,color=['blue','red'],label=['Diffim #Outliers '+str(anum)+
				', '+str(round(ap*100.,3))+'$\%$'
				,'SMP #Outliers '+str(bnum)+', '+str(round(bp*100.,3))+'$\%$'])
axs[3].set_title('Fake Mag > 23')
axs[3].set_xlabel('(Fit Mag - Fake Mag)',size=10)

for ax in axs:
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(8) 
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(10) 

axs[3].legend(fontsize = 'xx-small')
plt.subplots_adjust(bottom=.1, left=.1, right=.9, top=.90, hspace=.4)
plt.savefig(out+'fakemag_4hist_'+filt+'.pdf')


print 'here24'



a = (mcmc_mags-diffim_mags)

bin_edges = [0.,19.,21.,23.,99.]

f, axt = plt.subplots(2,2)
axs = axt.ravel()

if len(a[fakemags < 19.]) > 0:
	aa = a[fakemags < 19.]
	anum = len(aa[aa-np.mean(aa)>5.*np.std(aa)])
	atot = len(aa)
	ap = float(anum)/float(atot)
	axs[0].hist(a[fakemags < 19.],bins=bins,color='red',label='#Outliers '+str(anum)+
				', '+str(round(ap*100.,3))+'$\%$')
#else:
#	axs[0].hist([[0],[0]],color=['blue','green'],label=['Diffim','SMP'])
axs[0].set_title('Fake Mag < 19')
axs[0].set_xlabel('(SMP Mag - Diffim Mag)',size=10)
axs[0].legend(fontsize = 'xx-small')


aa = a[(fakemags < 21.) & (fakemags > 19.)]
anum = len(aa[aa-np.mean(aa)>5.*np.std(aa)])
atot = len(aa)
ap = float(anum)/float(atot)

axs[1].hist(a[(fakemags < 21.) & (fakemags > 19.)],bins=bins,color='red',label='#Outliers '+str(anum)+
				', '+str(round(ap*100.,3))+'$\%$')
axs[1].set_title('19 < Fake Mag < 21')
axs[1].set_xlabel('(SMP Mag - Diffim Mag)',size=10)
axs[1].legend(fontsize = 'xx-small')


aa = a[(fakemags < 23.) & (fakemags > 21.)]
anum = len(aa[aa-np.mean(aa)>5.*np.std(aa)])
atot = len(aa)
ap = float(anum)/float(atot)


axs[2].hist(a[(fakemags < 23.) & (fakemags > 21.)],bins=bins,color='red',label='#Outliers '+str(anum)+
				', '+str(round(ap*100.,3))+'$\%$')
axs[2].set_title('21 < Fake Mag < 23')
axs[2].set_xlabel('(SMP Mag - Diffim Mag)',size=10)
axs[2].legend(fontsize = 'xx-small')

print 'here25'


if len(a[fakemags > 23.]) > 0:
	aa = a[fakemags > 23.]
	anum = len(aa[aa-np.mean(aa)>5.*np.std(aa)])
	atot = len(aa)
	ap = float(anum)/float(atot)
	axs[3].hist(a[fakemags > 23.],bins=bins,color='red',label='#Outliers '+str(anum)+
				', '+str(round(ap*100.,3))+'$\%$')
axs[3].set_title('Fake Mag > 23')
axs[3].set_xlabel('(SMP Mag - Diffim Mag)',size=10)

for ax in axs:
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(8) 
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(10) 

axs[3].legend(fontsize = 'xx-small')
plt.subplots_adjust(bottom=.1, left=.1, right=.9, top=.90, hspace=.4)
plt.savefig(out+'smpdiff_4hist_'+filt+'.pdf')







print 'DONE'









