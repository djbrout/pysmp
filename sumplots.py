import os
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(threshold=np.nan)
import sys
from copy import copy
import addcoltoDESlightcurve as lc
from scipy import stats
#import plotlightcurve as lcurve

def bindata(x,y,bins,returnn=False):
        medians = np.zeros(len(bins)-1)
        mads = np.zeros(len(bins)-1)
        nums = np.zeros(len(bins)-1)

        for i in np.arange(len(bins)-1):
                bs = bins[i]
                bf = bins[i+1]
                ww = [(x>bs)&(x<bf)]
                yhere = y[ww]
		yhere = yhere[np.isfinite(yhere)]
                ss = [abs(yhere) < 3.*np.std(yhere)]
                try:
			nums[i] = len(yhere[ss])
                        medians[i] = np.median(yhere[ss])
                        mads[i] = 1.48*np.median(abs(yhere[ss]-medians[i]))*1/np.sqrt(len(yhere[ss]))
                except IndexError:
			print 'excepted'
			nums[i] = 0.
                        medians[i] = np.nan
                        mads[i] = np.nan
        xvals = (bins[1:] + bins[:-1])/2.
	if returnn:
		return xvals,medians,mads,nums
        return xvals,medians,mads



#def lightcurve(mjd,fitmag,fitmagerr,fakemag,fitflux,fitfluxerr,fakeflux,filename):

save = '/global/u1/d/dbrout/PySMP/paperplots/'

#filtsu = ['g','r','i']
filtsu = ['r']
outfolder = '/global/cscratch1/sd/dbrout/v3/smp_y1y2_shallow_v3_64newsky_floatpos_galsim_skyerr'

radec = False
snradec = False
snradecoff = True
dolc = False

fitmag = []
magerr = []
fakemag = []
fitfluxs = []
zps = []
fakefluxs = []
fluxerrs = []
globalraoffsets = []
globaldecoffsets = []
globalmagresid = []
globalfakemag = []
globalsnra = []
globalsndec = []
snnames = []
filtsvec = []
mjdvec = []
hostgal = []
peakmags = []
pfm = []
pfitmag = []
pm = []
so = []
pff = []
#files = os.listdir(os.path.join(outfolder,'np_data/'+filt+'/'))
for filt in filtsu:
    print filt
    if snradec:
	    fff = np.load(outfolder+'/np_data/v4_'+filt+'_withoff.npz')
    else:
	    fff = np.load(outfolder+'/np_data/v4_'+filt+'.npz')
    #print fff.keys()
    #print fff['sns'][0:10]
    #raw_input()
    flux = fff['flux']
    fluxerr = fff['fluxerr']
    fakemags = fff['fakemags']
    #print len(fff['hostgal_bandmags'])
    #print len(flux)
    #print fff['hostgal_bandmags'][0:1000]
    hostmag = fff['hostgal_bandmags']
    #raw_input()
    sns = fff['sns']
    if snradecoff:
	    snraoff = fff['snraoff']
	    sndecoff = fff['sndecoff']
	    print snraoff
    else:
	    snraoff = flux*0.-9
	    sndecoff = flux*0.-9
    for sn in sns:
	    #print sn
	    #print sn[:-2]
	    #raw_input()
	    snnames.append(sn[:-2])
	    filtsvec.append(sn[-1])
    for sn in np.unique(sns):
	    ww = (sns == sn)
	    peakmags.extend(fakemags[ww]*0. + min(fakemags[ww][fakemags[ww] > 15.]))
	    pfm.extend(fakemags[ww])
	    pfitmag.extend(31.-2.5*np.log10(flux[ww]))
	    ff = 10**(.4*(31.-fakemags[ww]))
	    ff[ff<1.] = 1.
	    #pff.extend((flux[ww]-ff)/ff)
	    if snradec:
		    pm.append(min(fakemags[ww][fakemags[ww] > 15.]))
	    if snradec:
		    pff.extend((flux[ww]-ff)/ff)
		    so.extend(np.sqrt(fff['globalsnra'][ww]**2+fff['globalsndec'][ww]**2))
    fitfluxs.extend(flux)
    mjdvec.extend(fff['mjd'])
    fakefluxs.extend(10**(.4*(31. - fff['fakemags'])))
    fluxerrs.extend(fluxerr)
    fm = 31.-2.5*np.log10(flux)
    fme = -2.5*np.log10(flux) + 2.5*np.log10(flux+fluxerr)
    magerr.extend(fme)
    fitmag.extend(fm)
    fakemag.extend(fakemags)
    hostgal.extend(hostmag)
    if radec:
	    graoff = fff['globalraoffsets']
	    gdecoff = fff['globaldecoffsets']
    #print len(graoff)
    #print len(sns)
    #print len(np.unique(sns))
    if snradec:
	    snra = fff['globalsnra']
	    sndec = fff['globalsndec']
	    globalsnra.extend(snra)
	    globalsndec.extend(sndec)
    #raw_input()
    if radec:
	    for i,sn in enumerate(np.unique(sns)):
	        #globalraoffsets.extend(graoff[i])
	        #globaldecoffsets.extend(gdecoff[i])
		    for f,m in zip(fm[sns == sn],fakemags[sns == sn]):
			    globalmagresid.extend(graoff[i]*0. + f-m)
			    globalfakemag.extend(graoff[i]*0. + m)
			    globalraoffsets.extend(graoff[i])
			    globaldecoffsets.extend(gdecoff[i])

globalraoffsets = np.array(globalraoffsets)
globaldecoffsets = np.array(globaldecoffsets)
globalmagresid = np.array(globalmagresid)
globalfakemag = np.array(globalfakemag)
peakmags = np.array(peakmags)
pfm = np.array(pfm)
pff = np.array(pff)
pfitmag = np.array(pfitmag)
pm = np.array(pm)
so = np.array(so)
print pff.shape
print so.shape
#raw_input()
'''
if radec:
	ww = (globalfakemag<22)
	globalraoffsets = globalraoffsets[ww]
	globaldecoffsets = globaldecoffsets[ww]
	globalmagresidh = globalmagresid[ww]
	globalfakemag = globalfakemag[ww]
'''
fitmag = np.array(fitmag)
fakemag = np.array(fakemag)
magerr= np.array(magerr)
fitfluxs = np.array(fitfluxs)
fitfluxso = copy(fitfluxs)
fakemago = copy(fakemag)
fakefluxs = np.array(fakefluxs)
fakefluxs[fakefluxs < 1] = 1. 
fluxerrs = np.array(fluxerrs)
fluxerrso = copy(fluxerrs)
mjdvec = np.array(mjdvec)
hostgal = np.array(hostgal)
filtsvec = np.array(filtsvec,dtype='str')
snnames = np.array(snnames,dtype='str')

globalsnra = np.array(globalsnra)
globalsndec = np.array(globalsndec)
#print globalsnra
#print 'globalsnra'
#raw_input()

globalsnoff = np.sqrt(globalsnra**2+globalsndec**2)

out = outfolder+'/lightcurves/all/'
if not os.path.exists(out):
	os.makedirs(out)
if dolc:
	for sn in np.unique(snnames):
		ww = (snnames == sn)
	        #print sn
		#lcurve.lightcurve(mjdvec[ww],fitmag[ww],magerr[ww],fakemag[ww],fitfluxs[ww],fluxerrs[ww],fakefluxs[ww],filtsvec[ww],out+sn+'.png',title=sn)
	#raw_input()

for sn in np.unique(snnames):
	ww = (snnames == sn)
	if (np.min(fakemag[ww & (fakemag > 15.)]) < 20.):
		print fakemag[ww & (fakemag > 15.) & (fakemag < 21.)]
		print (fitmag-fakemag)[ww & (fakemag > 15.) & (fakemag < 21.)]
		print magerr[ww & (fakemag > 15.) & (fakemag < 21.)]
		print mjdvec[ww & (fakemag > 15.) & (fakemag < 21.)]
		print sn
		#raw_input()
#print np.unique(snnames)
#raw_input()
#fs = ['g','r','i','all']
print filtsu
filtsu.append('all')

fs = filtsu
print fs
for filt in fs:
	plt.clf()
	f = filt
        if filt == 'all':
		ww = (filtsvec != filt) & ((abs(fitfluxs -fakefluxs)/fakefluxs) <.7)
        else:
                ww = (filtsvec == filt) & ((abs(fitfluxs -fakefluxs)/fakefluxs) <.7)
	plt.scatter(peakmags[ww],(fitfluxs[ww]-fakefluxs[ww])/fakefluxs[ww],alpha=.1,color='green',label=str(len(fakemag[ww])))
	ax,ay,aystd = bindata(peakmags[ww],(fitfluxs[ww]-fakefluxs[ww])/fakefluxs[ww],np.arange(min(peakmags),max(peakmags),.5))
	plt.errorbar(ax,ay,aystd,markersize=10,color='green',fmt='o',label='SMP')

	plt.xlabel('Peak Fake Mag')
	plt.ylim(-.08,.08)
	plt.plot([min(fakemag),max(fakemag)],[0,0],color='black')
	plt.legend(fontsize=11)
	plt.xlim(19.8,24.2)
	plt.ylabel('Percentage Flux Difference')
	plt.title(f+' band')
	out = outfolder+'/plots/'
	plt.savefig(out+'Pflux_vs_PeakFakeMag_'+f+'band.png')
	print out+'Pflux_vs_PeakFakeMag_'+f+'band.png'

fitmag[fitfluxs <= 0.] = 99
for filt in fs:
        plt.clf()
        f = filt
        if filt == 'all':
                ww = (filtsvec != filt)  & ((abs(fitfluxs -fakefluxs)/fakefluxs) <.7)
        else:
                ww = (filtsvec == filt)  & ((abs(fitfluxs -fakefluxs)/fakefluxs) <.7)
        plt.scatter(fakemag[ww],fitmag[ww]-fakemag[ww],alpha=.1,color='green',label=str(len(fakemag[ww])))
        ax,ay,aystd = bindata(fakemag[ww],fitmag[ww]-fakemag[ww],np.arange(min(fakemag),max(fakemag),.3))
        plt.errorbar(ax,ay,aystd,markersize=10,color='green',fmt='o',label='SMP')

	plt.xlabel('Fake Mag')
        plt.ylim(-.08,.08)
        plt.plot([min(fakemag),max(fakemag)],[0,0],color='black')
        plt.legend(fontsize=11)
        plt.xlim(19.8,24.2)
        plt.ylabel('Fit - Fake Mag')
        plt.title(f+' band')
        out = outfolder+'/plots/'
        plt.savefig(out+'MagResid_vs_FakeMag_'+f+'band.png')
        print out+'MagResid_vs_FakeMag_'+f+'band.png'



for filt in fs:
        plt.clf()
        f = filt
        if filt == 'all':
                ww = (filtsvec != filt)  & ((abs(fitfluxs -fakefluxs)/fakefluxs) <.7)
        else:
                ww = (filtsvec == filt)  & ((abs(fitfluxs -fakefluxs)/fakefluxs) <.7)
        #plt.scatter(fakemag[ww],fitmag[ww]-fakemag[ww],alpha=.1,color='green',label=str(len(fakemag[ww])))
        #ax,ay,aystd = bindata(fakemag[ww],fitmag[ww]-fakemag[ww],np.arange(min(fakemag),max(fakemag),.3))
        #plt.errorbar(ax,ay,aystd,markersize=10,color='green',fmt='o',label='SMP')
	plt.hist(fluxerrs[ww]/fakefluxs[ww],normed=True,bins=np.arange(0.,2,.05))
	plt.xlabel('Fit Flux Err / Fake Flux')
        #plt.ylim(-.08,.08)
        plt.plot([min(fakemag),max(fakemag)],[0,0],color='black')
        plt.legend(fontsize=11)
        plt.xlim(0,2.)
        plt.ylabel('Count')
        plt.title(f+' band')
        out = outfolder+'/plots/'
        plt.savefig(out+'ErrHist_'+f+'band.png')
        print out+'ErrHist_'+f+'band.png'




plt.clf()
for filt in fs:
	plt.clf()
	f = filt
	if filt == 'all':
		ww = (filtsvec != filt)
	else:
		ww = (filtsvec == filt)
	plt.scatter(fakemag[ww],(fitfluxs[ww]-fakefluxs[ww])/fakefluxs[ww],alpha=.1,color='green')#,label=str(len(fakemag[ww])))
	ax,ay,aystd = bindata(fakemag[ww],(fitfluxs[ww]-fakefluxs[ww])/fakefluxs[ww],np.arange(min(fakemag),max(fakemag),.35))
	plt.errorbar(ax,ay,aystd,markersize=10,color='green',fmt='o',label='SMP Exact Fake Pos')

	plt.xlabel('Fake Mag')
	plt.ylim(-.1,.1)
	plt.plot([min(fakemag),max(fakemag)],[0,0],color='black')
	plt.legend(fontsize=11)
	plt.xlim(19.5,24.5)
	plt.ylabel('Percentage Flux Difference')
	plt.title(f+' band')
	out = outfolder+'/plots/'
	plt.savefig(out+'Pflux_vs_FakeMag_'+f+'band.png')
	plt.savefig(save+'Pflux_vs_FakeMag_'+f+'band.png')

	print out+'Pflux_vs_FakeMag_'+f+'band.png'


	plt.clf()
	plt.scatter(fakemag[ww],(fitfluxs[ww]-fakefluxs[ww])/fluxerrs[ww],alpha=.1,color='green')
	ax,ay,aystd,l = bindata(fakemag[ww],(fitfluxs[ww]-fakefluxs[ww])/fluxerrs[ww],np.arange(min(fakemag),max(fakemag),.35),returnn=True)
	plt.errorbar(ax,ay,aystd*l**.5,markersize=10,color='green',fmt='o',label='SMP Astrometry')

	plt.xlabel('Fake Mag')
	plt.ylim(-2,2)
	plt.plot([min(fakemag),max(fakemag)],[0,0],color='black')
	plt.legend(fontsize=11)
	plt.xlim(19.5,28)
	plt.ylabel('Stdev Fit - Fake')
	plt.title(f+' Band')
	out = outfolder+'/plots/'
	plt.savefig(out+'StdResid_vs_TrueMag_'+f+'band.png')
	plt.savefig(save+'StdResid_vs_FakeMag_'+f+'band.png')
	print out+'StdResid_vs_TrueMag_'+f+'band.png'

'''
plt.clf()
ww = ((abs(fitfluxs -fakefluxs)/fakefluxs) < .7)
globsnoff = (snraoff**2 + sndecoff**2)**.5
plt.scatter(globsnoff[ww],(fitfluxs[ww]-fakefluxs[ww])/fakefluxs[ww],alpha=.1,color='green')
ax,ay,aystd = bindata(globalsnoff[ww],(fitfluxs[ww]-fakefluxs[ww])/fakefluxs[ww],np.arange(-.0011,.2,.01))
plt.errorbar(ax,ay,aystd,markersize=10,color='black',fmt='o',label='SMP Astrometry')

plt.xlabel('Global SN Offset')
plt.ylim(-.08,.08)
plt.xlim(-.01,.2)
plt.axhline(0)
plt.legend(fontsize=11)
plt.ylabel('Percentage Flux Residual')
plt.title('all bands')
out = outfolder+'/plots/'
plt.savefig(out+'Pflux_vs_Globaloffset_all.png')
print out+'Pflux_vs_Globaloffset_all.png'


plt.clf()
ww = ((abs(fitfluxs -fakefluxs)/fakefluxs) < .7)

plt.scatter(peakmags[ww],globalsnoff[ww],alpha=.1,color='green')
ax,ay,aystd = bindata(peakmags[ww],globalsnoff[ww],np.arange(19.8,24.2,.3))
plt.errorbar(ax,ay,aystd,markersize=10,color='black',fmt='o',label='SMP Astrometry')

plt.ylabel('Global SN Offset')
plt.ylim(-0.01,.2)
plt.xlim(19.8,24.2)
plt.xlabel('Peak Magnitude')
plt.legend(fontsize=11)
plt.title('all bands')
out = outfolder+'/plots/'
plt.savefig(out+'Globaloffset_vs_peakmag_all.png')
print out+'Globaloffset_vs_peakmag_all.png'
'''



plt.clf()
ww = (abs(pfitmag-pfm) < 1.)
plt.scatter(pfitmag[ww],pfitmag[ww]-pfm[ww],alpha=.1,color='green')
ax,ay,aystd = bindata(pfitmag[ww],pfitmag[ww]-pfm[ww],np.arange(20,25,.35))
plt.errorbar(ax,ay,aystd,markersize=10,color='green',fmt='o',label='SMP Astrometry')

plt.xlabel('Peak Fit Mag')
plt.ylim(-.08,.08)
plt.plot([19.8,25],[0,0],color='black')
plt.legend(fontsize=11)
plt.xlim(19.8,25)
plt.ylabel('Fit Mag - Fake Mag')
plt.title('r band')
out = outfolder+'/plots/'
plt.savefig(out+'MagResid_vs_PeakFitmag_withoffset_r.png')
print out+'MagResid_vs_PeakFitmag_withoffset_r.png'

if snradec:
	plt.clf()
	ww = (peakmags < 90.) & (abs(pfitmag-pfm) < 1.)
	plt.scatter(pfitmag[ww],pff[ww],alpha=.1,color='green')
	ax,ay,aystd = bindata(pfitmag[ww],pff[ww],np.arange(20,25,.35))
	plt.errorbar(ax,ay,aystd,markersize=10,color='green',fmt='o',label='SMP Astrometry')

	plt.xlabel('Peak Fit Mag')
	plt.ylim(-.05,.05)
	plt.plot([19.8,25],[0,0],color='black')
	plt.legend(fontsize=11)
	plt.xlim(19.8,25)
	plt.ylabel('Percentage Flux Difference')
	plt.title('r bands')
	out = outfolder+'/plots/'
	plt.savefig(out+'Pflux_vs_PeakFitmag_withoffset_r.png')
	plt.savefig(save+'Pflux_vs_PeakFitmag_rband.png')

	print out+'Pflux_vs_PeakFitmag_withoffset_r.png'

'''
plt.clf()
plt.scatter(fitmag,fitmag-fakemag,alpha=.1,color='green')
ax,ay,aystd = bindata(fitmag,fitmag-fakemag,np.arange(20,25,.35))
plt.errorbar(ax,ay,aystd,markersize=10,color='green',fmt='o',label='SMP Astrometry')

plt.xlabel('Fit Mag')
plt.ylim(-.08,.08)
plt.plot([19.8,25],[0,0],color='black')
plt.legend(fontsize=11)
plt.xlim(19.8,25)
plt.ylabel('Fit Mag - Fake Mag')
plt.title('r band')
out = outfolder+'/plots/'
plt.savefig(out+'MagResid_vs_Fitmag_withoffset_r.png')
print out+'MagResid_vs_Fitmag_withoffset_r.png'
'''

if snradec:
	plt.clf()
	so = so[np.isfinite(pff)]
	pff = pff[np.isfinite(pff)]
	plt.scatter(so,pff,alpha=.1,color='green')
	ax,ay,aystd = bindata(so,pff,np.arange(0.,.1,.005))
	plt.errorbar(ax,ay,aystd,markersize=10,color='black',fmt='o')
	print ay
	plt.ylabel('Percentage Flux Difference')
	plt.ylim(-.5,.5)
#plt.plot([19.8,25],[0,0],color='black')
	plt.legend(fontsize=11)
	plt.xlim(0,.1)
	plt.xlabel('Global SN RA DEC Offset (arcsec)')
#plt.title('allbands band')
	out = outfolder+'/plots/'
	plt.savefig(out+'SNOff_vs_Pflux.png')
	plt.savefig(save+'SNOff_vs_Pflux.png')

	print out+'SNOff_vs_Pflux.png'

#raw_input()
plt.clf()
if snradec:
	ww = abs(fitmag-fakemag) < 1.

        #plt.scatter(globalsnoff[ww],fitmag[ww]-fakemag[ww],alpha=.1,color='green')
	ax,ay,aystd,n = bindata(globalsnoff[ww],fitmag[ww]-fakemag[ww],np.arange(-.0525,.1,.005),returnn=True)
	#plt.errorbar(ax,ay,aystd,markersize=10,color='green',fmt='o',label='SMP Astrometry')
	plt.scatter(ax,n)
	plt.xlabel('Global SN Offset')
        #plt.ylim(-.02,.02)
        #plt.plot([0,.1],[0,0],color='black')
	plt.legend(fontsize=11)
	plt.xlim(-.01,.06)
	plt.ylabel('Bin count')
	plt.title('all bands')
	out = outfolder+'/plots/'
	plt.savefig(out+'MagResid_vs_GlobaloffsetN_all.png')
	print out+'MagResid_vs_GlobaloffsetN_all.png'


plt.clf()

fig,axs = plt.subplots(1,2)
ne = (fitfluxs-fakefluxs)/fluxerrs
ww = (fitfluxs!=0) & (fluxerrs > 0) & (hostgal != 0) & (hostgal < 23.)
d = ne[ww][abs(ne[ww]) < 5]
axs[0].hist(ne[ww],bins=np.arange(-4.,4,.25),normed=True,label='RMS: '+str(round(np.sqrt(np.nanmean(np.square(d))),2)))
ne = (fitfluxs-fakefluxs)/fluxerrs
ww = (fitfluxs!=0) & (fluxerrs > 0) & (hostgal != 0) & (hostgal > 26.)
d = ne[ww][abs(ne[ww]) < 5]
axs[1].hist(ne[ww],bins=np.arange(-4.,4,.25),normed=True,label='RMS: '+str(round(np.sqrt(np.nanmean(np.square(d))),2)))

import matplotlib.mlab as mlab
import math
mean = 0
variance = 1
sigma = math.sqrt(variance)
x = np.arange(-5,5,.1)
axs[0].plot(x,mlab.normpdf(x,mean,sigma),color='black',label='Gaussian Normal')
axs[1].plot(x,mlab.normpdf(x,mean,sigma),color='black',label='Gaussian Normal')

axs[0].set_xlim(-4.,4.)
axs[1].set_xlim(-4.,4.)
axs[0].set_ylim(0,.5)
axs[1].set_ylim(0,.5)

axs[0].set_xlabel('Stdev (fitflux-fakeflux)/err')
axs[1].set_xlabel('Stdev (fitflux-fakeflux)/err')

axs[0].set_ylabel('Counts')
axs[1].set_ylabel('Counts')

axs[0].set_title('Host SB < 23.')
axs[1].set_title('Host SB > 26.')

axs[0].legend(fontsize=8)
axs[1].legend(fontsize=8)
plt.tight_layout()
plt.savefig(out+'StdHistSBMag.png')
print out+'StdHistSBMag.png'


plt.clf()
fig,axs = plt.subplots(2,2,sharex=True, sharey=True)
ax = axs.ravel()
ne = (fitfluxs-fakefluxs)/fluxerrs
snmag = copy(fitmag)
step = .01
ran = .25
sbmagrange = np.arange(20,29,step)
for i,f in enumerate(filtsu):
	prms = []
	for sbm in sbmagrange:
		hh = (fitfluxs!=0) & (fluxerrs > 0) & (hostgal != 0) & (hostgal > sbm-ran) & (hostgal < sbm+ran) & (abs(ne) < 5.) & (filtsvec == f)
		try:
			pne = ne[hh]
			prms.append(np.sqrt(np.nanmean(np.square(pne))))
		except:
			prms.append(np.nan)
	ax[i].plot(sbmagrange,prms,label=f+' band')
	ax[i].plot(sbmagrange,sbmagrange*0.+1.,linestyle='--',color='black')
	ax[i].legend(fontsize=8)
	ax[i].set_ylim(.8,3)

ax[2].set_xlabel('SB MAG BIN')
ax[3].set_xlabel('SB MAG BIN')
ax[0].set_ylabel('RMS')
ax[2].set_ylabel('RMS')

fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
plt.setp([a.get_xticklabels() for a in fig.axes[0:2]], visible=False)
plt.setp([a.get_yticklabels() for a in fig.axes[2:0]], visible=False)
plt.savefig(out+'RMSvsSBMag.png')
print out+'RMSvsSBMag.png'



plt.clf()
tt = (fakemag>50) & (fitfluxs!=0) & (fluxerrs > 0)
#fitmag = fitmag[tt]
#magerr = magerr[tt]
tfitfluxs = fitfluxs[tt]
tfakefluxs = fakefluxs[tt]
tfluxerrs = fluxerrs[tt]

#fakemag = fakemag[tt]

d = (tfitfluxs-tfakefluxs)/tfluxerrs
rms = np.sqrt(np.nanmean(np.square(d)))

chisq = (tfitfluxs-tfakefluxs)**2/tfluxerrs**2
chisq = np.mean(chisq[abs(d)<3])

d = d[abs(d)<3]
rms = np.sqrt(np.nanmean(np.square(d)))



print 'mean',np.mean(d)
print 'std',np.std(d)
print 'rms',rms
print 'lend',len(d)
print 'mode',stats.mode(d)
plt.hist((tfitfluxs-tfakefluxs)/tfluxerrs,bins=np.arange(-10,10,.25),normed=True,label='RMS: '+str(round(rms,3))+'\nChiSq (3sig cut) '+str(round(chisq,3))+'\nMedian '+str(round(np.median(d),3))+' +- '+str(round(np.std(d),3)) )
import matplotlib.mlab as mlab                                                                                                                                                               
import math                                                                                                                                                                                      
mean = 0                                                                                                                                                                                        
variance = 1                                                                                                                                                                                   
sigma = math.sqrt(variance)                                                                                                                                                                    
x = np.arange(-5,5,.1)                                                                                                                                                                        
plt.plot(x,mlab.normpdf(x,mean,sigma),color='black',label='Gaussian Normal')                                                                                                               
plt.xlim(-4.,4.)
plt.xlabel('Stdev (fitflux-fakeflux)/err')
plt.ylabel('Counts')
plt.title('Epochs WITHOUT Flux')
plt.legend(fontsize=8)
plt.savefig(out+'StdHistNoFlux.png')
print out+'StdHistNoFlux.png'


plt.clf()
plt.hist(tfitfluxs*10**(.4*(23.9-31.)),bins=np.arange(-2.05,2,.1),normed=True,label='Median '+str(round(np.median(tfitfluxs)*10**(.4*(23.9-31.)),3))+' +- '+str(round(np.std(tfitfluxs*10**(.4*(23.9-31.)))/np.sqrt(len(tfitfluxs)),3)))
plt.xlabel('fitflux')
plt.ylabel('Counts')
plt.xlim(-2,2)
plt.title('Epochs WITHOUT Flux')
plt.legend(fontsize=8)
plt.savefig(out+'FluxHistNoFlux.png')
print out+'FluxHistNoFlux.png'


plt.clf()
tt = (fakemag<50) & (fitfluxs!=0) & (fluxerrs > 0)
#fitmag = fitmag[tt]
#magerr = magerr[tt]
tfitfluxs = fitfluxs[tt]
tfakefluxs = fakefluxs[tt]
tfluxerrs = fluxerrs[tt]

#fakemag = fakemag[tt]

d = (tfitfluxs-tfakefluxs)/tfluxerrs
rms = np.sqrt(np.nanmean(np.square(d)))

d = d[abs(d)<5]
rms = np.sqrt(np.nanmean(np.square(d)))
print 'mean',np.mean(d)
print 'std',np.std(d)
print 'rms',rms
print 'lend',len(d)
print 'mode',stats.mode(d)
plt.hist((tfitfluxs-tfakefluxs)/tfluxerrs,bins=np.arange(-10,10,.25),normed=True,label='RMS: '+str(round(rms,3)))
import matplotlib.mlab as mlab
import math
mean = 0
variance = 1
sigma = math.sqrt(variance)
x = np.arange(-5,5,.1)
plt.plot(x,mlab.normpdf(x,mean,sigma),color='black',label='Gaussian Normal')
plt.xlim(-4.,4.)
plt.xlabel('Stdev (fitflux-fakeflux)/err')
plt.ylabel('Counts')
plt.title('Epochs WITH Flux')
plt.legend(fontsize=8)
plt.savefig(out+'StdHistWithFlux.png')
print out+'StdHistWithFlux.png'


raw_input('stopped')



if not radec:
	sys.exit()

from matplotlib.colors import LogNorm

plt.clf()
x = np.sqrt(globalmagresidh**2)
y = np.sqrt(globalraoffsets**2+globaldecoffsets**2)
wh = (np.isfinite(x)&np.isfinite(y))
y = y*3600
#plt.scatter(globalmagresidh,np.sqrt(globalraoffsets**2+globaldecoffsets**2),alpha=.05,color='green',label='fakemag<22')
#plt.hexbin(x[wh],y[wh],gridsize=400,bins='log')
plt.hist2d(x[wh],y[wh],bins=[3000,20000],norm=LogNorm())
plt.axis([0, .1, 0, .13])

plt.colorbar(label='Counts')
plt.ylabel('abs(Global - SCAMP STARCAT RA/DEC) arcsec')
plt.xlabel('abs(Fit SN Mag - Fake SN Mag)')
plt.axhline(0,color='black')
plt.axvline(0,color='black')
plt.xlim(0.,.1)
plt.ylim(0,.13)
plt.legend()
out = outfolder+'/plots/'
plt.savefig(out+'GlobalStarcatResids.png')
print out+'GlobalStarcatResids.png'

scr = np.sqrt(globalraoffsets**2+globaldecoffsets**2)
scr = scr[np.isfinite(scr)]
scr = scr*3600#convert to arcsec
plt.clf()
plt.hist(scr,bins=np.arange(0,.13,.0025))
plt.gca().set_yscale("log")
plt.xlabel('Global Starcat Offset (arcsec)')
plt.ylabel('Counts')
plt.savefig(out+'GlobalStarcatHist.png')
print out+'GlobalStarcatHist.png'


globalsnra = np.array(globalsnra)
globalsndec = np.array(globalsndec)
#print len(globalmagresid)
#print len(globalsnra)

ww= fakemag < 50.
fitmag = fitmag[ww]
fakemag = fakemag[ww]
globalsnra = globalsnra[ww]
globalsndec = globalsndec[ww]
from matplotlib.colors import LogNorm

rrr = np.sqrt((fitmag-fakemag)**2)
ggg = np.sqrt(globalsnra**2+globalsndec**2)
www = (np.isfinite(rrr) & np.isfinite(ggg))
plt.clf()
#plt.scatter(fitmag-fakemag,np.sqrt(globalsnra**2+globalsndec**2),alpha=.05,color='green')
#plt.hist2d(rrr[www],ggg[www], bins=40, norm=LogNorm())
#plt.hexbin(rrr[www],ggg[www],gridsize=00,bins='log')
plt.hist2d(rrr[www],ggg[www],bins=[2000,400], norm=LogNorm())
plt.axis([0, .1, 0, .1])
plt.colorbar(label='Counts')
plt.ylabel('Global SN RA DEC Resid (pixels)')
plt.xlabel('abs(Fit SN Mag - Fake SN Mag)')
plt.axhline(0,color='black')
plt.axvline(0,color='black')
plt.xlim(0.,.1)
plt.ylim(0,.1)
plt.legend()
out = outfolder+'/plots/'
plt.savefig(out+'GlobalSNResids.png')
print out+'GlobalSNResids.png'


plt.clf()
weights = np.ones_like(ggg[www])/float(len(ggg[www]))
plt.hist(ggg[www],bins=np.arange(0,.5,.01))#,weights=weights)
plt.xlabel(' Global SN Position Offset (pixels)')
plt.ylabel('Counts')
plt.savefig(out+'GlobalSNHist.png')
print out+'GlobalSNHist.png'






########NOW ADD TO LIGHTCURVE FILE #############
lcfpath = '/global/cscratch1/sd/dbrout/DESY1_imgList_fake/'
newpath = '/global/cscratch1/sd/dbrout/v3/updated_lcf_v1/'
print fitfluxso
raw_input()
for sn in np.unique(snnames):
	snloc = lcfpath+sn+'.dat'
	print snloc
	ww = (snnames == sn)
	#raw_input()
	lc.addtolightcurve(snloc,newpath,'SMP',filtsvec[ww],mjdvec[ww],fitfluxso[ww],fluxerrso[ww],fakemago[ww])
	#raw_input()



