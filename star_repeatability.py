

import os
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(threshold=np.nan)
import sys

out = '/global/u1/d/dbrout/PySMP/paperplots/'
f = np.load('/global/cscratch1/sd/dbrout/v3/smp_y1y2_shallow_v3_40des/np_data/r/des_fake_00229567_r_mcmc_input.npz')
#print f.keys()

nominalzpt = 31.

tmjds = f['mjd']
zptfiles = f['zpt_files']
#print zptfiles
#raw_input()
mags = []
catmags = []
ras = []
decs = []
mjds = []
dict = {}
fitzpts = []
for m,z in zip(tmjds,zptfiles):
    if m > 0:
        a = np.load(z)
        print a.keys()
        #print a['mpfit_zpt']
        #print len(a['mpfit_zpt'])
        #print len(a['mpfit_mag'])
        raw_input()
        fitzpt = a['mpfit_zpt']
        mag = a['mpfit_mag']
        catmag = a['cat_mag']
        ra = a['ra']
        dec = a['dec']
        for r,d,mm in zip(ra,dec,mags):
            try:
                dict[(r,d)].append(mm)
            except:
                dict[(r,d)] = []
                dict[(r,d)].append(mm)

        mags.extend(mag)
        catmags.extend(catmag)
        ras.extend(ra)
        decs.extend(dec)
        mjds.extend(mag*0. + m)
        fitzpts.extend(ra*0. + fitzpt)
mags = np.array(mags)
ras = np.array(ras)
decs = np.array(decs)
fitzpts = np.array(fitzpts)
catmags = np.array(catmags)
mags = mags - nominalzpt + fitzpts



print dict.keys()
print len(dict.keys())
mag_minus_mean = []
fig,ax = plt.subplots(10,10,figsize=(50,50))
axs = ax.ravel()
cntr=-1
tcatmags = []
goodcounts = 0.
badcounts = 0.
totalcounts = 0.
for k in dict.keys():
    cntr+=1
    tra = k[0]
    tdec = k[1]
    maghere = mags[(ras == tra) & (decs == tdec)]
    catmaghere = catmags[(ras == tra) & (decs == tdec)]
    tmean = np.mean(maghere)

    if np.std(maghere-tmean) < .04:
        mag_minus_mean.extend(maghere-tmean)
        tcatmags.extend(catmaghere)
        goodcounts += 1.
        totalcounts += 1.
    else:
        badcounts += 1.
        totalcounts += 1.
    try:
        axs[cntr].hist(maghere-tmean,bins=np.arange(-.1,.1,.01))
    except:
        continue

print 'Fraction of stars that dont make repeatability cut of std<0.04:',badcounts/totalcounts
raw_input()
mag_minus_mean = np.array(mag_minus_mean)
tcatmags = np.array(tcatmags)
plt.tight_layout()
plt.savefig(out+'star_repeatability.png')

plt.clf()
fig, ax = plt.subplots(1,1,figsize=(10,10))
std = np.std(mag_minus_mean)
ax.hist(mag_minus_mean,bins=np.arange(-.1,.1,.005),label='Std: '+str(round(std,3)))
ax.legend()
ax.set_xlabel('r-r_mean')
ax.set_ylabel('Count')
plt.savefig(out+'star_repeatabilityall.png')
print '/global/cscratch1/sd/dbrout/v3/smp_y1y2_shallow_v3_40globalstars/plots/r/star_repeatabilityall.png'

print min(tcatmags)
print max(tcatmags)
raw_input()
mag_lims = [17.,18.,19.,20.,21.]
plt.clf()
fig, axs = plt.subplots(2,2,figsize=(10,10))
for i,ax in enumerate(axs.ravel()):
    ppp = mag_minus_mean[(tcatmags>mag_lims[i])&(tcatmags<=mag_lims[i+1])]
    std = np.std(ppp)
    print std
    print mag_lims[i]
    print mag_lims[i+1]
    ax.hist(ppp,bins=np.arange(-.1,.1,.003),label='Star Mag '+str(mag_lims[i])+'-'+str(mag_lims[i+1])+' \nStd: '+str(round(std,3)))
    ax.legend(fontsize=8)
    ax.set_xlim(-.05,.05)
    ax.set_xlabel('r-r_mean')
    ax.set_ylabel('Count')

plt.savefig(out+'star_repeatabilitybins.png')
print '/global/cscratch1/sd/dbrout/v3/smp_y1y2_shallow_v3_40globalstars/plots/r/star_repeatabilitybins.png'

print 'Overall: ',np.std(mag_minus_mean)
