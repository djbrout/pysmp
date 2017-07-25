import numpy as np
import os
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt

filts = ['g','r','i','z']
pth = '/project/projectdirs/dessn/dbrout/simdummytest3/np_data/'
myskyerr = []
sexrms = []
hostmag = []
aps = []
cntr = 1
for filt in filts:
    for f in os.listdir(pth+'/'+filt+'/'):
        if cntr > 1000: continue
        if 'smpDict' in f:
            try:
                d = np.load(pth + '/' + filt + '/' + f)
                mse = d['skyerr']
                srms = d['sexrms']
                hm = d['hostgal_sbmag']
                ap = d['aper_skyerr']
            except:
                continue
            cntr += 1
            print cntr

            myskyerr.extend(mse)
            sexrms.extend(srms)
            hostmag.extend(hm)
            aps.extend(ap)

myskyerr = np.array(myskyerr)
sexrms = np.array(sexrms)
hostmag = np.array(hostmag)
aps = np.array(aps)

plt.scatter(hostmag,myskyerr-sexrms,alpha=.4,color='black',label='My Skyerr - SexRMS')
plt.scatter(hostmag,aps-sexrms,alpha=.4,color='orange',label='Aper Skyerr - SexRMS')
plt.axhline(0)
plt.legend()
plt.xlim(20,30)
plt.ylim(-50,150)
plt.xlabel('Hostmag')
plt.ylabel('Fit Skyerr - SexRMS')
plt.savefig('skyerrtest.png')