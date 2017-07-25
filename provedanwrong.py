import numpy as np
import os
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt

filts = ['g','r','i','z']
pth = '/project/projectdirs/des/djbrout/simtest/np_data/'
myskyerr = []
sexrms = []
hostmag = []
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
            except:
                continue
            cntr += 1
            print cntr

            myskyerr.extend(mse)
            sexrms.extend(srms)
            hostmag.extend(hm)


myskyerr = np.array(myskyerr)
sexrms = np.array(sexrms)
hostmag = np.array(hostmag)

plt.scatter(hostmag,myskyerr-sexrms,alpha=.4,color='black')
plt.axhline(0)
plt.xlabel('Hostmag')
plt.ylabel('My Sky - SexRMS')
plt.savefig('skyerrtest.png')