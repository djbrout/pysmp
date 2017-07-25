import numpy as np
import os
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt

filt = 'g'
pth = '/project/projectdirs/des/djbrout/simtest/np_data/g/'
myskyerr = []
sexrms = []
hostmag = []
cntr = 1
for f in os.listdir(pth):
    if cntr > 1000: continue
    if 'smpDict' in f:

        print f
        d = np.load(pth+'/'+f)
        try:
            mse = d['skyerr']
            srms = d['sexrms']
            hm = d['hostgal_sbmag']
        except:
            continue
        cntr += 1

        myskyerr.extend(mse)
        sexrms.extend(srms)
        hostmag.extend(hm)


myskyerr = np.array(myskyerr)
sexrms = np.array(sexrms)
hostmag = np.array(hostmag)

plt.scatter(hostmag,myskyerr-sexrms,alph=.4,color=black)
plt.axhline(0)
plt.xlabel('Hostmag')
plt.ylabel('My Sky - SexRMS')
plt.savefig('skyerrtest.png')