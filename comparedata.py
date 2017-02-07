import numpy as np
import os

import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt

v6dir = '/pnfs/des/scratch/pysmp/smp_v622/np_data/r/'
v4dir = '/pnfs/des/scratch/pysmp/smp_v42/np_data/r/'

f6 = os.listdir(v6dir)
ff6 = []
for f in f6:
    if 'withSn' in f:
        ff6.append(f)

f4 = os.listdir(v4dir)
ff4 = []
for f in f4:
    if 'withSn' in f:
        ff4.append(f)


commonfiles = []
for f in ff6:
    if f in ff4:
        print f
        commonfiles.append(f)

bigv6stamps = []
bigv4stamps = []
bigv6fakemags = []
bigv4fakemags = []
bigv6mjd = []
bigv4mjd = []
bigv6fitflux = []
bigv4fitflux = []
resid = []
residstamp = []

for i,f in enumerate(commonfiles):
    v6dat = np.load(v6dir+f)
    v4dat = np.load(v4dir + f)
    print i

    try:
        for j,m in enumerate(v6dat['mjd']):

            if m in v4dat['mjd']:
                if m != 0 :
                    v4dat['sky']
                    v6dat['sky']
                    ww = v4dat['mjd'] == m

                    bigv6fakemags.append(v6dat['fakemag'][j])
                    bigv4fakemags.append(v4dat['fakemag'][ww][0])

                    bigv6stamps.append(v6dat['data'][j])
                    bigv4stamps.append(v4dat['data'][ww])

                    bigv6mjd.append(v6dat['mjd'][j])
                    bigv4mjd.append(v4dat['mjd'][ww][0])

                    bigv6fitflux.append(v6dat['modelvec'][j])
                    bigv4fitflux.append(v4dat['modelvec'][ww][0])

                    resid.append(np.sum((v6dat['data'][j,:,:] - v4dat['data'][ww,:,:]  + v4dat['sky'][ww] ).ravel()))
    except:
        print 'column not in file'

bigv6mjd = np.array(bigv6mjd)
bigv4mjd = np.array(bigv4mjd)
resid = np.array(resid)
fakemag = np.array(bigv6fakemags)
print bigv4mjd.shape,bigv6mjd.shape,resid.shape
for r in resid:
    print r

plt.scatter(10**(.4*(31.-fakemag)),resid)
plt.xlabel('Fake Flux')
plt.ylabel('DATA RESIDUAL FLUX')
plt.savefig('resid.png')
