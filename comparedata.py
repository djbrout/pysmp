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

substamp = 30
mask = np.zeros((substamp,substamp))
skyerr_radius = 7.
for x in np.arange(substamp):
    for y in np.arange(substamp):
        if np.sqrt((substamp / 2. - x) ** 2 + (substamp / 2. - y) ** 2) < skyerr_radius:

            mask[int(x), int(y)] = 1.

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
                    #print v6dat['fakemag'][j],v4dat['fakemag'][ww][0]
                    bigv6fakemags.append(v6dat['fakemag'][j])
                    bigv4fakemags.append(v4dat['fakemag'][ww][0])

                    bigv6stamps.append(v6dat['data'][j])
                    bigv4stamps.append(v4dat['data'][ww])

                    bigv6mjd.append(v6dat['mjd'][j])
                    bigv4mjd.append(v4dat['mjd'][ww][0])

                    bigv6fitflux.append(v6dat['modelvec'][j])
                    bigv4fitflux.append(v4dat['modelvec'][ww][0])


                    resid.append(-1*np.sum(((v6dat['data'][j,:,:] - v6dat['sky'][j] - v4dat['data'][ww,:,:]  + v4dat['sky'][ww] )*mask).ravel()))

                    residstamp.append(v6dat['data'][j,:,:] - v6dat['sky'][j] - v4dat['data'][ww,:,:]  + v4dat['sky'][ww])
    except:
        print 'column not in file'

bigv6mjd = np.array(bigv6mjd)
bigv4mjd = np.array(bigv4mjd)
resid = np.array(resid)
fakemag = np.array(bigv6fakemags)
bigv6fitflux = np.array(bigv6fitflux)
bigv4fitflux = np.array(bigv4fitflux)
print bigv4mjd.shape,bigv6mjd.shape,resid.shape
#for r in resid:
#    print r

plt.scatter(10**(.4*(31.-fakemag)),resid)
plt.xlabel('Fake Flux')
plt.ylabel('V4-V6 RESIDUAL DATA FLUX')
plt.savefig('resid.png')

plt.clf()
plt.scatter(resid,bigv4fitflux-bigv6fitflux)
plt.xlabel('DATA Residual Flux (V4-V6)')
plt.ylabel('MODEL Residual Flux (V4-V6)')
plt.savefig('residc.png')

from matplotlib.backends.backend_pdf import PdfPages
pdf_pages = PdfPages('v4v6_resid.pdf')
fig = plt.figure()
cntr = 0
plt.clf()
for i,r in enumerate(residstamp):
    #print np.array(r).shape
    if fakemag[i] < 24.:
        fig = plt.figure()
        plt.clf()
        ax = plt.subplot(111)
        print np.array(r[0,:,:]).shape
        print max(np.array(r[0,:,:]).ravel()), min(np.array(r[0,:,:]).ravel())


        ax.imshow(np.array(r[0,:,:]), cmap='gray', interpolation='nearest')

        try:
            cbar = fig.colorbar(ax)
        except:
            print 'could not produce a color bar'
        ax.set_title('Fakemag: '+str(round(fakemag[i],2))+' Residual Flux: '+str(round(resid[i])))

        #if cntr%1 == 0:
        pdf_pages.savefig(fig)

        cntr += 1

#if cntr%9 > 0:
#    pdf_pages.savefig(fig)
pdf_pages.close()