import numpy as np
import os

import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import pyfits as pf

v6dir = '/pnfs/des/scratch/pysmp/smp_v627/np_data/r/'
v4dir = '/pnfs/des/scratch/pysmp/smp_v42/np_data/r/'

v6root = '/pnfs/des/persistent/smp/v62/'
v4root = '/pnfs/des/persistent/smp/v4/'

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
bigv6ostamps = []
bigv4ostamps = []
bigv6fakemags = []
bigv4fakemags = []
bigv6mjd = []
bigv4mjd = []
bigv6fitflux = []
bigv4fitflux = []
bigv6sky = []
bigv4sky = []
bigv6zpt = []
bigv4zpt = []
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
    print v6dat.keys()
    #if i > 0:
    #    continue
    #if True:
    #try:
        #print round(np.sum((v6dat['galmodel_params']*mask).ravel())),round(np.sum((v4dat['galmodel_params']*mask).ravel())), round(np.sum((v6dat['galmodel_params']*mask).ravel())-np.sum((v4dat['galmodel_params']*mask).ravel()))
    if True:
        for j,m in enumerate(v6dat['mjd']):

            if m in v4dat['mjd']:
                if m != 0 :
                    v4dat['sky']
                    v6dat['sky']
                    v6dat['x']
                    v4dat['y']
                    ww = v4dat['mjd'] == m
                    #print v6dat['fakemag'][j],v4dat['fakemag'][ww][0]
                    bigv6fakemags.append(v6dat['fakemag'][j])
                    bigv4fakemags.append(v4dat['fakemag'][ww][0])

                    bigv6stamps.append(v6dat['data'][j,:,:])#*10**(.4*(31-v6dat['fitzpt'][j])))
                    bigv4stamps.append(v4dat['data'][ww][0,:,:])#*10**(.4*(31-v4dat['fitzpt'][ww][0])))

                    v6data = pf.getdata(v6dat['datafilenames'][j])
                    v4data = pf.getdata(v4dat['datafilenames'][ww][0])

                    x = int(v6dat['x'][j])
                    y = int(v6dat['y'][j])
                    v6data = v6data[y-15:y+15,x-15:x+15]
                    v4data = v4data[y-15:y+15,x-15:x+15]
                    #print v6data.shape
                    #print v4data.shape

                    bigv6ostamps.append(v6data)

                    bigv6mjd.append(v6dat['mjd'][j])
                    bigv4mjd.append(v4dat['mjd'][ww][0])

                    bigv6fitflux.append(v6dat['modelvec'][j])
                    bigv4fitflux.append(v4dat['modelvec'][ww][0])

                    bigv6sky.append(v6dat['sky'][j])
                    bigv4sky.append(v4dat['sky'][ww][0])

                    v6scalefactor = 10**(.4*(31.-v6dat['fitzpt'][j]))
                    v4scalefactor = 10**(.4*(31.-v4dat['fitzpt'][ww][0]))

                    if v6dat['modelvec'][j] == 0:
                        v6gal = round(np.sum((v6dat['sims'][j,:,:]*mask).ravel()))
                        v4gal = round(np.sum((v4dat['sims'][ww][0, :, :]*mask).ravel()))
                        print v6gal,v4gal,v6gal-v4gal

                    #print v6dat['fitzpt'][j]-v4dat['fitzpt'][ww][0]

                    #print
                    #resid.append(np.sum(((v6dat['data'][j,:,:] - v6dat['sky'][j] - v4dat['data'][ww,:,:]  + v4dat['sky'][ww] )*mask).ravel()))
                    resid.append(np.sum(((v6data*v6scalefactor - v6dat['sky'][j] - v4data*v4scalefactor  + v4dat['sky'][ww] )*mask).ravel()))

                    #residstamp.append(v6dat['data'][j,:,:] - v6dat['sky'][j] - v4dat['data'][ww,:,:]  + v4dat['sky'][ww])
                    residstamp.append((v6data*v6scalefactor - v6dat['sky'][j] - v4data*v4scalefactor  + v4dat['sky'][ww] ))

        #raw_input()

    # except:
    #     print 'column not in file'

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
plt.ylabel('V6-V4 RESIDUAL DATA FLUX')
plt.savefig('resid.png')

plt.clf()
plt.scatter(resid,bigv6fitflux-bigv4fitflux,alpha=.5)
plt.xlabel('DATA Residual Flux (V6-V4)')
plt.ylabel('MODEL Residual Flux (V6-V4)')
plt.xlim(-1000,1000)
plt.ylim(-1000,1000)
plt.plot([-1000,1000],[-1000,1000],color='black')
plt.savefig('residc.png')

from matplotlib.backends.backend_pdf import PdfPages
pdf_pages = PdfPages('v4v6_resid.pdf')
fig = plt.figure()
cntr = 0
plt.clf()
for i,r in enumerate(residstamp):
    #print np.array(r).shape
    if fakemag[i] < 24.:
        fig = plt.figure(figsize=(20,20))
        plt.clf()
        plt.subplot(131)
        plt.title('Fakemag: '+str(round(fakemag[i],2))+'\n Residual Flux: '+str(round(resid[i])))

        plt.imshow(np.array(r[:,:]), cmap='gray', interpolation='nearest')
        plt.colorbar()
        plt.subplot(132)
        plt.title('V6')

        mn = min([min((bigv6stamps[i]-bigv6sky[i]).ravel()),min((bigv4stamps[i]-bigv4sky[i]).ravel())])
        mx = max([max((bigv6stamps[i]-bigv6sky[i]).ravel()),max((bigv4stamps[i]-bigv4sky[i]).ravel())])

        plt.imshow(np.array(bigv6stamps[i]-bigv6sky[i]), cmap='gray', interpolation='nearest',vmin=mn,vmax=mx)
        plt.colorbar()
        plt.subplot(133)
        plt.title('V4')
        plt.imshow(np.array(bigv4stamps[i]-bigv4sky[i]), cmap='gray', interpolation='nearest',vmin=mn,vmax=mx)
        plt.colorbar()

        # try:
        #     cbar = fig.colorbar(ax)
        # except:
        #     print 'could not produce a color bar'

        #if cntr%1 == 0:
        pdf_pages.savefig(fig)

        cntr += 1

#if cntr%9 > 0:
#    pdf_pages.savefig(fig)
pdf_pages.close()