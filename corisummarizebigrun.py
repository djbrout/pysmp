
import numpy as np
import exceptions
import os
import sys
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import NullFormatter
import iterstat
import pyfits as pf
import dilltools as dt
from copy import copy
from scipy.stats import sigmaclip


def go(fakedir,resultsdir,cacheddata,cd,filter,tfield,dostars,deep_or_shallow,isfermigrid=False,real=False):

    if isfermigrid:
        useifdh = True
    else:
        useifdh = False
    tmpwriter = dt.tmpwriter(useifdh=useifdh)

    #getparametriczpt("/global/cscratch1/sd/dbrout/v6/","/global/cscratch1/sd/dbrout/v6/stardata_"+filter)

    if not cacheddata:
        #dostars = True
        cd = cd[0]
        if dostars:
            grabstardata("/project/projectdirs/dessn/dbrout/zptfiles/","/global/cscratch1/sd/dbrout/v7/stardata",tfield,filter)
            stardata = np.load('/global/cscratch1/sd/dbrout/v7/stardata.npz')
            plotstarrms(stardata['starflux'], np.sqrt(stardata['starfluxerr'] ** 2), stardata['starzpt'],
                    stardata['catmag'], stardata['chisq'], stardata['rmsaddin'], stardata['sky'], stardata['skyerr'],
                    stardata['poisson'],stardata['ids'],stardata['centroidedras'],stardata['centroideddecs'],
                        stardata['fwhm'],stardata['zptscat'],stardata['mjd'],stardata['field'],stardata['ccd'],stardata['band'],
                    title=tfield+'_'+filter+'_',outdir='/global/cscratch1/sd/dbrout/v6/')
        data = grabdata(tmpwriter,resultsdir,cd,tfield,filter=filter,real=real)
    else:

        #dostars = Trueasdf
        if dostars:
            stardata = np.load('/global/cscratch1/sd/dbrout/v7/stardata.npz')
            #plotstarlc(stardata['starflux'],stardata['starfluxerr'],stardata['starzpt'],stardata['ids'],stardata['mjd'],stardata['catmag'])

            plotstarrms(stardata['starflux'], np.sqrt(stardata['starfluxerr'] ** 2), stardata['starzpt'],
                        stardata['catmag'], stardata['chisq'], stardata['rmsaddin'], stardata['sky'], stardata['skyerr'],
                        stardata['poisson'],stardata['ids'],stardata['centroidedras'],stardata['centroideddecs'],
                        stardata['fwhm'],stardata['zptscat'],stardata['mjd'],stardata['field'],stardata['ccd'],stardata['band'],
                        title=tfield+'_'+filter+'_',outdir='/global/cscratch1/sd/dbrout/v6/')
        data = {'filter':[]}
        for k in np.load(cd[0]).keys():
            data[k] = []
        for c in cd:
            print c
            td = np.load(c)
            for key in td.keys():
                print key, len(td[key])
                data[key].extend(td[key])
            filt = c.split('.npz')[0].split('_')[-1]
            filts =  [ filt for x in range(len(td['Flux']))]
            data['filter'].extend(filts)
    #print len(data['Flux'])
    #raw_input()
    #for key in data.keys():
    #    data[key] = np.array(data[key])
    if not os.path.exists(resultsdir+'/Summary/'):
        os.mkdir(resultsdir+'/Summary/')
    if not os.path.exists(resultsdir+'/Summary/'+filter+'/'):
        os.mkdir(resultsdir+'/Summary/'+filter+'/')
    #raw_input()
    plotpercentageresid(data['Flux'],data['Fluxerr'],data['FakeMag'],data['FitZPT'],data['FakeZPT'], data['diffimflux'],data['diffimfluxerr'],
                        data['sky'],data['skyerr'],data['DPMJD'],data['Chisq'],data['imfiles'],data['ra'],data['dec'],
                        data['image_stamp'],resultsdir+'/Summary/'+filter+'/',data['fakefiles'],data['HostMag'],
                        filter,data['FakeZPT'],data['rmsaddin'],data['filter'],data['flag'],data['fwhm'],real=real)

    #data['FakeMag'] = np.array(data['FakeMag'])
    #print len(data['Flux']),len(data['FakeMag'][data['FakeMag']<26.])
    #raw_input('ppp')
    plotsigmaresid(data['Flux'],data['Fluxerr'],data['FakeMag'], data['FitZPT'], data['FakeZPT'],data['HostMag'],
                   data['Chisq'],data['rmsaddin'],data['field'],resultsdir+'/Summary/'+filter+'/',data['rmsaddin'],
                   data['diffimflux'], data['diffimfluxerr'],filter,data['filter'],deep_or_shallow,data['sky'],data['skyerr'],
                   data['fwhm'],data['imfiles'],data['DPMJD'],data['flag'],data['snid'],data['mjd'],real=real)#resultsdir)


def lookup_rms_addin(smpfile,obsid):
    pass



def getparametriczpt(imagedir,outfile):
    bigdata = {'starflux': [], 'starfluxerr': [], 'starzpt': [], 'diffimzpt': [], 'catmag': [], 'chisq': [],
               'rmsaddin': [],'resid':[],
               'sky': [], 'skyerr': [], 'psf': [], 'poisson': [], 'ids': [], 'centroidedras': [], 'centroideddecs': [],
               'numzptstars': [], 'fwhm': []}
    zptfiles = []
    cntr = 0
    goodbigdata = copy(bigdata)
    for dirName, subdirList, fileList in os.walk(imagedir):
        if cntr > 10000.: break
        # print('Found directory: %s' % dirName)
        for fname in fileList:
            # print fname
            if 'globalstar.npz' in fname:
                # print('\t%s' % fname)
                # print os.path.join(imagedir,dirName,fname)
                #if not 'SN-S2' in fname: continue
                #    if not 'SN-S1' in fname: continue
                # try:
                #     os.system('cp ' + os.path.join(imagedir, dirName, fname) + ' test.npz')
                #     zptdata = np.load('test.npz')
                # except:
                #     print 'could not load'
                #     continue
                # zptdata = np.load('/pnfs/des/persistent/smp/v2/20131119_SN-S2/r_21/SNp1_256166_SN-S2_tile20_r_21+fakeSN_rband_dillonzptinfo_globalstar.npz')
                # print zptdata.keys()
                # raw_input()
                if not fname in zptfiles:
                    try:
                        # if True:
                        test = zptdata['chisqu']
                        test = zptdata['fwhm']
                        try:
                            if len(zptdata['flux_star_std']) != len(zptdata['flux_starh']):
                                print 'skippeddddd'
                                continue
                        except:
                            print 'FAIL'
                            continue

                        try:
                            cm = zptdata['cat_magsmp']
                        except:
                            print 'FAO:'
                            continue
                        if len(zptdata['cat_magsmp']) != len(zptdata['flux_starh']):
                            print 'skipperrrrr'
                            continue
                        psfs = zptdata['psfs']

                        if len(psfs) != len(zptdata['cat_magsmp']):
                            print 'skippeppppppp'
                            continue
                        # if True:

                        # if max(zptdata['cat_mag'])>21.1:
                        #    continue
                        # if True:
                        bigdata = copy(goodbigdata)

                        bigdata['ids'].extend(zptdata['ids'])
                        #bigdata['centroidedras'].extend(zptdata['centroidedras'])
                        #bigdata['centroideddecs'].extend(zptdata['centroideddecs'])

                        #bigdata['skyerr'].extend(zptdata['skyerr'])
                        #bigdata['sky'].extend(zptdata['sky'])
                        #bigdata['starflux'].extend(zptdata['flux_starh'])
                        #bigdata['starfluxerr'].extend(zptdata['flux_star_std'])
                        #bigdata['starzpt'].extend(zptdata['flux_starh'] * 0. + zptdata['fit_zpt'])
                        #bigdata['fwhm'].extend(zptdata['flux_starh'] * 0. + zptdata['fwhm'])

                        # bigdata['diffimzpt'].extend(zptdata['fakezpt'])
                        #psfs = zptdata['psfs']
                        #for i in range(len(psfs)):
                        #    bigdata['psf'].append(psfs[i, :, :])
                        #    bigdata['poisson'].append(
                        #        np.sqrt(np.sum(psfs[i, :, :].ravel() ** 2 * zptdata['flux_starh'][i])))
                        #    # print zptdata['flux_starnormm'][i],zptdata['flux_star_std'][i],bigdata['poisson'][-1]
                        #    # raw_input()

                        fs = zptdata['flux_starh']
                        zp = zptdata['fit_zpt']
                        ww = (cm < 19.) & (cm > 17.)

                        # plt.scatter(cm[ww],float(zp) - cm[ww] - 2.5*np.log10(fs[ww]))
                        # plt.scatter(cm[ww],- 2.5*np.log10(fs[ww]))
                        # plt.savefig('testzpt.png')
                        md, std = iterstat.iterstat(float(zp) - cm[ww] - 2.5 * np.log10(fs[ww]),
                                                    startMedian=True, sigmaclip=3, iter=10)

                        if std == 0: continue

                        ww = (cm > 15.5)

                        bigdata['catmag'].extend(cm[ww])

                        bigdata['resid'].extend((float(zp) - cm[ww] - 2.5 * np.log10(fs[ww]))/(float(std)/.01))


                        print 'worked',cntr, std
                        cntr += 1
                        #bigdata['numzptstars'].extend(zptdata['flux_starh'] * 0. + len(cm[ww]))
                        #bigdata['rmsaddin'].extend(zptdata['flux_starh'] * 0. + std / np.sqrt(len(cm[ww])))
                        # print 'read in ',fname
                        #zptfiles.append(fname)

                    except:
                        print 'failed'

    np.savez('zptparam.npz', **bigdata)
    plt.clf()
    plt.scatter(bigdata['catmag'],bigdata['resid'],alpha=.05,color='blue')
    ax, ay, aystd, n = bindata(np.array(bigdata['catmag']),np.array(bigdata['resid']),np.arange(15.5, 21.5, .5),returnn=True)
    aystd *= n**.5
    plt.plot(ax, ay, linewidth=3, color='orange', label='SMP')
    plt.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP')
    plt.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP')
    plt.axhline(0,color='black')
    plt.xlabel('catmag')
    plt.ylim(-.2,.2)
    plt.xlim(15.5,21)
    plt.ylabel('scaled residuals')
    plt.savefig('zptresids.png')
    print 'saved'
    #sys.exit()

def grabstardata(imagedir,outfile,tfield,filt):
    bigdata = {'starflux': [], 'starfluxerr': [], 'starzpt': [], 'diffimzpt':[], 'catmag': [], 'chisq': [], 'rmsaddin': [], 'zptscat':[],
               'sky':[], 'skyerr': [],'psf':[],'poisson':[],'ids':[],'centroidedras':[],'centroideddecs':[],
               'numzptstars':[], 'fwhm':[],'mjd':[],'band':[],'field':[],'ccd':[]}
    zptfiles = []
    cntr = 0
    goodbigdata = copy(bigdata)
    imglist = np.array(os.listdir(imagedir),dtype='str')
    numimages = len(imglist)

    for fname in imglist[np.array(np.random.uniform(0,numimages-1,size=100000),dtype='int')]:
        #if cntr > 5000.: break
        #print('Found directory: %s' % dirName)
        if True:
        #for fname in fileList:
            #print fname
            if 'globalstar.npz' in fname:
                #print('\t%s' % fname)
                #print os.path.join(imagedir,dirName,fname)
                #if not tfield in fname: continue
                #if not '_'+filt+'_' in fname: continue

                field = fname.split('-')[1].split('_')[0]
                ccd = fname.split('_')[5]
                band = fname.split('_')[4]

                zptdata = np.load(os.path.join(imagedir,fname))
                #print zptdata.keys()
                #raw_input()
                #    if not 'SN-S1' in fname: continue
                # try:
                #     os.system('cp ' + os.path.join(imagedir,dirName,fname) + ' test.npz')
                #     zptdata = np.load('test.npz')
                # except:
                #     print 'could not load'
                #     continue
                #zptdata = np.load('/pnfs/des/persistent/smp/v2/20131119_SN-S2/r_21/SNp1_256166_SN-S2_tile20_r_21+fakeSN_rband_dillonzptinfo_globalstar.npz')
                #print zptdata.keys()
                #raw_input()
                if not fname in zptfiles:
                    try:
                    #if True:
                        test = zptdata['chisqu']
                        test = zptdata['fwhm']
                        test = zptdata['zptscat']
                        test = zptdata['flux_starmp']
                        test = zptdata['fit_magmy']
                        #print zptdata['flux_star_stdlm']
                        #raw_input()
                        try:
                            if len(zptdata['flux_star_stdmp']) != len(zptdata['flux_starmp']):
                                print 'skippeddddd'
                                continue
                        except:
                            print 'FAIL'
                            continue

                        try:
                            cm = zptdata['cat_magsmp']
                        except:
                            print 'FAO:'
                            continue
                        if len(zptdata['cat_magsmp']) != len(zptdata['flux_starmp']):
                            print 'skipperrrrr'
                            continue
                        psfs = zptdata['psfs']

                        if len(psfs) != len(zptdata['cat_magsmp']):
                            print 'skippeppppppp'
                            continue
                        #if True:

                        #if max(zptdata['cat_mag'])>21.1:
                        #    continue
                        #if True:
                        bigdata = copy(goodbigdata)
                        bigdata['field'].extend([field for x in zptdata['ids'] ])
                        bigdata['ccd'].extend([ccd for x in zptdata['ids']])
                        bigdata['band'].extend([band for x in zptdata['ids']])

                        bigdata['ids'].extend(zptdata['ids'])
                        bigdata['centroidedras'].extend(zptdata['centroidedras'])
                        bigdata['centroideddecs'].extend(zptdata['centroideddecs'])

                        bigdata['skyerr'].extend(zptdata['skyerr'])
                        bigdata['sky'].extend(zptdata['sky'])
                        bigdata['starflux'].extend(zptdata['flux_starmp'])
                        bigdata['starfluxerr'].extend(zptdata['flux_star_stdmp'])
                        bigdata['starzpt'].extend(zptdata['flux_starmp']*0. + zptdata['fit_zpt'])
                        bigdata['fwhm'].extend(zptdata['flux_starmp']*0. + zptdata['fwhm'])
                        bigdata['mjd'].extend(zptdata['flux_starmp']*0. + zptdata['mjd'])
                        bigdata['catmag'].extend(zptdata['cat_magsmp'])
                        #bigdata['diffimzpt'].extend(zptdata['fakezpt'])
                        psfs = zptdata['psfs']

                        #for i in range(len(psfs)):
                        #    bigdata['psf'].append(psfs[i,:,:])
                        #    bigdata['poisson'].append(np.sqrt(np.sum(psfs[i,:,:].ravel()**2*zptdata['flux_starmp'][i])))
                            #print zptdata['flux_starnormm'][i],zptdata['flux_star_std'][i],bigdata['poisson'][-1]
                            #raw_input()

                        fs = zptdata['flux_starmp']
                        zp = zptdata['fit_zpt']
                        ww = (cm < 19.) & (cm > 16.)

                        #plt.scatter(cm[ww],float(zp) - cm[ww] - 2.5*np.log10(fs[ww]))
                        # plt.scatter(cm[ww],- 2.5*np.log10(fs[ww]))
                        # plt.savefig('testzpt.png')
                        md, std = iterstat.iterstat(float(zp) - cm[ww] - 2.5*np.log10(fs[ww]),
                                                     startMedian=True, sigmaclip=3, iter=10)

                        #print zptdata['zptscat'],np.std(float(zp) - cm[ww] - 2.5*np.log10(fs[ww]))
                        #raw_input()
                        #print 'worked now std',std/np.sqrt(len(cm[ww]))
                        bigdata['numzptstars'].extend(zptdata['flux_starmp']*0. + len(cm[ww]))
                        bigdata['rmsaddin'].extend(zptdata['flux_starmp']*0. + std)
                        bigdata['zptscat'].extend(zptdata['flux_starmp']*0.+np.std(float(zp) - cm[ww] - 2.5*np.log10(fs[ww])))
                        #print 'read in ',fname
                        zptfiles.append(fname)

                        bigdata['chisq'].extend(zptdata['chisqu'])

                        cntr += 1
                        print 'CNTR',cntr

                        goodbigdata = copy(bigdata)
                    except:
                        print 'FAILED', fname
                        pass

    # try:
    #     bigdata['centroidedras'] = np.array(bigdata['centroidedras'])
    #     bigdata['ids'] = np.array(bigdata['ids'])
    #     stds = []
    #     for i in zptdata['ids']:
    #         #print i
    #         #print bigdata['centroidedras'].shape
    #         #print bigdata['ids'] == i
    #         print i,np.mean(bigdata['centroidedras'][bigdata['ids'] == i]),np.std(bigdata['centroidedras'][bigdata['ids'] == i])
    #         if np.std(bigdata['centroidedras'][bigdata['ids'] == i]) > 0.:
    #             stds.append(np.std(bigdata['centroidedras'][bigdata['ids'] == i]))
    #
    #     print 'std is', np.mean(stds)
    # except:
    #     print 'ids not in archiv'
    #sys.exit()

    np.savez(outfile, **bigdata)
    print np.unique(bigdata['mjd'])
    print np.unique(bigdata['field'])
    #sys.exit()
    #os.system('ifdh rm ' + outfile)
    #os.system('ifdh cp ' + 'dat.dat' + ' ' + outfile)
    #os.system('rm dat.dat')
    #sys.exit()

def grabdata(tmpwriter,resultsdir,cd,tfield,filter = 'g',oldformat=False,real=False):

    # dofakefilt,dofakemjd,dofakemag,dofakera,dofakedec = np.loadtxt('data/grepalldofake_'+filter+'.txt',usecols=(3, 9, 10, 14, 15), unpack=True, dtype='string', skiprows=0)
    # dofakemjd = np.array(dofakemjd,dtype='float')
    # dofakemag = np.array(dofakemag,dtype='float')
    # dofakera = np.array(dofakera,dtype='float')
    # dofakedec = np.array(dofakedec,dtype='float')


    # diffzpts = dt.readcol('ZP.out', delim=' ')
    # dz = np.array(diffzpts['zpt'],dtype='float')
    # dccd = np.array(diffzpts['ccd'],dtype='float')
    # dexp = np.array(diffzpts['expnum'],dtype='float')
    # dfield = np.array(diffzpts['field'],dtype='str')



    # expnum,dofakeccds,dofakefilt2,dofakeid,dofakemjd2,dofakemag2,dofaketflux,dofakeflux,dofakera2,dofakedec2 = np.loadtxt('data/doFake.out',
    #                                         usecols=(1,2,3,5, 9, 10,11,12, 14,15), unpack=True, dtype='string', skiprows=1)
    # #dofakemjd2 = np.array(dofakemjd2, dtype='float')
    # dofakemag2 = np.array(dofakemag2, dtype='float')
    # #dofakera2 = np.array(dofakera2, dtype='float')
    # #dofakedec2 = np.array(dofakedec2, dtype='float')
    # dofaketflux = np.array(dofaketflux, dtype='float')
    # dofakeflux = np.array(dofakeflux, dtype='float')
    # #dofakerand = dofakeflux/dofaketflux
    # dofakeexpnum = np.array(expnum, dtype='float')
    # #dofakeccds = np.array(dofakeccds,dtype='float')
    # #dofakefilt2 = np.array(dofakefilt2,dtype='string')
    # dofakeid = np.array(dofakeid,dtype='float')
    # ##print dofakemjd
    # dofakezpt = dofakemag2 + 2.5*np.log10(dofaketflux)

    import pandas as pd

    dofakedata = pd.read_csv('data/bigdoFake.out', delim_whitespace=True, header=0)
    # print dofakedata
    dofakeexpnum = dofakedata['EXPNUM'].values
    dofakemag2 = dofakedata['TRUEMAG'].values
    dofaketflux = dofakedata['TRUEFLUXCNT'].values
    dofakeflux = dofakedata['FLUXCNT'].values
    dofakeid = dofakedata['FAKEID'].values
    print dofakemag2,''
    dofakezpt = dofakemag2 + 2.5 * np.log10(dofaketflux)
    print 'done reading dofake'


    #raw_input('dofakemjd')
    files = os.listdir(os.path.join(resultsdir, 'lightcurves'))
    smpfiles = []
    for f in files:
        if filter+'.smp' in f:
            print f
            smpfiles.append(os.path.join(resultsdir, 'lightcurves', f))

    #print "Found " + str(len(smpfiles)) + " .smp files"

    if not os.path.exists(os.path.join(resultsdir,'Summary')):
        os.makedirs(os.path.join(resultsdir,'Summary'))
    #outfile = os.path.join(resultsdir,'Summary','sumdata.npz')
    outfile = cd
    bigdata = {'Flux':[],'Fluxerr':[],'FakeMag':[],'FitZPT':[],'FakeZPT':[],'HostMag':[],'Chisq':[],'DPMJD':[],
               'starflux':[],'starfluxerr':[],'starzpt':[],'catmag':[],'rmsaddin':[],'field':[],'sky':[],'imfiles':[],
               'mjd':[],'fakefile':[],'ra':[],'dec':[],'image_stamp':[],'fakefiles':[],'diffzpt':[],'diffimflux':[],
               'diffimfluxerr':[],'fakeid':[],'skyerr':[],'acceptance':[],'filter':[],'fwhm':[],'flag':[],'snid':[]}
    zptfiles = []
    #deep = 0
    print smpfiles
    tot = len(smpfiles)
    cntr = 0
    for f in smpfiles[:]:
        #if f == '/project/projectdirs/des/djbrout/106x3/lightcurves/des_fake_00212070_g.smp':
        #    raw_input('filee')

        #print f
        #raw_input()
        fakeid = f.split('_')[-2]
        snid = f.split('_')[-2]
        cntr += 1
        if cntr > 5000: continue
        #if cntr == 34: continue
        #if cntr == 53: continue
        #if not '_r.smp' in f: continue
        print cntr, 'of',tot
        deep = 0
        #os.system('cp '+f+' test.npz')
        data = dt.readcol(f)
        #print data.keys()
        #print data['PSF_FILE']
        #raw_input('stop')

        #print tra
        try:
            if len(data['RA']) == 0:
                print 'empty'
                continue
        except:
            print 'empty'
            continue
        tra = data['RA']
        #print data.keys()
        #raw_input()
        #print tra[0]
        #dra = np.zeros(len(dofakera))+tra[0]
        #cra = np.isclose(dra,dofakera,atol=1.e-0)
        #tdec = data['DEC']
        #print tra[0],tdec[0],
        #ddec = np.zeros(len(dofakedec))+tdec[0]
        #cdec = np.isclose(ddec,dofakedec,atol=1.e-0)
        #print fakeid
        # fw = np.where(dofakeid == fakeid)
        # print np.unique(dofakeid), int(fakeid)
        # print fw,len(fw)
        # raw_input()
        #if fakeid in dofakeid:
        #    raw_input(fakeid)
        #    continue
        #    #print fakeid
        #print dofakeid[dofakeid == fakeid]
        #raw_input()

        #print len(np.unique(dofakeid))
        #raw_input('dofakeid')
        if not real:
            if not oldformat:
                if not int(fakeid) in dofakeid:
                    print 'dddd'
                    fakemag = data['FAKEMAG']*0. + 99.
                else:
                    fakemag = data['FAKEMAG']
            else:
                fakemag = data['FAKEMAG']
        else:
            pass
        #print fakemag
        #raw_input()
        #print data['DPMJD']
        #raw_input()
        if not real:
            if len(data['DPMJD'][data['DPMJD'] > 300.]) < 2:
                print 'bad DPMJD'*20
                continue
        if np.min(data['FLUX']) < -100000:
            print 'BAD FLUX'*20
            continue

        skipnewfakemag = False
        if real:
            skipnewfakemag = True
        #if not skipnewfakemag:
        if not skipnewfakemag:
            newfakemag = []
            for imf,fm,x,y in zip(data['IMAGE_FILE'],fakemag,data['XPOS'],data['YPOS']):

                #print imf
                try:
                    exn = imf.split('/')[-1].split('_')[1]
                    #ccd = float(imf.split('_')[7].split('+')[0])
                except:
                    #print 'skipping'
                    newfakemag.append(99)
                    bigdata['fakeid'].append(99.)
                    bigdata['FakeZPT'].append(31.)
                    continue
                    #continue
                # if fm == 99.:
                #     print 'fakemag already is 99'
                #     newfakemag.append(99)
                #     continue
                #dra = np.zeros(len(dofakera2)) + tra[0]
                #cra = np.isclose( dofakera2,dra, atol=1.e-3)
                #tdec = data['DEC']
                #ddec = np.zeros(len(dofakedec2)) + tdec[0]
                #cdec = np.isclose(dofakedec2, ddec, atol=1.e-3)



                expn = (dofakeexpnum == float(exn))

                #ccdw = (dofakeccds == ccd)
                #filtw = (dofakefilt2 == 'g')
                #print exn, fm

                # print dofakemag2[]
                #www = (expnum == float(exn)) & (np.isclose(float(fm), dofakemag2, atol=1.e-3))
                dfw = dofakeid == int(fakeid)
                www = expn & dfw
                #print dofakeid[www]
                #raw_input()
                # ifm = (dofakemag2 == fm)
                #print exn, fm
                #print fm,dofakemag2[expn],exn
                #raw_input()
                #print  len(dofakemag2[expn]),exn
                if not len(dofakemag2[www]) == 1:
                    print '99',
                    newfakemag.append(99.)
                    bigdata['fakeid'].append(99.)
                    bigdata['FakeZPT'].append(31.)
                else:
                    #nfm = float(fm) + 2.5*np.log10(dofaketflux[www][0]) - 2.5*np.log10(dofakeflux[www][0])
                    print fm,
                    newfakemag.append(fm)

                    bigdata['fakeid'].append(dofakeid[www][0])
                    bigdata['FakeZPT'].append(dofakezpt[www][0])

                #raw_input()

            #print np.array(newfakemag)-fakemag
            #raw_input()
            fakemag = np.array(newfakemag,dtype='float')

        '''
        sn = f.split('/')[-1][0:17]+'.dat'
        snd = open('/pnfs/des/scratch/pysmp/DESY1_imgList_fake/'+sn,'r').read()
        if not '-S1' in snd:
            if not '-S2' in snd:
                deep = 1
                print 'deep'
                continue
        '''
                #raw_input()
        #print data.keys()
        #raw_input()

        if True:

            bigdata['rmsaddin'].extend(data['ZPTERR'])

            #print data.keys()
            #raw_input()

            bigdata['Flux'].extend(data['FLUX'])
            bigdata['Fluxerr'].extend(data['FLUXERR'])
            if not real:
                if not skipnewfakemag:
                    bigdata['FakeMag'].extend(fakemag)
                else:
                    bigdata['FakeMag'].extend(data['FAKEMAG'])
            bigdata['FitZPT'].extend(data['ZPT'])
            #print fakemag.shape
            #print data['CHI2'].shape
            #raw_input()
            #print data['ZPT'],data['FAKEZPT']
            #raw_input('aaa')
            #bigdata['FakeZPT'].extend(data['FAKEZPT'])
            bigdata['Chisq'].extend(data['CHI2'])
            bigdata['sky'].extend(data['SKY'])
            bigdata['DPMJD'].extend(data['DPMJD'])
            bigdata['mjd'].extend(data['MJD'])
            bigdata['imfiles'].extend(data['IMAGE_FILE'])
            bigdata['fakefiles'].extend([f for i in range(len(data['FLUX']))])
            for p in data['PSF_FILE']:
                try:
                    #print pf.open(p)[1].header['PSF_FWHM'],pf.open(p)[1].header['PSF_FWHM']*  2.235 * 0.27
                    #print p
                    bigdata['fwhm'].append(pf.open(p)[1].header['PSF_FWHM'] *  2.235 * 0.27)
                except:
                    print 'excepted',p
                    bigdata['fwhm'].append(-999.)
            bigdata['diffimflux'].extend(data['DIFFIM_FLUX'])
            bigdata['diffimfluxerr'].extend(data['DIFFIM_FLUXERR'])
            bigdata['skyerr'].extend(data['SKYERR'])
            bigdata['filter'].extend([filter for i in range(len(data['SKYERR']))])
            #print data.keys()
            #raw_input()
            bigdata['snid'].extend(data['DESCRIPTIVE_FLAG']*0+int(snid))
            bigdata['flag'].extend(data['DESCRIPTIVE_FLAG'])
            #bigdata['fwhm'].extend(data['FWHM'])

            #print data['FLUX']
            #print np.mean(data['CHI2']), bigdata['fakeid'][-10:]

            # #print data['IMAGE_FILE']
            # for e in data['IMAGE_FILE']:
            #     try:
            #         #print e
            #         expnum = float(e.split('_')[3])
            #         ccd = float(e.split('_')[7].split('+')[0])
            #         #bigdata['expnums'].append(expnum)
            #         #bigdata['ccds'].append(ccd)
            #         #if dz[(dccd == ccd) & (dexp == expnum)]:
            #         #    print 'pass'
            #         #else:
            #         #    print 'fail'
            #         diffzpt = dz[(dccd == ccd) & (dexp == expnum)]
            #         #print ccd, expnum,diffzpt
            #         #raw_input()
            #         #print diffzpt[0]
            #         bigdata['diffzpt'].append(diffzpt[0])
            #     except:
            #         bigdata['diffzpt'].append(0)
            #         #print 'nanana'
            # #raw_input()

            bigdata['ra'].extend(data['RA'])
            bigdata['dec'].extend(data['DEC'])
            bigdata['image_stamp'].extend(data['IMAGE_STAMP'])

            if real:
                fakemag = data['FLUX']*0. + 99
                bigdata['FakeZPT'].extend(data['FLUX']*0. + 31.)
                bigdata['FakeMag'].extend(fakemag)
            fakeflux = 10 ** (.4 * (31. - fakemag))

            # www = (fakemag < 21.5) & (data['FLUX']-fakeflux < -600.) & (data['FLUX']-fakeflux > -1000.)
            # if len(fakemag[www]) > 0:
            #     #print f
            #
            #     print 'stopped because has a bad outlier'
                #raw_input()

            #for m, faz, fiz in zip(data['MJD'],data['FAKEZPT'], data['ZPT']):
            #    if abs(faz - fiz) > 1:
                    #print f
                    #print m
                    #raw_input('STOPPED')

                    # try:
            #bigdata['rmsaddin'].extend(data['RMSADDIN'])
            #     #print data['RMSADDIN']
            #     #print np.mean(data['RMSADDIN'])
            #     #raw_input()
            # except:
            #     data2 = dt.readcol('./working/lightcurves/'+f.split('/')[-1])
            #     rms = np.mean(data2['RMSADDIN'][data2['RMSADDIN'] > 0.0])
            #     print rms
            #     raw_input()
            #     bigdata['rmsaddin'].extend(data['CHI2']*0. + rms)
            bigdata['field'].extend(data['CHI2']*0 + np.float(deep))
            #print f,'read in'
        # except:
        #     print 'Columns missing in file '+f

        # for sf in data['ZPTFILE']:
        #     zptdata = np.load(sf)
        #     if not sf in zptfiles:
        #         try:
        #             bigdata['starfluxerr'].extend(zptdata['flux_star_std'])
        #             bigdata['starflux'].extend(zptdata['flux_star'])
        #             bigdata['starzpt'].extend(zptdata['fit_zpt'])
        #             bigdata['catmag'].extend(zptdata['cat_mag'])
        #             print 'read in ',sf
        #             zptfiles.append(sf)
        #         except:
        #             print 'Missing flux_star_std'

        fakef = f.split('/')[-1][:17]
        filt = f.split('/')[-1][18]
        fakefile = os.path.join(fakedir,fakef+'.dat')
        ff = open(fakefile,'r').readlines()
        #print 'fileter',filt
        hostmag = -999
        #raw_input()
        filters = np.array(['u','g','r','i','z'],dtype='str')
        arg = np.argwhere(filters == filter)
        for l in ff:
            key = l.split(':')[0]
            if key == 'HOSTGAL_SB_FLUXCAL':
                if filt == filter:
                    if float(l.split()[int(arg)]) <= 0.:
                        hgf = 1.
                    else:
                        hgf = float(l.split()[int(arg)])
                    hostmag = 27.5 - 2.5 * np.log10(hgf)
        #print 'hostmag',hostmag
        #raw_input()
        bigdata['HostMag'].extend(data['FLUX']*0 + hostmag)

        if cntr%200 == 0:
            np.savez(outfile, **bigdata)

        #raw_input()
    #print bigdata['diffzpt']
    #raw_input()
    #print np.unique(np.array(bigdata['FakeMag']))
    #raw_input('here')

    os.system('rm '+cd+' -f')

    print 'saving to cachfile'
    np.savez(outfile,**bigdata)
    print 'saved'
    #tmpwriter.savez(outfile,*bigdata)
    return bigdata


def plotpercentageresid(flux,fluxerr,fakemag,fitzpt,fakezpt,diffimflux,diffimfluxerr,sky,skyerr,dpmjd,chisq,imfiles,ra,dec,imstamp,outdir,fakefiles,hostmag,filter,oldfakezpt,zptstd,filterarr,flag,fwhm,real=False):
    flux = np.asarray(flux)
    fakemag = np.asarray(fakemag,dtype='float')
    dpmjd = np.asarray(dpmjd,dtype='float')
    sky = np.asarray(sky)
    fluxerr = np.asarray(fluxerr)
    fakezpt= np.asarray(fakezpt)
    fitzpt = np.asarray(fitzpt)
    fwhm = np.asarray(fwhm)
    chisq = np.asarray(chisq)
    imfiles = np.asarray(imfiles,dtype='str')
    ra = np.asarray(ra)
    dec = np.asarray(dec)
    imstamp = np.asarray(imstamp)
    fakefiles = np.asarray(fakefiles,dtype='str')
    filterarr = np.asarray(filterarr,dtype='str')

    hostmag = np.asarray(hostmag)
    print np.unique(np.sort(hostmag))[:25]
    print fakezpt
    #raw_input('hostmags')

    if real:
        ww = (dpmjd > 220.) | (dpmjd < -40.)
        flux = flux[ww]
        fluxerr = fluxerr[ww]
        fakemag = fakemag[ww]
        sky = sky[ww]
        chisq = chisq[ww]
        fakezpt = fakezpt[ww]
        fitzpt = fitzpt[ww]
        hostmag = hostmag[ww]
        dpmjd = dpmjd[ww]
        fwhm = fwhm[ww]


    skyerr = np.asarray(skyerr)
    zptstd=np.asarray(zptstd)
    # print hostmag.shape
    # raw_input()
    #print np.unique(fakemag)
    #print np.unique(flux)
    #raw_input('asdf#')

    #for ft,fa in zip(fitzpt,fakezpt):
    #    print ft,fa
    #raw_input('asdffff')
    fluxerr = np.asarray(fluxerr)
    fluxerr = np.sqrt(fluxerr**2.+flux + 10**(.4*(31.-hostmag)))
    fakezpt = np.asarray(fakezpt)



    #dmag  = -2.5*np.log10(diffimflux) + oldfakezpt
    #diffimflux = 10**(.4*(31-dmag))

    oldfakezpt = np.array(oldfakezpt)

    #print fakefiles[(fluxerr<30.) & (fluxerr>0.)]
    #raw_input()
    #print fakefiles[(chisq < 0.1)]
    #raw_input()
    # for fl , ff in zip(fakefiles[(hostmag<21.4) & (fluxerr>0.) & (fakemag > 90.)],flux[(hostmag<21.4) & (fluxerr>0.) & (fakemag > 90.)]):
    #     print ff, fl
    #for fl , ff in zip(fakefiles[(fluxerr>0.) & (fakemag < 90.)],flux[(fluxerr>0.) & (fakemag < 90.)]):
    #    print ff, fl
    #raw_input()
    # print fakezpt
    #
    # print fakemag.shape
    # print flux.shape
    #raw_input()
    #print fakemag[0].shape
    #sys.exit()
    #diffimflux *= 10 ** (.4 * (31. - fakezpt))

    # plt.clf()
    # plt.scatter(diffimflux,flux)
    # plt.plot([min(diffimflux),max(diffimflux)],[min(diffimflux),max(diffimflux)],color='black')
    # plt.xlabel('Diffim Flux')
    # plt.ylabel('SMP Flux')
    #
    # plt.savefig(outdir+'/fluxvs.png')
    #
    # plt.clf()
    # plt.scatter(dpmjd, fluxerr, alpha=.1)
    # plt.xlabel('DPMJD')
    # plt.ylabel('SMP Fluxerr')
    # plt.ylim(50,550)
    # plt.xlim(-40,300)
    # plt.savefig(outdir + '/fluxerrvsDPMJD.png')
    #



    fitzpt = np.asarray(fitzpt)
    fakezpt = np.asarray(fakezpt)

    fakeflux = 10**(.4*(31. - fakemag))


    #diffimflux *= 10**(.4*(31. - oldfakezpt))

    # for fm,ff,fl in zip(fakemag,fakeflux,flux):
    #
    #     print fm,ff,fl
    #raw_input()
    fakefluxo = copy(fakeflux)
    fakeflux *= 10**(-1*.4*(fitzpt - fakezpt))
    fluxerr = (fluxerr**2+flux)**.5

    #fluxerr *= 10**(1*.4*(fitzpt - fakezpt))
    sky *= 10**(.4*(31-fitzpt))

    # plt.clf()
    # plt.scatter(flux, fluxerr, alpha=.1)
    # plt.xlabel('SMP Flux')
    # plt.ylabel('SMP Fluxerr')
    # plt.ylim(50, 550)
    # plt.xlim(-1000, 7000)
    # plt.savefig(outdir + '/fluxerrvsflux.png')
    #
    # plt.clf()
    # plt.scatter(hostmag, fluxerr, alpha=.1)
    # plt.xlabel('Hostmag')
    # plt.ylabel('SMP Fluxerr')
    # plt.ylim(50, 550)
    # plt.xlim(21, 30)
    # plt.savefig(outdir + '/fluxerrvshostmag.png')
    #
    # plt.clf()
    # plt.scatter(sky, fluxerr, alpha=.1)
    # plt.xlabel('Sky')
    # plt.ylabel('SMP Fluxerr')
    # plt.ylim(50, 550)
    # plt.xlim(-2000, 10000)
    # plt.savefig(outdir + '/fluxerrvssky.png')
    #
    # plt.clf()
    # plt.scatter(fitzpt, fluxerr, alpha=.1)
    # plt.xlabel('SMP ZPT')
    # plt.ylabel('SMP Fluxerr')
    # plt.ylim(50, 550)
    # plt.xlim(31.5, 32.5)
    # plt.savefig(outdir + '/fluxerrvszpt.png')
    #
    # ww = fakemag == 99
    # plt.clf()
    # plt.scatter(skyerr[ww], fluxerr[ww], alpha=.1)
    # plt.xlabel('SMP Skyerr')
    # plt.ylabel('SMP Fluxerr')
    # plt.ylim(50, 550)
    # plt.xlim(0, 200)
    # plt.savefig(outdir + '/fluxerrvsskyerr.png')
    #
    # plt.clf()
    # plt.hist(skyerr, bins=np.arange(0,100,2.),normed=True)
    # plt.xlabel('SMP Skyerr')
    # plt.xlim(0, 100)
    # plt.savefig(outdir + '/skyerrhist.png')


    d = (flux - fakeflux) / (fluxerr-np.sqrt(10.**(.4*(31.-hostmag))))
    ww = (flux != 0.)
    d = d[ww]
    fm = np.array(fakemag,dtype='float')[ww]
    dc99 = d[fm > 90.]
    rms99 = np.sqrt(np.nanmean(np.square(dc99[abs(dc99) < 10.])))

    ww = (flux != 0.) & (fakemag != 0) & (fakemag < 28.5) & (flux != 0.) & (abs((flux-fakeflux)/fluxerr) < 10.)
    plt.clf()
    fig = plt.figure(figsize=(15, 10))


    # plt.scatter(fakemag[ww], (-2.5*np.log10(diffimflux[ww]) + oldfakezpt[ww] - fakemag[ww]), alpha=.1, color='red')
    # ax, ay, aystd = bindata(fakemag[ww],(-2.5*np.log10(diffimflux[ww]) + oldfakezpt[ww] - fakemag[ww]),
    #                         np.arange(min(fakemag[ww]), max(fakemag[ww]), .5))
    # plt.errorbar(ax, ay, aystd, markersize=10, color='red', fmt='o', label='Diffim',alpha=.4)

    #print fakemag[ww].shape,flux[ww].shape,fakeflux[ww].shape

    #for x in np.unique(fakeflux):
    #    print x
    #print np.unique(fakeflux)
    #raw_input('fakeflux')

    # plt.scatter(fakemag[ww], (diffimflux[ww] - fakeflux[ww]) / fakeflux[ww], alpha=.2,color='red')
    # ax, ay, aystd = bindata(fakemag[ww], (diffimflux[ww] - fakeflux[ww]) / fakeflux[ww],
    #                         np.arange(19, 28, .5))
    # plt.errorbar(ax, ay, aystd, markersize=15, color='red', fmt='o', label='DIFFIMG')

    plt.scatter(fakemag[ww],(flux[ww]-fakeflux[ww])/fakeflux[ww],alpha=.2,color='green')
    ax, ay, aystd = bindata(fakemag[ww],(flux[ww]-fakeflux[ww])/fakeflux[ww],
                            np.arange(19,28, .5))
    plt.errorbar(ax, ay, aystd, markersize=15, color='green', fmt='o', label='SMP')



    plt.axhline(0)
    plt.xlim(18,29)
    #plt.ylim(-.1,.1)
    plt.ylim(-.099,.099)
    plt.xlabel('Fake Mag',fontsize=30.)
    plt.ylabel('Fractional Flux Difference',fontsize=30.)
    plt.title(filter+' band',fontsize=30.)
    plt.savefig(outdir+'/percentagefluxdiff.png')

    plt.clf()

    plt.scatter(fwhm[ww]/2.235,(flux[ww]-fakeflux[ww])/fluxerr[ww],alpha=.2,color='green')
    ax, ay, aystd = bindata(fwhm[ww]/2.235,(flux[ww]-fakeflux[ww])/fluxerr[ww],
                            np.arange(.5,3, .2))
    plt.errorbar(ax, ay, aystd, markersize=15, color='green', fmt='o', label='SMP')



    plt.axhline(0)
    plt.xlim(.5,3.)
    #plt.ylim(-.1,.1)
    plt.ylim(-6,6)
    plt.xlabel('FWHM (arcsec)',fontsize=30.)
    plt.ylabel('Fit-Fake/Err',fontsize=30.)
    plt.title(filter+' band',fontsize=30.)
    plt.savefig(outdir+'/stdfluxdifffwhm.png')

    plt.clf()

    plt.hist(fwhm[ww]/2.235,bins=100)
    plt.xlabel('PSF_FWHM')
    plt.savefig(outdir+'/fwhmhist.png')
    print 'saved',outdir+'/fwhmhist.png'
    #raw_input()
    plt.clf()

    ax,ay = dt.binrms(fakemag[ww],(flux[ww]-fakeflux[ww])/fluxerr[ww],np.arange(19,28, .1),.5)
    plt.scatter(fakemag[ww],(flux[ww]-fakeflux[ww])/fluxerr[ww],alpha=.2,color='green')
    ax, ay, aystd,n = bindata(fakemag[ww],(flux[ww]-fakeflux[ww])/fluxerr[ww],
                            np.arange(19,28, .5),returnn=True)
    plt.errorbar(ax, ay,yerr=aystd*n**.5, color='green', label='SMP')



    plt.axhline(1,linestyle='--')
    plt.xlim(18,29)
    #plt.ylim(.6,1.4)
    plt.ylim(-2.2,2.2)
    plt.axhline(-1)
    plt.axhline(1)
    plt.xlabel('Fake Mag',fontsize=30.)
    plt.ylabel('RMS',fontsize=30.)
    plt.title(filter+' band',fontsize=30.)
    plt.savefig(outdir+'/rms.png')

    print 'saved', outdir+'/rms.png'
    # print min(hostmag[ww]),max(hostmag[ww])
    plt.clf()
    #raw_input()

    try:
        fig = plt.figure(figsize=(15, 10))
        plt.scatter(hostmag[ww],(flux[ww]-fakeflux[ww]),alpha=.5)
        ax, ay, aystd = bindata(hostmag[ww],(flux[ww]-fakeflux[ww]),
                                np.arange(min(hostmag[ww]),max(hostmag[ww]), .2))
        plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')

        plt.axhline(0)
        plt.xlim(19,29)
        #plt.ylim(-.1,.1)
        plt.ylim(-500,500)
        plt.xlabel('Host Mag')
        plt.ylabel('Flux Difference')
        plt.title(filter + ' band')

        plt.savefig(outdir+'/fluxdiff_hostmag.png')
    except:
        print 'bad hostmags'

    ww = (flux != 0.) & (fakemag != 0) & (fakemag < 25.)


    try:
        plt.clf()
        fig = plt.figure(figsize=(15, 10))
        plt.scatter(hostmag[ww], (flux[ww] - fakeflux[ww]), alpha=.5)
        ax, ay, aystd = bindata(hostmag[ww], (flux[ww] - fakeflux[ww]),
                                np.arange(min(hostmag[ww]), max(hostmag[ww]), .2))
        plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')

        plt.axhline(0)
        plt.xlim(19, 29)
        # plt.ylim(-.1,.1)
        plt.ylim(-1000, 1000)
        plt.xlabel('Host Mag')
        plt.ylabel('Flux Difference')
        plt.title(filter + ' band')

        plt.savefig(outdir + '/fluxdiffgt23_hostmag.png')
    except:
        print 'bad hostmags'

    #print fakefiles[ww][((flux[ww]-fakeflux[ww])/fakeflux[ww] < -.04) & ((flux[ww]-fakeflux[ww])/fakeflux[ww] > -.4) & (fakemag[ww]<22.)]
    #raw_input()
    #for k in np.unique(fakefiles):
    #    print k
    #raw_input()

    ww = (flux != 0.) & (np.array(fakemag,dtype='float') > 90.)  # (fakemag < 28.5) & (flux != 0.)
    plt.clf()
    # fig = plt.figure(figsize=(15, 10))
    # plt.hist(flux[ww],bins=np.arange(-650,600,100))
    # #ax, ay, aystd = bindata(fakeflux[ww], (flux[ww] - fakeflux[ww]),
    # #                        np.arange(-100, 1000, 200))
    # #plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')
    #
    # #plt.axhline(0)
    # plt.xlim(-500, 500)
    # # plt.ylim(-.1,.1)
    # #plt.ylim(-600, 600)
    # plt.xlabel('Fake Flux')
    # plt.ylabel('Flux Difference ')
    # plt.savefig(outdir + '/efluxdiff.png')

    # plt.clf()
    # fig = plt.figure(figsize=(15, 10))
    # plt.hist(flux[ww]/fluxerr[ww], bins=np.arange(-6.2, 6, .4),normed=True,label='RMS Fakemag = 99: ' + str(round(rms99, 3)))
    # # ax, ay, aystd = bindata(fakeflux[ww], (flux[ww] - fakeflux[ww]),
    # #                        np.arange(-100, 1000, 200))
    # # plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')
    #
    # # plt.axhline(0)
    # plt.xlim(-5, 5)
    #
    # import matplotlib.mlab as mlab
    # import math
    # mean = 0
    # variance = 1
    # sigma = math.sqrt(variance)
    # x = np.arange(-5, 5, .1)
    # plt.plot(x, mlab.normpdf(x, mean, sigma), color='black', label='Gaussian Normal')
    # plt.legend()
    # # plt.ylim(-.1,.1)
    # # plt.ylim(-600, 600)
    # plt.xlabel('flux/fluxerr')
    # plt.ylabel('Count')
    # plt.title(filter+' band')
    #
    #
    # plt.savefig(outdir + '/efluxdiffstd.png')

    ww = (flux != 0.) & (fakemag != 0)  # (fakemag < 28.5) & (flux != 0.)
    plt.clf()
    fig = plt.figure(figsize=(15, 10))
    plt.scatter(fakeflux[ww], (flux[ww] - fakeflux[ww]), alpha=.5)
    ax, ay, aystd = bindata(fakeflux[ww], (flux[ww] - fakeflux[ww]),
                            np.arange(-100,5000, 200))
    plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')

    plt.axhline(0)
    plt.xlim(0, 5000)
    # plt.ylim(-.1,.1)
    plt.ylim(-2000, 600)
    plt.xlabel('Fake Flux')
    plt.ylabel('Fit - Fake Flux')
    plt.title(filter+' band')

    plt.savefig(outdir + '/fluxdiff.png')

    plt.clf()
    fig = plt.figure(figsize=(15, 10))
    try:
        plt.hist((flux[ww] - fakeflux[ww])/fakeflux[ww],bins=np.arange(-.51,.5,.02), histtype='step',color='blue',label='SMP STD '+
                                    str(round(1.48 * np.median(abs((flux[ww] - fakeflux[ww])/fakeflux[ww] - np.median((flux[ww] - fakeflux[ww])/fakeflux[ww]))), 3)))
    except:
        print 'empty'
    # plt.hist((diffimflux[ww] - fakeflux[ww])/fakeflux[ww],bins=np.arange(-.51,.5,.02), histtype='step',color='red',label='DIFFIM STD '+
    #                                 str(round(1.48 * np.median(abs((diffimflux[ww] - fakeflux[ww])/fakeflux[ww] - np.median((diffimflux[ww] - fakeflux[ww])/fakeflux[ww]))), 3)))


    plt.axhline(0)
    plt.xlim(-.5, .5)
    # plt.ylim(-.1,.1)
    #plt.ylim(-600, 600)
    #plt.xlabel('Fake Flux')
    #plt.ylabel('Flux Difference ')
    plt.legend()
    plt.title(filter + ' band')

    plt.savefig(outdir + '/pfluxdiffhist.png')



    ww = (flux != 0.) & (fakemag != 0)  # (fakemag < 28.5) & (flux != 0.)
    plt.clf()
    # fig = plt.figure(figsize=(15, 10))
    # plt.scatter(fakeflux[ww], (flux[ww] - fakeflux[ww])/fakeflux[ww], alpha=.5)
    # ax, ay, aystd = bindata(fakeflux[ww], (flux[ww] - fakeflux[ww])/fakeflux[ww],
    #                         np.arange(-100,5000, 200))
    # plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')
    #
    # plt.axhline(0)
    # plt.xlim(0, 5000)
    # # plt.ylim(-.1,.1)
    # plt.ylim(-.1, .1)
    # plt.xlabel('Fake Flux')
    # plt.ylabel('Flux Difference ')
    # plt.savefig(outdir + '/pfluxdiff.png')


    skyresid = 2000. - sky / 10 ** (.4 * (31 - fitzpt))
    # plt.clf()
    # fig = plt.figure(figsize=(15,10))
    # plt.scatter(skyresid[ww],(flux[ww]-fakeflux[ww]),alpha=.15)
    # ax, ay, aystd = bindata(skyresid[ww], (flux[ww] - fakeflux[ww]),
    #                         np.arange(-100,100,3.))
    # plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')
    #
    # plt.axhline(0)
    # plt.xlim(-25,25)
    # # plt.ylim(-.1,.1)
    # plt.ylim(-200, 200)
    # plt.xlabel('Sky Resid')
    # plt.ylabel('Flux Difference')
    # plt.savefig(outdir + '/fluxdiff_sky.png')

    plt.clf()
    # fig = plt.figure(figsize=(15, 10))
    # plt.scatter(skyresid[ww], (flux[ww] - fakeflux[ww])/fakeflux[ww], alpha=.15)
    # ax, ay, aystd = bindata(skyresid[ww], (flux[ww] - fakeflux[ww])/fakeflux[ww],
    #                         np.arange(-100, 100, 3.))
    # plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')
    #
    # plt.axhline(0)
    # plt.xlim(-25, 25)
    # # plt.ylim(-.1,.1)
    # plt.ylim(-.2, .2)
    # plt.xlabel('Sky Resid')
    # plt.ylabel('Percentage Flux Difference')
    # plt.savefig(outdir + '/pfluxdiff_sky.png')
    #
    # plt.clf()
    # fig = plt.figure(figsize=(15, 10))
    # plt.scatter(fakemag[ww], flux[ww], alpha=.5)
    # ax, ay, aystd = bindata(fakemag[ww], flux[ww] ,
    #                        np.arange(min(fakemag[ww]), max(fakemag[ww]), .15))
    # plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')
    #
    # plt.scatter(fakemag[ww],fakeflux[ww],alpha=.5,color='red')
    #
    # plt.axhline(0)
    # plt.xlim(22, 29)
    # # plt.ylim(-.1,.1)
    # plt.ylim(-500, 5000)
    # plt.xlabel('Fake Mag')
    # plt.ylabel('Fit Flux')
    # plt.savefig(outdir + '/fvf.png')

    #print imfiles[ww][(fakemag[ww]<21) & ((flux[ww]-fakeflux[ww])/fakeflux[ww] < -.9)]
    #print ra[ww][(fakemag[ww]<21) & ((flux[ww]-fakeflux[ww])/fakeflux[ww] < -.9)]
    #print dec[ww][(fakemag[ww]<21) & ((flux[ww]-fakeflux[ww])/fakeflux[ww] < -.9)]
    #print imstamp[ww][(fakemag[ww]<21) & ((flux[ww]-fakeflux[ww])/fakeflux[ww] < -.9)]

    #raw_input('imfiles bad')

    plt.clf()
    fig = plt.figure(figsize=(15, 10))

    ww = (flux != 0.) & (fakemag != 0) & ((flux - fakeflux) / fakeflux > -10.) & ((flux - fakeflux) / fakeflux < 10.)

    plt.scatter(chisq[ww], (flux[ww] - fakeflux[ww]) / fakeflux[ww], alpha=.5)
    ax, ay, aystd = bindata(chisq[ww], (flux[ww] - fakeflux[ww]) / fakeflux[ww],
                            np.arange(0,10, .05))
    plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')

    plt.axhline(0)
    plt.xlim(.6, 1.2)
    plt.ylim(-.2, .2)
    plt.xlabel('Chi Sq')
    plt.ylabel('Percentage Flux Difference')
    plt.title(filter+' band')

    plt.savefig(outdir + '/percentagefluxdiffchi.png')





    ww = (flux != 0.) & (fakemag != 0) & (fakemag ==99.)

    # plt.clf()
    # plt.scatter(sky[ww],(flux[ww]-fakeflux[ww]),alpha=.1)
    # ax, ay, aystd = dt.bindata(sky[ww],(flux[ww]-fakeflux[ww]),
    #                         np.arange(-10,10, 1),window=2.)
    # plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')
    # plt.axhline(0)
    # plt.xlim(-10,10)
    # plt.ylim(-2000,2000)
    # plt.xlabel('sky')
    # plt.title('Fakemag < 22')
    # plt.ylabel('Flux Difference')
    # plt.savefig(outdir + '/skyfluxdifflt22.png')

    ww = (flux != 0.) & (fakemag != 0) & (fakemag == 99.)

    # plt.clf()
    # plt.scatter(sky[ww],(flux[ww]-fakeflux[ww]),alpha=.1)
    # ax, ay, aystd = dt.bindata(sky[ww],(flux[ww]-fakeflux[ww]),
    #                         np.arange(-10,10, 1),window=2.)
    # plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')
    # plt.axhline(0)
    # plt.xlim(-10,10)
    # plt.ylim(-2000,2000)
    # plt.xlabel('sky')
    # plt.title('Fakemag > 22')
    # plt.ylabel('Flux Difference')
    # plt.savefig(outdir + '/skyfluxdiffgt22.png')

    ww = (dpmjd > 320.) & (flux != 0)
    # plt.clf()
    # plt.hist(flux[ww],bins=np.arange(-1050,1000,100))
    # plt.xlabel('Flux for epochs without fake SN Flux')
    # plt.xlim(-1000,1000)
    # plt.savefig(outdir + '/emptyflux.png')

    import matplotlib.mlab as mlab
    import math
    mean = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.arange(-5, 5, .1)

    ww = (fakemag == 99) & (flux != 0)
    plt.clf()
    # plt.hist(flux[ww]/fluxerr[ww], bins=np.arange(-4.1, 4, .2),normed=True)
    # plt.xlim(-4, 4)
    # plt.plot(x,mlab.normpdf(x, mean, sigma), color='black', label='Gaussian Normal')
    # plt.xlabel('flux/fluxerr (epochs without fake SN Flux)')
    # plt.savefig(outdir + '/emptyfluxstd.png')

    plt.clf()

    #try:
    if True:
        plt.scatter(hostmag[ww],flux[ww] / fluxerr[ww])
        plt.xlim(19, 29)
        plt.ylim(-3,3)
        ax, ay, aystd = dt.bindata(hostmag[ww], (flux[ww] / fluxerr[ww]),
                                   np.arange(19, 28, .2))
        plt.errorbar(ax, ay, aystd, markersize=10, color='gold', fmt='o', label='SMP: 3 sigma %.4f \n     '
                                                                                '5 sigma %.4f'
                                                                                ''%(float(len(flux[ww][abs(flux[ww]/fluxerr[ww]   > 3.)]))/len(flux[ww]),
                                                                                    float(len(flux[ww][abs(flux[ww]/fluxerr[ww]  > 5.)]))/len(flux[ww])))
        plt.xlabel('Hostmag')
        plt.ylabel('Flux/Fluxerr (epochs without fake SN Flux)')
        plt.legend()
        plt.axhline(0)
        plt.title(filter + ' band')

        plt.savefig(outdir + '/emptyfluxstdvshostmag.png')

        plt.clf()
        plt.scatter(fwhm[ww], flux[ww] / fluxerr[ww])
        plt.xlim(1.5, 4.)
        plt.ylim(-6, 6)
        ax, ay, aystd = dt.bindata(fwhm[ww], (flux[ww] / fluxerr[ww]),
                                   np.arange(1, 4.5, .2))
        plt.errorbar(ax, ay, aystd, markersize=10, color='gold', fmt='o', label='SMP: 3 sigma %.4f \n     '
                                                                                '5 sigma %.4f'
                                                                                '' % (float(
            len(flux[ww][abs(flux[ww] / fluxerr[ww] > 3.)])) / len(flux[ww]),
                                                                                      float(len(flux[ww][abs(
                                                                                          flux[ww] / fluxerr[
                                                                                              ww] > 5.)])) / len(
                                                                                          flux[ww])))
        plt.xlabel('Hostmag')
        plt.ylabel('Flux/Fluxerr (epochs without fake SN Flux)')
        plt.legend()
        plt.axhline(0)
        plt.title(filter + ' band')

        plt.savefig(outdir + '/emptyfluxstdvsfwhm.png')

        plt.clf()
        plt.scatter(hostmag[ww], flux[ww])
        plt.xlim(21, 30)
        plt.ylim(-700, 700)
        ax, ay, aystd = dt.bindata(hostmag[ww], (flux[ww]),
                                   np.arange(19, 30, .2))
        plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')
        plt.xlabel('Hostmag')
        plt.ylabel('Flux (epochs without fake SN Flux)')
        plt.axhline(0)
        plt.title(filter + ' band')

        plt.savefig(outdir + '/emptyfluxvshostmag.png')
    #except:
    #    print 'bad hostmags'
    # plt.clf()
    # plt.scatter(fakemag[ww], flux[ww])
    # plt.xlim(19, 28)
    # plt.ylim(-2500, 2500)
    # ax, ay, aystd = dt.bindata(fakemag[ww], (flux[ww]),
    #                            np.arange(19, 28, .1))
    # plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')
    # plt.xlabel('Hostmag')
    # plt.ylabel('Flux (epochs without fake SN Flux)')
    # plt.axhline(0)
    # plt.savefig(outdir + '/emptyfluxvsfm.png')

    print 'saved png'

def plotsigmaresid(flux,fluxerr,fakemag,fitzpt,fakezpt,hostmag,chisqarr,rmsaddin,deep,outdir,zptstd,diffimflux,diffimfluxerr,filter,filterarr,deep_or_shallow,sky,skyerr,fwhm,imfiles,dpmjd,flag,snid,mjd,real=False):



    flux = np.asarray(flux)
    fakemag = np.asarray(fakemag)
    #print len(fakemag[fakemag<26.])
    #raw_input('fffm111')
    fluxerr = np.asarray(fluxerr)
    hostmag = np.asarray(hostmag)
    dpmjd = np.asarray(dpmjd)
    #fluxerr = (fluxerr**2+flux+10**(.4*(31.-hostmag)))**.5
    fluxadd = copy(flux)
    fluxadd[flux<0] = 0.
    fluxerr = (fluxerr**2+abs(flux)+10**(.4*(31.-hostmag)))**.5
    fitzpt = np.asarray(fitzpt)
    fakezpt = np.asarray(fakezpt)
    rmsaddin = np.asarray(rmsaddin)
    #print rmsaddin[:500]
    #raw_input('rmsaddin')
    fakeflux = 10 ** (.4 * (31. - fakemag))
    fakeflux *= 10**(-1*.4*(fitzpt - fakezpt))
    fluxerrz = (fluxerr**2 + abs(flux) + 10**(.4*(31.-hostmag))+ (rmsaddin*flux)**2)**.5
    #fluxerr = (fluxerr**2 + abs(flux) + 10**(.4*(31.-hostmag))+ (rmsaddin*flux)**2)**.5

    #fluxerrz = fluxerr
    diffimgflux = np.array(diffimflux)
    diffimgfluxerr = np.array(diffimfluxerr)
    diffimflux = 10 ** (.4 * (31 - 27.5)) * np.array(diffimflux)
    diffimfluxerr = 10 ** (.4 * (31 - 27.5)) * np.array(diffimfluxerr)
    sky = np.array(sky)
    skyerr = np.array(skyerr)
    filterarr = np.array(filterarr,dtype='str')
    imfiles = np.array(imfiles,dtype='str')
    fwhm = np.array(fwhm)
    flag = np.array(flag)
    zptstd = np.array(zptstd)
    snid = np.array(snid)

    sky = sky*10**(.4*(fitzpt-31.))-10000
    skyerr = skyerr*10**(-.4*(fitzpt-31.))

    #print sky[:1000]
    #print flux[:100]
    #raw_input()

    #fluxerr *= 10**(.4*(fitzpt - fakezpt))

    #fakeflux *= 10 ** (-1 * .4 * (fitzpt - fakezpt))
    chisqarr = np.asarray(chisqarr)

    mjd = np.array(mjd)
    if real:
        ww = (dpmjd > 250.) | (dpmjd < -40.)
        flux = flux[ww]
        fluxerr = fluxerr[ww]
        fakemag = fakemag[ww]
        sky = sky[ww]
        chisqarr = chisqarr[ww]
        fakezpt = fakezpt[ww]
        fitzpt = fitzpt[ww]
        hostmag = hostmag[ww]
        skyerr = skyerr[ww]
        fluxerrz = fluxerrz[ww]
        filterarr = filterarr[ww]
        diffimflux = diffimflux[ww]
        diffimfluxerr = diffimfluxerr[ww]
        rmsaddin = rmsaddin[ww]
        dpmjd = dpmjd[ww]
        fakeflux = fakeflux[ww]
        fwhm = fwhm[ww]
        #flag = flag[ww]
        #zptstd=zptstd[ww]
    # dww = (dpmjd < 250.) | (dpmjd > -50.)
    # flag = flag[dww]
    # fwhm = fwhm[dww]
    # imfiles = imfiles[dww]
    # snid = snid[dww]
    # bfilt = filterarr[dww]
    # mjd = mjd[dww]
    #
    # plt.clf()
    # fig = plt.figure(figsize=(15, 10))
    # plt.hist(fwhm[flag>0],bins=np.arange(0,8,.1))
    # plt.xlim(0,8)
    # plt.xlabel('FWHM[flag>0]')
    # plt.savefig('fwhmflag.png')
    # plt.clf()
    # fig = plt.figure(figsize=(15, 10))
    # plt.hist(sky[flag > 0], bins=50)
    # plt.xlabel('Sky[flag>0]')
    # plt.savefig('skyflag.png')
    # plt.clf()
    # #fig = plt.figure(figsize=(15, 10))
    # #plt.hist(flag[flag > 0],bins=np.arange(0,66000,2))
    # #plt.xlabel('FLAGS')
    # #plt.savefig('flagflag.png')
    # #plt.clf()
    # print 'FLAGS'
    # for fl in np.unique(flag):
    #     print fl,len(flag[flag==fl])
    # for mj,sn in zip(mjd[flag==1],snid[flag==1]):
    #     print mj,sn
    # print '-'*15
    # #print np.unique(bfilt)
    # for sn in np.unique(snid):
    #     print len(flag[(flag == 4096) & (snid == sn) & (bfilt=='g') ]), len(flag[(flag == 4096) & (snid == sn) & (bfilt=='r') ]), len(flag[(flag == 4096) & (snid == sn) & (bfilt=='i') ]), len(flag[(flag == 4096) & (snid == sn) & (bfilt=='z') ]), sn
    #     #print len(flag[(flag==4096) & (snid==sn) ]),sn
    # print len(np.unique(snid[flag==4096]))
    # print '-'*15

    # fig = plt.figure(figsize=(15, 10))
    # plt.hist(zptstd[flag > 0], bins=50)
    # plt.xlim(0,8)
    # plt.xlabel('Zptstd[flag>0]')
    # plt.savefig('zptstdflag.png')
    # plt.clf()

    #print np.sqrt(10**(.4*(fitzpt - hostmag))/3.)
    #am = np.argmax(np.sqrt(10**(.4*(fitzpt - hostmag))/3.))
    #for a,f,fe in zip(np.sqrt(10**(.4*(fitzpt - hostmag))/3.),flux,fluxerr):
    #    print a,f,fe
    #print ,flux[am],fluxerr[am]
    #raw_input()

    #fluxerr = np.sqrt(np.asarray(fluxerr)**2+(abs(flux)/3.) + 10**(.4*(fitzpt - hostmag))/3.)
    #fitmag = 31-2.5*np.log10(flux)


    #MY CALCULATION TO GET FLUX RMS ADDIN
    #.1 = -2.5*np.log10(flux) + 2.5*np.log10(flux+rmsfluxerr)
    #rmsaddin + 2.5*log(flux) = 2.5*np.log10(flux+rmsfluxerr)
    #(rmsaddin + 2.5*log(flux))/2.5 = np.log10(flux+rmsfluxerr)
    #10**(.4*(rmsaddin + 2.5*log(flux))) = flux + rmsfluxerr

    #fluxrmsaddin = 10 ** (.4 * (rmsaddin + 2.5 * np.log10(flux))) - flux
    #fluxerr = (fluxerr**2 + fluxrmsaddin**2)**.5

    #frms = 10**(.4*(2.5*np.log10(abs(flux))+rmsaddin)) - abs(flux)

    #print frms[0:100]
    #print flux[0:100]*rmsaddin[0:100]
    #raw_input()

    fitmag = 31-2.5*np.log10(flux)
    fitmagerr = ((-2.5*np.log10(flux)+2.5*np.log10(flux+fluxerr))**2)**.5
    fakemag = fakemag + fakezpt - fitzpt

    hostmag = np.array(hostmag)

    fifx = 10**(.4*(31-fitmag))
    fafx = 10**(.4*(31-fakemag))
    fime = 10**(.4*(31-fitmag+fitmagerr)) - 10**(.4*(31-fitmag))

    #d = (fifx-fafx)/fime
    #d = (fitmag - fakemag)/(fitmagerr*1.08)


    ff = copy(fakeflux)
    # ff[ff < 1.] = 1.
    # plt.hist(fluxerr/ff,bins=20,normed=True)
    # plt.title(filter+' band')
    #
    # plt.savefig('sinosn.png')
    plt.clf()
    #np.savez('simnosn.npz',flux=flux,fakeflux=ff,fluxerr=np.sqrt(fluxerr**2 + abs(flux)/3.8))


    #for u in range(len(flux)):
    #    print fluxerr[u],np.sqrt(10.**(.4*(31.-hostmag[u]))),np.sqrt(fluxerr[u]**2-10.**(.4*(31.-hostmag[u])))
    #raw_input()

    d = (flux - fakeflux) / (fluxerr)
    #print d[:100]
    #print 'nnn'
    #dz = (flux - fakeflux) / ((fluxerrz**2 )**.5)
    df = (diffimflux - fakeflux) / ((diffimfluxerr**2 )**.5)
    ww = (flux != 0.) & (np.array(fakemag, dtype='float') > 0.) & (fluxerr > 0.) & (abs((flux-fakeflux)/fluxerr) < 4.)

    #fakemag[fakemag==99] = 29.5
    flux = flux[ww]
    fakemag = fakemag[ww]
    fitzpt = fitzpt[ww]
    fakezpt = fakezpt[ww]
    fakeflux= fakeflux[ww]
    fluxerr=fluxerr[ww]
    fluxerrz = fluxerrz[ww]
    diffimflux = diffimflux[ww]
    diffimfluxerr = diffimfluxerr[ww]
    sky = sky[ww]
    skyerr = skyerr[ww]
    d = d[ww]
    fwhm = fwhm[ww]
    diffimd = diffimflux/diffimfluxerr
    #print d[:100]
    #raw_input()
    #dz = dz[ww]
    #df = df[ww]
    hostmag = hostmag[ww]
    chisqarr = chisqarr[ww]
    filterarr = filterarr[ww]

    #print flux[0:10]
    #print fakeflux[0:10]
    #print flux.shape
    #print fakeflux.shape
    #d = (flux - fakeflux) / fluxerr

    chisq = (flux - fakeflux) ** 2 / fluxerr ** 2
    chisq = np.nanmean(chisq[abs(d) < 3])

    plt.clf()

    dc = d[abs(d) < 3.]

    dc99 = d[(np.array(fakemag,dtype='float') > 90.) & (chisqarr > .05) & (chisqarr < 2.5)]
    rms99 = np.sqrt(np.nanmean(np.square(dc99[abs(dc99) < 3.])))

    diffimrms = np.sqrt(np.nanmean(np.square(diffimd[(abs(diffimd) < 3.) & (chisqarr < 2.5) & (np.array(fakemag,dtype='float') > 90.) ])))

    rms992 = np.std(dc99[abs(dc99) < 3.])
    print rms99,rms992
    #raw_input()

    dcr = d[fakemag < 28.]
    rmsr = np.sqrt(np.nanmean(np.square(dcr[abs(dcr) < 3.])))

    #f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

    # fig = plt.figure(figsize=(16, 12))
    # gs = gridspec.GridSpec(1, 3, width_ratios=[4, 1])
    # ax1 = plt.subplot(gs[0])
    # ax2 = plt.subplot(gs[1])

    ww = (flux != 0.) & (np.array(fakemag,dtype='float') > 90.) & (fluxerr >0.) & (np.isfinite(flux)) & \
         (np.isfinite(fluxerr)) &(~np.isnan(flux)) &(~np.isnan(fluxerr)) & (chisqarr < 2.5) & (chisqarr > .05)

    #ww = (flux != 0) & (fakeflux < 1.)
    #print rms99
    fig = plt.figure(figsize=(15, 10))
    plt.hist(flux[ww] / fluxerr[ww], bins=np.arange(-6.1, 6., .2), normed=True,
             label='SMP RMS Fakemag = 99: ' + str(round(rms99, 3)),color='blue',alpha=.4)
    plt.hist(diffimflux[ww] / diffimfluxerr[ww], bins=np.arange(-6.1, 6., .2), normed=True,
             label='DIFFIM RMS Fakemag = 99: ' + str(round(diffimrms, 3)),color='red',alpha=.4)
    # ax, ay, aystd = bindata(fakeflux[ww], (flux[ww] - fakeflux[ww]),
    #                        np.arange(-100, 1000, 200))
    # plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')

    # plt.axhline(0)
    plt.xlim(-5, 5)

    import matplotlib.mlab as mlab
    import math
    mean = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.arange(-5, 5, .1)
    plt.plot(x, mlab.normpdf(x, mean, sigma), color='black', label='Gaussian Normal')
    plt.legend()
    # plt.ylim(-.1,.1)
    # plt.ylim(-600, 600)
    plt.xlabel('flux/fluxerr')
    plt.ylabel('Count')
    #plt.title(filter + ' band')

    plt.savefig(outdir + '/efluxdiffstd.png')

    plt.clf()
    ww = (flux != 0) & (np.array(fakemag, dtype='float') > 90.) & (fluxerr > 0.) & (np.isfinite(flux)) & \
         (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .05) & (chisqarr < 2.5)

    # ww = (flux != 0) & (fakeflux < 1.)
    print rms99
    fig = plt.figure(figsize=(15, 10))
    plt.hist(flux[ww], bins=np.arange(-1025.,1000., 50.), normed=True,color='blue',histtype='step',
             label='SMP STD Fakemag = 99: ' + str(round(1.48 * np.median(abs(flux[ww] - np.median(flux[ww]))), 3)))
    plt.hist(diffimflux[ww], bins=np.arange(-1025,1000., 50.), normed=True,color='red',histtype='step',
             label='DIFFIM STD Fakemag = 99: ' + str(round(1.48 * np.median(abs(diffimflux[ww] - np.median(diffimflux[ww]))), 3)))
    plt.legend()
    # plt.ylim(-.1,.1)
    # plt.ylim(-600, 600)
    plt.xlabel('flux')
    plt.ylabel('Count')
    #print filter.shape
    plt.title(filter + ' band')
    plt.xlim(-1000, 1000)

    plt.savefig(outdir + '/efluxdiff.png')

    plt.clf()
    ww = (flux != 0) & (np.array(fakemag, dtype='float') > 90.) & (fluxerr > 0.) & (np.isfinite(flux)) & \
         (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .05) & (chisqarr < 2.5)

    # ww = (flux != 0) & (fakeflux < 1.)
    print rms99
    fig = plt.figure(figsize=(15, 10))
    plt.hist(fluxerr[ww], bins=np.arange(-102.5, 400., 5.), normed=True,
             label='Mean Err Fakemag = 99: ' + str(round(np.mean(fluxerr[ww]), 3)))
    plt.legend()
    # plt.ylim(-.1,.1)
    # plt.ylim(-600, 600)
    plt.xlabel('fluxerr')
    plt.ylabel('Count')
    plt.xlim(0, 400)

    plt.title(filter + ' band')

    plt.savefig(outdir + '/efluxerrdiff.png')

    # for fmm,fffl,fffle,fafl in zip(fakemag[(flux != 0.) & ((np.abs(flux-fakeflux) / fluxerr) < .1) & (np.array(fakemag, dtype='float') < 27.)& (np.array(fakemag, dtype='float') > 0.)],
    #                     flux[(flux != 0.) & ((np.abs(flux-fakeflux) / fluxerr) < .1) & (np.array(fakemag, dtype='float') < 27.)& (np.array(fakemag, dtype='float') > 0.)],
    #                          fluxerr[(flux != 0.) & ((np.abs(flux - fakeflux) / fluxerr) < .1) & (
    #                              np.array(fakemag, dtype='float') < 27.) & (np.array(fakemag, dtype='float') > 0.)],
    #                     fakeflux[(flux != 0.) & ((np.abs(flux-fakeflux) / fluxerr) < .1) & (np.array(fakemag, dtype='float') < 27.)& (np.array(fakemag, dtype='float') > 0.)]):
    #
    #     print fmm,fffl,fafl, fffl-fafl, fffle
    #raw_input()
    ww = (flux != 0) & (np.array(fakemag, dtype='float') > 0.) & (np.array(fakemag, dtype='float') < 30.) & (fluxerr > 0.) & (np.isfinite(flux)) & \
         (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .05) & (chisqarr < 2.5) #& \
         #(abs(flux - fakeflux) > 0.1)

    #ww = (flux != 0) & (fakeflux < 1.)
    print rms99
    fig = plt.figure(figsize=(15, 10))
    plt.hist((flux[ww]-fakeflux[ww]) / fluxerr[ww], bins=np.arange(-6.1, 6., .2), normed=True,
             label='RMS Fakemag < 99: ' + str(round(rmsr, 3)))
    # ax, ay, aystd = bindata(fakeflux[ww], (flux[ww] - fakeflux[ww]),
    #                        np.arange(-100, 1000, 200))
    # plt.errorbar(ax, ay, aystd, markersize=10, color='green', fmt='o', label='SMP')

    # plt.axhline(0)
    plt.xlim(-5, 5)

    import matplotlib.mlab as mlab
    import math
    mean = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.arange(-5., 5, .1)
    plt.plot(x, mlab.normpdf(x, mean, sigma), color='black', label='Gaussian Normal')
    plt.legend()
    # plt.ylim(-.1,.1)
    # plt.ylim(-600, 600)
    plt.xlabel('flux/fluxerr')
    plt.ylabel('Count')
    plt.title(filter + ' band')

    plt.savefig(outdir + '/fluxdiffstd.png')

    nullfmt = NullFormatter()  # no labels

    # definitions for the axes
    left, width = 0.13, 0.61
    bottom, height = 0.1, 0.6
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom + height / 2., width+.2, height / 2.]
    rect_scatterflux = [left, bottom, width+.2, height]
    rect_histx = [left, bottom_h - .04, width+.2, .18]
    rect_histy = [left_h, bottom + height / 2., 0.2, height / 2.]
    rect_histyflux = [left_h, bottom, 0.2, height / 2.]

    # start with a rectangular Figure
    plt.figure(1, figsize=(45, 40))

    #ax1 = plt.axes(rect_scatter)
    ax3 = plt.axes(rect_histx)
    #ax2 = plt.axes(rect_histy)
    ax4 = plt.axes(rect_scatterflux)
    #ax5 = plt.axes(rect_histyflux)

    # no labels
    #ax2.yaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    #ax5.yaxis.set_major_formatter(nullfmt)

    #outliers3 = float(len(d[(abs(d)>3.) & (chisqarr > .05) & (chisqarr < 2.5) & (np.array(fakemag, dtype='float') > 0.)]))/float(len(d))
    #outliers5 = float(len(d[(abs(d)>5.) & (chisqarr > .05) & (chisqarr < 2.5) & (np.array(fakemag, dtype='float') > 0.)]))/float(len(d))

    #ax2.hist(d[np.isfinite(d)], bins=np.arange(-10, 10, .25), normed=True,label='RMS Fakemag = 99: ' + str(round(rms99, 3))+
    #                                                            '\nRMS Fakemag < 99: '+ str(round(rmsr, 3))+'\n3sig Outlier'#+
    #                                                            ' Fraction: '+str(round(outliers3,3))+'\n5sig Outlier'+
    #                                                            ' Fraction: '+str(round(outliers5,3))
    #                                                            ,orientation='horizontal')
             #label='RMS: ' + str(round(rms, 3)) + '\nChiSq (3sig cut) ' + str(round(chisq, 3)) + '\nMedian ' + str(
             #   round(np.median(d), 3)) + ' +- ' + str(round(np.std(d), 3)),

    import matplotlib.mlab as mlab
    import math
    mean = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.arange(-5, 5, .1)
    #ax2.plot(mlab.normpdf(x, mean, sigma),x, color='black', label='Gaussian Normal')

    #ax2.set_ylim(-4, 4)
    #ax2.set_xlim(0,.5)
    #.xlabel('STDEV')
    #plt.ylabel('Normalized Count')
    #ax2.legend(fontsize='xx-small',loc=(0.,1.25))
    #plt.savefig('stdresid.png')

    #plt.clf()
    #fakemag[fakemag == 99] = 28.5




    #ax1.scatter(fakemag,d,alpha=.3,color='blue')
    #ax, ay, aystd = dt.bindata(fakemag[(d<3.)& (np.array(fakemag, dtype='float') > 0.)], d[(d<3.)& (np.array(fakemag, dtype='float') > 0.)], np.arange(19., 28, .1),window=.5)
    #ax1.plot([19, 28.7], [0, 0],color='grey')
    #ax1.plot(ax, ay, linewidth=3, color='black', label='SMP')
    #ax1.plot(ax, ay+aystd, linewidth=2, color='black',linestyle='--', label='SMP')
    #ax1.plot(ax, ay-aystd, linewidth=2, color='black',linestyle='--', label='SMP')

    #ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

    #ax1.set_xlim(19, 28.7)
    #ax1.set_ylim(-3., 3.)
    #ax1.set_xlabel('Fake Mag')
    #ax1.set_ylabel('STD')

    #ax, ayrms= dt.binrms(fakemag, d, np.arange(19.5, max(fakemag), .1),.5)
    #ax3.plot(ax, ayrms, color='blue',label='RMS',linewidth=3)


    ax3.plot([0,100],[1.,1.],linestyle='--',color='black')
    ax3.set_ylim(.5,2.0)
    ax3.legend(fontsize='x-small')

    fresid = np.zeros(flux.shape)
    for i,f,ff in zip(range(len(flux)),flux,fakeflux):
        if f == 0.:
            fresid[i] = np.nan
        else:
            fresid[i] = (f - ff) / max([abs(ff),1.])
    #fresid[abs(fakeflux) < 1.] = flux[abs(fakeflux) < 1.] - fakeflux[abs(fakeflux) < 1.]

    #ax5.hist(fresid, bins=np.arange(-.155,.15,.01),color='blue', orientation='horizontal')

    filts = ['g','r','i','z']
    colors = ['green','red','indigo','black']
    if filter == 'all':
    #ax4.scatter(fakemag,fresid,alpha=.03,color='black')
        for filt,col in zip(filts,colors):
            ww = (filterarr == filt) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.)\
                 & (fluxerr > 0.) & (np.isfinite(flux)) & \
                 (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .05) \
                 & (chisqarr < 2.5)
            axa, aya, aystd = dt.bindata(fakemag[ww],fresid[ww],
                                    np.arange(20., 26., .1),window=2.,dontrootn=True)
            ax4.plot([19, 28.7], [0, 0],color='grey')

            ax, ayrms = dt.binrms(fakemag[ww][d[ww]<10.], d[ww][d[ww]<10.], np.arange(20., 28, .1), 1.5)
            ax3.plot(ax, ayrms, color=col, label=filt+' band', linewidth=3,alpha=.8)
            ax3.plot(ax, ax * 0 + 1., linestyle='--', color='black',alpha=.8)
            ax4.plot(axa, aya, linewidth=3, color=col,label=filt+' band',alpha=.8)
            ax4.plot(axa, aya + aystd, linewidth=2, color=col, linestyle='--',alpha=.8)
            ax4.plot(axa, aya - aystd, linewidth=2, color=col, linestyle='--',alpha=.8)
    else:
        axa, aya, aystd = dt.bindata(fakemag, fresid,
                                     np.arange(18., 26., .5), window=.5)
        #ax4.plot([19, 28.7], [0, 0], color='grey')
        #print len(fakemag[fakemag<26])
        #raw_input('fmmmmm')

        ax, ayrms,num = dt.binrms(fakemag[d < 10.], d[d < 10.], np.arange(18., 28, .5), .5,returnn=True)
        ax3.plot(ax, ayrms, color='blue', label='ALL SNe', linewidth=3)
        # ax, ayrms = dt.binrms(fakemag, dz, np.arange(20., 28, .1), 1.5)
        # ax3.plot(ax, ayrms, color='blue',linestyle='--', label='ALL SNe', linewidth=3)
        # ax, ayrms = dt.binrms(fakemag, df, np.arange(20., 28, .1), 1.5)
        # ax3.plot(ax, ayrms, color='red', linestyle='--', label='DIFFIMG', linewidth=3)
        #ax3.plot(ax, ax * 0 + 1., linestyle='--', color='black')
        ax3.axhline(0,c='k',linestyle='--')
        ax4.axhline(0,c='k')

        # ww = hostmag > 25.
        # ax, ayrms = dt.binrms(fakemag[ww], d[ww], np.arange(19.5, max(fakemag), .1), .5)
        # ax3.plot(ax, ayrms, color='red', label='HostMag > 25.', linewidth=3)
        #
        # ww = hostmag < 23.
        # ax, ayrms = dt.binrms(fakemag[ww], d[ww], np.arange(19.5, max(fakemag), .1), .5)
        # ax3.plot(ax, ayrms, color='green', label='HostMag < 23', linewidth=3)
        #ax3.legend(fontsize='x-small',location='upper right')

        ax4.plot(axa, aya, linewidth=3, color='blue')
        ax4.plot(axa, aya+(aystd), linewidth=2, color='blue',linestyle='--')
        ax4.plot(axa, aya-aystd, linewidth=2, color='blue',linestyle='--')
    ax4.set_xlim(18.5,26)
    ax4.set_ylim(-.05,.05)
    ax4.set_xlabel('Fake Mag')
    #ax5.set_xlabel('Counts')
    ax3.set_ylabel('RMS')
    #ax3.set_title(filter+' band')
    ax4.set_ylabel('(fitflux - fakeflux)/fakeflux')

    if not filter == 'all':
        ax3.set_title(filter+' band')
    else:
        ax3.set_title(deep_or_shallow.upper()+' Fields')


    ax3.set_xlim(ax4.get_xlim())
    ax4.legend(fontsize='x-small',loc='upper right')
    #ax2.set_ylim(ax1.get_ylim())
    #ax5.set_ylim(ax4.get_ylim())
    #ax2.xaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    #ax1.xaxis.set_major_formatter(nullfmt)

    plt.subplots_adjust(wspace=0.01,hspace=0.01)

    plt.savefig(outdir+'/'+deep_or_shallow+'std.png')
    print 'saved' , outdir+'/'+deep_or_shallow+'std.png'
    #raw_input('press to continue')

    #--------------------------------------------------------------------------------------








    plt.clf()


    #nullfmt = NullFormatter()  # no labels

    # definitions for the axes
    left, width = 0.13, 0.61
    bottom, height = 0.1, 0.6
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom + height / 2., width + .2, height / 2.]
    #rect_scatterflux = [left, bottom, width + .2, height]
    rect_histx = [left, bottom, width + .2, height+.18]
    rect_histy = [left_h, bottom + height / 2., 0.2, height / 2.]
    rect_histyflux = [left_h, bottom, 0.2, height / 2.]

    # start with a rectangular Figure
    plt.figure(1, figsize=(45, 40))

    # ax1 = plt.axes(rect_scatter)
    ax3 = plt.axes(rect_histx)

    import math
    variance = 1


    ax3.plot([0, 100], [1., 1.], linestyle='--', color='black')
    ax3.set_ylim(.5, 3.0)
    ax3.legend(fontsize='x-small')

    fresid = np.zeros(flux.shape)
    for i, f, ff in zip(range(len(flux)), flux, fakeflux):
        if f == 0.:
            fresid[i] = np.nan
        else:
            fresid[i] = (f - ff) / max([abs(ff), 1.])

    filts = ['g', 'r', 'i', 'z']
    colors = ['green', 'red', 'indigo', 'black']
    if filter == 'all':
        # ax4.scatter(fakemag,fresid,alpha=.03,color='black')
        for filt, col in zip(filts, colors):
            ww = filterarr == filt
            axa, aya, aystd = dt.bindata(hostmag[ww], fresid[ww],
                                         np.arange(20., 26., .1), window=2., dontrootn=True)
            #ax4.plot([19, 28.7], [0, 0], color='grey')

            ax, ayrms = dt.binrms(hostmag[ww][d[ww] < 3.], d[ww][d[ww] < 3.], np.arange(18., 35., .3), 2.5)
            ax3.plot(ax, ayrms, color=col, label=filt + ' band', linewidth=3, alpha=.8)
            ax3.plot(ax, ax * 0 + 1., linestyle='--', color='black', alpha=.8)
            #ax4.plot(axa, aya, linewidth=3, color=col, label=filt + ' band', alpha=.8)
            #ax4.plot(axa, aya + aystd, linewidth=2, color=col, linestyle='--', alpha=.8)
            #ax4.plot(axa, aya - aystd, linewidth=2, color=col, linestyle='--', alpha=.8)
    else:
        axa, aya, aystd = dt.bindata(hostmag, fresid,
                                     np.arange(20., 26., .1), window=2., dontrootn=True)
        #ax4.plot([19, 28.7], [0, 0], color='grey')

        ax, ayrms = dt.binrms(hostmag[d < 3.], d[d < 3.], np.arange(20., 28, .5), .5)
        ax3.plot(ax, ayrms, color='blue', label='ALL SNe', linewidth=3)

        ax3.plot(ax, ax * 0 + 1., linestyle='--', color='black')


    ax3.set_ylabel('RMS')
    # ax3.set_title(filter+' band')
    #ax4.set_ylabel('(fitflux - fakeflux)/fakeflux')

    if not filter == 'all':
        ax3.set_title(filter + ' band')
    else:
        ax3.set_title(deep_or_shallow.upper() + ' Fields')

    ax3.set_xlim(20,28)

    ax3.set_xlabel('Host Surface Brightness Mag')
    plt.subplots_adjust(wspace=0.01, hspace=0.01)

    plt.savefig(outdir + '/' + deep_or_shallow + 'hoststd.png')
    print 'saved', outdir + '/' + deep_or_shallow + 'hoststd.png'

    plt.clf()
    ax, ayrms = dt.binrms(fwhm[fakemag>30]*2.35*.27,d[fakemag>30], np.arange(1.5, 4., .1), .1)
    plt.plot(ax, ayrms, color='blue', label='ALL SNe', linewidth=3)
    #plt.scatter(fwhm,d,color='grey')
    plt.xlabel('FWHM arcsec')
    #plt.plot(ax, ax * 0 + 1., linestyle='--', color='black')
    plt.axhline(1.)
    plt.savefig(outdir + '/' + deep_or_shallow + 'fwhmstd.png')




    plt.clf()

    sig = 1.48 * np.median(np.abs(d[(abs(d) < 3.)&(fakemag>30)&(hostmag >26.)]))
    plt.hist(d[(fakemag>30)&(hostmag >26.)],bins=np.arange(-5.1,5,.2),label='Sigma %.2f'%sig,normed=True)
    import math
    mean = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.arange(-5, 5, .1)
    plt.plot(x,mlab.normpdf(x, mean, sigma), color='black', label='Gaussian Normal')
    plt.title(filter + ' band')
    plt.savefig(outdir + '/' + 'emtpysig.png')











    #--------------------------------------------------------------------------------------


    plt.clf()

    f,axes = plt.subplots(2,sharex=True)
    axes = axes.ravel()

    ww = (fakemag >= 50) & (flux != 0.) & (hostmag < 22.) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.)\
                 & (fluxerr > 0.) & (np.isfinite(flux)) & \
                 (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .05) \
                 & (chisqarr < 2.5)
    ax, ayrms = dt.binrms(fitzpt[ww], (flux/fluxerr)[ww], np.arange(28., 35., .1), 1.)

    axes[1].plot(ax, ayrms, color='blue', label='Hostmag lt 21', linewidth=3, alpha=.8)
    axes[1].plot(ax, ax * 0 + 1., linestyle='--', color='black', alpha=.8)

    ww = (fakemag >= 50) & (flux != 0.) & (hostmag > 22.) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.)\
                 & (fluxerr > 0.) & (np.isfinite(flux)) & \
                 (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .05) \
                 & (chisqarr < 2.5)
    ax, ayrms = dt.binrms(fitzpt[ww], (flux / fluxerr)[ww], np.arange(28., 35., .1), 1.)

    axes[1].plot(ax, ayrms, color='green', label='Hostmag gt 21', linewidth=3, alpha=.8)

    axes[1].set_xlabel('Fit ZPT')
    axes[1].set_ylabel('RMS')
    axes[0].hist(fitzpt[ww],bins=np.arange(27,35,.2))
    axes[0].set_xlim(27,35)
    axes[1].set_xlim(27,35)
    plt.title('NO SN FLUX IN IMAGE')
    plt.legend()
    plt.savefig(outdir + '/' + deep_or_shallow +  '_'+filter+'_zptcorr.png')
    print 'saved',outdir + '/' + deep_or_shallow +  '_'+filter+'_zptcorr.png'








    plt.clf()

    f, axes = plt.subplots(2, sharex=True)
    axes = axes.ravel()

    ww = (fakemag >= 50) & (flux != 0.) & (hostmag < 299999.) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.) \
         & (fluxerr > 0.) & (np.isfinite(flux)) & \
         (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .05) \
         & (chisqarr < 2.5)
    ax, ayrms = dt.binrms(sky[ww], (flux / fluxerr)[ww], np.arange(-1005,1000,10.), 50.)

    axes[1].plot(ax, ayrms, color='blue', label='Hostmag lt 21', linewidth=3, alpha=.8)
    axes[1].plot(ax, ax * 0 + 1., linestyle='--', color='black', alpha=.8)

    # ww = (fakemag >= 50) & (flux != 0.) & (hostmag > 22.) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.) \
    #      & (fluxerr > 0.) & (np.isfinite(flux)) & \
    #      (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .2) \
    #      & (chisqarr < 1.2)
    #
    #
    # ax, ayrms = dt.binrms(sky[ww], (flux / fluxerr)[ww], np.arange(100., 25000., 100), 500.)
    #
    # plt.plot(ax, ayrms, color='green', label='Hostmag gt 21', linewidth=3, alpha=.8)

    axes[1].set_xlabel('Sky')
    axes[1].set_ylabel('RMS')
    axes[0].axvline(-50)
    axes[0].axvline(50)
    axes[1].axvline(-50)
    axes[1].axvline(50)
    axes[0].set_title('NO SN FLUX IN IMAGE')
    axes[0].hist(sky[ww], bins=np.arange(-1005., 1000, 1.))
    axes[0].set_xlim(-50, 50)
    axes[1].set_xlim(-50, 50)
    # plt.legend()
    plt.savefig(outdir + '/' + deep_or_shallow +  '_'+filter+'_skycorr.png')
    print 'saved', outdir + '/' + deep_or_shallow +  '_'+filter+'_skycorr.png'

    ww1 = (abs(sky)<50.) & (fakemag >= 50) & (flux != 0.) & (hostmag < 299999.) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.) \
         & (fluxerr > 0.) & (np.isfinite(flux)) & \
         (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .05) \
         & (chisqarr < 2.5)
    ww2 = (fakemag >= 50) & (flux != 0.) & (hostmag < 299999.) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.) \
         & (fluxerr > 0.) & (np.isfinite(flux)) & \
         (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .05) \
         & (chisqarr < 2.5)
    #print 'percentage loss from sky cut of 50',float(len(sky[ww2])-len(sky[ww1]))/len(sky[ww2])
    #print len(sky[ww2])
    #print len(sky[ww1])
    #raw_input()
    # ww = (abs(sky)<50.)
    #
    # print 'percentage loss from sky cut of 50 in good conditions',float(len(sky)-len(sky[ww]))/len(sky)



    plt.clf()

    f, axes = plt.subplots(2, sharex=True)
    axes = axes.ravel()

    ww = (fakemag >= 50) & (flux != 0.) & (hostmag < 299999.) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.) \
         & (fluxerr > 0.) & (np.isfinite(flux)) & \
         (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .05) \
         & (chisqarr < 2.5)
    ax, ayrms = dt.binrms(sky[ww]/skyerr[ww], (flux / fluxerr)[ww], np.arange(-30.25,30,.5), .5)

    axes[1].plot(ax, ayrms, color='blue', label='Hostmag lt 21', linewidth=3, alpha=.8)
    axes[1].plot(ax, ax * 0 + 1., linestyle='--', color='black', alpha=.8)

    # ww = (fakemag >= 50) & (flux != 0.) & (hostmag > 22.) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.) \
    #      & (fluxerr > 0.) & (np.isfinite(flux)) & \
    #      (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .2) \
    #      & (chisqarr < 1.2)
    #
    #
    # ax, ayrms = dt.binrms(sky[ww], (flux / fluxerr)[ww], np.arange(100., 25000., 100), 500.)
    #
    # plt.plot(ax, ayrms, color='green', label='Hostmag gt 21', linewidth=3, alpha=.8)

    axes[1].set_xlabel('Sky/Skyerr')
    axes[1].set_ylabel('RMS')
    axes[0].set_title('NO SN FLUX IN IMAGE')
    axes[0].hist(sky[ww]//skyerr[ww], bins=np.arange(-30.25, 30, .5))
    axes[0].set_xlim(-30, 30)
    axes[1].set_xlim(-30, 30)
    # plt.legend()
    plt.savefig(outdir + '/' + deep_or_shallow +  '_'+filter+'_skyskyerrcorr.png')
    print 'saved', outdir + '/' + deep_or_shallow +  '_'+filter+'_skyskyerrcorr.png'


    plt.clf()

    f, axes = plt.subplots(2, sharex=True)
    axes = axes.ravel()

    ww = (fakemag >= 50) & (flux != 0.) & (hostmag < 2999999.) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.) \
         & (fluxerr > 0.) & (np.isfinite(flux)) & \
         (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .05) \
         & (chisqarr < 2.5)
    ax, ayrms = dt.binrms(skyerr[ww], (flux / fluxerr)[ww], np.arange(0, 600, 20), 40.)

    axes[1].plot(ax, ayrms, color='blue', label='Hostmag lt 21', linewidth=3, alpha=.8)
    axes[1].plot(ax, ax * 0 + 1., linestyle='--', color='black', alpha=.8)

    # ww = (fakemag >= 50) & (flux != 0.) & (hostmag > 22.) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.) \
    #      & (fluxerr > 0.) & (np.isfinite(flux)) & \
    #      (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .2) \
    #      & (chisqarr < 1.2)
    #
    # ax, ayrms = dt.binrms(skyerr[ww], (flux / fluxerr)[ww], np.arange(0., 600., 20.), 40.)
    #
    # plt.plot(ax, ayrms, color='green', label='Hostmag gt 21', linewidth=3, alpha=.8)

    axes[1].set_xlabel('Skyerr')
    axes[1].set_ylabel('RMS')
    axes[0].set_title('NO SN FLUX IN IMAGE')
    axes[0].hist(skyerr[ww], bins=np.arange(-1005, 1000, 10))
    axes[0].set_xlim(0, 500)
    axes[1].set_xlim(0, 500)
    # plt.legend()
    plt.savefig(outdir + '/' + deep_or_shallow + '_'+filter+'_skyerrcorr.png')
    print 'saved', outdir + '/' + deep_or_shallow +  '_'+filter+'_skyerrcorr.png'






    plt.clf()

    f, axes = plt.subplots(2, sharex=True)
    axes = axes.ravel()

    ww = (fakemag >= 50) & (flux != 0.) & (hostmag < 299999.) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.) \
         & (fluxerr > 0.) & (np.isfinite(flux)) & \
         (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .05) \
         & (chisqarr < 2.5)
    ax, ayrms = dt.binrms((skyerr/sky)[ww], (flux / fluxerr)[ww], np.arange(-10.,10.,.1), .25)

    axes[1].plot(ax, ayrms, color='blue', label='Hostmag lt 21', linewidth=3, alpha=.8)
    axes[1].plot(ax, ax * 0 + 1., linestyle='--', color='black', alpha=.8)

    # ww = (fakemag >= 50) & (flux != 0.) & (hostmag > 22.) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.) \
    #      & (fluxerr > 0.) & (np.isfinite(flux)) & \
    #      (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .2) \
    #      & (chisqarr < 1.2)
    #
    #
    # ax, ayrms = dt.binrms(sky[ww], (flux / fluxerr)[ww], np.arange(100., 25000., 100), 500.)
    #
    # plt.plot(ax, ayrms, color='green', label='Hostmag gt 21', linewidth=3, alpha=.8)

    axes[1].set_xlabel('skyerr/Sky')
    axes[1].set_ylabel('RMS')
    #axes[0].axvline(-50)
    #axes[0].axvline(50)
    #axes[1].axvline(-50)
    #axes[1].axvline(50)
    axes[0].set_title('NO SN FLUX IN IMAGE')
    axes[0].hist(sky[ww], bins=np.arange(-10., 10., .1))
    axes[0].set_xlim(-10, 10.)
    axes[1].set_xlim(-10., 10.)
    # plt.legend()
    plt.savefig(outdir + '/' + deep_or_shallow +  '_'+filter+'_pskycorr.png')
    print 'saved',outdir + '/' + deep_or_shallow +  '_'+filter+'_pskycorr.png'


    # plt.clf()
    # ww = (fakemag >= 50) & (flux != 0.) & (hostmag < 2999999.) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.) \
    #      & (fluxerr > 0.) & (np.isfinite(flux)) & \
    #      (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .2) \
    #      & (chisqarr < 1.2)
    # ax, ayrms = dt.binrms((skyerr/sky)[ww], (flux / fluxerr)[ww], np.arange(0, .1, .001), .005)
    #
    # plt.plot(ax, ayrms, color='blue', label='Hostmag lt 21', linewidth=3, alpha=.8)
    # plt.plot(ax, ax * 0 + 1., linestyle='--', color='black', alpha=.8)
    #
    # # ww = (fakemag >= 50) & (flux != 0.) & (hostmag > 22.) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.) \
    # #      & (fluxerr > 0.) & (np.isfinite(flux)) & \
    # #      (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .2) \
    # #      & (chisqarr < 1.2)
    # #
    # # ax, ayrms = dt.binrms(skyerr[ww], (flux / fluxerr)[ww], np.arange(0., 600., 20.), 40.)
    # #
    # # plt.plot(ax, ayrms, color='green', label='Hostmag gt 21', linewidth=3, alpha=.8)
    #
    # plt.xlabel('Skyerr/Sky')
    # plt.ylabel('RMS')
    # plt.title('NO SN FLUX IN IMAGE')
    # # plt.legend()
    # plt.savefig(outdir + '/' + deep_or_shallow + '_'+filter+'_pskyerrcorr.png')
    # print 'saved', outdir + '/' + deep_or_shallow + '_'+filter+'_pskyerrcorr.png'



    plt.clf()
    plt.scatter(sky,fresid,alpha=.2)
    ax, ay, aystd = dt.bindata(sky, fresid, np.arange(-500, 500, 50), window=100)
    plt.plot(ax, ay, linewidth=3, color='orange', label='SMP')
    plt.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP')
    plt.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP')

    plt.xlabel('Sky')
    plt.xlim(-500,500)
    plt.ylabel('Fractional flux diff')
    plt.ylim(-.5,.5)
    #plt.title('NO SN FLUX IN IMAGE')
    # plt.legend()
    plt.savefig(outdir + '/' + deep_or_shallow + 'zptvssky.png')

    #raw_input('press to continue')

    plt.clf()
    fig = plt.figure(figsize=(25, 20))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    nullfmt = NullFormatter()  # no labels

    # definitions for the axes
    # definitions for the axes
    left, width = 0.13, 0.61
    bottom, height = 0.1, 0.6
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom + height / 2., width, height / 2.]
    rect_scatterflux = [left, bottom, width, height / 2.]
    rect_histx = [left, bottom_h - .04, width, .18]
    rect_histy = [left_h, bottom + height / 2., 0.2, height / 2.]
    rect_histyflux = [left_h, bottom, 0.2, height / 2.]
    # start with a rectangular Figure
    try:
        plt.figure(1, figsize=(32, 24))

        ax1 = plt.axes(rect_scatter)
        ax3 = plt.axes(rect_histx)
        ax2 = plt.axes(rect_histy)
        ax4 = plt.axes(rect_scatterflux)
        ax5 = plt.axes(rect_histyflux)

        # no labels
        ax2.yaxis.set_major_formatter(nullfmt)
        ax3.xaxis.set_major_formatter(nullfmt)
        ax5.yaxis.set_major_formatter(nullfmt)


        ax2.hist(d, bins=np.arange(-10, 10, .25), normed=True,label='RMS Fakemag = 99: ' + str(round(rms99, 3))+
                                                                    '\nRMS Fakemag < 99: '+ str(round(rmsr, 3))
                 ,orientation='horizontal')
                 #label='RMS: ' + str(round(rms, 3)) + '\nChiSq (3sig cut) ' + str(round(chisq, 3)) + '\nMedian ' + str(
                 #   round(np.median(d), 3)) + ' +- ' + str(round(np.std(d), 3)),

        import matplotlib.mlab as mlab
        import math
        mean = 0
        variance = 1
        sigma = math.sqrt(variance)
        x = np.arange(-5, 5, .1)
        ax2.plot(mlab.normpdf(x, mean, sigma),x, color='black', label='Gaussian Normal')

        ax2.set_ylim(-4, 4)
        ax2.set_xlim(0,.5)
        #.xlabel('STDEV')
        #plt.ylabel('Normalized Count')
        ax2.legend(fontsize='xx-small',loc=(0.,1.25))
        #plt.savefig('stdresid.png')

        #plt.clf()
        ww = fakemag < 30000.
        ax1.scatter(hostmag[ww],d[ww],alpha=.3,color='blue')
        ax, ay, aystd = dt.bindata(hostmag[ww], d[ww], np.arange(19.5, 28.5, .1),window=1.5)
        ax1.plot([min(hostmag), max(hostmag)], [0, 0],color='grey')
        ax1.plot(ax, ay, linewidth=3, color='orange', label='SMP')
        ax1.plot(ax, ay+aystd, linewidth=2, color='orange',linestyle='--', label='SMP')
        ax1.plot(ax, ay-aystd, linewidth=2, color='orange',linestyle='--', label='SMP')

        #ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

        ax1.set_xlim(19.5,29)
        ax1.set_ylim(-3., 3.)
        ax1.set_xlabel('Host Mag')
        ax1.set_ylabel('STD')

        #ax3.plot(ax, ayrms, color='blue',label='RMS',linewidth=3)


        ax3.plot([0,100],[1.,1.],linestyle='--',color='black')
        ax3.set_ylim(.9,1.25)

        fresid = np.zeros(flux.shape)
        for i,f,ff in zip(range(len(flux)),flux,fakeflux):
            if f == 0.:
                fresid[i] = np.nan
            else:
                fresid[i] = (f - ff) / max([abs(ff),1.])
        #fresid[abs(fakeflux) < 1.] = flux[abs(fakeflux) < 1.] - fakeflux[abs(fakeflux) < 1.]

        ax5.hist(fresid[ww], bins=np.arange(-.155,.15,.01),color='blue', orientation='horizontal')

        ax4.scatter(hostmag[ww],fresid[ww],alpha=.3,color='blue')
        ax, ay, aystd = dt.bindata(hostmag[ww],fresid[ww],
                                np.arange(19.5, 28.5, .1),window=1.)
        ax4.plot([min(hostmag), max(hostmag)], [0, 0],color='grey')

        ax4.plot(ax, ay, linewidth=3, color='orange')
        ax4.plot(ax, ay+aystd, linewidth=2, color='orange',linestyle='--')
        ax4.plot(ax, ay-aystd, linewidth=2, color='orange',linestyle='--')
        ax4.set_xlim(ax1.get_xlim())
        ax4.set_ylim(-.2,.2)
        ax4.set_xlabel('Host Mag')
        ax5.set_xlabel('Counts')
        ax3.set_ylabel('RMS')
        ax3.set_title(filter + ' band')

        ax4.set_ylabel('(fitflux - fakeflux)/fakeflux')

        ax3.set_xlim(ax1.get_xlim())
        ax3.set_ylim(.9,1.25)
        ax2.set_ylim(ax1.get_ylim())
        ax5.set_ylim(ax4.get_ylim())
        ax2.xaxis.set_major_formatter(nullfmt)
        ax3.xaxis.set_major_formatter(nullfmt)
        ax1.xaxis.set_major_formatter(nullfmt)
        plt.subplots_adjust(wspace=0.001,hspace=0.001)

        ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag[ww]), max(hostmag[ww]), .1), .5)
        ax3.plot(ax, ayrms, color='blue', label='ALL SNe w/ Light', linewidth=3)
        ax3.plot(ax, ax * 0 + 1., linestyle='--',color='black')

        # ww = fakemag > 28.
        # ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag), max(hostmag), .1), 1.5)
        # ax3.plot(ax, ayrms, color='red', label='FakeMag = 99', linewidth=3)
        #
        # ww = fakemag < 22.
        # ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag), max(hostmag), .1), 1.5)
        # ax3.plot(ax, ayrms, color='green', label='FakeMag < 22', linewidth=3)
        # ax3.legend(fontsize='small')
        #
        # ww = (fakemag > 22.) & (fakemag < 28)
        # ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag), max(hostmag), .1), 1.5)
        # ax3.plot(ax, ayrms, color='purple', label='FakeMag < 22', linewidth=3)
        # ax3.legend(fontsize='small')
        plt.title(filter + ' band')

        plt.savefig(outdir+'/hostmagstd.png')

    except:
        print 'bad hostmags'

    plt.clf()
    fig = plt.figure(figsize=(25, 20))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    nullfmt = NullFormatter()  # no labels

    # definitions for the axes
    # definitions for the axes
    left, width = 0.13, 0.61
    bottom, height = 0.1, 0.6
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom + height / 2., width, height / 2.]
    rect_scatterflux = [left, bottom, width, height / 2.]
    rect_histx = [left, bottom_h - .04, width, .18]
    rect_histy = [left_h, bottom + height / 2., 0.2, height / 2.]
    rect_histyflux = [left_h, bottom, 0.2, height / 2.]
    # start with a rectangular Figure
    if True:
        plt.figure(1, figsize=(32, 24))

        ax1 = plt.axes(rect_scatter)
        ax3 = plt.axes(rect_histx)
        ax2 = plt.axes(rect_histy)
        #ax4 = plt.axes(rect_scatterflux)
        #ax5 = plt.axes(rect_histyflux)

        # no labels
        ax2.yaxis.set_major_formatter(nullfmt)
        ax3.xaxis.set_major_formatter(nullfmt)
        #ax5.yaxis.set_major_formatter(nullfmt)

        ax2.hist(d, bins=np.arange(-10, 10, .25), normed=True, label='RMS Fakemag = 99: ' + str(round(rms99, 3)) +
                                                                     '\nRMS Fakemag < 99: ' + str(round(rmsr, 3))
                 , orientation='horizontal')
        # label='RMS: ' + str(round(rms, 3)) + '\nChiSq (3sig cut) ' + str(round(chisq, 3)) + '\nMedian ' + str(
        #   round(np.median(d), 3)) + ' +- ' + str(round(np.std(d), 3)),

        import matplotlib.mlab as mlab
        import math
        mean = 0
        variance = 1
        sigma = math.sqrt(variance)
        x = np.arange(-5, 5, .1)
        ax2.plot(mlab.normpdf(x, mean, sigma), x, color='black', label='Gaussian Normal')

        ax2.set_ylim(-4, 4)
        ax2.set_xlim(0, .5)
        # .xlabel('STDEV')
        # plt.ylabel('Normalized Count')
        ax2.legend(fontsize='xx-small', loc=(0., 1.25))
        # plt.savefig('stdresid.png')



        # plt.clf()
        ww = (fakemag >= 50) & (flux != 0.) & (hostmag < 299999.) & (flux != 0) & (np.array(fakemag, dtype='float') > 0.) \
         & (fluxerr > 0.) & (np.isfinite(flux)) & \
         (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .05) \
         & (chisqarr < 2.5)

        #print 'fakemag',fakemag

        ax1.scatter(hostmag[ww], d[ww], alpha=.3, color='blue')
        ax, ay, aystd = dt.bindata(hostmag[ww], d[ww], np.arange(min(hostmag), 27.5, .1), window=1.5)
        ax1.plot([min(hostmag), max(hostmag)], [0, 0], color='grey')
        ax1.plot(ax, ay, linewidth=3, color='orange', label='SMP')
        ax1.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP')
        ax1.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP')

        # ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

        ax1.set_xlim(20, 29)
        ax1.set_ylim(-3., 3.)
        ax1.set_xlabel('Host Mag')
        ax1.set_ylabel('STD')

        # ax3.plot(ax, ayrms, color='blue',label='RMS',linewidth=3)


        ax3.plot([0, 100], [1., 1.], linestyle='--', color='black')
        ax3.set_ylim(.9, 1.25)

        fresid = np.zeros(flux.shape)
        for i, f, ff in zip(range(len(flux)), flux, fakeflux):
            if f == 0.:
                fresid[i] = np.nan
            else:
                fresid[i] = (f - ff) / max([abs(ff), 1.])
        # fresid[abs(fakeflux) < 1.] = flux[abs(fakeflux) < 1.] - fakeflux[abs(fakeflux) < 1.]

        #ax5.hist(fresid[ww], bins=np.arange(-.155, .15, .01), color='blue', orientation='horizontal')

        #ax4.scatter(hostmag[ww], fresid[ww], alpha=.3, color='blue')
        #ax, ay, aystd = dt.bindata(hostmag[ww], fresid[ww],
        #                           np.arange(min(hostmag), 27.5, .1), window=1.)
        #ax4.plot([min(hostmag), max(hostmag)], [0, 0], color='grey')

        #ax4.plot(ax, ay, linewidth=3, color='orange')
        #ax4.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--')
        #ax4.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--')
        #ax4.set_xlim(ax1.get_xlim())
        #ax4.set_ylim(-.2, .2)
        #ax3.set_xlabel('Host Mag')
        ax5.set_xlabel('Counts')
        ax3.set_ylabel('RMS')
        ax3.set_title(filter + ' band')

        #ax4.set_ylabel('(fitflux - fakeflux)/fakeflux')

        ax3.set_xlim(ax1.get_xlim())
        ax3.set_ylim(.9, 1.25)
        ax2.set_ylim(ax1.get_ylim())
        #ax5.set_ylim(ax4.get_ylim())
        ax2.xaxis.set_major_formatter(nullfmt)
        ax3.xaxis.set_major_formatter(nullfmt)
        #ax1.xaxis.set_major_formatter(nullfmt)
        plt.subplots_adjust(wspace=0.001, hspace=0.001)

        #ww = fakemag > 90.
        ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag[ww]), max(hostmag[ww]), .1), .5)
        ax3.plot(ax, ayrms, color='blue', label='ALL SNe', linewidth=3)
        ax3.plot(ax, ax * 0 + 1., linestyle='--', color='black')

        # ww = fakemag > 28.
        # ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag), max(hostmag), .1), 1.5)
        # ax3.plot(ax, ayrms, color='red', label='FakeMag = 99', linewidth=3)
        #
        # ww = fakemag < 22.
        # ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag), max(hostmag), .1), 1.5)
        # ax3.plot(ax, ayrms, color='green', label='FakeMag < 22', linewidth=3)
        # ax3.legend(fontsize='small')
        #
        # ww = (fakemag > 22.) & (fakemag < 28)
        # ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag), max(hostmag), .1), 1.5)
        # ax3.plot(ax, ayrms, color='purple', label='FakeMag < 22', linewidth=3)
        # ax3.legend(fontsize='small')
        plt.title(filter + ' band')

        plt.savefig(outdir + '/emptyhostmagstd.png')

    #except:
    #    print 'bad hostmags'



    # ----------------------------------------------------------------------------------------------------
    plt.clf()
    fig = plt.figure(figsize=(25, 20))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    nullfmt = NullFormatter()  # no labels

    # definitions for the axes
    # definitions for the axes
    left, width = 0.13, 0.61
    bottom, height = 0.1, 0.6
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom + height / 2., width, height / 2.]
    rect_scatterflux = [left, bottom, width, height / 2.]
    rect_histx = [left, bottom_h - .04, width, .18]
    rect_histy = [left_h, bottom + height / 2., 0.2, height / 2.]
    rect_histyflux = [left_h, bottom, 0.2, height / 2.]
    # start with a rectangular Figure
    plt.figure(1, figsize=(25, 20))

    ax1 = plt.axes(rect_scatter)
    ax3 = plt.axes(rect_histx)
    ax2 = plt.axes(rect_histy)
    ax4 = plt.axes(rect_scatterflux)
    ax5 = plt.axes(rect_histyflux)

    # no labels
    ax2.yaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax5.yaxis.set_major_formatter(nullfmt)

    ax2.hist(d[np.isfinite(d)], bins=np.arange(-10, 10, .25), normed=True, label='RMS Fakemag = 99: ' + str(round(rms99, 3)) +
                                                                 '\nRMS Fakemag < 99: ' + str(round(rmsr, 3))
             , orientation='horizontal')
    # label='RMS: ' + str(round(rms, 3)) + '\nChiSq (3sig cut) ' + str(round(chisq, 3)) + '\nMedian ' + str(
    #   round(np.median(d), 3)) + ' +- ' + str(round(np.std(d), 3)),

    import matplotlib.mlab as mlab
    import math
    mean = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.arange(-5, 5, .1)
    ax2.plot(mlab.normpdf(x, mean, sigma), x, color='black', label='Gaussian Normal')

    ax2.set_ylim(-4, 4)
    ax2.set_xlim(0, .5)
    # .xlabel('STDEV')
    # plt.ylabel('Normalized Count')
    ax2.legend(fontsize='xx-small',loc=(0.,1.25))
    # plt.savefig('stdresid.png')

    # plt.clf()
    ax1.scatter(chisqarr, d, alpha=.3, color='blue')
    ax, ay, aystd = dt.bindata(chisqarr, d, np.arange(0.6, 1.2, .001), window=.01)
    ax1.plot([min(chisqarr), max(chisqarr)], [0, 0], color='grey')
    ax1.plot(ax, ay, linewidth=3, color='orange', label='SMP')
    ax1.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP')
    ax1.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP')

    # ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

    ax1.set_xlim(0.6, 1.2)
    ax1.set_ylim(-3., 3.)
    ax1.set_xlabel('Chi Sq')
    ax1.set_ylabel('STD')

    ax, ayrms = dt.binrms(chisqarr, d, np.arange(0.6, 1.2, .001), .01)
    # ax3.plot(ax, ayrms, color='blue',label='RMS',linewidth=3)


    #ax3.plot([0, 100], [1., 1.], linestyle='--', color='black')
    ax3.set_ylim(.9, 1.25)

    fresid = np.zeros(flux.shape)
    for i, f, ff in zip(range(len(flux)), flux, fakeflux):
        if f == 0.:
            fresid[i] = np.nan
        else:
            fresid[i] = (f - ff) / max([abs(ff), 1.])
    # fresid[abs(fakeflux) < 1.] = flux[abs(fakeflux) < 1.] - fakeflux[abs(fakeflux) < 1.]

    ax5.hist(fresid, bins=np.arange(-.155, .15, .01), color='blue', orientation='horizontal')

    ax4.scatter(chisqarr, fresid, alpha=.3, color='blue')
    ax, ay, aystd = dt.bindata(chisqarr, fresid,
                               np.arange(0.6, 1.2, .001), window=.01)
    ax4.plot([min(chisqarr), max(chisqarr)], [0, 0], color='grey')

    ax4.plot(ax, ay, linewidth=3, color='orange')
    ax4.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--')
    ax4.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--')
    ax4.set_xlim(ax1.get_xlim())
    ax4.set_ylim(-.15, .15)
    ax4.set_xlabel('Chi Sq')
    ax5.set_xlabel('Counts')
    ax3.set_ylabel('RMS')
    ax3.set_title(filter+' band')

    ax4.set_ylabel('(fitflux - fakeflux)/fakeflux')

    ax3.set_xlim(ax1.get_xlim())
    ax2.set_ylim(ax1.get_ylim())
    ax5.set_ylim(ax4.get_ylim())
    ax2.xaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax1.xaxis.set_major_formatter(nullfmt)
    plt.subplots_adjust(wspace=0.001, hspace=0.001)

    ax, ayrms = dt.binrms(chisqarr, d, np.arange(0.6, 1.2, .005), .02)
    ax3.plot(ax, ayrms, color='blue', label='ALL SNe', linewidth=3)
    ax3.plot(ax, ax * 0 + 1., linestyle='--', color='black')

    # try:
    #     ww = hostmag > 25
    #     ax, ayrms = dt.binrms(chisqarr[ww], d[ww], np.arange(0.8, 1.2, .005), .02)
    #     ax3.plot(ax, ayrms, color='red', label='Hostmag > 25', linewidth=3)
    #
    #     ww = hostmag < 23.
    #     ax, ayrms = dt.binrms(chisqarr[ww], d[ww], np.arange(0.8, 1.2, .005), .02)
    #     ax3.plot(ax, ayrms, color='green', label='HostMag < 23', linewidth=3)
    #     ax3.legend(fontsize='x-small',loc=2)
    # except:
    #     print 'bad hostmags'

    plt.title(filter+' band')

    plt.savefig(outdir+'/chisqstd.png')

    plt.clf()
    plt.hist(chisqarr,bins=np.arange(0.,1.5,.01))
    plt.xlim(.6,1.4)
    plt.xlabel('chi sq')
    plt.title(filter+' band')

    plt.savefig(outdir+'/chisqhist.png')

    plt.clf()
    #plt.hist(chisqarr, bins=np.arange(0., 1.5, .01))
    plt.scatter(chisqarr,d)
    plt.xlim(.0, 1.4)
    plt.ylim(-4,4)
    plt.xlabel('chi sq')
    plt.ylabel('rms')
    plt.title(filter + ' band')

    plt.savefig(outdir + '/chisqvsrms.png')



    print 'saved stdresid.png'




def plotstarrms(flux,fluxerr,zpt,catmag,chisq,rmsaddin,sky,skyerr,poisson,indices,ras,decs,fwhm,zptscat,mjd,field,ccd,band,title='',outdir='.'):
    catflux = 10 ** (.4 * (zpt - catmag))
    ff = (flux - catflux) / catflux
    st = np.std(ff)
    #luxerr *= 1.0857
    fluxerrorig = fluxerr
    fluxerr = np.sqrt(fluxerr**2+ (zptscat*flux)**2 )




    # plt.clf()
    # plt.scatter(catmag[(np.floor(mjd) == 56575)],fluxerr[(np.floor(mjd) == 56575)], color='black')
    # plt.xlabel('SkyERR')
    # plt.ylabel('#')
    # plt.savefig(outdir + '/' + title + 'photerr.png')
    # print outdir + '/' + title + 'photerr.png'
    # #raw_input('skyerrdone')
    #
    # plt.clf()
    # plt.hist(skyerr[np.floor(mjd)==56575],bins=np.arange(10,70,1.),color='black')
    # plt.xlabel('SkyERR')
    # plt.ylabel('#')
    # plt.savefig(outdir + '/' + title + 'skyerrhist.png')
    # print outdir + '/' + title + 'skyerrhist.png'
    #raw_input('skyerrdone')

    #print max(catmag)
    #print catmag.shape,rmsaddin.shape,ff.shape,indices.shape,flux.shape,fluxerr.shape,zpt.shape
    #raw_input()
    # ww = (catmag < 29.) & (rmsaddin < 1.) & (abs(ff) < 5*st)
    #
    # flux = flux[ww]
    # fluxerr = fluxerr[ww]
    # zpt = zpt[ww]
    # catmag = catmag[ww]
    # skyerr= skyerr[ww]
    # sky = sky[ww]
    # #rmsaddin = rmsaddin[ww]
    # poisson = poisson[ww]
    # indices = indices[ww]
    # ras = ras[ww]
    # decs = decs[ww]



    #print -2.5*np.log10(skyerr)+zpt
    #raw_input('skyerr in mags')
    #starmagerr = np.sqrt((-2.5*np.log10(sky)+2.5*np.log10(skyerr))**2+rmsaddin**2)
    #starmagerr = rmsaddin
    fluxerro = copy(fluxerr)
    #fluxerr = np.sqrt(fluxerr**2)
    catflux = 10**(.4*(np.array(zpt)-np.array(catmag)))

    #chisq = chisq[ww]
    #chisq = chisq[ww]*1/np.sqrt((abs(flux) / 3.))
    # plt.clf()
    # plt.scatter(catmag,(flux-catflux)/catflux)
    # plt.ylim(-.5,.5)
    # plt.savefig('starresid.png')
    # plt.clf()
    # plt.scatter(catmag,(flux-catflux)/fluxerr)
    # plt.ylim(-5,5)
    # plt.savefig('starstd.png')

    #fluxerr = np.sqrt(fluxerr**2 + )

    starmag = -2.5*np.log10(flux) + zpt






    #print 'fluxerr vs rmsadding' ,np.median((-2.5*np.log10(flux) + 2.5*np.log10(flux+fluxerr))), np.median(rmsaddin)
    #raw_input()
    #starmagerr2 = ((-2.5*np.log10(flux) + 2.5*np.log10(flux+fluxerr))**2 + rmsaddin**2 + (-2.5*np.log10(flux) + 2.5*np.log10(flux+poisson))**2 )**.5
    #starmagerr3 = ((-2.5*np.log10(sky) + 2.5*np.log10(sky+skyerr))**2 + rmsaddin[ww]**2)**.5
    skymagerr = -2.5*np.log10(sky) + 2.5*np.log10(sky+skyerr)
    #starmagerr = (-2.5*np.log10(flux) + 2.5*np.log10(flux+fluxerr+poisson)) + rmsaddin# + (-2.5*np.log10(flux) + 2.5*np.log10(flux+poisson))#+ rmsaddin #+ (-2.5*np.log10(flux) + 2.5*np.log10(flux+poisson))**2)**.5

    #starmagerrr = 1.0857*fluxerr/flux
    starmagerr = 1.0857*np.sqrt(fluxerr**2)/flux #+ rmsaddin
    starmagerr = np.sqrt(fluxerr**2)/flux #+ rmsaddin
    starmagerrorig = np.sqrt(fluxerrorig**2)/flux
    starmagerrzpt = np.sqrt(fluxerr**2)/flux #+ zptscat #+ rmsaddin

    rv = (flux-catflux)/(fluxerr**2)**.5
    rvorig = (flux-catflux)/((fluxerrorig)**2)**.5
    rvz = (flux-catflux)/((zptscat*flux)**2)**.5
    #print rv.shape,catmag.shape
    #raw_input()

    plt.clf()
    # repeatability = []
    # uindices = []
    # for ind in np.unique(indices):
    #     #print starmag[indices == ind]
    #     #raw_input()
    #     starobs = starmag[indices == ind]
    #     starobs = starobs[starobs-np.mean(starobs) <.1]
    #     starobs = sigmaclip(starobs)
    #     starmeanmag = np.mean(starobs)
    #     if len(starobs)>5.:
    #         repeatability.append(np.std(starobs))
    #         #print np.std(starobs)
    #         uindices.append(ind)

    # repeatability = np.array(repeatability)
    # uindices = np.array(uindices)

    plt.clf()
    plt.hist(fwhm,bins=np.arange(0,10,.1),color='black')
    plt.xlabel('FWHM (arcsec)')
    plt.ylabel('# of Obs')
    plt.savefig(outdir + '/' + title + 'fwhmhist.png')

    print outdir + '/' + title + 'fwhmhist.png'
    plt.clf()

    maxpoints = 5000000

    load = False
    if load:
        a = np.load(outdir + '/pltstarvec.npz')
        pltvecy = a['pltvecy']
        pltvecfield = a['pltvecfield']
        pltvecfwhm = a['pltvecfwhm']
        pltvecmjd = a['pltvecmjd']
        pltvecbigfield = a['pltvecbigfield']
        pltvecfieldr = a['pltvecfieldr']
        pltvecfieldb = a['pltvecfieldb']
        pltvecband = a['pltvecband']
        pltvecccd = a['pltvecccd']
        pltvecccdr = a['pltvecccdr']
        pltvecccdb = a['pltvecccdb']
    else:
        cntr = 0
        pltvecx = []
        pltvecy = []
        pltvecband = []
        pltvecfwhm = []
        pltvecmjd = []
        pltvecfield = []
        pltvecbigfield = []
        pltvecfieldr = []
        pltvecfieldb = []
        pltvecccd =[]
        pltvecccdr =[]
        pltvecccdb =[]


        stardictras = np.array([11,22])
        stardictdecs =  np.array([11,22])
        stardictcatmags =  np.array([11,22])
        stardictmeans = np.array([11,22])
        stardictlens = np.array([11,22])
        stardictww = np.array([])
        for sme, sm, ind, r, d, cm, f, fe, fh,tfield,tccd,tband,tmjd in zip(starmagerr[::-1], starmag[::-1],
                indices[::-1], ras[::-1], decs[::-1], catmag[::-1], flux[::-1], fluxerr[::-1], fwhm[::-1],
                field[::-1],ccd[::-1],band[::-1],mjd[::-1]):
            cntr += 1
            if cntr > maxpoints: continue
            if cntr > 100000: continue
            if cntr % 1 == 0: print cntr,'of',len(starmagerr[::-1])

            # print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
            # print starmag[indices == ind]
            # raw_input()

            firstww = np.where(np.isclose(stardictras, r, rtol=1.e-5) & np.isclose(stardictdecs, d, rtol=1.e-5) & (stardictcatmags == cm))
            if firstww[0] != []:
                starmean = stardictmeans[firstww[0]]
                starlen = stardictlens[firstww[0]]
                starww = stardictww[firstww[0]]
            else:
                starww = np.isclose(ras, r, rtol=1.e-5) & np.isclose(decs, d, rtol=1.e-5) & (catmag == cm)
                starwwmag = starmag[starww]
                starmean = np.mean(starwwmag)
                starlen = len(starww)
                stardictras = np.append(stardictras,r)
                stardictdecs = np.append(stardictdecs,d)
                stardictcatmags = np.append(stardictcatmags,cm)
                stardictmeans = np.append(stardictmeans,starmean)
                stardictlens = np.append(stardictlens,starlen)
                stardictww = np.append(stardictww,starww)

            #starmean = np.mean(starww)
            #repeatability = np.std(starww)
            # repeatability = np.std(starmag[indices == ind])
            if starlen > 5.:
                pltvecy.append(sm - starmean)
                pltvecband.append(tband)
                pltvecmjd.append(tmjd)
                pltvecfwhm.append(fh)
                pltvecbigfield.append(tfield)

                if len(np.unique(field[starww])) > 1:
                    #print tfield
                    pltvecfield.append(tfield)
                    pltvecfieldr.append(sm-starmean)
                    pltvecfieldb.append(tband)
                if len(np.unique(ccd[starww])) > 1:
                    pltvecccd.append(tccd)
                    pltvecccdr.append(sm-starmean)
                    pltvecccdb.append(tband)

        pltvecy = np.array(pltvecy)
        pltvecband = np.array(pltvecband,dtype='str')
        pltvecmjd = np.array(pltvecmjd,dtype='str')
        pltvecfwhm = np.array(pltvecfwhm,dtype='str')
        pltvecbigfield = np.array(pltvecbigfield,dtype='str')
        pltvecfield = np.array(pltvecfield,dtype='str')
        #print pltvecfield.shape
        #raw_input()
        pltvecfieldr = np.array(pltvecfieldr)

        pltvecfieldb = np.array(pltvecfieldb,dtype='str')
        pltvecccd = np.array(pltvecccd,dtype='str')
        pltvecccdr = np.array(pltvecccdr)
        pltvecccdb = np.array(pltvecccdb,dtype='str')

        np.savez(outdir +'/pltstarvec',pltvecy=pltvecy,pltvecfield=pltvecfield,pltvecfwhm=pltvecfwhm,pltvecmjd=pltvecmjd,
                 pltvecbigfield=pltvecbigfield,pltvecccd=pltvecccd,pltvecband=pltvecband,
                 pltvecfieldr=pltvecfieldr, pltvecccdr=pltvecccdr,pltvecfieldb=pltvecfieldb,pltvecccdb=pltvecccdb)

    plt.clf()
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 9))
    axes = axes.flatten()
    for i, b in enumerate(np.unique(pltvecband)):
        ax = axes[i]
        prms = np.sqrt(np.nanmean(np.square(pltvecy[(abs(pltvecy)<.05)&(pltvecband==b)])))
        ax.hist(pltvecy[pltvecband==b], alpha=.99, color='black',histtype='step',
                 bins=np.arange(-.05025,.05,.0005),label=b+' RMS:'+str(round(prms,4)))
        ax.set_xlabel(r'$'+b+' - '+b+'_{mean}$')

        ax.legend()
        ax.set_xlim(-.05,.05)
    #plt.ylim(-.02, .01)

    # ax, ay, aystd = dt.bindata(np.array(pltvecx), np.array(pltvecy), np.arange(2, 10, .1), window=.1,
    #                            dontrootn=True)
    # plt.plot(ax, ay, linewidth=3, color='orange', label='SMP', alpha=.6)
    # plt.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP', alpha=.6)
    # plt.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP', alpha=.6)

    #plt.title(title.split('_')[0]+ ' '+title.split('_')[1]+' band')

    print 'finished snls'
    #plt.plot([2, 10], [0, 0], color='black')
    plt.savefig(outdir + '/repeatabilitylikesnls.png')
    print outdir + '/' + title + '_repeatabilitylikesnls.png'

    plt.clf()

    fig, axes = plt.subplots(figsize=(12/1.5, 9/1.5))
    for i, b, c in zip(np.arange(4),['g','r','i','z'],['green','red','indigo','black']):
        ww = pltvecccdb==b
        for j,chip in enumerate(np.sort(np.unique(pltvecccd))):
            yy = (pltvecccd == chip)
            try:
                mean = np.mean(pltvecccdr[(abs(pltvecccdr)<.05)& (ww) & (yy)])
                rms = np.sqrt(np.nanmean(np.square(pltvecccdr[(abs(pltvecccdr)<.05) & ww & yy])))
                print chip,b,mean,len(pltvecccdr[(abs(pltvecccdr)<.05)& (ww) & (yy)])
                if len(pltvecccdr[(abs(pltvecccdr)<.05)& (ww) & (yy)]) > 10:
                    if j == 0:
                        plt.errorbar([int(chip)], [mean], yerr=[rms], fmt='o', mew=0, c=c,alpha=.7,
                                    label=b+' band')
                    else:
                        plt.errorbar([int(chip)], [mean], yerr=[rms], fmt='o', mew=0, c=c,alpha=.7)
            except:
                pass
    plt.ylabel('Chip Mean - All Mean')
    plt.xlabel('Chip')
    plt.legend(fontsize='small')
    plt.axhline(0,color='grey',linestyle='--')
    plt.savefig(outdir+'/chipdependence.png')
    plt.clf()
    fielddict = {}
    cnt = 0
    for f in np.sort(np.unique(pltvecfield)):
        fielddict[f] = cnt
        cnt += 1

    print np.unique(pltvecbigfield)

    fig, ax = plt.subplots(figsize=(12/1.5, 9/1.5))
    for i, b, c in zip(np.arange(4),['g','r','i','z'],['green','red','indigo','black']):
        ww = pltvecfieldb==b
        for j,tf in enumerate(np.sort(np.unique(pltvecfield))):
            yy = pltvecfield == tf
            #print ww.shape,yy.shape,pltvecfield.shape
            if True:
            #try:
                mean = np.mean(pltvecfieldr[(abs(pltvecfieldr)<.05)& (ww) & (yy)])
                rms = np.sqrt(np.nanmean(np.square(pltvecfieldr[(abs(pltvecfieldr)<.05) & ww & yy])))
                #print b, field

                if j == 0:
                    ax.errorbar([fielddict[tf]], [mean], yerr=[rms], fmt='o', mew=0, c=c,alpha=.7,
                                label=b+' band')
                else:
                    ax.errorbar([fielddict[tf]], [mean], yerr=[rms], fmt='o', mew=0, c=c,alpha=.7)
            #except:
            #    pass
    #ax.set_xlabel('Field')
    ax.set_ylabel('Field Mean - All Mean')
    ax.legend(fontsize='small')
    ax.axhline(0,color='grey',linestyle='--')
    ax.set_xlim(-1,10)
    plt.xticks(np.arange(0,10))
    ax.set_xticklabels(np.sort(np.unique(pltvecfield)))
    plt.savefig(outdir+'/fielddependence.png')

    sys.exit()
    plt.clf()
    cntr = 0
    pltvecx = []
    pltvecy = []
    for sme, sm, ind, r, d, cm, f, fe,fh in zip(starmagerr, starmag, indices, ras, decs, catmag, flux, fluxerr, fwhm):
        cntr += 1
        if cntr > maxpoints: continue
        if cntr > 500: continue

        # print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        # print starmag[indices == ind]
        # raw_input()
        starww = starmag[np.isclose(ras, r, rtol=1.e-5) & np.isclose(decs, d, rtol=1.e-5) & (catmag == cm)]
        repeatability = np.std(starww)
        # repeatability = np.std(starmag[indices == ind])
        if len(starww) > 5.:
            # if repeatability < .3:

            pltvecy.append(sme - repeatability)
            pltvecx.append(fh)

    #plt.xscale('log')
    plt.scatter(pltvecx,pltvecy, alpha=.3, color='black')
    plt.xlabel('PSF FWHM')
    plt.ylabel('PhotErr - Repeatability')
    plt.xlim(3.5, 7)
    plt.ylim(-.02, .01)

    ax, ay, aystd = dt.bindata(np.array(pltvecx), np.array(pltvecy), np.arange(2, 10, .1), window=.1,
                               dontrootn=True)
    plt.plot(ax, ay, linewidth=3, color='orange', label='SMP', alpha=.6)
    plt.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP', alpha=.6)
    plt.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP', alpha=.6)

    plt.title(title + 'BAND')
    print 'finished fwhm'
    plt.plot([2, 10], [0, 0], color='black')
    plt.savefig(outdir + '/' + title + '_repeatability_vs_fwhm.png')

    plt.clf()
    cntr = 0
    pltvecx = []
    pltvecy = []
    for sme,sm,ind,r,d,cm,f,fe in zip(starmagerr,starmag,indices,ras,decs,catmag,flux,fluxerr):
        cntr+=1
        if cntr > maxpoints: continue
        if cntr > 500: continue

        #print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        #print starmag[indices == ind]
        #raw_input()
        starww = starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        repeatability = np.std(starww)
        #repeatability = np.std(starmag[indices == ind])
        if len(starww) > 5.:
            #if repeatability < .3:

            pltvecy.append(sme-repeatability)
            pltvecx.append(sme)

    plt.scatter(sme, sme - repeatability, alpha=.3, color='black')
    plt.xscale('log')
    plt.xlabel('Photometric Error')
    plt.ylabel('Repeatability - PhotErr')
    plt.xlim(.0005,.05)
    plt.ylim(-.02,.01)

    ax, ay, aystd = dt.bindata(np.array(pltvecx),np.array(pltvecy), np.arange(.0005,.05, .00001), window=.002,dontrootn=True)
    plt.plot(ax, ay, linewidth=3, color='orange', label='SMP',alpha=.6)
    plt.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP',alpha=.6)
    plt.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP',alpha=.6)

    plt.title(title+'BAND')
    print 'finished resid'

    plt.plot([.0003,.02],[0,0],color='black')
    plt.savefig(outdir+'/'+title+'_repeatability-photerr_vs_photerrr.png')

    plt.clf()
    # cntr = 0
    # pltvecx = []
    # pltvecy = []
    # for sme, sm, ind, r, d, cm, cs, se, rai in zip(starmagerr, starmag, indices, ras, decs, catmag, chisq, skyerr,
    #                                                rmsaddin):
    #     cntr += 1
    #     if cntr > maxpoints: continue
    #     if cntr > 10000: continue
    #
    #     # print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
    #     # print starmag[indices == ind]
    #     # raw_input()
    #     starww = starmag[np.isclose(ras, r, rtol=1.e-5) & np.isclose(decs, d, rtol=1.e-5) & (catmag == cm)]
    #     repeatability = np.std(starww)
    #     # repeatability = np.std(starmag[indices == ind])
    #     if len(starww) > 4.:
    #         # if repeatability < .3:
    #         plt.scatter(rai, float(sme) - repeatability, alpha=.3, color='black')
    #         pltvecy.append(sme - repeatability)
    #         pltvecx.append(rai)
    #
    # plt.xscale('log')
    # plt.xlim(0.003, .05)
    # plt.ylim(-.02, .01)
    #
    # ax, ay, aystd = dt.bindata(np.array(pltvecx), np.array(pltvecy), np.arange(0., .1, .00005), window=.0001,
    #                            dontrootn=True)
    # plt.plot(ax, ay, linewidth=3, color='orange', label='SMP', alpha=.6)
    # plt.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP', alpha=.6)
    # plt.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP', alpha=.6)
    #
    # plt.xlabel('ZPT Uncertainty')
    # plt.ylabel('PhotErr - Repeatability')
    # plt.title(title+'BAND')
    # plt.plot([0., 1000], [0, 0], color='black')
    # plt.plot([0., 1000], [0, 1000], color='black', linestyle='--')
    #
    # plt.savefig(outdir + '/' + title + '_repeatability_vs_zptuncertainty.png')

    plt.clf()
    print 'saved zptu'

    cntr = 0
    pltvecx = []
    pltvecy = []
    for sme,sm,ind,r,d,cm,f,fe in zip(starmagerr,starmag,indices,ras,decs,catmag,flux,fluxerr):
        cntr+=1
        if cntr > maxpoints: continue
        if cntr > 500: continue

        #print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        #print starmag[indices == ind]
        #raw_input()
        starww = starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        repeatability = np.std(starww)#/np.sqrt(len(starww))
        #repeatability = np.std(starmag[indices == ind])
        if len(starww) > 5.:
            #if repeatability < .3:
            plt.scatter(sme,repeatability,alpha=.3,color='black')
            pltvecy.append(repeatability)
            pltvecx.append(sme)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Photometric Error')
    plt.ylabel('Repeatability')
    plt.xlim(.002,.05)
    plt.ylim(.002,.05)

    ax, ay, aystd = dt.bindata(np.array(pltvecx),np.array(pltvecy), np.arange(.003,.05, .0001), window=.0005,dontrootn=True)
    photerr = copy(ax)
    repeaterr = copy(ay)
    plt.plot(ax, ay, linewidth=3, color='orange', label='SMP',alpha=.6)
    plt.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP',alpha=.6)
    plt.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP',alpha=.6)

    plt.title(title+'BAND')

    plt.plot([min(starmagerr),max(starmagerr)],[min(starmagerr),max(starmagerr)],color='grey')
    plt.savefig(outdir+'/'+title+'_repeatability_vs_photerr.png')

    #sys.exit()

    print 'saved repeat vs photerr1'
    plt.clf()

    cntr = 0
    pltvecx = []
    pltvecy = []
    for sme, sm, ind, r, d, cm, f, fe,cf in zip(starmagerrzpt, starmag, indices, ras, decs, catmag, flux, fluxerr,catflux):
        cntr += 1
        if cntr > maxpoints: continue
        if cntr > 500: continue

        # print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        # print starmag[indices == ind]
        # raw_input()
        starww = starmag[np.isclose(ras, r, rtol=1.e-5) & np.isclose(decs, d, rtol=1.e-5) & (catmag == cm)]
        repeatability = np.std(starww)
        # repeatability = np.std(starmag[indices == ind])
        if len(starww) > 5.:
            # if repeatability < .3:
            #plt.scatter(sme, repeatability, alpha=.3, color='black')
            pltvecy.append(repeatability)
            pltvecx.append(sme)
    plt.scatter(pltvecx,pltvecy, alpha=.3, color='black')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Photometric Error + Zpt Scatter')
    plt.ylabel('Repeatability')
    plt.xlim(.0005, .05)
    plt.ylim(.0005, .05)

    ax, ay, aystd = dt.bindata(np.array(pltvecx), np.array(pltvecy), np.arange(.0005, .05, .0001), window=.0002,
                               dontrootn=True)
    photerr = copy(ax)
    repeaterr = copy(ay)
    plt.plot(ax, ay, linewidth=3, color='orange', label='SMP', alpha=.6)
    plt.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP', alpha=.6)
    plt.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP', alpha=.6)

    plt.title(title + 'BAND')

    plt.plot([min(starmagerr), max(starmagerr)], [min(starmagerr), max(starmagerr)], color='black')
    plt.savefig(outdir + '/' + title + '_repeatability_vs_photerrzpt.png')
    print 'saved',outdir + '/' + title + '_repeatability_vs_photerrzpt.png'
    plt.clf()



    cntr = 0
    pltvecx = []
    pltvecy = []
    rep = copy(starmagerr)
    for i, sme, sm, ind, r, d, cm in zip(range(len(starmagerr)),starmagerr, starmag, indices, ras, decs, catmag):
        cntr += 1
        if cntr > maxpoints: continue
        if cntr > 500: continue
        # print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        # print starmag[indices == ind]
        # raw_input()
        starww = starmag[np.isclose(ras, r, rtol=1.e-5) & np.isclose(decs, d, rtol=1.e-5) & (catmag == cm)]
        repeatability = np.std(starww)
        rep[i] = repeatability#*np.sqrt(len(starww))
        # repeatability = np.std(starmag[indices == ind])
        if len(starww) > 4.:
            # if repeatability < .3:
            plt.scatter(cm, float(sme) - repeatability, alpha=.3, color='black')
            pltvecy.append(sme - repeatability)
            pltvecx.append(cm)
    #plt.yscale('log')
    plt.xlim(15.,21.5)
    plt.ylim(-.02, .01)
    plt.title(title+'BAND')
    plt.xlabel('Catalog Magnitude')
    plt.ylabel('PhotErr - Repeatability')

    ax, ay, aystd = dt.bindata(np.array(pltvecx), np.array(pltvecy), np.arange(15, 22, .1), window=.5,dontrootn=True)
    plt.plot(ax, ay, linewidth=3, color='orange', label='SMP',alpha=.6)
    plt.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP',alpha=.6)
    plt.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP',alpha=.6)

    plt.plot([15., 21.5], [0,0], color='black')
    plt.savefig(outdir + '/' + title + '_repeatability_vs_catmag.png')

    plt.clf()

    cntr = 0
    pltvecx = []
    pltvecy = []
    for sme, sm, ind, r, d, cm,cs in zip(starmagerr, starmag, indices, ras, decs, catmag,chisq):
        cntr += 1
        if cntr > maxpoints: continue
        if cntr > 500: continue
        # print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        # print starmag[indices == ind]
        # raw_input()
        starww = starmag[np.isclose(ras, r, rtol=1.e-5) & np.isclose(decs, d, rtol=1.e-5) & (catmag == cm)]
        repeatability = np.std(starww)
        # repeatability = np.std(starmag[indices == ind])
        if len(starww) > 4.:
            # if repeatability < .3:
            plt.scatter(cs, float(sme) - repeatability, alpha=.3, color='black')
            pltvecy.append(sme - repeatability)
            pltvecx.append(cs)

    plt.xscale('log')
    plt.xlim(0.5, 1000.)
    plt.ylim(-.02, .01)

    ax, ay, aystd = dt.bindata(np.array(pltvecx), np.array(pltvecy), np.arange(.5, 100, .1), window=1.,dontrootn=True)
    plt.plot(ax, ay, linewidth=3, color='orange', label='SMP',alpha=.6)
    plt.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP',alpha=.6)
    plt.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP',alpha=.6)

    plt.xlabel('Chi Squared')
    plt.ylabel('PhotErr - Repeatability')
    plt.plot([0., 1000], [0, 0], color='black')
    plt.title(title+'BAND')
    plt.savefig(outdir + '/' + title + '_repeatability_vs_chisq.png')

    plt.clf()
    cntr = 0
    pltvecx = []
    pltvecy = []
    for sme, sm, ind, r, d, cm, cs,se in zip(starmagerr, starmag, indices, ras, decs, catmag, chisq,skyerr):
        cntr += 1
        if cntr > maxpoints: continue
        if cntr > 500: continue
        # print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        # print starmag[indices == ind]
        # raw_input()
        starww = starmag[np.isclose(ras, r, rtol=1.e-5) & np.isclose(decs, d, rtol=1.e-5) & (catmag == cm)]
        repeatability = np.std(starww)
        # repeatability = np.std(starmag[indices == ind])
        if len(starww) > 4.:
            # if repeatability < .3:
            plt.scatter(se, float(sme) - repeatability, alpha=.3, color='black')
            pltvecy.append(sme - repeatability)
            pltvecx.append(se)

    plt.xscale('log')
    plt.xlim(10, 1000)
    plt.ylim(-.02, .01)

    ax, ay, aystd = dt.bindata(np.array(pltvecx), np.array(pltvecy), np.arange(5, 500, .5), window=1., dontrootn=True)
    plt.plot(ax, ay, linewidth=3, color='orange', label='SMP',alpha=.6)
    plt.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP',alpha=.6)
    plt.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP',alpha=.6)

    plt.title(title+'BAND')
    plt.xlabel('Sky Error')
    plt.ylabel('PhotErr - Repeatability')
    plt.plot([0., 1000], [0, 0], color='black')
    plt.savefig(outdir + '/' + title + '_repeatability_vs_skyerr.png')




    print starmag[0:10]
    print catmag[0:10]

    # from scipy.interpolate import interp1d
    # f = interp1d(photerr, repeaterr)
    # starmagerrinterp = copy(starmagerr)
    # for i,sme in enumerate(starmagerr):
    #     if sme <= min(photerr):
    #         sme = min(photerr)+.0001
    #     if sme >= max(photerr):
    #         sme = max(photerr)-.0001
    #     starmagerrinterp[i] = f(sme)
    # #starmagerrinterp = f(starmagerr)

    dmz = (starmag - catmag) / starmagerr
    print starmagerr[:50]
    print rep[:50]
    #raw_input('ccc')
    dmam = (starmag - catmag) / np.maximum(starmagerr,rep)
    #dmam = (starmag - catmag) / rep
    #dmas = (starmag - catmag) / starmagerr3
    dsss = (starmag - catmag) / skymagerr


    #raw_input('printing mags')
    #dmam = (flux - catflux) / np.sqrt(fluxerr**2 + flux)
    ds = (flux - catflux) / skyerr
    dp = (flux-catflux) / poisson

    #chisq = (flux - catflux) ** 2 / catflux

    plt.clf()
    plt.scatter(catmag,chisq)
    plt.ylim(0,15)
    plt.xlim(17,21)
    plt.axhline(1)
    plt.xlabel('Cat Mag')
    plt.ylabel('Chi Sq')
    plt.savefig(outdir + '/' + title + 'chivscat.png')
    #chisq = np.nanmean(chisq[abs(d) < 3])

    plt.clf()

    #dc = dmam[abs(dmam) < 3]
    #rmsscat = r

    rms = np.sqrt(np.nanmean(np.square(rv[abs(rv) < 3.])))


    plt.clf()

    nullfmt = NullFormatter()  # no labels

    # definitions for the axes
    left, width = 0.13, 0.61
    bottom, height = 0.1, 0.6
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom + height / 2., width+.2, height / 2.]
    rect_scatterflux = [left, bottom, width+.2, height / 2.]
    rect_histx = [left, bottom_h - .04, width+.2, .18]
    rect_histy = [left_h, bottom + height / 2., 0.2, height / 2.]
    rect_histyflux = [left_h, bottom, 0.2, height / 2.]

    # start with a rectangular Figure
    plt.figure(1, figsize=(25, 20))

    ax1 = plt.axes(rect_scatter)
    ax3 = plt.axes(rect_histx)
    #ax2 = plt.axes(rect_histy)
    ax4 = plt.axes(rect_scatterflux)
    #ax5 = plt.axes(rect_histyflux)

    # no labels
    #ax2.yaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    #ax5.yaxis.set_major_formatter(nullfmt)

    #ax2.hist(rv[np.isfinite(rv)], bins=np.arange(-10.1, 10, .2), normed=True, label='RMS: ' + str(round(rms, 3))
    #         , orientation='horizontal',color='black')
    # label='RMS: ' + str(round(rms, 3)) + '\nChiSq (3sig cut) ' + str(round(chisq, 3)) + '\nMedian ' + str(
    #   round(np.median(d), 3)) + ' +- ' + str(round(np.std(d), 3)),

    # import matplotlib.mlab as mlab
    # import math
    # mean = 0
    # variance = 1
    # sigma = math.sqrt(variance)
    # x = np.arange(-5, 5, .1)
    # ax2.plot(mlab.normpdf(x, mean, sigma), x, color='grey', label='Gaussian Normal')
    #
    # ax2.set_ylim(-4, 4)
    # ax2.set_xlim(0, .5)
    # .xlabel('STDEV')
    # plt.ylabel('Normalized Count')
    # ax2.legend(fontsize='x-small',loc = 'center right', bbox_to_anchor = (1.1, 1.15))
    # plt.savefig('stdresid.png')

    # plt.clf()

    #ax1.scatter(catmag, rv, alpha=.001, color='black')
    ax, ay, aystd = dt.bindata(catmag, rv, np.arange(min(catmag), max(catmag), .1), window=.3,dontrootn=True)
    ax1.plot([min(catmag), max(catmag)], [0, 0], color='grey')
    ax1.plot(ax, ay, linewidth=3, color='black', label='SMP',alpha=.8)
    ax1.plot(ax, ay + aystd, linewidth=2, color='black', linestyle='--', label='SMP',alpha=.8)
    ax1.plot(ax, ay - aystd, linewidth=2, color='black', linestyle='--', label='SMP',alpha=.8)

    # ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

    ax1.set_xlim(16.5, max(catmag))
    ax1.set_ylim(-4., 4.)
    ax1.set_xlabel('Cat Mag')
    ax1.set_ylabel('STD')

    # ax, ayrms= dt.binrms(fakemag, d, np.arange(19.5, max(fakemag), .1),.5)
    # ax3.plot(ax, ayrms, color='blue',label='RMS',linewidth=3)


    ax3.plot([0, 100], [1., 1.], linestyle='--', color='grey')
    ax3.set_ylim(0., 3.)
    ax3.legend(fontsize='small')

    fresid = np.zeros(flux.shape)
    for i, f, ff in zip(range(len(flux)), flux, catflux):
        if f == 0.:
            fresid[i] = np.nan
        else:
            fresid[i] = (f - ff) / max([abs(ff), 1.])
    # fresid[abs(fakeflux) < 1.] = flux[abs(fakeflux) < 1.] - fakeflux[abs(fakeflux) < 1.]

    #ax5.hist(fresid, bins=np.arange(-.155, .15, .001), color='black', orientation='horizontal')
    ww = abs(fresid) < .1
    ax4.scatter(catmag[ww], fresid[ww], alpha=.001, color='black')
    ax, ay, aystd = dt.bindata(catmag[ww], fresid[ww],
                               np.arange(16., max(catmag[ww]), .1), window=1.,dontrootn=True)
    ax4.plot([19, 28.7], [0, 0], color='grey')

    #ax, ayrms = dt.binrms(catmag, d, np.arange(16., max(catmag), .1), .5)
    #ax3.plot(ax, ayrms, color='blue', label='Chisq Min Err', linewidth=3,alpha=.7)
    # ax, ayrms = dt.binrms(catmag, dsss, np.arange(16., max(catmag), .1), .5)
    # ax3.plot(ax, ayrms, color='green', label='Skyerr', linewidth=3,alpha=.4)
    #ax, ayrms = dt.binrms(catmag, dmz, np.arange(16., max(catmag), .1), .5)
    # print ayrms
    # raw_input('zpt scatter err')
    #ax3.plot(ax, ayrms, color='red', label='ZPT Scatter Err', linewidth=3,alpha=.7)

    #ax, ayrms = dt.binrms(catmag, dp, np.arange(16., max(catmag), .1), .5)
    #ax3.plot(ax, ayrms, color='green', label='Poisson Err', linewidth=3,alpha=.7)

    # ax, ayrms = dt.binrms(catmag, dmas, np.arange(16., max(catmag), .1), .1)
    # ax3.plot(ax, ayrms, color='orange', label='ZPT Scatter Err and Sky Err', linewidth=3,alpha=.4)
    ax, ayrms = dt.binrms(catmag, rv, np.arange(16., max(catmag), .1), .1)
    ax3.plot(ax, ayrms, color='black', label='PSF Fit + ZPT Scatter', linewidth=3,alpha=.8)
    ax, ayrms = dt.binrms(catmag, rvorig, np.arange(16., max(catmag), .1), .1)
    ax3.plot(ax, ayrms, color='green', label='PSF Fit', linewidth=3, alpha=.8)
    ax, ayrms = dt.binrms(catmag, rvz, np.arange(16., max(catmag), .1), .1)
    ax3.plot(ax, ayrms, color='blue', label='ZPT Scatter', linewidth=3, alpha=.8)
    ax3.plot(ax, ax * 0 + 1., linestyle='--', color='grey')
    ax3.legend(fontsize='xx-small',loc = 'upper left')#, bbox_to_anchor = (1.4, 0.8))
    # ww = hostmag > 25.
    # ax, ayrms = dt.binrms(catmag[ww], d[ww], np.arange(19.5, max(fakemag), .1), .5)
    # ax3.plot(ax, ayrms, color='red', label='HostMag > 25.', linewidth=3)
    #
    # ww = hostmag < 23.
    # ax, ayrms = dt.binrms(fakemag[ww], d[ww], np.arange(19.5, max(fakemag), .1), .5)
    # ax3.plot(ax, ayrms, color='green', label='HostMag < 23', linewidth=3)
    # ax3.legend(fontsize='small')

    ax4.plot(ax, ay, linewidth=3, color='black',alpha=.8)
    ax4.plot(ax, ay + aystd, linewidth=2, color='black', linestyle='--',alpha=.8)
    ax4.plot(ax, ay - aystd, linewidth=2, color='black', linestyle='--',alpha=.8)
    ax4.set_xlim(ax1.get_xlim())
    ax4.set_ylim(-.025, .025)
    ax4.set_xlabel('Catalog Mag')
    #ax5.set_xlabel('Counts')
    ax3.set_ylabel('RMS')
    ax3.set_title(filter+' band')

    ax4.set_ylabel(r'$\frac{fitflux - catalog flux}{catalog flux}$')

    ax3.set_xlim(ax1.get_xlim())
    ax3.set_ylim(.5,2.)
    #ax2.set_ylim(ax1.get_ylim())
    #ax5.set_ylim(ax4.get_ylim())
    #ax2.xaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax1.xaxis.set_major_formatter(nullfmt)

    plt.subplots_adjust(wspace=0.001, hspace=0.001)
    plt.savefig(outdir + '/' + title +'starstd.png')
    print 'saved',outdir + '/' + title +'starstd.png'


    #------------------------------------------------------------------------------------------------





    plt.clf()

    dc = d[abs(d) < 3]

    rms = np.sqrt(np.nanmean(np.square(dc[abs(dc) < 3.])))

    plt.clf()

    nullfmt = NullFormatter()  # no labels

    # definitions for the axes
    # definitions for the axes
    left, width = 0.13, 0.61
    bottom, height = 0.1, 0.6
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom + height / 2., width, height / 2.]
    rect_scatterflux = [left, bottom, width, height / 2.]
    rect_histx = [left, bottom_h - .04, width, .18]
    rect_histy = [left_h, bottom + height / 2., 0.2, height / 2.]
    rect_histyflux = [left_h, bottom, 0.2, height / 2.]
    # start with a rectangular Figure
    plt.figure(1, figsize=(25, 20))

    ax1 = plt.axes(rect_scatter)
    ax3 = plt.axes(rect_histx)
    ax2 = plt.axes(rect_histy)
    ax4 = plt.axes(rect_scatterflux)
    ax5 = plt.axes(rect_histyflux)

    # no labels
    ax2.yaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax5.yaxis.set_major_formatter(nullfmt)

    ax2.hist(d, bins=np.arange(-10, 10, .25), normed=True, label='RMS: ' + str(round(rms, 3))
             , orientation='horizontal',color='black')
    # label='RMS: ' + str(round(rms, 3)) + '\nChiSq (3sig cut) ' + str(round(chisq, 3)) + '\nMedian ' + str(
    #   round(np.median(d), 3)) + ' +- ' + str(round(np.std(d), 3)),

    import matplotlib.mlab as mlab
    import math
    mean = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.arange(-5, 5, .1)
    ax2.plot(mlab.normpdf(x, mean, sigma), x, color='red', label='Gaussian Normal')

    ax2.set_ylim(-4, 4)
    ax2.set_xlim(0, .5)
    # .xlabel('STDEV')
    # plt.ylabel('Normalized Count')
    ax2.legend(fontsize='small')
    # plt.savefig('stdresid.png')

    # plt.clf()

    ax1.scatter(chisq, d, alpha=.3, color='blue')
    ax, ay, aystd = dt.bindata(chisq, d, np.arange(.1, 5, .01), window=.1)
    ax1.plot([.1, 5.], [0, 0], color='grey')
    ax1.plot(ax, ay, linewidth=3, color='orange', label='SMP')
    ax1.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP')
    ax1.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP')

    # ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

    ax1.set_xlim(.1, 5.)
    ax1.set_ylim(-3., 3.)
    ax1.set_xlabel('Chisq')
    ax1.set_ylabel('STD')

    # ax, ayrms= dt.binrms(fakemag, d, np.arange(19.5, max(fakemag), .1),.5)
    # ax3.plot(ax, ayrms, color='blue',label='RMS',linewidth=3)


    ax3.plot([0, 100], [1., 1.], linestyle='--', color='black')
    ax3.set_ylim(.7, 1.7)
    ax3.legend(fontsize='small')

    fresid = np.zeros(flux.shape)
    for i, f, ff in zip(range(len(flux)), flux, catflux):
        if f == 0.:
            fresid[i] = np.nan
        else:
            fresid[i] = (f - ff) / max([abs(ff), 1.])
    # fresid[abs(fakeflux) < 1.] = flux[abs(fakeflux) < 1.] - fakeflux[abs(fakeflux) < 1.]

    ax5.hist(fresid, bins=np.arange(-.155, .15, .01), color='blue', orientation='horizontal')

    ax4.scatter(chisq, fresid, alpha=.3, color='blue')
    ax, ay, aystd = dt.bindata(chisq, fresid,
                               np.arange(.1, 5, .01), window=.1)
    ax4.plot([.1, 5], [0, 0], color='grey')

    ax, ayrms = dt.binrms(chisq, d, np.arange(.1, 5., .01), .1)
    ax3.plot(ax, ayrms, color='blue', label='ALL Standard Stars in r Band', linewidth=3)
    ax3.plot(ax, ax * 0 + 1., linestyle='--', color='black')

    # ww = hostmag > 25.
    # ax, ayrms = dt.binrms(catmag[ww], d[ww], np.arange(19.5, max(fakemag), .1), .5)
    # ax3.plot(ax, ayrms, color='red', label='HostMag > 25.', linewidth=3)
    #
    # ww = hostmag < 23.
    # ax, ayrms = dt.binrms(fakemag[ww], d[ww], np.arange(19.5, max(fakemag), .1), .5)
    # ax3.plot(ax, ayrms, color='green', label='HostMag < 23', linewidth=3)
    # ax3.legend(fontsize='small')

    ax4.plot(ax, ay, linewidth=3, color='orange')
    ax4.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--')
    ax4.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--')
    ax4.set_xlim(ax1.get_xlim())
    ax4.set_ylim(-.025, .025)
    ax4.set_xlabel('Chisq')
    ax5.set_xlabel('Counts')
    ax3.set_ylabel('RMS')
    ax4.set_ylabel('(fitflux - catflux)/catflux')

    ax3.set_xlim(ax1.get_xlim())
    ax2.set_ylim(ax1.get_ylim())
    ax5.set_ylim(ax4.get_ylim())
    ax2.xaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax1.xaxis.set_major_formatter(nullfmt)

    plt.subplots_adjust(wspace=0.001, hspace=0.001)
    plt.savefig(outdir + '/' + title +'starchi.png')

    print 'saved starstd.png'

def plotstarlc(flux,fluxerr,zpt,ids,mjd,catmag):
    from matplotlib.backends.backend_pdf import PdfPages
    pdf_pages = PdfPages('allstarlc.pdf')

    plt.clf()
    #plt.figure(figsize=(20,20))
    fig, axs = plt.subplots(nrows=4, ncols=4, figsize=(30,25))
    pagescounter = 0
    icntr = 0
    justpaged = False
    for i,id in enumerate(np.unique(ids)[:]):
        if pagescounter >50: break
        print i,

        ww = ids == id

        cm = catmag[ww]
        for c in np.unique(cm):
            ww = (ids == id) & (catmag == c)
            if len(flux[ww]) < 10:
                #print 'doesnt pass'
                continue

            if icntr % 16 == 0:
                if not justpaged:
                    print 'saving page',pagescounter
                    if pagescounter > 0:
                        pdf_pages.savefig(fig)
                    plt.clf()
                    fig, axs = plt.subplots(nrows=4, ncols=4, figsize=(30,25))
                    pagescounter += 1
                    justpaged = True


            justpaged = False
            #print flux[ww]*10**(.4*(31-zpt[ww]))
            tm = zpt[ww] - 2.5*np.log10(flux[ww])
            axs.ravel()[int(icntr%16)].scatter(np.array(mjd[ww],dtype='float'),tm - np.mean(tm),color='black')
            axs.ravel()[int(icntr%16)].set_ylabel('Fit Mag - Mean')
            axs.ravel()[int(icntr%16)].set_xlabel('mjd')
            icntr += 1
            #print 'saved'
        #axs.ravel()[int(i%16)].errorbar(np.array(mjd[ww],dtype='float'),flux[ww]*10**(-.4*(31.-zpt[ww])),yerr=fluxerr[ww]*10**(-.4*(31-zpt[ww])),fmt='o',color='black')
    pdf_pages.close()

    #plt.savefig('allstarlc.png')
    print 'saved allstarlc.pdf'

def bindata(x, y, bins, returnn=False):
    medians = np.zeros(len(bins) - 1)
    mads = np.zeros(len(bins) - 1)
    nums = np.zeros(len(bins) - 1)

    for i in np.arange(len(bins) - 1):
        bs = bins[i]
        bf = bins[i + 1]
        ww = [(x > bs) & (x < bf)]
        yhere = y[ww]
        yhere = yhere[np.isfinite(yhere)]
        ss = [abs(yhere) < 3. * np.std(yhere)]
        try:
            nums[i] = len(yhere[ss])
            medians[i] = np.median(yhere[ss])
            mads[i] = 1.48 * np.median(abs(yhere[ss] - medians[i])) * 1 / np.sqrt(len(yhere[ss]))
        except IndexError:
            print 'excepted'
            nums[i] = 0.
            medians[i] = np.nan
            mads[i] = np.nan
    xvals = (bins[1:] + bins[:-1]) / 2.
    if returnn:
        return xvals, medians, mads, nums
    return xvals, medians, mads


if __name__ == "__main__":
    fakedir = '/project/projectdirs/des/djbrout/pysmp/imglist/all/'
    #fakedir = '/project/projectdirs/des/djbrout/pysmp/imglist/spec/'

    #resultsdir = '/pnfs/des/scratch/pysmp/smp_04_modelerrors'
    #resultsdir = '/pnfs/des/scratch/pysmp/smp_02_simnosnnoskyerr'
    resultsdir = '/project/projectdirs/dessn/dbrout/allsim/'
    #resultsdir = '/project/projectdirs/des/djbrout/simtest/'
    #resultsdir= './working/'
    #resultsdir= '/export/scratch0/ps1sn1/data/v10.0/GPC1v3/eventsv1/smpworkspace/PS_TEST1/'
    #resultsdir = './workingsimnosn'
    isfermigrid = False
    cacheddata = False

    deep_or_shallow = 'shallow'

    #cd = '/project/projectdirs/des/djbrout/v67pixshift//summary_results.npz'
    #cd = '/pnfs/des/scratch/pysmp/smp_02_simnosnnoskyerr/np_data/summary_results.npz'
    import sys, getopt

    try:
        args = sys.argv[1:]

        opt, arg = getopt.getopt(
            args, "fd:rd:cd:cdf:b",
            longopts=["fakedir=", "resultsdir=", "cacheddata", "cashedfile=","filter=",'dostars'])

    except getopt.GetoptError as err:
        print "No command line arguments"

    filter = 'r'
    dostars = False
    for o,a in opt:
        if o in ["-h","--help"]:
            print __doc__
            sys.exit(0)
        elif o in ["-fd","--fakedir"]:
            fakedir = a
        elif o in ["-rd","--resultsdir"]:
            resultsdir = a
        elif o in ["-cd","--cacheddata"]:
            cacheddata = True
        elif o in ["-cdf","--cachedfile"]:
            cd = a
        elif o in ["-b","--filter"]:
            filter = str(a)
        elif o in ['--dostars']:
            dostars = True
    print filter
    tfield = 'SN-S2'
    #cd = '/global/cscratch1/sd/dbrout/v7/summary_results_'+tfield+'_'+filter+'.npz'
    if filter == 'all':
        cd = []
        for filt in ['g','r','i','z']:
            cd.append('/global/cscratch1/sd/dbrout/summary_results_'+deep_or_shallow+'_' + filt + '.npz')
    else:
        cd = ['/global/cscratch1/sd/dbrout/summary_results_'+deep_or_shallow+'_' + filter + '.npz']

    #print cd
    #raw_input()
    go(fakedir,resultsdir,cacheddata,cd,filter,tfield,dostars,deep_or_shallow,real=False)