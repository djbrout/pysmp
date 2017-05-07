
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

import dilltools as dt
from copy import copy
from scipy.stats import sigmaclip


def go(fakedir,resultsdir,cacheddata,cd,filter,tfield,isfermigrid=False):

    if isfermigrid:
        useifdh = True
    else:
        useifdh = False
    tmpwriter = dt.tmpwriter(useifdh=useifdh)

    #getparametriczpt("/global/cscratch1/sd/dbrout/v6/","/global/cscratch1/sd/dbrout/v6/stardata_"+filter)

    if not cacheddata:
        dostars = True
        if dostars:
            grabstardata("/global/cscratch1/sd/dbrout/v6/","/global/cscratch1/sd/dbrout/v6/stardata_"+tfield+"_"+filter,tfield,filter)
            stardata = np.load('/global/cscratch1/sd/dbrout/v6/stardata_'+tfield+"_"+filter+'.npz')
            plotstarrms(stardata['starflux'], np.sqrt(stardata['starfluxerr'] ** 2), stardata['starzpt'],
                    stardata['catmag'], stardata['chisq'], stardata['rmsaddin'], stardata['sky'], stardata['skyerr'],
                    stardata['poisson'],stardata['ids'],stardata['centroidedras'],stardata['centroideddecs'],
                        stardata['fwhm'],stardata['zptscat'],stardata['mjd'],
                    title=tfield+'_'+filter+'_',outdir='/global/cscratch1/sd/dbrout/v6/')
        data = grabdata(tmpwriter,resultsdir,cd,tfield,filter=filter)
    else:

        dostars = True
        if dostars:
            stardata = np.load('/global/cscratch1/sd/dbrout/v6/stardata_'+tfield+"_"+filter+'.npz')
            plotstarrms(stardata['starflux'], np.sqrt(stardata['starfluxerr'] ** 2), stardata['starzpt'],
                        stardata['catmag'], stardata['chisq'], stardata['rmsaddin'], stardata['sky'], stardata['skyerr'],
                        stardata['poisson'],stardata['ids'],stardata['centroidedras'],stardata['centroideddecs'],
                        stardata['fwhm'],stardata['zptscat'],stardata['mjd'],
                        title=tfield+'_'+filter+'_',outdir='/global/cscratch1/sd/dbrout/v6/')
        data = np.load(cd)
    if not os.path.exists(resultsdir+'/Summary/'):
        os.mkdir(resultsdir+'/Summary/')
    if not os.path.exists(resultsdir+'/Summary/'+filter+'/'):
        os.mkdir(resultsdir+'/Summary/'+filter+'/')
    #raw_input()
    plotpercentageresid(data['Flux'],data['Fluxerr'],data['FakeMag'],data['FitZPT'],data['FakeZPT'], data['diffimflux'],data['diffimfluxerr'],
                        data['sky'],data['skyerr'],data['DPMJD'],data['Chisq'],data['imfiles'],data['ra'],data['dec'],
                        data['image_stamp'],resultsdir+'/Summary/'+filter+'/',data['fakefiles'],data['HostMag'],
                        filter,data['FakeZPT'],data['rmsaddin'])
    plotsigmaresid(data['Flux'],data['Fluxerr'],data['FakeMag'], data['FitZPT'], data['FakeZPT'],data['HostMag'],
                   data['Chisq'],data['rmsaddin'],data['field'],resultsdir+'/Summary/'+filter+'/',data['rmsaddin'],
                   data['diffimflux'], data['diffimfluxerr'])#resultsdir)


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
        if cntr > 1000.: break
        # print('Found directory: %s' % dirName)
        for fname in fileList:
            # print fname
            if 'globalstar.npz' in fname:
                # print('\t%s' % fname)
                # print os.path.join(imagedir,dirName,fname)
                if not 'SN-S2' in fname: continue
                #    if not 'SN-S1' in fname: continue
                try:
                    os.system('cp ' + os.path.join(imagedir, dirName, fname) + ' test.npz')
                    zptdata = np.load('test.npz')
                except:
                    print 'could not load'
                    continue
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
               'numzptstars':[], 'fwhm':[],'mjd':[]}
    zptfiles = []
    cntr = 0
    goodbigdata = copy(bigdata)
    for dirName, subdirList, fileList in os.walk(imagedir):
        if cntr > 1000.: break
        #print('Found directory: %s' % dirName)
        for fname in fileList:
            #print fname
            if 'globalstar.npz' in fname:
                #print('\t%s' % fname)
                #print os.path.join(imagedir,dirName,fname)
                if not tfield in fname: continue
                if not '_'+filt+'_' in fname: continue
                #    if not 'SN-S1' in fname: continue
                try:
                    os.system('cp ' + os.path.join(imagedir,dirName,fname) + ' test.npz')
                    zptdata = np.load('test.npz')
                except:
                    print 'could not load'
                    continue
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
                        test = zptdata['fit_mag']
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
                        for i in range(len(psfs)):
                            bigdata['psf'].append(psfs[i,:,:])
                            bigdata['poisson'].append(np.sqrt(np.sum(psfs[i,:,:].ravel()**2*zptdata['flux_starmp'][i])))
                            #print zptdata['flux_starnormm'][i],zptdata['flux_star_std'][i],bigdata['poisson'][-1]
                            #raw_input()

                        fs = zptdata['flux_starmp']
                        zp = zptdata['fit_zpt']
                        ww = (cm < 20.) & (cm > 17.)

                        #plt.scatter(cm[ww],float(zp) - cm[ww] - 2.5*np.log10(fs[ww]))
                        # plt.scatter(cm[ww],- 2.5*np.log10(fs[ww]))
                        # plt.savefig('testzpt.png')
                        md, std = iterstat.iterstat(float(zp) - cm[ww] - 2.5*np.log10(fs[ww]),
                                                     startMedian=True, sigmaclip=3, iter=10)

                        print 'worked now std',std/np.sqrt(len(cm[ww]))
                        bigdata['numzptstars'].extend(zptdata['flux_starmp']*0. + len(cm[ww]))
                        bigdata['rmsaddin'].extend(zptdata['flux_starmp']*0. + std/np.sqrt(len(cm[ww])))
                        bigdata['zptscat'].extend(zptdata['flux_starmp']*0.+zptdata['zptscat'])
                        #print 'read in ',fname
                        zptfiles.append(fname)

                        bigdata['chisq'].extend(zptdata['chisqu'])

                        cntr += 1
                        print 'CNTR',cntr

                        goodbigdata = copy(bigdata)
                    except:
                        print 'FAILED', fname
                        pass

    try:
        bigdata['centroidedras'] = np.array(bigdata['centroidedras'])
        bigdata['ids'] = np.array(bigdata['ids'])
        stds = []
        for i in zptdata['ids']:
            #print i
            #print bigdata['centroidedras'].shape
            #print bigdata['ids'] == i
            print i,np.mean(bigdata['centroidedras'][bigdata['ids'] == i]),np.std(bigdata['centroidedras'][bigdata['ids'] == i])
            if np.std(bigdata['centroidedras'][bigdata['ids'] == i]) > 0.:
                stds.append(np.std(bigdata['centroidedras'][bigdata['ids'] == i]))

        print 'std is', np.mean(stds)
    except:
        print 'ids not in archiv'
    #sys.exit()

    np.savez(outfile, **bigdata)
    print np.unique(bigdata['mjd'])
    sys.exit()
    #os.system('ifdh rm ' + outfile)
    #os.system('ifdh cp ' + 'dat.dat' + ' ' + outfile)
    #os.system('rm dat.dat')
    #sys.exit()

def grabdata(tmpwriter,resultsdir,cd,tfield,filter = 'g',oldformat=False):

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

    dofakedata = pd.read_csv('data/doFake.out', delim_whitespace=True, header=0)
    # print dofakedata
    dofakeexpnum = dofakedata['EXPNUM'].values
    dofakemag2 = dofakedata['TRUEMAG'].values
    dofaketflux = dofakedata['TRUEFLUXCNT'].values
    dofakeflux = dofakedata['FLUXCNT'].values
    dofakeid = dofakedata['FAKEID'].values
    dofakezpt = dofakemag2 + 2.5 * np.log10(dofaketflux)
    print 'done reading dofake'


    #raw_input('dofakemjd')
    files = os.listdir(os.path.join(resultsdir, 'lightcurves'))
    smpfiles = []
    for f in files:
        if filter+'.smp' in f:
            smpfiles.append(os.path.join(resultsdir, 'lightcurves', f))

    #print "Found " + str(len(smpfiles)) + " .smp files"

    if not os.path.exists(os.path.join(resultsdir,'Summary')):
        os.makedirs(os.path.join(resultsdir,'Summary'))
    #outfile = os.path.join(resultsdir,'Summary','sumdata.npz')
    outfile = cd
    bigdata = {'Flux':[],'Fluxerr':[],'FakeMag':[],'FitZPT':[],'FakeZPT':[],'HostMag':[],'Chisq':[],'DPMJD':[],
               'starflux':[],'starfluxerr':[],'starzpt':[],'catmag':[],'rmsaddin':[],'field':[],'sky':[],'imfiles':[],
               'mjd':[],'fakefile':[],'ra':[],'dec':[],'image_stamp':[],'fakefiles':[],'diffzpt':[],'diffimflux':[],
               'diffimfluxerr':[],'fakeid':[],'skyerr':[],'acceptance':[]}
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
        cntr += 1
        if cntr > 250: continue
        #if cntr == 34: continue
        #if cntr == 53: continue
        #if not '_r.smp' in f: continue
        print cntr, 'of',tot
        deep = 0
        #os.system('cp '+f+' test.npz')
        data = dt.readcol(f)
        tra = data['RA']
        #print tra
        if len(tra) == 0:
            print 'empty'
            continue
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

        if not oldformat:
            if not int(fakeid) in dofakeid:
                print 'dddd'
                fakemag = data['FAKEMAG']*0. + 99.
            else:
                fakemag = data['FAKEMAG']
        else:
            fakemag = data['FAKEMAG']
        print fakemag
        #raw_input()
        if len(data['DPMJD'][data['DPMJD'] > 300.]) < 2:
            print 'bad DPMJD'*20
            continue
        if np.min(data['FLUX']) < -10000:
            print 'BAD FLUX'*20
            continue

        skipnewfakemag = True
        #if not skipnewfakemag:
        if True:
            newfakemag = []
            for imf,fm,x,y in zip(data['IMAGE_FILE'],fakemag,data['XPOS'],data['YPOS']):

                #print imf
                try:
                    exn = imf.split('/')[-1].split('_')[1]
                    #ccd = float(imf.split('_')[7].split('+')[0])
                except:
                    print 'skipping'
                    newfakemag.append(99)
                    bigdata['fakeid'].append(99.)
                    bigdata['FakeZPT'].append(31.)

                    continue
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
                if not len(dofakemag2[www]) > 0:
                    newfakemag.append(99.)
                    bigdata['fakeid'].append(99.)
                    bigdata['FakeZPT'].append(31.)
                else:
                    nfm = float(fm) + 2.5*np.log10(dofaketflux[www][0]) - 2.5*np.log10(dofakeflux[www][0])
                    newfakemag.append(nfm)
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
            #try:
            #print data['ID_OBS']
            #raw_input()
            #print len(data['FLUX']),len(data['FLUXERR']),len(data['FAKEMAG']),len(data['ZPT']),(data['FAKEZPT'])
            #data2 = dt.readcol('./working/lightcurves/' + f.split('/')[-1])
            #rms = np.mean(data2['RMSADDIN'][data2['RMSADDIN'] > 0.0])
            #print data.keys()
            #raw_input( )
            bigdata['rmsaddin'].extend(data['ZPTERR'])

            #print data.keys()
            #raw_input()

            bigdata['Flux'].extend(data['FLUX'])
            bigdata['Fluxerr'].extend(data['FLUXERR'])
            if not skipnewfakemag:
                bigdata['FakeMag'].extend(fakemag)
            else:
                bigdata['FakeMag'].extend(data['FAKEMAG'])
            bigdata['FitZPT'].extend(data['ZPT'])
            #print data['ZPT'],data['FAKEZPT']
            #raw_input('aaa')
            #bigdata['FakeZPT'].extend(data['FAKEZPT'])
            bigdata['Chisq'].extend(data['CHI2'])
            bigdata['sky'].extend(data['SKY'])
            bigdata['DPMJD'].extend(data['DPMJD'])
            bigdata['mjd'].extend(data['MJD'])
            bigdata['imfiles'].extend(data['IMAGE_FILE'])
            bigdata['fakefiles'].extend([f for i in range(len(data['FLUX']))])
            bigdata['diffimflux'].extend(data['DIFFIM_FLUX'])
            bigdata['diffimfluxerr'].extend(data['DIFFIM_FLUXERR'])
            bigdata['skyerr'].extend(data['SKYERR'])

            print np.mean(data['CHI2']), bigdata['fakeid'][-10:]

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


def plotpercentageresid(flux,fluxerr,fakemag,fitzpt,fakezpt,diffimflux,diffimfluxerr,sky,skyerr,dpmjd,chisq,imfiles,ra,dec,imstamp,outdir,fakefiles,hostmag,filter,oldfakezpt,zptstd):
    flux = np.asarray(flux)
    fakemag = np.asarray(fakemag,dtype='float')
    dpmjd = np.asarray(dpmjd,dtype='float')
    sky = np.asarray(sky)
    chisq = np.asarray(chisq)
    imfiles = np.asarray(imfiles,dtype='str')
    ra = np.asarray(ra)
    dec = np.asarray(dec)
    imstamp = np.asarray(imstamp)
    fakefiles = np.asarray(fakefiles,dtype='str')
    hostmag = np.asarray(hostmag)
    skyerr = np.asarray(skyerr)
    zptstd=np.asarray(zptstd)
    # print hostmag.shape
    # raw_input()
    #print np.unique(fakemag)
    #print np.unique(flux)
    #raw_input('asdf#')

    for ft,fa in zip(fitzpt,fakezpt):
        print ft,fa
    #raw_input('asdffff')
    fluxerr = np.asarray(fluxerr)
    fakezpt = np.asarray(fakezpt)

    diffimflux = 10**(.4*(31-27.5))*np.array(diffimflux)
    diffimfluxerr = 10**(.4*(31-27.5))*np.array(diffimfluxerr)

    #dmag  = -2.5*np.log10(diffimflux) + oldfakezpt
    #diffimflux = 10**(.4*(31-dmag))

    oldfakezpt = np.array(oldfakezpt)

    #print fakefiles[(fluxerr<30.) & (fluxerr>0.)]
    #raw_input()
    #print fakefiles[(chisq < 0.1)]
    #raw_input()
    # for fl , ff in zip(fakefiles[(hostmag<21.4) & (fluxerr>0.) & (fakemag > 90.)],flux[(hostmag<21.4) & (fluxerr>0.) & (fakemag > 90.)]):
    #     print ff, fl
    for fl , ff in zip(fakefiles[(fluxerr>0.) & (fakemag < 90.)],flux[(fluxerr>0.) & (fakemag < 90.)]):
        print ff, fl
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
    fluxerr = (fluxerr**2 + (flux*zptstd)**2)**.5

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


    d = (flux - fakeflux) / ((fluxerr ** 2) ** .5)
    ww = (flux != 0.)
    d = d[ww]
    fm = np.array(fakemag,dtype='float')[ww]
    dc99 = d[fm > 90.]
    rms99 = np.sqrt(np.nanmean(np.square(dc99[abs(dc99) < 3.])))

    ww = (flux != 0.) & (fakemag != 0)#(fakemag < 28.5) & (flux != 0.)
    plt.clf()
    fig = plt.figure(figsize=(15, 10))


    # plt.scatter(fakemag[ww], (-2.5*np.log10(diffimflux[ww]) + oldfakezpt[ww] - fakemag[ww]), alpha=.1, color='red')
    # ax, ay, aystd = bindata(fakemag[ww],(-2.5*np.log10(diffimflux[ww]) + oldfakezpt[ww] - fakemag[ww]),
    #                         np.arange(min(fakemag[ww]), max(fakemag[ww]), .5))
    # plt.errorbar(ax, ay, aystd, markersize=10, color='red', fmt='o', label='Diffim',alpha=.4)

    #print fakemag[ww].shape,flux[ww].shape,fakeflux[ww].shape

    for x in np.unique(fakeflux):
        print x
    #print np.unique(fakeflux)
    #raw_input('fakeflux')

    plt.scatter(fakemag[ww], (diffimflux[ww] - fakeflux[ww]) / fakeflux[ww], alpha=.2,color='red')
    ax, ay, aystd = bindata(fakemag[ww], (diffimflux[ww] - fakeflux[ww]) / fakeflux[ww],
                            np.arange(19, 28, .5))
    plt.errorbar(ax, ay, aystd, markersize=15, color='red', fmt='o', label='DIFFIMG')

    plt.scatter(fakemag[ww],(flux[ww]-fakeflux[ww])/fakeflux[ww],alpha=.2,color='green')
    ax, ay, aystd = bindata(fakemag[ww],(flux[ww]-fakeflux[ww])/fakeflux[ww],
                            np.arange(19,28, .5))
    plt.errorbar(ax, ay, aystd, markersize=15, color='green', fmt='o', label='SMP')



    plt.axhline(0)
    plt.xlim(18,29)
    #plt.ylim(-.1,.1)
    plt.ylim(-.2,.2)
    plt.xlabel('Fake Mag',fontsize=30.)
    plt.ylabel('Fractional Flux Difference',fontsize=30.)
    plt.title(filter+' band',fontsize=30.)
    plt.savefig(outdir+'/percentagefluxdiff.png')

    print 'saved', outdir+'/percentagefluxdiff.png'
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
    plt.ylim(-600, 600)
    plt.xlabel('Fake Flux')
    plt.ylabel('Flux Difference ')
    plt.title(filter+' band')

    plt.savefig(outdir + '/fluxdiff.png')

    plt.clf()
    fig = plt.figure(figsize=(15, 10))
    plt.hist((flux[ww] - fakeflux[ww])/fakeflux[ww],bins=np.arange(-.51,.5,.02), histtype='step',color='blue',label='SMP STD '+
                                    str(round(1.48 * np.median(abs((flux[ww] - fakeflux[ww])/fakeflux[ww] - np.median((flux[ww] - fakeflux[ww])/fakeflux[ww]))), 3)))
    plt.hist((diffimflux[ww] - fakeflux[ww])/fakeflux[ww],bins=np.arange(-.51,.5,.02), histtype='step',color='red',label='DIFFIM STD '+
                                    str(round(1.48 * np.median(abs((diffimflux[ww] - fakeflux[ww])/fakeflux[ww] - np.median((diffimflux[ww] - fakeflux[ww])/fakeflux[ww]))), 3)))


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

    try:
        plt.scatter(hostmag[ww],flux[ww] / fluxerr[ww])
        plt.xlim(21, 29)
        plt.ylim(-3,3)
        ax, ay, aystd = dt.bindata(hostmag[ww], (flux[ww] / fluxerr[ww]),
                                   np.arange(19, 28, .2))
        plt.errorbar(ax, ay, aystd, markersize=10, color='gold', fmt='o', label='SMP')
        plt.xlabel('Hostmag')
        plt.ylabel('Flux/Fluxerr (epochs without fake SN Flux)')
        plt.axhline(0)
        plt.title(filter + ' band')

        plt.savefig(outdir + '/emptyfluxstdvshostmag.png')

        plt.clf()
        plt.scatter(hostmag[ww], flux[ww])
        plt.xlim(21, 30)
        plt.ylim(-700, 700)
        ax, ay, aystd = dt.bindata(hostmag[ww], (flux[ww]),
                                   np.arange(19, 30, .2))
        plt.errorbar(ax, ay, aystd, markersize=10, color='gold', fmt='o', label='SMP')
        plt.xlabel('Hostmag')
        plt.ylabel('Flux (epochs without fake SN Flux)')
        plt.axhline(0)
        plt.title(filter + ' band')

        plt.savefig(outdir + '/emptyfluxvshostmag.png')
    except:
        print 'bad hostmags'
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

def plotsigmaresid(flux,fluxerr,fakemag,fitzpt,fakezpt,hostmag,chisqarr,rmsaddin,deep,outdir,zptstd,diffimflux,diffimfluxerr):
    flux = np.asarray(flux)
    fakemag = np.asarray(fakemag)
    fluxerr = np.asarray(fluxerr)
    fitzpt = np.asarray(fitzpt)
    fakezpt = np.asarray(fakezpt)
    rmsaddin = np.asarray(rmsaddin)
    fakeflux = 10 ** (.4 * (31. - fakemag))
    fakeflux *= 10**(-1*.4*(fitzpt - fakezpt))
    fluxerrz = (fluxerr**2 + (flux*np.array(zptstd))**2 + flux)**.5
    diffimgflux = np.array(diffimflux)
    diffimgfluxerr = np.array(diffimfluxerr)
    diffimflux = 10 ** (.4 * (31 - 27.5)) * np.array(diffimflux)
    diffimfluxerr = 10 ** (.4 * (31 - 27.5)) * np.array(diffimfluxerr)

    #fluxerr *= 10**(.4*(fitzpt - fakezpt))

    #fakeflux *= 10 ** (-1 * .4 * (fitzpt - fakezpt))
    chisqarr = np.asarray(chisqarr)
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
    ff[ff < 1.] = 1.
    plt.hist(fluxerr/ff,bins=20,normed=True)
    plt.title(filter+' band')

    plt.savefig('sinosn.png')
    plt.clf()
    np.savez('simnosn.npz',flux=flux,fakeflux=ff,fluxerr=np.sqrt(fluxerr**2 + abs(flux)/3.8))

    d = (flux - fakeflux) / ((fluxerr**2 )**.5)
    dz = (flux - fakeflux) / ((fluxerrz**2 )**.5)
    df = (diffimflux - fakeflux) / ((diffimfluxerr**2 )**.5)
    ww = (flux != 0.) & (np.array(fakemag, dtype='float') > 0.) #& (deep == 0)

    #fakemag[fakemag==99] = 29.5
    flux = flux[ww]
    fakemag = fakemag[ww]
    fitzpt = fitzpt[ww]
    fakezpt = fakezpt[ww]
    fakeflux= fakeflux[ww]
    fluxerr=fluxerr[ww]
    fluxerrz = fluxerrz[ww]
    d = d[ww]
    dz = dz[ww]
    df = df[ww]
    hostmag = hostmag[ww]
    chisqarr = chisqarr[ww]

    #print flux[0:10]
    #print fakeflux[0:10]
    #print flux.shape
    #print fakeflux.shape
    #d = (flux - fakeflux) / fluxerr

    chisq = (flux - fakeflux) ** 2 / fluxerr ** 2
    chisq = np.nanmean(chisq[abs(d) < 3])

    plt.clf()

    dc = d[abs(d) < 3.]

    dc99 = d[(np.array(fakemag,dtype='float') > 90.) & (chisqarr > .5) & (chisqarr < 1.4)]
    rms99 = np.sqrt(np.nanmean(np.square(dc99[abs(dc99) < 3.])))

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
         (np.isfinite(fluxerr)) &(~np.isnan(flux)) &(~np.isnan(fluxerr)) & (chisqarr > .2) & (chisqarr < 2.)

    #ww = (flux != 0) & (fakeflux < 1.)
    #print rms99
    fig = plt.figure(figsize=(15, 10))
    plt.hist(flux[ww] / fluxerr[ww], bins=np.arange(-6.1, 6., .2), normed=True,
             label='RMS Fakemag = 99: ' + str(round(rms99, 3)))
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
    plt.title(filter + ' band')

    plt.savefig(outdir + '/efluxdiffstd.png')

    plt.clf()
    ww = (flux != 0) & (np.array(fakemag, dtype='float') > 90.) & (fluxerr > 0.) & (np.isfinite(flux)) & \
         (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .5) & (chisqarr < 1.7)

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
    plt.title(filter + ' band')
    plt.xlim(-1000, 1000)

    plt.savefig(outdir + '/efluxdiff.png')

    plt.clf()
    ww = (flux != 0) & (np.array(fakemag, dtype='float') > 90.) & (fluxerr > 0.) & (np.isfinite(flux)) & \
         (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .5) & (chisqarr < 1.4)

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

    for fmm,fffl,fffle,fafl in zip(fakemag[(flux != 0.) & ((np.abs(flux-fakeflux) / fluxerr) < .1) & (np.array(fakemag, dtype='float') < 27.)& (np.array(fakemag, dtype='float') > 0.)],
                        flux[(flux != 0.) & ((np.abs(flux-fakeflux) / fluxerr) < .1) & (np.array(fakemag, dtype='float') < 27.)& (np.array(fakemag, dtype='float') > 0.)],
                             fluxerr[(flux != 0.) & ((np.abs(flux - fakeflux) / fluxerr) < .1) & (
                                 np.array(fakemag, dtype='float') < 27.) & (np.array(fakemag, dtype='float') > 0.)],
                        fakeflux[(flux != 0.) & ((np.abs(flux-fakeflux) / fluxerr) < .1) & (np.array(fakemag, dtype='float') < 27.)& (np.array(fakemag, dtype='float') > 0.)]):

        print fmm,fffl,fafl, fffl-fafl, fffle
    #raw_input()
    ww = (flux != 0) & (np.array(fakemag, dtype='float') > 0.) & (np.array(fakemag, dtype='float') < 30.) & (fluxerr > 0.) & (np.isfinite(flux)) & \
         (np.isfinite(fluxerr)) & (~np.isnan(flux)) & (~np.isnan(fluxerr)) & (chisqarr > .2) & (chisqarr < 2.) & \
         (abs(flux - fakeflux) > 0.1)

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

    rect_scatter = [left, bottom + height / 2., width, height / 2.]
    rect_scatterflux = [left, bottom, width, height / 2.]
    rect_histx = [left, bottom_h - .04, width, .18]
    rect_histy = [left_h, bottom + height / 2., 0.2, height / 2.]
    rect_histyflux = [left_h, bottom, 0.2, height / 2.]

    # start with a rectangular Figure
    plt.figure(1, figsize=(45, 40))

    ax1 = plt.axes(rect_scatter)
    ax3 = plt.axes(rect_histx)
    ax2 = plt.axes(rect_histy)
    ax4 = plt.axes(rect_scatterflux)
    ax5 = plt.axes(rect_histyflux)

    # no labels
    ax2.yaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax5.yaxis.set_major_formatter(nullfmt)

    outliers3 = float(len(d[(abs(d)>3.) & (chisqarr > .5) & (chisqarr < 1.5) & (np.array(fakemag, dtype='float') > 0.)]))/float(len(d))
    outliers5 = float(len(d[(abs(d)>5.) & (chisqarr > .6) & (chisqarr < 1.5) & (np.array(fakemag, dtype='float') > 0.)]))/float(len(d))

    ax2.hist(d[np.isfinite(d)], bins=np.arange(-10, 10, .25), normed=True,label='RMS Fakemag = 99: ' + str(round(rms99, 3))+
                                                                '\nRMS Fakemag < 99: '+ str(round(rmsr, 3))+'\n3sig Outlier'#+
                                                                ' Fraction: '+str(round(outliers3,3))+'\n5sig Outlier'+
                                                                ' Fraction: '+str(round(outliers5,3))
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
    #fakemag[fakemag == 99] = 28.5




    ax1.scatter(fakemag,d,alpha=.3,color='blue')
    ax, ay, aystd = dt.bindata(fakemag[(d<3.)& (np.array(fakemag, dtype='float') > 0.)], d[(d<3.)& (np.array(fakemag, dtype='float') > 0.)], np.arange(19., 28, .1),window=.5)
    ax1.plot([19, 28.7], [0, 0],color='grey')
    ax1.plot(ax, ay, linewidth=3, color='orange', label='SMP')
    ax1.plot(ax, ay+aystd, linewidth=2, color='orange',linestyle='--', label='SMP')
    ax1.plot(ax, ay-aystd, linewidth=2, color='orange',linestyle='--', label='SMP')

    #ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

    ax1.set_xlim(19, 28.7)
    ax1.set_ylim(-3., 3.)
    ax1.set_xlabel('Fake Mag')
    ax1.set_ylabel('STD')

    #ax, ayrms= dt.binrms(fakemag, d, np.arange(19.5, max(fakemag), .1),.5)
    #ax3.plot(ax, ayrms, color='blue',label='RMS',linewidth=3)


    ax3.plot([0,100],[1.,1.],linestyle='--',color='black')
    ax3.set_ylim(.9,1.25)
    ax3.legend(fontsize='x-small')

    fresid = np.zeros(flux.shape)
    for i,f,ff in zip(range(len(flux)),flux,fakeflux):
        if f == 0.:
            fresid[i] = np.nan
        else:
            fresid[i] = (f - ff) / max([abs(ff),1.])
    #fresid[abs(fakeflux) < 1.] = flux[abs(fakeflux) < 1.] - fakeflux[abs(fakeflux) < 1.]

    ax5.hist(fresid, bins=np.arange(-.155,.15,.01),color='blue', orientation='horizontal')


    ax4.scatter(fakemag,fresid,alpha=.3,color='blue')
    axa, aya, aystd = dt.bindata(fakemag,fresid,
                            np.arange(20., 26., .1),window=2.)
    ax4.plot([19, 28.7], [0, 0],color='grey')

    ax, ayrms = dt.binrms(fakemag, d, np.arange(20., 28, .1), 1.5)
    ax3.plot(ax, ayrms, color='blue', label='ALL SNe', linewidth=3)
    ax, ayrms = dt.binrms(fakemag, dz, np.arange(20., 28, .1), 1.5)
    ax3.plot(ax, ayrms, color='blue',linestyle='--', label='ALL SNe', linewidth=3)
    ax, ayrms = dt.binrms(fakemag, df, np.arange(20., 28, .1), 1.5)
    ax3.plot(ax, ayrms, color='red', linestyle='--', label='DIFFIMG', linewidth=3)
    ax3.plot(ax, ax * 0 + 1., linestyle='--', color='black')

    # ww = hostmag > 25.
    # ax, ayrms = dt.binrms(fakemag[ww], d[ww], np.arange(19.5, max(fakemag), .1), .5)
    # ax3.plot(ax, ayrms, color='red', label='HostMag > 25.', linewidth=3)
    #
    # ww = hostmag < 23.
    # ax, ayrms = dt.binrms(fakemag[ww], d[ww], np.arange(19.5, max(fakemag), .1), .5)
    # ax3.plot(ax, ayrms, color='green', label='HostMag < 23', linewidth=3)
    # ax3.legend(fontsize='small')

    ax4.plot(axa, aya, linewidth=3, color='orange')
    ax4.plot(axa, aya+aystd, linewidth=2, color='orange',linestyle='--')
    ax4.plot(axa, aya-aystd, linewidth=2, color='orange',linestyle='--')
    ax4.set_xlim(ax1.get_xlim())
    ax4.set_ylim(-.1,.1)
    ax4.set_xlabel('Fake Mag')
    ax5.set_xlabel('Counts')
    ax3.set_ylabel('RMS')
    ax3.set_title(filter+' band')
    ax4.set_ylabel('(fitflux - fakeflux)/fakeflux')
    plt.title(filter+' band')

    ax3.set_xlim(ax1.get_xlim())
    ax2.set_ylim(ax1.get_ylim())
    ax5.set_ylim(ax4.get_ylim())
    ax2.xaxis.set_major_formatter(nullfmt)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax1.xaxis.set_major_formatter(nullfmt)

    plt.subplots_adjust(wspace=0.01,hspace=0.01)

    plt.savefig(outdir+'/std.png')


    #--------------------------------------------------------------------------------------
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
        ww = fakemag < 30.
        ax1.scatter(hostmag[ww],d[ww],alpha=.3,color='blue')
        ax, ay, aystd = dt.bindata(hostmag[ww], d[ww], np.arange(min(hostmag), 27.5, .1),window=1.5)
        ax1.plot([min(hostmag), max(hostmag)], [0, 0],color='grey')
        ax1.plot(ax, ay, linewidth=3, color='orange', label='SMP')
        ax1.plot(ax, ay+aystd, linewidth=2, color='orange',linestyle='--', label='SMP')
        ax1.plot(ax, ay-aystd, linewidth=2, color='orange',linestyle='--', label='SMP')

        #ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

        ax1.set_xlim(21,28)
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
                                np.arange(min(hostmag), 27.5, .1),window=1.)
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

        ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag), max(hostmag), .1), 1.5)
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
        ax4 = plt.axes(rect_scatterflux)
        ax5 = plt.axes(rect_histyflux)

        # no labels
        ax2.yaxis.set_major_formatter(nullfmt)
        ax3.xaxis.set_major_formatter(nullfmt)
        ax5.yaxis.set_major_formatter(nullfmt)

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
        ww = fakemag > 90.
        ax1.scatter(hostmag[ww], d[ww], alpha=.3, color='blue')
        ax, ay, aystd = dt.bindata(hostmag[ww], d[ww], np.arange(min(hostmag), 27.5, .1), window=1.5)
        ax1.plot([min(hostmag), max(hostmag)], [0, 0], color='grey')
        ax1.plot(ax, ay, linewidth=3, color='orange', label='SMP')
        ax1.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP')
        ax1.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP')

        # ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

        ax1.set_xlim(21, 28)
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

        ax5.hist(fresid[ww], bins=np.arange(-.155, .15, .01), color='blue', orientation='horizontal')

        ax4.scatter(hostmag[ww], fresid[ww], alpha=.3, color='blue')
        ax, ay, aystd = dt.bindata(hostmag[ww], fresid[ww],
                                   np.arange(min(hostmag), 27.5, .1), window=1.)
        ax4.plot([min(hostmag), max(hostmag)], [0, 0], color='grey')

        ax4.plot(ax, ay, linewidth=3, color='orange')
        ax4.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--')
        ax4.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--')
        ax4.set_xlim(ax1.get_xlim())
        ax4.set_ylim(-.2, .2)
        ax4.set_xlabel('Host Mag')
        ax5.set_xlabel('Counts')
        ax3.set_ylabel('RMS')
        ax3.set_title(filter + ' band')

        ax4.set_ylabel('(fitflux - fakeflux)/fakeflux')

        ax3.set_xlim(ax1.get_xlim())
        ax3.set_ylim(.9, 1.25)
        ax2.set_ylim(ax1.get_ylim())
        ax5.set_ylim(ax4.get_ylim())
        ax2.xaxis.set_major_formatter(nullfmt)
        ax3.xaxis.set_major_formatter(nullfmt)
        ax1.xaxis.set_major_formatter(nullfmt)
        plt.subplots_adjust(wspace=0.001, hspace=0.001)

        #ww = fakemag > 90.
        ax, ayrms = dt.binrms(hostmag[ww], d[ww], np.arange(min(hostmag), max(hostmag), .1), 1.5)
        ax3.plot(ax, ayrms, color='blue', label='ALL SNe w/ Light', linewidth=3)
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




def plotstarrms(flux,fluxerr,zpt,catmag,chisq,rmsaddin,sky,skyerr,poisson,indices,ras,decs,fwhm,zptscat,mjd,title='',outdir='.'):
    catflux = 10 ** (.4 * (zpt - catmag))
    ff = (flux - catflux) / catflux
    st = np.std(ff)
    fluxerr = np.sqrt(fluxerr**2)

    plt.clf()
    plt.scatter(catmag[(np.floor(mjd) == 56575)],fluxerr[(np.floor(mjd) == 56575)], color='black')
    plt.xlabel('SkyERR')
    plt.ylabel('#')
    plt.savefig(outdir + '/' + title + 'photerr.png')
    print outdir + '/' + title + 'photerr.png'
    #raw_input('skyerrdone')

    plt.clf()
    plt.hist(skyerr[np.floor(mjd)==56575],bins=np.arange(10,70,1.),color='black')
    plt.xlabel('SkyERR')
    plt.ylabel('#')
    plt.savefig(outdir + '/' + title + 'skyerrhist.png')
    print outdir + '/' + title + 'skyerrhist.png'
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
    starmagerr = 1.0857*np.sqrt(fluxerr**2+flux)/flux #+ rmsaddin
    starmagerr = np.sqrt(fluxerr**2)/flux #+ rmsaddin

    starmagerrzpt = np.sqrt(fluxerr**2+zptscat*flux)/flux #+ rmsaddin


    rv = (flux-catflux)/(fluxerr**2+zptscat*flux)**.5

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

    maxpoints = 10000

    cntr = 0
    pltvecx = []
    pltvecy = []
    for sme, sm, ind, r, d, cm, f, fe,fh in zip(starmagerr, starmag, indices, ras, decs, catmag, flux, fluxerr, fwhm):
        cntr += 1
        if cntr > maxpoints: continue
        if cntr > 100: continue

        # print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        # print starmag[indices == ind]
        # raw_input()
        starww = starmag[np.isclose(ras, r, rtol=1.e-5) & np.isclose(decs, d, rtol=1.e-5) & (catmag == cm)]
        repeatability = np.std(starww) / np.sqrt(len(starww))
        # repeatability = np.std(starmag[indices == ind])
        if len(starww) > 5.:
            # if repeatability < .3:
            plt.scatter(fh, sme - repeatability, alpha=.3, color='black')
            pltvecy.append(sme - repeatability)
            pltvecx.append(fh)

    #plt.xscale('log')
    plt.xlabel('PSF FWHM')
    plt.ylabel('Repeatability - PhotErr')
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
        if cntr > 100: continue

        #print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        #print starmag[indices == ind]
        #raw_input()
        starww = starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        repeatability = np.std(starww)/np.sqrt(len(starww))
        #repeatability = np.std(starmag[indices == ind])
        if len(starww) > 5.:
            #if repeatability < .3:
            plt.scatter(sme,sme-repeatability,alpha=.3,color='black')
            pltvecy.append(sme-repeatability)
            pltvecx.append(sme)

    plt.xscale('log')
    plt.xlabel('Photometric Error')
    plt.ylabel('Repeatability - PhotErr')
    plt.xlim(.0003,.02)
    plt.ylim(-.02,.01)

    ax, ay, aystd = dt.bindata(np.array(pltvecx),np.array(pltvecy), np.arange(.0003,.007, .00001), window=.00004,dontrootn=True)
    plt.plot(ax, ay, linewidth=3, color='orange', label='SMP',alpha=.6)
    plt.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP',alpha=.6)
    plt.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP',alpha=.6)

    plt.title(title+'BAND')
    print 'finished resid'

    plt.plot([.0003,.02],[0,0],color='black')
    plt.savefig(outdir+'/'+title+'_repeatability-photerr_vs_photerrr.png')

    plt.clf()
    cntr = 0
    pltvecx = []
    pltvecy = []
    for sme, sm, ind, r, d, cm, cs, se, rai in zip(starmagerr, starmag, indices, ras, decs, catmag, chisq, skyerr,
                                                   rmsaddin):
        cntr += 1
        if cntr > maxpoints: continue
        if cntr > 100: continue

        # print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        # print starmag[indices == ind]
        # raw_input()
        starww = starmag[np.isclose(ras, r, rtol=1.e-5) & np.isclose(decs, d, rtol=1.e-5) & (catmag == cm)]
        repeatability = np.std(starww) / np.sqrt(len(starww))
        # repeatability = np.std(starmag[indices == ind])
        if len(starww) > 4.:
            # if repeatability < .3:
            plt.scatter(rai, float(sme) - repeatability, alpha=.3, color='black')
            pltvecy.append(sme - repeatability)
            pltvecx.append(rai)

    plt.xscale('log')
    plt.xlim(0.003, .05)
    plt.ylim(-.02, .01)

    ax, ay, aystd = dt.bindata(np.array(pltvecx), np.array(pltvecy), np.arange(0., .1, .00005), window=.0001,
                               dontrootn=True)
    plt.plot(ax, ay, linewidth=3, color='orange', label='SMP', alpha=.6)
    plt.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP', alpha=.6)
    plt.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP', alpha=.6)

    plt.xlabel('ZPT Uncertainty')
    plt.ylabel('PhotErr - Repeatability')
    plt.title(title+'BAND')
    plt.plot([0., 1000], [0, 0], color='black')
    plt.plot([0., 1000], [0, 1000], color='black', linestyle='--')

    plt.savefig(outdir + '/' + title + '_repeatability_vs_zptuncertainty.png')

    plt.clf()
    print 'saved zptu'

    cntr = 0
    pltvecx = []
    pltvecy = []
    for sme,sm,ind,r,d,cm,f,fe in zip(starmagerr,starmag,indices,ras,decs,catmag,flux,fluxerr):
        cntr+=1
        if cntr > maxpoints: continue
        if cntr > 20000: continue

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
    plt.xlim(.0003,.01)
    plt.ylim(.0003,.02)

    ax, ay, aystd = dt.bindata(np.array(pltvecx),np.array(pltvecy), np.arange(.0005,.005, .0001), window=.0001,dontrootn=True)
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
    for sme, sm, ind, r, d, cm, f, fe in zip(starmagerrzpt, starmag, indices, ras, decs, catmag, flux, fluxerr):
        cntr += 1
        if cntr > maxpoints: continue
        if cntr > 100000: continue

        # print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        # print starmag[indices == ind]
        # raw_input()
        starww = starmag[np.isclose(ras, r, rtol=1.e-5) & np.isclose(decs, d, rtol=1.e-5) & (catmag == cm)]
        repeatability = np.std(starww) / np.sqrt(len(starww))
        # repeatability = np.std(starmag[indices == ind])
        if len(starww) > 5.:
            # if repeatability < .3:
            plt.scatter(sme, repeatability, alpha=.3, color='black')
            pltvecy.append(repeatability)
            pltvecx.append(sme)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Photometric Error + Zpt Scatter')
    plt.ylabel('Repeatability')
    plt.xlim(.0003, .02)
    plt.ylim(.0003, .02)

    ax, ay, aystd = dt.bindata(np.array(pltvecx), np.array(pltvecy), np.arange(.0003, .02, .0001), window=.0001,
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
        # print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        # print starmag[indices == ind]
        # raw_input()
        starww = starmag[np.isclose(ras, r, rtol=1.e-5) & np.isclose(decs, d, rtol=1.e-5) & (catmag == cm)]
        repeatability = np.std(starww)/np.sqrt(len(starww))
        rep[i] = repeatability*np.sqrt(len(starww))
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
        if cntr > 10000: continue
        # print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        # print starmag[indices == ind]
        # raw_input()
        starww = starmag[np.isclose(ras, r, rtol=1.e-5) & np.isclose(decs, d, rtol=1.e-5) & (catmag == cm)]
        repeatability = np.std(starww) / np.sqrt(len(starww))
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
        if cntr > 100: continue
        # print starmag[np.isclose(ras,r,rtol=1.e-5) & np.isclose(decs,d,rtol=1.e-5) & (catmag == cm)]
        # print starmag[indices == ind]
        # raw_input()
        starww = starmag[np.isclose(ras, r, rtol=1.e-5) & np.isclose(decs, d, rtol=1.e-5) & (catmag == cm)]
        repeatability = np.std(starww) / np.sqrt(len(starww))
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

    ax2.hist(rv[np.isfinite(rv)], bins=np.arange(-10.1, 10, .2), normed=True, label='RMS: ' + str(round(rms, 3))
             , orientation='horizontal',color='black')
    # label='RMS: ' + str(round(rms, 3)) + '\nChiSq (3sig cut) ' + str(round(chisq, 3)) + '\nMedian ' + str(
    #   round(np.median(d), 3)) + ' +- ' + str(round(np.std(d), 3)),

    import matplotlib.mlab as mlab
    import math
    mean = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.arange(-5, 5, .1)
    ax2.plot(mlab.normpdf(x, mean, sigma), x, color='grey', label='Gaussian Normal')

    ax2.set_ylim(-4, 4)
    ax2.set_xlim(0, .5)
    # .xlabel('STDEV')
    # plt.ylabel('Normalized Count')
    ax2.legend(fontsize='x-small',loc = 'center right', bbox_to_anchor = (1.2, 1.3))
    # plt.savefig('stdresid.png')

    # plt.clf()

    ax1.scatter(catmag, rv, alpha=.02, color='black')
    ax, ay, aystd = dt.bindata(catmag, rv, np.arange(min(catmag), max(catmag), .1), window=.3)
    ax1.plot([min(catmag), max(catmag)], [0, 0], color='grey')
    ax1.plot(ax, ay, linewidth=3, color='orange', label='SMP')
    ax1.plot(ax, ay + aystd, linewidth=2, color='orange', linestyle='--', label='SMP')
    ax1.plot(ax, ay - aystd, linewidth=2, color='orange', linestyle='--', label='SMP')

    # ax1.errorbar(ax, ay, aystd, markersize=20, color='green', fmt='o', label='SMP')

    ax1.set_xlim(16., max(catmag))
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

    ax5.hist(fresid, bins=np.arange(-.155, .15, .001), color='black', orientation='horizontal')
    ww = abs(fresid) < .1
    ax4.scatter(catmag[ww], fresid[ww], alpha=.02, color='black')
    ax, ay, aystd = dt.bindata(catmag[ww], fresid[ww],
                               np.arange(16., max(catmag[ww]), .1), window=1.)
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
    ax3.plot(ax, ayrms, color='orange', label='ZPT Scat + ChisMin', linewidth=3,alpha=.95)
    ax3.plot(ax, ax * 0 + 1., linestyle='--', color='grey')
    ax3.legend(fontsize='x-small',loc = 'center right', bbox_to_anchor = (1.4, 0.8))
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
    ax4.set_xlabel('Cat Mag')
    ax5.set_xlabel('Counts')
    ax3.set_ylabel('RMS')
    ax3.set_title(filter+' band')

    ax4.set_ylabel('(fitflux - catflux)/catflux')

    ax3.set_xlim(ax1.get_xlim())
    ax3.set_ylim(.5,2.)
    ax2.set_ylim(ax1.get_ylim())
    ax5.set_ylim(ax4.get_ylim())
    ax2.xaxis.set_major_formatter(nullfmt)
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
    resultsdir = '/project/projectdirs/des/djbrout/115sim/'
    #resultsdir= './working/'
    #resultsdir= '/export/scratch0/ps1sn1/data/v10.0/GPC1v3/eventsv1/smpworkspace/PS_TEST1/'
    #resultsdir = './workingsimnosn'
    isfermigrid = False
    cacheddata = False
    #cd = '/project/projectdirs/des/djbrout/v67pixshift//summary_results.npz'
    #cd = '/pnfs/des/scratch/pysmp/smp_02_simnosnnoskyerr/np_data/summary_results.npz'
    import sys, getopt

    try:
        args = sys.argv[1:]

        opt, arg = getopt.getopt(
            args, "fd:rd:cd:cdf:b",
            longopts=["fakedir=", "resultsdir=", "cacheddata", "cashedfile=","filter="])

    except getopt.GetoptError as err:
        print "No command line argument    s"

    filter = 'r'

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

    print filter
    tfield = 'SN-S2'
    cd = '/global/cscratch1/sd/dbrout/v6/summary_results_'+tfield+'_'+filter+'.npz'
    go(fakedir,resultsdir,cacheddata,cd,filter,tfield)