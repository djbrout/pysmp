
import numpy as np
import os
from copy import copy
import dilltools as dt
import sys

readmetext = '## Example output lightcurve file\n\
Contains the original forced photometry lightcurve file structure with the addition of several columns\n\
```\n\
SMP_FLUX - fit flux (at zpt of 31)\n\
SMP_FLUXERR - fit flux error (at zpt of 31)\n\
SMP_FLUX_ZPT - zeropoint of smp flux (always 31)\n\
SMP_FIT_ZPT - smp fit zeropoint which is used to scale all fluxes to 31\n\
SMP_FIT_ZPT_STD - uncertainty in smp fit zeropoint which is used to scale all fluxes to 31\n\
SMP_CHISQ - (sim_stamp - image_stamp)^2/(err)^2\n\
SMP_SKY - fit value for the sky\n\
SMP_SKYERR - uncertainty in the sky\n\
SMP_FIX - 1 means SN flux was fixed to zero in smp fit\n\
SMP_FLAG - 1 means this epoch was flagged for some reason inside smp pipeline\n\
```\n\
(just to be clear: ALL FLUXES, SKYS, SKYERRS, ETC... ARE REPORTED AT A ZEROPOINT OF 31)\n\
-999 means missing data\n\
\n\
##Flag Bit Definitions\
> CHISQ_FLAG = 2 #none have been set\
> BADSKY_FLAG = 4 # FITSKY > 999.\
> BADSKYERR_FLAG = 8 # FIT_SKYERR > .5\
> BADZPT_FLAG = 16 # 27 > FITZPT or FITZPT > 35\
> BADSCATTER_FLAG = 32# the scatter about the zeropoint greater that .2 mags\
> SMP_PIPELINE_FLAG = 32768 # A flag was set inside the SMP pipeline (ex. bad PSF, too much masking, etc...)\
> DIDNT_ENTER_SMP_FLAG = 65536 #These epochs never even made it into the scene modeling pipeline (for example the image never existed on NERSC)\
\
> 1 = Error : SN Coordinates %s,%s are not within image\
> 2 = Error: Bad PSF file\
> 4 = Error: Bad ZPT\
> 8 = Error: SKY > SKYMAX_band\
> 16 = Error: SKY < SKYMIN_band\
> 32 = Error: SKYERR > SKYERRMAX_band\
> 64 = Error: xsn > 25 and ysn > 25 and xsn < snparams.nxpix-25 and ysn < snparams.nypix-25\
> 128 = Error: Image is empty\
> 256 = Error: did not pass FWHM cut\
> 512 = Error: N zpt fit stars too small \
> 1024 = Error: ZPT_SCATTER/SQRT(N) LARGER THAN .01 MAGS\
> 2048 = Error: Scale factor is not finite\
```\n'
#

CHISQ_FLAG = 2
PIPELINE_FLAG = 1
BADSKY_FLAG = 2
BADSKYERR_FLAG = 4
BADZPT_FLAG = 16
BADZPTERR_FLAG = 32
DONTFIT_FLAG = 32768
FAILED_SMP_FLAG = 65536

#dofakefilt,dofakemjd,dofakemag,dofakera,dofakedec = np.loadtxt('data/grepalldofake_'+filter+'.txt',usecols=(3, 9, 10, 14, 15), unpack=True, dtype='string', skiprows=0)
print 'reading dofake'
#expnum, dofakeccds, dofakefilt2, dofakeid, dofakemjd2, dofakemag2, dofaketflux, dofakeflux, dofakera2, dofakedec2 = np.loadtxt(
#    'data/doFake.out', usecols=(1, 2, 3, 5, 9, 10, 11, 12, 14, 15), unpack=True, dtype='string', skiprows=1)


def addtolightcurve(lightcurvefile,saveloc,mjd,flux,fluxerr,zpt,zptrms,chisq,sky,skyerr,flag,zptfiles,idobs,pkmjd,imfiles,dflag,gain,
                    hostsbfluxcals,zpterr,
                    dofakes=False,faketrueflux=False,filt=None,saveinplace=False,mjdplus=80.,mjdminus=25.):
    pkmjd = float(pkmjd)
    idobs=np.array(idobs,dtype='int')
    #print 'inside'
    if not os.path.exists(os.path.basename(saveloc)):
        #print 'making'
        os.makedirs(os.path.basename(saveloc))

    if saveinplace:
        try:
            origfile = open(saveloc, 'r')
        except:
            print 'doesnt exist'
            return False
        lines = origfile.readlines()
        origfile.close()
        #print lines
        #raw_input()
    else:
        #print 'else'
        origfile = open(lightcurvefile, 'r')
        lines = origfile.readlines()
        origfile.close()


    savefile = open(saveloc,'w')

    if dofakes:
        fakeid = lightcurvefile.split('_')[-1].split('.')[0]
        #print fakeid
        fakeisthere = True
        #print dofakeid[0:100]
        if not int(fakeid) in dofakeid:
            fakeisthere = False
        #else:
        #    print 'FOUNDIT'*10


    mjd = np.array(mjd)
    flux = np.array(flux)
    fluxerr = np.array(fluxerr)
    zpt = np.array(zpt)
    zpterr = np.array(zpterr)
    zptrms = np.array(zptrms)
    chisq = np.array(chisq)
    sky = np.array(sky)
    skyerr = np.array(skyerr)
    zptfiles = np.array(zptfiles,dtype='str')
    idobs = np.array(idobs)
    imfiles = np.array(imfiles,dtype='str')
    dflag = np.array(dflag,dtype='int')
    flag = np.array(flag,dtype='int')
    fix = copy(flux)
    gain = np.array(gain)
    hostsbfluxcals = np.array(hostsbfluxcals)
    fix[fix == 0.] = int(1)
    fix[fix != 1] = int(0)
    fix = np.array(fix,dtype='int')


    writelines = ''
    #zp = np.array(zp)
    for line in lines:
        #print line
        #raw_input()
        #print len(line.replace('#', '').split()),line
        wline = line
        #if saveinplace:
        #    print len(line.replace('#', '').split()),line
        #raw_input()
        try:
        #if True:
            if len(line.replace('#','').split()) == 28:
                wline = line
                print 'here'*100
                #pass
            elif line.split(' ')[0] == 'VARNAMES:':
                wline = line.strip()+' SMP_FLUX SMP_FLUXERR SMP_FLUX_ZPT SMP_FIT_ZPT SMP_FIT_ZPT_STD SMP_CHISQ SMP_SKY SMP_SKYERR SMP_FIX SMP_FLAG SMP_GAIN \n'
            elif line.split(' ')[0] == 'OBS:':
                #print line
                #raw_input()
                #print line.split(' ')[3]
                #raw_input()
                # if float(line.split(' ')[3]) < pkmjd - mjdminus:
                #     wline = ''
                # elif float(line.split(' ')[3]) > pkmjd + mjdplus:
                #     wline = ''
                # else:
                if True:

                    #print len(line.replace('#', '').split())
                    #band = line.split()[4]

                    #if filt is None:
                    #    wline = line.strip() + ' -999 -999 -999 -999 -999 -999 -999 -999 -999 '+str(int(FAILED_SMP_FLAG))+'\n'
                    #elif band != filt:
                    #    wline = line
                    #    #continue
                    #print line
                    #raw_input()
                    #else:
                    #tidobs = int(line.split()[1])
                    #tmjd = float(line.split()[3])
                    tim = '/global/cscratch1/sd/dbrout/v7/'+line.split()[13]#for fakes do 13, spec 12
                    tim = tim.split('/')[-1]
                    #print tim
                    # print tidobs
                    # if float(tidobs) == 414.:
                    #     print idobs
                    #
                    #     print 'stopped at this night'
                    #     raw_input()
                    #print tidobs,idobs
                    #raw_input()
                    #print imfiles,tim
                    #raw_input()
                    good = False
                    #print imfiles[:10]
                    #print tim
                    # # print line
                    # print dflag
                    #raw_input()
                    if tim.strip() in imfiles:
                        good = True
                        #print 'good'
                    elif tim.strip().replace('p1','Y1') in imfiles:
                        tim = tim.strip().replace('p1','Y1')
                        good = True
                    elif tim.strip().replace('p1','Y2') in imfiles:
                        tim = tim.strip().replace('p1','Y2')
                        good = True
                    elif tim.strip().replace('p1','Y3') in imfiles:
                        tim = tim.strip().replace('p1','Y3')
                        good = True
                    # print good
                    # raw_input()
                    if good:
                        #print 'inside'
                        #raw_input()
                        #ww = np.isclose(idobs,tidobs,atol=0.005)# & (filt == band)
                        ww = imfiles == tim
                        #print fluxerr[ww]
                        #raw_input()
                        # if float(tidobs) == 414.:
                        #     print 'stopped at this night'
                        #     raw_input()
                        keepgoing = True
                       # print len(fluxerr[ww])
                        if len(fluxerr[ww]) == 1:
                            #print 'here',dofakes
                            # try:
                            #     #print zptfiles[ww][0]
                            #     zptdata = np.load(zptfiles[ww][0].replace('v6','v7'))
                            # except:
                            #     print 'excepted',line
                            #     wline = line.strip() + ' -999 -999 -999 -999 -999 -999 -999 -999 -999 -999 ' + str(
                            #             int(dflag[ww][0])) + '\n'
                            #     #print len(wline.split())
                            #     keepgoing = False
                            if (faketrueflux) & (keepgoing):
                                if fakeisthere:
                                    tmag = float(line.split()[12])
                                    exn = line.split()[13].split('/')[-1].split('_')[1]

                                    expn = (dofakeexpnum == float(exn))
                                    #print tmag,exn,
                                    dfw = dofakeid == int(fakeid)
                                    www = expn & dfw
                                    #print dofakemag2[www]
                                    if not len(dofakemag2[www]) > 0:
                                        tmag = 99.
                                        tzpt = 31.
                                        flux_zpt = 31.
                                    else:
                                        tzpt = dofakezpt[www][0]
                                        flux_zpt = dofakezpt[www][0]
                                    #tzpt = float(line.split()[7])

                                else:
                                    tmag = 99.
                                    tzpt = 31.
                                    flux_zpt = 31.

                                fit_zpt = zptdata['fit_zpt']
                                fit_zpt_std = 0.
                                #tflux = 10 ** (.4 * (tzpt - tmag ))
                                tflux = 10**(.4*(tzpt - tmag))
                                tflux *= 10 ** (-1. * .4 * (fit_zpt - tzpt))
                                #print tflux,fit_zpt,tzpt
                                #print 'hereeee'
                                #raw_input()
                                tfluxerr = tflux**.5
                                fit_zpt = tzpt
                                #tfluxerr *= 10 ** (-1 * .4 * (fitzpt - fakezpt))

                                #print exn, tzpt, tmag, tflux
                                #raw_input()
                            elif (dofakes) & (keepgoing):
                                sys.exit('dofakes no longer allowed')
                            #     if fakeisthere:
                            #         #tmag = float(line.split()[12])
                            #         exn = line.split()[13].split('/')[-1].split('_')[1]
                            #
                            #         expn = (dofakeexpnum == float(exn))
                            #         #print tmag,exn,
                            #         dfw = dofakeid == int(fakeid)
                            #         www = expn & dfw
                            #         #print dofakemag2[www]
                            #         if not len(dofakemag2[www]) > 0:
                            #             #tmag = 99.
                            #             tzpt = 31.
                            #             #flux_zpt = 31.
                            #         else:
                            #             tzpt = dofakezpt[www][0]
                            #             #flux_zpt = dofakezpt[www][0]
                            #         #tzpt = float(line.split()[7])
                            #
                            #     else:
                            #         #tmag = 99.
                            #         tzpt = 31.
                            #         #flux_zpt = 31.
                            #     tflux = flux[ww][0]
                            #     fit_zpt = zptdata['fit_zpt']
                            #     fit_zpt_std = zptdata['fit_zpt_std']
                            #     flux_zpt = 31.
                            #     tfluxerr = fluxerr[ww][0]
                            #     tflux *= 10 ** (-1 * .4 * (fit_zpt - tzpt))
                            #     tfluxerr *= 10 ** (-1 * .4 * (fit_zpt - tzpt))

                            elif keepgoing:
                                #print 'here'
                                tflux = flux[ww][0]
                                #fit_zpt = zptdata['fit_zpt']
                                #fit_zpt_std = zptdata['fit_zpt_std']
                                flux_zpt = 31.
                                #raw_input()
                                hostsbflux = hostsbfluxcals[ww][0]*10**(.4*(31.-27.5))

                                tfluxerr = np.sqrt(fluxerr[ww][0]**2 + max(0,flux[ww][0]) + hostsbflux + (zpterr[ww][0]*max(0,flux[ww][0]))**2)
                                #print tflux

                            if keepgoing:
                                tsky = sky[ww][0] - 10000.*10**(.4*(31.-zpt[ww][0]))
                                tskyerr = skyerr[ww][0]
                                thisflag = 0
                                # if flag[ww][0] == 1:
                                #     thisflag = DONTFIT_FLAG
                                # #if chisq[ww][0] > 1.0:
                                # #    thisflag = DONTFIT_FLAG
                                # if abs(tsky) > 1000:
                                #     thisflag = BADSKY_FLAG
                                # if abs(tskyerr) < .5:
                                #     thisflag = BADSKYERR_FLAG
                                # if (fit_zpt < 27.) | (fit_zpt > 35.):
                                #     thisflag = BADZPT_FLAG
                                # if (fit_zpt_std > 0.2):
                                #     thisflag = BADZPTERR_FLAG
                                #print thisflag,chisq[ww][0]
                                if fix[ww][0] == 1:
                                    tdflag = 1
                                else:
                                    tdflag = dflag[ww][0]
                                wline = line.strip() + ' ' + str(round(tflux, 5)) + ' ' + str(round(tfluxerr, 5)) + \
                                       ' '+str(round(flux_zpt, 5))+' '+str(round(zpt[ww][0], 5))+' '+str(round(zpterr[ww][0], 5))+ \
                                       ' '+str(round(chisq[ww][0], 5))+ \
                                       ' ' + str(round(tsky, 5)) + ' ' + str(round(tskyerr, 5)) + \
                                       ' ' + str(fix[ww][0]) + ' ' + str(int(tdflag))+ ' ' + str(float(round(gain[ww][0],5))) + '\n'
                                #print len(wline.split())
                            #raw_input()
                            #print wline
                            #raw_input()
                            #print line
                        #raw_input()
                        else:
                            wline = line.strip() + ' -999 -999 -999 -999 -999 -999 -999 -999 -999 -999 ' + str(
                                int(FAILED_SMP_FLAG)) + '\n'
                            print mjd,tmjd

                            print 'baddddddd'*100
                            raw_input()
                    else:
                        if float(line.split()[3]) < 57524.371:
                            print 'NOTHERE'*5
                            print line
                        wline = line.strip() + ' -999 -999 -999 -999 -999 -999 -999 -999 -999 -999 ' + str(int(FAILED_SMP_FLAG)) + '\n'

                   # print len(wline.split())
                    #raw_input()
        except:
            e = sys.exc_info()[0]
            #print e
            wline = line.strip() + ' -999 -999 -999 -999 -999 -999 -999 -999 -999 -999 ' + str(
                int(FAILED_SMP_FLAG)) + '\n'
            #print len(wline.split())
        #print len(wline.split())
        #print wline
        writelines += wline
        #savefile.write(wline)
    savefile.write(writelines)
    savefile.close()
    print ''
    return True
    #raw_input()



if __name__ == "__main__":
    print 'started'

    resultsdir = '/project/projectdirs/des/djbrout/redospec/'
    #resultsdir = '/project/projectdirs/des/djbrout/allsim/'
    fakeheader = False
    #if isfake:
    dodiffim = False

    lcdir = '/project/projectdirs/dessn/dbrout/imgList/all/'
    #else:
    #    lcdir = '/project/projectdirs/des/djbrout/pysmp/imglist/spec/'

    savelcdir = None
    fakes = False
    faketrueflux = False

    filts = ['g','r','i','z']

    import sys, getopt

    try:
    #if True:
        args = sys.argv[1:]

        opt, arg = getopt.getopt(
            args, "fd:rd:cd:cdf:b",
            longopts=["index=","lcdir=", "resultsdir=", "savelcdir=","dofakes","faketrueflux","fakeheader","dodiffim"])

    except getopt.GetoptError as err:
        print "No command line arguments"

    index = None

    for o, a in opt:
        #print o
        if o in ["-h", "--help"]:
            print __doc__
            sys.exit(0)
        elif o in ["--lcdir"]:
            lcdir = a
        elif o in ["-rd", "--resultsdir"]:
            resultsdir = a
        elif o in ["--savelcdir"]:
            savelcdir = a
        elif o in ["--dofakes"]:
            fakes = True
        elif o in ["--faketrueflux"]:
            faketrueflux = True
            fakes = True
        elif o in ["--index"]:
            index = int(a)
        elif o in ["--fakeheader"]:
            fakeheader = True
        elif o in ["--dodiffim"]:
            dodiffim = True
    #print fakes
    #raw_input()

    if savelcdir is None:
        savelcdir = resultsdir + '/SMP_SPECv2_2'
        if fakeheader: savelcdir = resultsdir + '/SMP_SIMv2_1'

    if not os.path.exists(os.path.basename(savelcdir)):
        os.mkdir(os.path.basename(savelcdir))
    if not os.path.exists(savelcdir):
        os.mkdir(savelcdir)


    missingfile = resultsdir+'/missing.txt'

    # readme = open(savelcdir + '/' + savelcdir.split('/')[-1] + '.README', 'w')
    # readme.write(readmetext)
    # readme.close()
    if fakeheader:
        os.system('cp fakes.HEADER '+savelcdir + '/' + savelcdir.split('/')[-1] + '.HEADER')
    else:
        os.system('cp specHD.HEADER '+ savelcdir + '/' + savelcdir.split('/')[-1] + '.HEADER')
    #raw_input()
    #snlist = open(savelcdir + '/' + savelcdir.split('/')[-1] + '.LIST', 'w')



    #for i, filt in enumerate(filts):
    #sne = os.listdir(resultsdir + '/lightcurves')
    #tsne = []
    #for sn in sne:
    #    if '.smp' in sn:
    #        print 'here ',sn
    #        tsne.append('_'.join(sn.split('_')[:-1]))
    #sne = np.array(tsne,dtype='str')
    #sne = np.unique(sne)

    #sne = open('data/allslightcurves.txt').readlines()
    sne = open('data/speclist.txt').readlines()
    if fakeheader: sne = open('data/alllightcurves.txt').readlines()


    tsne = []
    for sn in sne:
        print sn
        tsne.append(sn.split('.')[0].split('/')[-1])
    sne = np.array(tsne, dtype='str')
    sne = np.unique(sne)
    if index is None:
        sne = sne
    else:
        index = int(os.environ['$SLURM_JOBID'])
        sne = [sne[index-1]]
    #print tsne
    tsneold = tsne
    if fakes:
        import pandas as pd

        dofakedata = pd.read_csv('data/doFake.out', delim_whitespace=True, header=0)
        # print dofakedata
        dofakeexpnum = dofakedata['EXPNUM'].values
        dofakemag2 = dofakedata['TRUEMAG'].values
        dofaketflux = dofakedata['TRUEFLUXCNT'].values
        dofakeid = dofakedata['FAKEID'].values
        dofakezpt = dofakemag2 + 2.5 * np.log10(dofaketflux)
        print 'done reading dofake'

    a = open(missingfile,'w')
    tsne = []
    numbad = 0
    for sn in sne[:]:
        snbad = False
        tbad = 0
        badfilts = []
        for i, filt in enumerate(filts):
            if 'starfits' in sn:
                continue

            lcfile = lcdir+'/'+sn+'.dat'
            #if not filt is None:
            smpfile = resultsdir+'/lightcurves/'+sn+'_'+filt+'.smp'
            #else:
            #    continue
            savelcfile = savelcdir+'/'+sn+'.dat'
            if not os.path.exists(smpfile):
                #print 'SMP RESULTS DO NOT EXIST FOR ',smpfile
                badfilts.append(filt)
                #a.write('_'.join(smpfile.split('/')[-1].split('.')[0].split('_')[:-1])+' '+filt+' \n')
                #os.system('echo '+sn+' '+filt+' >> '+missingfile)
                snbad = True
                tbad += 1
        if not snbad:
            tsne.append(sn)
        else:
            if tbad < 1:
                for filt in badfilts:
                    a.write(smpfile.split('/')[-1]+' \n')
                numbad += 1

    sne = tsne


    #sne = tsneold
    a.close()
    print 'TOTAL SNe:',len(sne),'Missing SNe:',numbad
    cntr = 0
    donesne = []
    for sn in sne[::-1]:
        print sn
        print '-'*100
        #sn = 'des_real_01248907'
        if dodiffim:
            os.popen('cp '+lcdir+'/'+sn.split('.')[0]+'.dat '+savelcdir+'/')
            donesne.append(sn.split('.')[0].split('_')[-1])
            continue
        mjd = []
        flux = []
        fluxerr = []
        zpt = []
        zpterr = []
        rmsaddin = []
        chi2 = []
        sky = []
        skyerr = []
        smpflag = []
        zptfile = []
        imfiles = []
        idobs = []
        band = []
        dflag = []
        gain = []
        hostsbfluxcals = []
        cntr += 1

        for i, filt in enumerate(filts):
            #for sn in sne[:]:
            #cntr += 1
            #print sn
            if cntr > 10000:
                continue
            if 'starfits' in sn:
                continue

            lcfile = lcdir+'/'+sn+'.dat'
            #if not filt is None:
            smpfile = resultsdir+'/lightcurves/'+sn+'_'+filt+'.smp'
            print smpfile
            #else:
            #    continue
            savelcfile = savelcdir+'/'+sn+'.dat'
            # if not os.path.exists(smpfile):
            #     print 'SMP RESULTS DO NOT EXIST FORR ',smpfile
            #     #os.system('echo '+sn+' '+filt+' >> '+missingfile)
            #     continue
            # else:
            #     pass

                #print 'DOES EXIST',smpfile
            #print lcfile
            pkmjd = open(lcfile).readlines()[10].split()[1]
            if filt == 'g': fi = 1
            if filt == 'r': fi = 2
            if filt == 'i': fi = 3
            if filt == 'z': fi = 4

            hostsbfluxcal = open(lcfile).readlines()[20].split()[fi]
            #print pkmjd
            #raw_input()
            inplace = False
            #if i > 0: inplace = True
            sndata = dt.readcol(smpfile,1,2)
            try:
                a = sndata['MJD']
            except:
                print 'could not grab mjd '*100
                continue

            #print sndata['FLUX'].shape
            mjd.extend(sndata['MJD'])
            flux.extend(sndata['FLUX'])
            fluxerr.extend(sndata['FLUXERR'])
            zpt.extend(sndata['ZPT'])
            zpterr.extend(sndata['ZPTERR'])
            rmsaddin.extend(sndata['RMSADDIN'])
            chi2.extend(sndata['CHI2'])
            sky.extend(sndata['SKY'])
            skyerr.extend(sndata['SKYERR'])
            smpflag.extend(sndata['SMP_FLAG'])
            zptfile.extend(sndata['ZPTFILE'])
            idobs.extend(sndata['ID_OBS'])
            band.extend(sndata['BAND'])
            imfiles.extend(sndata['IMAGE_FILE'])
            dflag.extend(sndata['DESCRIPTIVE_FLAG'])
            gain.extend(sndata['GAIN'])
            # print sndata['DESCRIPTIVE_FLAG']
            hostsbfluxcals.extend(np.array(sndata['GAIN'],dtype='float')*0. + float(hostsbfluxcal))
            #print sndata.keys()
            #raw_input()
        #print len(band),len(idobs),len(sky),len(flux),len(mjd)
        #print np.sort(idobs)
        #print band
        #raw_input()
        #print 'lcfile',lcfile
        #print flux
        if not os.path.exists(savelcfile):

            imfiles = np.array([imf.split('/')[-1].replace('+fakeSN', '') for imf in imfiles], dtype='str')
            #print imfiles
            #raw_input()
            successful = addtolightcurve(lcfile,savelcfile,mjd,flux,fluxerr,
                     zpt, rmsaddin,
                     chi2,sky,skyerr,smpflag,zptfile,
                     idobs,pkmjd,imfiles,dflag,gain,hostsbfluxcals,zpterr, dofakes=fakes, saveinplace=False,faketrueflux=faketrueflux)

        print int(cntr),'SAVED SUCCESSFULLY',savelcfile,'\n'
        donesne.append(sn+'.dat')#.split('.')[0].split('_')[-1])

        #raw_input('stoppppp')
    a = open(savelcdir+'/'+savelcdir.split('/')[-1]+'.README','w')
    a.write(readmetext)
    a.close()

    a = open(savelcdir+'/'+savelcdir.split('/')[-1]+'.LIST','w')
    for sn in donesne:
        a.write(sn+' \n')
    a.close()

        #raw_input()
        #if filt == None and successful:
        #snlist.write(sn + '_smp.dat\n')
            #except:
            #    print 'SMP RESULTS DO NOT EXIST FOR ', smpfile
            #raw_input()
    #snlist.close()

    #print 'cd '+resultsdir+'\n tar -zcf '+savelcdir.split('/')[-1]+'.tar.gz '+savelcdir.split('/')[-1]+'/'

    #os.popen('cd '+resultsdir+'\n tar -zcf '+savelcdir.split('/')[-1]+'.tar.gz '+savelcdir.split('/')[-1]+'/')


#addtolightcurve('testlc.dat','./testdats/','testcol',,[888,777,000,111],[8,8,8,8],[31.,31.,31.,31.])

       
