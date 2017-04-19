
import numpy as np
import os
from copy import copy
import dilltools as dt


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
## Flag Bit Definitions\n\
```\n\
CHISQ_FLAG = 2\n\
PIPELINE_FLAG = 1\n\
BADSKY_FLAG = 4\n\
BADSKYERR_FLAG = 8\n\
BADZPT_FLAG = 16\n\
BADZPTERR_FLAG = 32\n\
DONTFIT_FLAG = 65536\n\
```\n'
#
# CHISQ_FLAG = 8192
# PIPELINE_FLAG =4096
# BADSKY_FLAG = 16384
# BADSKYERR_FLAG = 32768
# BADZPT_FLAG = 65536
# BADZPTERR_FLAG = 32
DONTFIT_FLAG= 32768
FAILED_SMP_FLAG = 65536

#dofakefilt,dofakemjd,dofakemag,dofakera,dofakedec = np.loadtxt('data/grepalldofake_'+filter+'.txt',usecols=(3, 9, 10, 14, 15), unpack=True, dtype='string', skiprows=0)
print 'reading dofake'
#expnum, dofakeccds, dofakefilt2, dofakeid, dofakemjd2, dofakemag2, dofaketflux, dofakeflux, dofakera2, dofakedec2 = np.loadtxt(
#    'data/doFake.out', usecols=(1, 2, 3, 5, 9, 10, 11, 12, 14, 15), unpack=True, dtype='string', skiprows=1)
import pandas as pd
dofakedata = pd.read_csv('data/doFake.out', delim_whitespace=True, header=0)
#print dofakedata
dofakeexpnum = dofakedata['EXPNUM'].values
dofakemag2 = dofakedata['TRUEMAG'].values
dofaketflux = dofakedata['TRUEFLUXCNT'].values
dofakeid = dofakedata['FAKEID'].values

#dofakeexpnum = np.array(expnum, dtype='float')
#dofakemag2 = np.array(dofakemag2, dtype='float')
#dofaketflux = np.array(dofaketflux, dtype='float')
dofakezpt = dofakemag2 + 2.5 * np.log10(dofaketflux)
#dofakeid = np.array(dofakeid, dtype='float')
print 'done reading dofake'

def addtolightcurve(lightcurvefile,saveloc,mjd,flux,fluxerr,zpt,zptrms,chisq,sky,skyerr,flag,zptfiles,idobs,
                    dofakes=False,faketrueflux=False,filt=None,saveinplace=False):
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
    zptrms = np.array(zptrms)
    chisq = np.array(chisq)
    sky = np.array(sky)
    skyerr = np.array(skyerr)
    zptfiles = np.array(zptfiles,dtype='str')
    idobs = np.array(idobs)

    flag = np.array(flag,dtype='int')
    fix = copy(flux)

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
        if len(line.replace('#','').split()) == 28:
            wline = line
            #pass
        elif line.split(' ')[0] == 'VARNAMES:':
            wline = line.strip()+' SMP_FLUX SMP_FLUXERR SMP_FLUX_ZPT SMP_FIT_ZPT SMP_FIT_ZPT_STD SMP_CHISQ SMP_SKY SMP_SKYERR SMP_FIX SMP_FLAG\n'
        elif line.split(' ')[0] == 'OBS:':
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
            tidobs = float(line.split()[1])
            print tidobs,idobs
            raw_input()
            if int(tidobs) in idobs:
                ww = np.isclose(idobs,tidobs,atol=0.1)# & (filt == band)
                #print fluxerr[ww]
                #raw_input()
                keepgoing = True
                if len(fluxerr[ww]) == 1:
                    #print 'here',dofakes
                    try:
                        zptdata = np.load(zptfiles[ww][0])
                    except:
                        wline = line.strip() + ' -999 -999 -999 -999 -999 -999 -999 -999 -999 ' + str(
                                int(FAILED_SMP_FLAG)) + '\n'
                        keepgoing = False
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

                        fit_zpt = tzpt
                        fit_zpt_std = 0.
                        #tflux = 10 ** (.4 * (tzpt - tmag ))
                        tflux = 10**(.4*(tzpt - tmag))
                        #tflux *= 10 ** (-1 * .4 * (fitzpt - fakezpt))
                        tfluxerr = tflux**.5
                        #tfluxerr *= 10 ** (-1 * .4 * (fitzpt - fakezpt))

                        #print exn, tzpt, tmag, tflux
                        #raw_input()
                    elif (dofakes) & (keepgoing):
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
                        tflux = flux[ww][0]
                        fit_zpt = zptdata['fit_zpt']
                        fit_zpt_std = zptdata['fit_zpt_std']
                        flux_zpt = 31.
                        tfluxerr = fluxerr[ww][0]
                        tflux *= 10 ** (-1 * .4 * (fit_zpt - tzpt))
                        tfluxerr *= 10 ** (-1 * .4 * (fit_zpt - tzpt))

                    elif keepgoing:
                        tflux = flux[ww][0]
                        fit_zpt = zptdata['fit_zpt']
                        fit_zpt_std = zptdata['fit_zpt_std']
                        flux_zpt = 31.
                        tfluxerr = fluxerr[ww][0]
                        #print tflux

                    if keepgoing:
                        tsky = sky[ww][0] - 10000.*10**(.4*(31.-fit_zpt))
                        tskyerr = skyerr[ww][0]
                        thisflag = 0
                        if flag[ww][0] == 1:
                            thisflag = DONTFIT_FLAG
                        if chisq[ww][0] > 1.25:
                            thisflag = DONTFIT_FLAG
                        if abs(tsky) > 1000:
                            thisflag = DONTFIT_FLAG
                        if abs(tskyerr) < .5:
                            thisflag = DONTFIT_FLAG
                        if (fit_zpt < 27.) | (fit_zpt > 35.):
                            thisflag = DONTFIT_FLAG
                        if (fit_zpt_std > 0.2):
                            thisflag = DONTFIT_FLAG
                        #print thisflag,chisq[ww][0]
                        wline = line.strip() + ' ' + str(round(tflux, 3)) + ' ' + str(round(tfluxerr, 3)) + \
                               ' '+str(round(flux_zpt, 3))+' '+str(round(fit_zpt, 3))+' '+str(round(fit_zpt_std, 3))+ \
                               ' '+str(round(chisq[ww][0], 3))+ \
                               ' ' + str(round(tsky, 3)) + ' ' + str(round(tskyerr, 3)) + \
                               ' ' + str(fix[ww][0]) + ' ' + str(int(thisflag)) + '\n'

                    #print wline
                    #raw_input()
                    #print line
                #raw_input()
            else:
                wline = line.strip() + ' -999 -999 -999 -999 -999 -999 -999 -999 -999 ' + str(
                    int(FAILED_SMP_FLAG)) + '\n'
        writelines += wline
        #savefile.write(wline)
    savefile.write(writelines)
    savefile.close()
    return True
    #raw_input()



if __name__ == "__main__":
    print 'started'
    lcdir = '/project/projectdirs/des/djbrout/pysmp/imglist/all/'
    resultsdir = '/project/projectdirs/des/djbrout/114sim/'

    savelcdir = resultsdir+'/SMP_RAW_SIM_v1_3'
    fakes = False
    faketrueflux = False


    filts = ['g','r','i','z']

    import sys, getopt

    try:
        args = sys.argv[1:]

        opt, arg = getopt.getopt(
            args, "fd:rd:cd:cdf:b",
            longopts=["lcdir=", "resultsdir=", "savelcdir=","dofakes","faketrueflux"])

    except getopt.GetoptError as err:
        print "No command line argument    s"

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

    #print fakes
    #raw_input()

    if not os.path.exists(os.path.basename(savelcdir)):
        os.mkdir(os.path.basename(savelcdir))
    if not os.path.exists(savelcdir):
        os.mkdir(savelcdir)

    readme = open(savelcdir + '/' + savelcdir.split('/')[-1] + '.README', 'w')
    readme.write(readmetext)
    readme.close()
    #raw_input()
    snlist = open(savelcdir + '/' + savelcdir.split('/')[-1] + '.LIST', 'w')



    #for i, filt in enumerate(filts):
    sne = os.listdir(resultsdir + '/SNe')
    cntr = 0
    for sn in sne[:]:
        mjd = []
        flux = []
        fluxerr = []
        zpt = []
        rmsaddin = []
        chi2 = []
        sky = []
        skyerr = []
        smpflag = []
        zptfile = []
        idobs = []
        cntr += 1
        for i, filt in enumerate(filts):
            #for sn in sne[:]:
            #cntr += 1
            #print sn
            #if cntr > 50:
            #    continue
            if 'starfits' in sn:
                continue

            lcfile = lcdir+'/'+sn+'.dat'
            #if not filt is None:
            smpfile = resultsdir+'/lightcurves/'+sn+'_'+filt+'.smp'
            #else:
            #    continue
            savelcfile = savelcdir+'/'+sn+'_smp.dat'
            if not os.path.exists(smpfile):
                print 'SMP RESULTS DO NOT EXIST FOR ',smpfile
                continue

            inplace = False
            #if i > 0: inplace = True
            sndata = dt.readcol(smpfile,1,2)

            mjd.extend(sndata['MJD'])
            flux.extend(sndata['FLUX'])
            fluxerr.extend(sndata['FLUXERR'])
            zpt.extend(sndata['ZPT'])
            rmsaddin.extend(sndata['RMSADDIN'])
            chi2.extend(sndata['CHI2'])
            sky.extend(sndata['SKY'])
            skyerr.extend(sndata['SKYERR'])
            smpflag.extend(sndata['SMP_FLAG'])
            zptfile.extend(sndata['ZPTFILE'])
            print sndata['ID_OBS']
            idobs.extend(sndata['ID_OBS'])


        print idobs
        raw_input()
        successful = addtolightcurve(lcfile,savelcfile,mjd,flux,fluxerr,
                     zpt, rmsaddin,
                     chi2,sky,skyerr,smpflag,zptfile,
                     sndata['ID_OBS'], dofakes=fakes, saveinplace=False,faketrueflux=faketrueflux)

        print int(cntr),'SAVED SUCCESSFULLY',savelcfile,'\n'
        #if filt == None and successful:
        snlist.write(sn + '_smp.dat\n')
            #except:
            #    print 'SMP RESULTS DO NOT EXIST FOR ', smpfile
            #raw_input()
    snlist.close()
    print 'cd '+resultsdir+'\n tar -zcf '+savelcdir.split('/')[-1]+'.tar.gz '+savelcdir.split('/')[-1]+'/'
    os.popen('cd '+resultsdir+'\n tar -zcf '+savelcdir.split('/')[-1]+'.tar.gz '+savelcdir.split('/')[-1]+'/')


#addtolightcurve('testlc.dat','./testdats/','testcol',,[888,777,000,111],[8,8,8,8],[31.,31.,31.,31.])

       
