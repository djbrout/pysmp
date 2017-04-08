
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

def addtolightcurve(lightcurvefile,saveloc,mjd,flux,fluxerr,zpt,zptrms,chisq,sky,skyerr,flag,zptfiles,idobs,filt=None,saveinplace=False):

    if not os.path.exists(os.path.basename(saveloc)):
        #print 'making'
        os.makedirs(os.path.basename(saveloc))

    if saveinplace:
        try:
            origfile = open(saveloc, 'r')
        except:
            #print 'doesnt exist'
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



    #zp = np.array(zp)
    for line in lines:
        #print line
        #print len(line.replace('#', '').split()),line
        wline = line
        #if saveinplace:
        #    print len(line.replace('#', '').split()),line
        #    #raw_input()
        if len(line.replace('#','').split()) == 27:
            pass
        elif line.split(' ')[0] == 'VARNAMES:':
            wline = line.strip()+' SMP_FLUX SMP_FLUXERR SMP_FLUX_ZPT SMP_FIT_ZPT SMP_FIT_ZPT_STD SMP_CHISQ SMP_SKY SMP_SKYERR SMP_FIX SMP_FLAG\n'
        elif line.split(' ')[0] == 'OBS:':
            #print len(line.replace('#', '').split())
            if filt is None:
                wline = line.strip() + ' -999 -999 -999 -999 -999 -999 -999 -999 -999 '+str(int(FAILED_SMP_FLAG))+'\n'
            #id = int(line.split()[1])
            #tmjd = float(line.split()[3])
            band = line.split()[4]
            #texpnum = line.split()[12].split('/')[2].split('_')[1]
            tidobs = float(line.split()[1])
            ww = np.isclose(idobs,tidobs,atol=0.1) & (filt == band)
            #ww = np.core.defchararray.find(expnums, texpnum) != -1
            #print mjd,tmjd
            #raw_input()
            #print len(fluxerr[ww])
            if len(fluxerr[ww]) == 1:
                zptdata = np.load(zptfiles[ww][0])
                fit_zpt = zptdata['fit_zpt']
                fit_zpt_std = zptdata['fit_zpt_std']
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
                wline = line.strip() + ' ' + str(round(flux[ww][0], 3)) + ' ' + str(round(fluxerr[ww][0], 3)) + \
                       ' 31. '+str(round(fit_zpt, 3))+' '+str(round(fit_zpt_std, 3))+ \
                       ' '+str(round(chisq[ww][0], 3))+ \
                       ' ' + str(round(tsky, 3)) + ' ' + str(round(tskyerr, 3)) + \
                       ' ' + str(fix[ww][0]) + ' ' + str(int(thisflag)) + '\n'

                #print line
        savefile.write(wline)
    savefile.close()
    return True
    #raw_input()



if __name__ == "__main__":
    lcdir = '/project/projectdirs/des/djbrout/pysmp/imglist/all/'
    resultsdir = '/project/projectdirs/des/djbrout/109sim/'

    savelcdir = resultsdir+'/SMP_RAW_SIM_v1_1'



    filts = ['g','r','i','z',None]

    import sys, getopt

    try:
        args = sys.argv[1:]

        opt, arg = getopt.getopt(
            args, "fd:rd:cd:cdf:b",
            longopts=["lcdir=", "resultsdir=", "savelcdir"])

    except getopt.GetoptError as err:
        print "No command line argument    s"

    for o, a in opt:
        if o in ["-h", "--help"]:
            print __doc__
            sys.exit(0)
        elif o in ["--lcdir"]:
            lcdir = a
        elif o in ["-rd", "--resultsdir"]:
            resultsdir = a
        elif o in ["--savelcdir"]:
            savelcdir = a

    if not os.path.exists(os.path.basename(savelcdir)):
        os.mkdir(os.path.basename(savelcdir))
    if not os.path.exists(savelcdir):
        os.mkdir(savelcdir)

    readme = open(savelcdir + '/' + savelcdir.split('/')[-1] + '.README', 'w')
    readme.write(readmetext)
    readme.close()
    #raw_input()
    snlist = open(savelcdir + '/' + savelcdir.split('/')[-1] + '.LIST', 'w')

    for i, filt in enumerate(filts):

        sne = os.listdir(resultsdir+'/SNe')

        for sn in sne[:]:
            #print sn
            if 'starfits' in sn:
                continue

            lcfile = lcdir+'/'+sn+'.dat'
            if not filt is None:
                smpfile = resultsdir+'/lightcurves/'+sn+'_'+filt+'.smp'
            else:
                smpfile = resultsdir+'/lightcurves/'+sn+'_g.smp'
            savelcfile = savelcdir+'/'+sn+'_smp.dat'
            if not os.path.exists(smpfile):
                print 'SMP RESULTS DO NOT EXIST FOR ',smpfile
                continue

            inplace = False
            if i > 0: inplace = True
            sndata = dt.readcol(smpfile,1,2)
            #print sndata.keys()

            if True:
                #print 'adding'
                successful = addtolightcurve(lcfile,savelcfile,sndata['MJD'],sndata['FLUX'],sndata['FLUXERR'],
                            sndata['ZPT'], sndata['RMSADDIN'],
                            sndata['CHI2'],sndata['SKY'],sndata['SKYERR'],sndata['SMP_FLAG'],sndata['ZPTFILE'],
                                sndata['ID_OBS'],filt=filt,saveinplace=inplace)
                print 'SAVED SUCCESSFULLY',filt,savelcfile,'\n'
                if filt == None and successful:
                    snlist.write(sn + '_smp.dat\n')
            #except:
            #    print 'SMP RESULTS DO NOT EXIST FOR ', smpfile
            #raw_input()
    snlist.close()
    print 'cd '+resultsdir+'\n tar -zcf '+savelcdir.split('/')[-1]+'.tar.gz '+savelcdir.split('/')[-1]+'/'
    os.popen('cd '+resultsdir+'\n tar -zcf '+savelcdir.split('/')[-1]+'.tar.gz '+savelcdir.split('/')[-1]+'/')


#addtolightcurve('testlc.dat','./testdats/','testcol',,[888,777,000,111],[8,8,8,8],[31.,31.,31.,31.])

       
