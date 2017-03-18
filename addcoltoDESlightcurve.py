
import numpy as np
import os
from copy import copy
import dilltools as dt

def addtolightcurve(lightcurvefile,saveloc,mjd,flux,fluxerr,zpt,zptrms,chisq,sky,skyerr,flag,filt=None,saveinplace=False):

    if not os.path.exists(os.path.basename(saveloc)):
        os.makedirs(os.path.basename(saveloc))

    if saveinplace:
        origfile = open(saveloc, 'r')
        lines = origfile.readlines()
        origfile.close()
    else:
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

    flag = np.array(flag,dtype='int')
    fix = copy(flux)

    fix[fix == 0.] = int(1)
    fix[fix != 1] = int(0)
    fix = np.array(fix,dtype='int')



    #zp = np.array(zp)
    for line in lines:
        if len(line.replace('#','').split(' ')) == 44:
            continue
        elif line.split(' ')[0] == 'VARNAMES:':
            line = line.strip()+' SMP_FLUX SMP_FLUXERR SMP_ZPT SMP_ZPT_RMS SMP_CHISQ SMP_SKY SMP_SKYERR SMP_FIX SMP_FLAG\n'
        elif line.split(' ')[0] == 'OBS:':
            print line

            if filt is None:
                line = line.strip() + ' -999 -999 -999 -999 -999 -999 -999\n'
            id = int(line.split()[1])
            tmjd = round(float(line.split()[3]),3)
            band = line.split()[4]
            ww = (np.round(mjd,3) == tmjd) & (filt == band)
            if fluxerr[ww] > 0:
                line = line.strip() + ' ' + str(round(flux[ww][0], 3)) + ' ' + str(round(fluxerr[ww][0], 3)) + \
                       ' 31. '+str(round(zptrms[ww][0], 3))+\
                       ' '+str(round(chisq[ww][0], 3))+ \
                       ' ' + str(round(sky[ww][0], 3)) + ' ' + str(round(skyerr[ww][0], 3)) + \
                       ' ' + str(fix[ww][0]) + ' ' + str(flag[ww][0]) + '\n'

        savefile.write(line)
    savefile.close()



if __name__ == "__main__":
    lcdir = '/project/projectdirs/des/djbrout/pysmp/imglist/spec/'
    resultsdir = '/project/projectdirs/des/djbrout/spec4/'

    savelcdir = resultsdir+'/newlc'



    filts = ['g','r','i','z']

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

    for i, filt in enumerate(filts):
        sne = os.listdir(resultsdir+'/SNe')

        for sn in sne:
            if 'starfits' in sn:
                continue
            lcfile = lcdir+'/'+sn+'.dat'
            smpfile = resultsdir+'/lightcurves/'+sn+'_'+filt+'.smp'
            savelcfile = savelcdir+'/'+sn+'_smp.dat'
            if not os.path.exists(smpfile):
                print 'SMP RESULTS DO NOT EXIST FOR ',smpfile
                continue

            inplace = False
            if i > 0: inplace = True
            sndata = dt.readcol(smpfile,1,2)

            addtolightcurve(lcfile,savelcfile,sndata['MJD'],sndata['FLUX'],sndata['FLUXERR'],
                            sndata['ZPT'], sndata['RMSADDIN'],
                            sndata['CHI2'],sndata['SKY'],sndata['SKYERR'],sndata['SMP_FLAG'],filt=filt,saveinplace=inplace)
            print 'SAVED SUCCESSFULLY',filt,savelcfile,
            raw_input()




#addtolightcurve('testlc.dat','./testdats/','testcol',,[888,777,000,111],[8,8,8,8],[31.,31.,31.,31.])

       
