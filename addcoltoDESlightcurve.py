
import numpy as np
import os
from copy import copy

def addtolightcurve(lightcurvefile,saveloc,mjd,flux,fluxerr,zpt,chisq,sky,skyerr,filt=None,saveinplace=False):

    if not os.path.exists(saveloc):
        os.makedirs(saveloc)

    origfile = open(lightcurvefile, 'r')
    lines = origfile.readlines()
    origfile.close()

    savefile = open(saveloc+'/'+lightcurvefile.split('/')[-1],'w')




    mjd = np.array(mjd)
    flux = np.array(flux)
    fluxerr = np.array(fluxerr)
    zpt = np.array(zpt)
    chisq = np.array(chisq)
    sky = np.array(sky)
    skyerr = np.array(skyerr)
    fix = copy(flux)

    fix[fix == 0.] = int(1)
    fix[fix != 1] = int(0)



    #zp = np.array(zp)
    for line in lines:
        if len(line) == 20:
            continue
        elif line.split(' ')[0] == 'VARNAMES:':
            line = line.strip()+' SMP_FLUX SMP_FLUXERR SMP_ZPT SMP_CHISQ SMP_SKY SMP_SKYERR SMP_FIX\n'
        elif line.split(' ')[0] == 'OBS:':
            if filt is None:
                line = line.strip() + ' -999 -999 -999 -999 -999 -999 -999\n'
            id = int(line.split()[1])
            tmjd = float(line.split()[3])
            band = line.split()[4]
            ww = (mjd == tmjd) & (filt == band)
            if fluxerr[ww] > 0:
                line = line.strip() + ' ' + str(round(flux[ww][0], 3)) + ' ' + str(round(fluxerr[ww][0], 3)) + \
                       ' '+str(round(zpt[ww][0], 3))+' '+str(round(chisq[ww][0], 3))+ \
                       ' ' + str(round(sky[ww][0], 3)) + ' ' + str(round(skyerr[ww][0], 3)) + \
                       ' ' + str(round(fix[ww][0], 3)) + '\n'

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

    if not os.path.exists(savelcdir):
        os.mkdir(savelcdir)

    for i, filt in enumerate(filts):
        sne = os.listdir(resultsdir+'/SNe')
        for sn in sne:
            lcfile = lcdir+'/'+sn+'.dat'
            print lcfile
        raw_input()




#addtolightcurve('testlc.dat','./testdats/','testcol',,[888,777,000,111],[8,8,8,8],[31.,31.,31.,31.])

       
