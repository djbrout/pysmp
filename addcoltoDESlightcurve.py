
import numpy as np
import os

def addtolightcurve(lightcurvefile,saveloc,column_name,mjd,flux,fluxerr,fakemag,filt=None):
    if not os.path.exists(saveloc):
        os.makedirs(saveloc)
    savefile = open(saveloc+'/'+lightcurvefile.split('/')[-1],'w')
    origfile = open(lightcurvefile,'r')
    lines = origfile.readlines()
    mjd = np.array(mjd)
    if not filt is None:
        filt = np.array(filt,dtype='str')
    flux = np.array(flux)
    fluxerr = np.array(fluxerr)
    fakemag = np.array(fakemag)
    #zp = np.array(zp)
    for line in lines:
        if line.split(' ')[0] == 'VARNAMES:':
            line = line.strip()+' SMPFLUX_'+column_name.upper()+' SMPFLUXERR_'+column_name.upper()+' SMPZPT_'+column_name.upper()+'\n'
        elif line.split(' ')[0] == 'OBS:':
            id = int(line.split()[1])
            tmjd = float(line.split()[3])
            band = line.split()[4]
            ww = (mjd == tmjd) & (filt == band)
            if fluxerr[ww] > 0:
                line = line.strip()+' '+str(round(flux[ww][0],2))+' '+str(round(fluxerr[ww][0],2))+' '+str(round(fakemag[ww][0],3))+'\n'
            #else:
            #    line = line.strip()+' -999 -999 -999\n'
        savefile.write(line)
        #print line
    savefile.close()
    origfile.close()
#addtolightcurve('testlc.dat','./testdats/','testcol',,[888,777,000,111],[8,8,8,8],[31.,31.,31.,31.])

       
