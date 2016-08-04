
import numpy as np
import os

def wraplightcurves(listfile,filedir,npzdir,lightcurveoutdir,filt=None):
    files = open(listfile).readlines()
    for fl in files:
        f = fl.split("/")[-1].strip()
        print f
        print filedir
        print os.path.join(npzdir,f.strip(".psmp")+'_'+filt+'_withSN.npz')
        print os.path.exists(os.path.join(npzdir,f.strip(".psmp")+'_'+filt+'_withSN.npz'))
        #raw_input()
        if os.path.exists(os.path.join(npzdir,f.strip(".psmp")+'_'+filt+'_withSN.npz')):
            print "it existssssssssss"
            fin = os.path.join(filedir,f)
            fout = os.path.join(filedir,f+'_dillon')
            lcout = os.path.join(lightcurveoutdir,f.strip(".psmp")+'.png')
            print fin
            print fout
            print lcout
            raw_input()

def addtolightcurve(lightcurvef ile,saveloc,column_name,filt,mjd,flux,fluxerr,fakemag):
    #if not os.path.exists(saveloc):
    #    os.makedirs(saveloc)
    savefile = open(saveloc+'/'+lightcurvefile.split('/')[-1],'w')
    origfile = open(lightcurvefile,'r')
    lines = origfile.readlines()
    mjd = np.array(mjd)
    filt = np.array(filt,dtype='str')
    flux = np.array(flux)
    fluxerr = np.array(fluxerr)
    fakemag = np.array(fakemag)
    #zp = np.array(zp)
    for line in lines:
        if line.split(' ')[0] == 'VARNAMES:':
            line = line.strip()+' FLUX_'+column_name.upper()+' FLUXERR_'+column_name.upper()+' FLUX_ZPT_'+column_name.upper()+' FIT_ZPT_'+column_name.upper()+'\n'
        elif line.split(' ')[0] == 'OBS:':
            id = int(line.split()[1])
            tmjd = float(line.split()[3])
            band = line.split()[4]
            ww = (mjd == tmjd) & (filt == band)
            if fluxerr[ww] > 0:
                line = line.strip()+' '+str(round(flux[ww][0],2))+' '+str(round(fluxerr[ww][0],2))+' '+str(round(fakemag[ww][0],3))+'\n'
            else:
                line = line.strip()+' -999 -999 -999\n'
        savefile.write(line)
        #print line
    savefile.close()
    origfile.close()
#addtolightcurve('testlc.dat','./testdats/','testcol',,[888,777,000,111],[8,8,8,8],[31.,31.,31.,31.])

       


if __name__ == "__main__":
    wraplightcurves("data/snfiles_ps.txt","/home/dscolnic/",
                    "/export/scratch0/ps1sn1/data/v10.0/GPC1v3/eventsv1/smpworkspace/PS_TEST1/np_data/g/",
                    "/export/scratch0/ps1sn1/data/v10.0/GPC1v3/eventsv1/smpworkspace/PS_TEST1/lightcurves/g",
                    filt='g')