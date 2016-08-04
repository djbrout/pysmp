
import numpy as np
import os
import plotlightcurve as lc

def wraplightcurves(listfile,filedir,npzdir,lightcurveoutdir,filt=None):
    files = open(listfile).readlines()
    for fl in files:
        f = fl.split("/")[-1].strip()
        npzfile = os.path.join(npzdir,f.strip(".psmp")+'_'+filt+'_withSn.npz')
        inputfile = os.path.join(npzdir,f.strip(".psmp")+'_'+filt+'_mcmc_input.npz')
        if os.path.exists(npzfile) & os.path.exists(inputfile):
            fin = os.path.join(filedir,f)
            fout = os.path.join(filedir,f+'_dillon')
            lcout = os.path.join(lightcurveoutdir,f.strip(".psmp")+'.png')
            data = np.load(npzfile)
            input = np.load(inputfile)

            fitmag = 31.-2.5*np.log10(data['modelvec'])
            fitmagerr = -2.5*np.log10(data['modelvec'])+2.5*np.log10(data['modelvec']+data['modelvec_uncertainty'])
            fitflux = data['modelvec']
            fitfluxerr = data['modelvec_uncertainty']


            diffimmag = 31.-2.5*np.log10(input['diffim_flux'])
            diffimmagerr = -2.5 * np.log10(input['diffim_flux']) + 2.5 * np.log10(
                input['diffim_flux'] + input['diffim_fluxerr'])
            diffimflux = input['diffim_flux']
            diffimfluxerr = input['diffim_fluxerr']

            fitmag[fitmag > 50.] = 99; fitmagerr[fitmagerr > 50]=0.;fitmag[np.isnan(fitmag)] = 99; fitmagerr[np.isnan(fitmag)] = 0

            print diffimflux,fitflux
            print diffimmag,fitmag

            lc.pslightcurve(input['mjd'], fitmag, fitmagerr, diffimmag, diffimmagerr, fitflux, fitfluxerr,
                                                                                   diffimflux, diffimfluxerr,
                                                                                   filt, lcout,
                                                                                   title='Dillon Test (Preliminary)')
            #addtolightcurve(fin,fout,'DILLON_SMP','g',mjd,flux,fluxerr,usezpt,fitzpt)
            raw_input()

def addtolightcurve(lightcurvefile,saveloc,column_name,filt,mjd,flux,fluxerr,fakemag):
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
    wraplightcurves("data/snfiles_psdill.txt","/home/dscolnic/",
                    "/export/scratch0/ps1sn1/data/v10.0/GPC1v3/eventsv1/smpworkspace/PS_TEST1/np_data/g/",
                    "/export/scratch0/ps1sn1/data/v10.0/GPC1v3/eventsv1/smpworkspace/PS_TEST1/lightcurves/g",
                    filt='g')