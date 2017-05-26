import numpy as np
import os






if __name__ == "__main__":
    print 'started'
    #lclist = 'data/allspec.txt'
    #lclist = 'data/specHD.HEADER'
    lclist = 'data/speclist.txt'
    isreal = True

    resultsdir = '/project/projectdirs/des/djbrout/spec_v7/'

    savelcdir = resultsdir+'/SMP_SIMdeep_v3'
    fakes = False
    faketrueflux = False

    filts = ['g','r','i','z']

    import sys, getopt

    try:
        args = sys.argv[1:]

        opt, arg = getopt.getopt(
            args, "fd:rd:cd:cdf:b",
            longopts=["index=","lclist=", "resultsdir=", "savelcdir=","dofakes","faketrueflux"])

    except getopt.GetoptError as err:
        print "No command line arguments"

    index = None

    for o, a in opt:
        #print o
        if o in ["-h", "--help"]:
            print __doc__
            sys.exit(0)
        elif o in ["--lclist"]:
            lclist = a
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
    #print fakes
    #raw_input()

    if not os.path.exists(os.path.basename(savelcdir)):
        os.mkdir(os.path.basename(savelcdir))
    if not os.path.exists(savelcdir):
        os.mkdir(savelcdir)


    missingfile = resultsdir+'/missing.txt'

    #readme = open(savelcdir + '/' + savelcdir.split('/')[-1] + '.README', 'w')
    #readme.write(readmetext)
    #readme.close()
    #raw_input()
    #snlist = open(savelcdir + '/' + savelcdir.split('/')[-1] + '.LIST', 'w')



    #for i, filt in enumerate(filts):
    #sne = os.listdir(resultsdir + '/SNe')
    if '.HEADER' in lclist:
        lines = open(lclist,'r').readlines()
        sne = []
        for l in lines:
            if 'SN:' in l:
                if isreal:
                    sne.append('des_real_0'+l.split()[1]+'.dat')
                else:
                    sne.append('des_fake_0'+l.split()[1]+'.dat')

    else:
        sne = open(lclist,'r').readlines()


    if index is None:
        sne = sne
    else:
        sne = [sne[index]]


    # if fakes:
    #     import pandas as pd
    #
    #     dofakedata = pd.read_csv('data/doFake.out', delim_whitespace=True, header=0)
    #     # print dofakedata
    #     dofakeexpnum = dofakedata['EXPNUM'].values
    #     dofakemag2 = dofakedata['TRUEMAG'].values
    #     dofaketflux = dofakedata['TRUEFLUXCNT'].values
    #     dofakeid = dofakedata['FAKEID'].values
    #     dofakezpt = dofakemag2 + 2.5 * np.log10(dofaketflux)
    #     print 'done reading dofake'

    a = open(missingfile,'w')
    tsne = []
    for sn in sne[:]:
        snbad = False

        for i, filt in enumerate(filts):
            if 'starfits' in sn:
                continue

            #lcfile = sn.split('.')[0]+'_'+filt+'.dat'
            #print lcfile
            #raw_input()
            #if not filt is None:
            smpfile = resultsdir+'/lightcurves/'+sn.split('.')[0]+'_'+filt+'.smp'
            #else:
            #    continue
            #savelcfile = savelcdir+'/'+sn+'_smp.dat'
            if not os.path.exists(smpfile):
                print 'SMP RESULTS DO NOT EXIST FOR ',smpfile
                a.write('_'.join(smpfile.split('/')[-1].split('.')[0].split('_')[:-1])+' '+filt+' \n')
                #os.system('echo '+sn+' '+filt+' >> '+missingfile)
                snbad = True
                continue
        if not snbad:
            tsne.append(sn)
    sne = tsne
    a.close()
    print 'TOTAL SNE:',len(sne)
    print 'Written',missingfile
    # cntr = 0
    # for sn in sne[::-1]:
    #     mjd = []
    #     flux = []
    #     fluxerr = []
    #     zpt = []
    #     rmsaddin = []
    #     chi2 = []
    #     sky = []
    #     skyerr = []
    #     smpflag = []
    #     zptfile = []
    #     idobs = []
    #     cntr += 1
    #     for i, filt in enumerate(filts):
    #         #for sn in sne[:]:
    #         #cntr += 1
    #         #print sn
    #         if cntr > 10000:
    #             continue
    #         if 'starfits' in sn:
    #             continue
    #
    #         lcfile = lcdir+'/'+sn+'.dat'
    #         #if not filt is None:
    #         smpfile = resultsdir+'/lightcurves/'+sn+'_'+filt+'.smp'
    #         #else:
    #         #    continue
    #         savelcfile = savelcdir+'/'+sn+'_smp.dat'
    #         if not os.path.exists(smpfile):
    #             print 'SMP RESULTS DO NOT EXIST FOR ',smpfile
    #             os.system('echo '+sn+' '+filt+' >> '+missingfile)
    #             continue
    #         #print lcfile
    #         pkmjd = open(lcfile).readlines()[10].split()[1]
    #         #print pkmjd
    #         #raw_input()
    #         inplace = False
    #         #if i > 0: inplace = True
    #         sndata = dt.readcol(smpfile,1,2)
    #
    #         mjd.extend(sndata['MJD'])
    #         flux.extend(sndata['FLUX'])
    #         fluxerr.extend(sndata['FLUXERR'])
    #         zpt.extend(sndata['ZPT'])
    #         rmsaddin.extend(sndata['RMSADDIN'])
    #         chi2.extend(sndata['CHI2'])
    #         sky.extend(sndata['SKY'])
    #         skyerr.extend(sndata['SKYERR'])
    #         smpflag.extend(sndata['SMP_FLAG'])
    #         zptfile.extend(sndata['ZPTFILE'])
    #         idobs.extend(sndata['ID_OBS'])
    #
    #     successful = addtolightcurve(lcfile,savelcfile,mjd,flux,fluxerr,
    #                  zpt, rmsaddin,
    #                  chi2,sky,skyerr,smpflag,zptfile,
    #                  idobs,pkmjd, dofakes=fakes, saveinplace=False,faketrueflux=faketrueflux)
    #
    #     print int(cntr),'SAVED SUCCESSFULLY',savelcfile,'\n'