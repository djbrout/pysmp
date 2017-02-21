import numpy as np
import os


def run(fakefiles):
    for ff in fakefiles:
        os.popen('python smptar.py -f r -r "." -s '+ff.strip())

def getfakefiles(ccd,field,num):
    #print 'data/'+field.lower()+'lightcurves.txt'
    allfieldfakes = open('data/'+field.lower()+'lightcurves.txt','r').readlines()
    goodfakes = []
    print allfieldfakes


    for fk in allfieldfakes:

        print "awk 'FNR > 50 { nextfile }; /g_"+ccd+".LIST/ { print FILENAME }' "+ fk.strip()
        out = os.popen("awk 'FNR > 50 { nextfile }; /g_"+ccd+".LIST/ { print FILENAME }' "+ fk.strip()).read()

        print out

        if len(out) > 0:
            goodfakes.append(fk)

        if len(goodfakes) == num:
            break

    return goodfakes


def ifdhtarballs(field,ccdnums):
    print ''
    # for ccd in ccdnums:
    #    print os.popen('ifdh cp -D --force=xrootd /pnfs/des/persistent/smp/v62/SN-'+field.upper()+'_CCD'+str(ccd)+'_v6.tar .').read()
    #    print os.popen('tar -xf SN-'+field.upper()+'_CCD'+str(ccd)+'_v6.tar').read()
    #
    # os.popen('ifdh log "Job $CLUSTER.$PROCESS has finished copying its input file at `date`"')

def setup(ccd,field,num):

    fakelist = getfakefiles(ccd,field,num)

    ccdnums = []
    for l in fakelist:
        ll = open(l).readlines()[80,120]
        for li in ll:
            imfile = li.split(' ')[13]
            ccdnums.append(int(imfile.split('/')[-1].split('_')[-1][:2]))


    ccdnums = np.unique(np.array(ccdnums,dtype='int'))

    print ccdnums

    ifdhtarballs(field,ccdnums)

    run(fakelist)


if __name__ == "__main__":
    import sys,getopt
    try:
        args = sys.argv[1:]
        opt, arg = getopt.getopt(
            args, "c:f:n",
            longopts=["ccd=","field=","num="])
    except getopt.GetoptError as err:
        print str(err)
        print "Error : incorrect option or missing argument."
        print __doc__
        sys.exit(1)


    ccd, field, num = None, None, 5



    for o, a in opt:
        if o in ["-c", "--ccd"]:
            ccd = a
        elif o in ["-f","--field"]:
            field = a

    print ccd,field,num
    setup(ccd,field,num)
