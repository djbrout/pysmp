import numpy as np



def getfakefiles(ccd,field,num):
    print 'data/'+field.lower()+'lightcurves.txt'
    allfieldfakes = open('data/'+field.lower()+'lightcurves.txt').readlines()
    goodfakes = []
    for fk in allfieldfakes:

        if len(goodfakes) == num:
            break

        stopreading = False
        for l in open(fk.strip()).readlines():
            if 'g_'+ccd+'.LIST' in l:
                print fk
                goodfakes.append(fk)
                stopreading = True
            if stopreading:
                break


    return goodfakes

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