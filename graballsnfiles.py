import os
import sys

f = open('sntarfilesv4.txt','r').read()
files = f.split()

#sys.exit()

#filename = 'SN-S1_CCD01_v3.tar'

#NEED TO DO THIS EVERY DAY
#'kx509; grid-proxy-init; voms-proxy-init -rfc -noregen -voms des:/des/Role=Analysis; '
cntr = 0
for filename in files[::-1]:
    #print filename
    if not '-X1' in filename:
        continue
    if '.idx' in filename: continue
    cntr += 1
    if cntr > 40: continue
    #if not 'SN-E1' in filename: continue
    if os.path.isfile("/pnfs/des/persistent/smp/v4/"+filename):
        if os.stat("/pnfs/des/persistent/smp/v4/"+filename).st_size > 0.:
            print filename,'already exists'
            continue
    print 'globus copying',filename
    out = os.popen( 'globus-url-copy -nodcau -cred /tmp/x509up_u48121 '
                '-ss "/DC=org/DC=opensciencegrid/O=Open Science Grid/OU=Services/CN=dtn03-garchive.nersc.gov"'
                ' gsiftp://garchive.nersc.gov:2811/home/projects/dessn/diffim/FinalPhoto/v4/' +
                filename+' gsiftp://fndca1.fnal.gov:2811/des/persistent/smp/v4/'+filename).read()
    print out
    print filename

