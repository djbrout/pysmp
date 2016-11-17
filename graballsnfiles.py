import os
import sys

f = open('sntarfilesv2.txt','r').read()
files = f.split()

#sys.exit()

#filename = 'SN-S1_CCD01_v3.tar'

#NEED TO DO THIS EVERY DAY
#'kx509; grid-proxy-init; voms-proxy-init -rfc -noregen -voms des:/des/Role=Analysis; '
cntr = 0
for filename in files:
    #print filename
    #if not '-X1' in filename:
    #    continue
    #filename = filename.replace('v4','v5')
    if '.idx' in filename: continue
    cntr += 1
    #if not '-' in filename: continue
    #if not '-S2' in filename: continue

    #if cntr < 30: continue
    #if not 'SN-E1' in filename: continue
    if os.path.isfile("/pnfs/des/persistent/smp/v2/"+filename):
        #if os.stat("/pnfs/des/persistent/smp/v2/"+filename).st_size > 0.:
        print filename,'already exists'
        continue
    print 'globus copying',filename
    out = os.popen( 'globus-url-copy -nodcau -cred /tmp/x509up_u48121 '
                '-ss "/DC=org/DC=opensciencegrid/O=Open Science Grid/OU=Services/CN=dtn03-garchive.nersc.gov"'
                ' gsiftp://garchive.nersc.gov:2811/home/projects/dessn/diffim/FinalPhoto/v2/' +
                filename+' gsiftp://fndca1.fnal.gov:2811/des/persistent/smp/v2/'+filename).read()
    print out
    print filename

