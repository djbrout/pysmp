import os

filename = 'SN-S1_CCD01_v3.tar'
'kx509; grid-proxy-init; voms-proxy-init -rfc -noregen -voms des:/des/Role=Analysis; '
out = os.popen( 'globus-url-copy -nodcau -cred /tmp/x509up_u48121 '
                '-ss "/DC=org/DC=opensciencegrid/O=Open Science Grid/OU=Services/CN=dtn03-garchive.nersc.gov"'
                ' gsiftp://garchive.nersc.gov:2811/home/projects/dessn/diffim/FinalPhoto/v3/' +
                filename+' gsiftp://fndca1.fnal.gov:2811/des/scratch/pysmp/v3/'+filename).read()
print out