


source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh
#source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
#setup jobsub_client
kx509
voms-proxy-init -rfc -noregen -valid 24:00 -voms des:/des/Role=Analysis

#alias sqs="jobsub_q -G des --role=DESGW --user=desgw"

# setup perl 5.18.1+6 || exit 134
setup Y2Nstack 1.0.6+18
setup diffimg gwdevel8
setup ftools v6.17
setup autoscan v3.1+0
setup easyaccess
setup extralibs 1.0
setup ifdhc
echo "EUPS setup complete"



#source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
#PYTHONPATH=$PYTHONPATH:/cvmfs/fermilab.opensciencegrid.org/products/common/prd/jobsub_client/v1_2_2/NULL:/cvmfs/fermilab.opensciencegrid.org/products/common/prd/pycurl/v7_16_4/Linux64bit-2-6-2-12/pycurl


#setup jobsub_client
echo 'done setting up jobsub_client'