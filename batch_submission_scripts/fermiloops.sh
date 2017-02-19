#source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh

#. /grid/fermiapp/products/common/etc/setups.sh
#source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups

#setup jobsub_client
rm /tmp/*x509up_u`id -u`*
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
setup jobsub_client
jobsub_submit -G des --memory=1500MB --maxConcurrent=50 --disk=59GB --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --generate-email-summary --email-to=bdrizzle@yahoo.com --mail_never --expected-lifetime=15h -M --verbose --OS=SL6\
 -N 1800 file://runfermirun.csh
echo done
