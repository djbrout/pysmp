#source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh

#. /grid/fermiapp/products/common/etc/setups.sh
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups

setup jobsub_client


jobsub_submit -G des --memory=1500MB --disk=59GB --force=xrootd --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --expected-lifetime=12h -M --verbose --OS=SL6\
 --email-to=djbrout@gmail.com -N 50 file://runfermirun.csh
echo done
