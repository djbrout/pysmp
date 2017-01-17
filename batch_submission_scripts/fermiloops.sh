#source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh

#. /grid/fermiapp/products/common/etc/setups.sh
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups

setup jobsub_client


jobsub_submit -G des --memory=2000MB --disk=59GB --cpu=4 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --expected-lifetime=8h -M --verbose --OS=SL6\
 --email-to=djbrout@gmail.com -N 1000 file://runfermirun.csh
echo done
