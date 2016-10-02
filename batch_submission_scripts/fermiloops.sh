source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh

. /grid/fermiapp/products/common/etc/setups.sh

setup jobsub_client


jobsub_submit -G des --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC -M --verbose --OS=SL6\
  --generate-email-summary -N 10 file://runfermirun.csh
echo 10 jobs submitted
echo done
