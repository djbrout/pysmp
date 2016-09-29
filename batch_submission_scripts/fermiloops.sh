source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh

. /grid/fermiapp/products/common/etc/setups.sh

setup jobsub_client


jobsub_submit -G des --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC -M --verbose --OS=SL6\
  --log_file=/pnfs/des/persistent/smp/logs/r1.log --generate-email-summary -N 2\
  file://$HOME/runfermirun.csh
echo 2 jobs submitted
echo done
