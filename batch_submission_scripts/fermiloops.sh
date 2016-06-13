source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh
. /grid/fermiapp/products/common/etc/setups.sh
setup jobsub_client

for i in {0..0}
do
    jobsub_submit -G des --resource-provides=usage_model=DEDICATED -M --verbose --OS=SL6 file:///data/des41.a/data/djbrout/pysmp/batch_submission_scripts/runfermirun.csh $i g --log_file=/pnfs/des/scratch/pysmp/logfiles/log_$i.log  
    echo job $i submitted
done
