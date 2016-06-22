source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh

. /grid/fermiapp/products/common/etc/setups.sh

setup jobsub_client

num_jobs = 1

jobsub_submit -G des --resource-provides=usage_model=DEDICATED -M\
     --verbose --OS=SL6 file:///data/des41.a/data/djbrout/pysmp/batch_submission_scripts/runfermirun.csh $i g\
      --log_file=/data/des41.a/data/djbrout/pysmp/logs/log_$i.log\
      --generate-email-summary -N $num_jobs
echo $num_jobs jobs submitted
done
