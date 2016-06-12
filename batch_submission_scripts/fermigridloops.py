import os
os.system('source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh')
os.system('. /grid/fermiapp/products/common/etc/setups.sh')
os.system('setup jobsub_client')
for i in range(0,2):
	os.system('jobsub_submit -G des --resource-provides=usage_model=DEDICATED -M --OS=SL6 file://batch_submission_scripts/runfermirun.csh '+str(i)+' g')
	print i,'submitted'
