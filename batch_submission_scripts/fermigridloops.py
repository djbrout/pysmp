import os
for i in range(0,2):
	os.system('jobsub_submit -G des --resource-provides=usage_model=DEDICATED -M --OS=SL6 file://batch_submission_scripts/runfermirun.csh '+str(i)+' g')
	print i,'submitted'
